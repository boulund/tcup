#!/usr/bin/env python3.5
# Wrapper for TCUP
# Fredrik Boulund 2016

from sys import argv, exit
import os
import logging
import argparse
import shlex
import subprocess
import platform

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("TCUP-wrapper")


def parse_args():
    """
    Parse command line arguments.
    """
    
    desc = """TCUP wrapper; align peptides to reference databases in parallel
    and run TCUP on alignment results. Fredrik Boulund 2016"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FASTA", 
            help="FASTA file with peptides from tandem MS.")
    parser.add_argument("GENOME_DB",
            help="Path to reference bacterial genome db (FASTA or blastdb format depending on OS).")
    parser.add_argument("RESISTANCE_DB",
            help="Path to antibiotic resistance gene db (FASTA or blastdb format depending on OS.")
    
    taxonomic_composition_parser = parser.add_argument_group("Taxonomic composition")
    taxonomic_composition_parser.add_argument("-t", "--taxref-db", dest="taxref_db",
            required=True,
            help="Path to taxref db (sqlite3).")
    taxonomic_composition_parser.add_argument("-a", "--annotation-db", dest="annotation_db",
            required=True,
            help="Path to annotation db (sqlite3).")

    ar_parser = parser.add_argument_group("Antibiotic resistance")
    ar_parser.add_argument("-r", "--resistance-db", dest="resistance_db",
            required=True,
            help="Path to resistance db (sqlite3).")
    
    if len(argv) < 2:
        parser.print_help()
        exit()

    options = parser.parse_args()

    return options
        


def run_blast(fasta, db, outfilename, task="blastx"):
    """
    Run BLAST using FASTA file as query to DB.
    Assumes running on Windows.
    """

    cmd_string = "{task}.exe -query {query} -db {db} -out {outfile} -outfmt 6"
    blast_cmd = shlex.split(cmd_string.format(task=task, 
        query=fasta,
        db=db,
        outfile=outfilename))
    logger.debug("Running %s", " ".join(blast_cmd))
    blast_process = subprocess.Popen(blast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return blast_process


def run_blat(fasta, db, outfilename, query_type, target_type, high_sens=False):
    """
    Run BLAT using FASTA as query to DB.
    """

    if high_sens:
        cmd_string = "blat {db} {fasta} -q={query_type} -t={target_type} -minIdentity=90 -stepSize=5 -tileSize=5 -minScore=10 -out=blast8 {outfile}"
    else:
        cmd_string = "blat {db} {fasta} -q={query_type} -t={target_type} -minIdentity=90 -out=blast8 {outfile}"
    blat_cmd = shlex.split(cmd_string.format(fasta=fasta,
        db=db,
        outfile=outfilename,
        query_type=query_type,
        target_type=target_type))
    logger.debug("Running %s", " ".join(blat_cmd))
    blat_process = subprocess.Popen(blat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return blat_process


def run_ar_detection(blast8, resistance_db, outfilename):
    """
    Run TCUP antibiotic resistance detection.
    """

    cmd_string = "antibiotic_resistance -r {resdb} {blast8} -o {outfile}"
    ar_cmd = shlex.split(cmd_string.format(resdb=resistance_db,
        blast8=blast8,
        outfile=outfilename))
    logger.debug("Running %s", " ".join(ar_cmd))
    ar_process = subprocess.Popen(ar_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return ar_process


def run_taxonomic_composition(blast8, taxref_db, annotation_db, outfilename):
    """
    Run TCUP taxonomic composition estimation.
    """

    cmd_string = "taxonomic_composition --taxref-db {taxref} --annotation-db {annotation} --output {outfile}.txt --write-xlsx {outfile}.xlsx {blast8}"
    tax_cmd = shlex.split(cmd_string.format(taxref=taxref_db,
        annotation=annotation_db,
        outfile=outfilename,
        blast8=blast8))
    logger.debug("Running %s", " ".join(tax_cmd))
    tax_process = subprocess.Popen(tax_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return tax_process




def main():
    """
    Run the entire TCUP workflow.
    """

    options = parse_args()

    genome_mapping_output = os.path.basename(options.FASTA)+".genomes.blast8"
    ar_mapping_output = os.path.basename(options.FASTA)+".ar.blast8"
    ar_output = os.path.basename(options.FASTA)+".antibiotic_resistance.txt"
    taxcomp_output = os.path.basename(options.FASTA)+".taxonomic_composition"


    print("Running genome mapping...")
    print("Running antibiotic resistance gene mapping...")
    if platform.system().startswith("Linux"):
        genome_mapping_process = run_blat(options.FASTA, 
                options.GENOME_DB, 
                genome_mapping_output, 
                query_type="prot",
                target_type="dnax",
                high_sens=True)
        ar_mapping_process = run_blat(options.FASTA, 
                options.RESISTANCE_DB, 
                ar_mapping_output,
                query_type="prot",
                target_type="prot")
    elif platform.system().startswith("Windows"):
        genome_mapping_process = run_blast(options.FASTA, 
                options.GENOME_DB, 
                genome_mapping_output,
                task="tblastn")
        ar_mapping_process = run_blast(options.FASTA, 
                options.RESISTANCE_DB, 
                ar_mapping_output,
                task="blastp")

    ar_mapping_out = ar_mapping_process.communicate()
    print("Antibiotic resistance gene mapping completed.")
    logger.debug("Finished AR mapping: %s", ar_mapping_out)
    ar_process = run_ar_detection(ar_mapping_output, options.resistance_db, ar_output)

    ar_process_out = ar_process.communicate()
    if ar_process.returncode == 2:
        print(ar_process_out[1].decode("utf-8").split("\n")[-2].rsplit(":",1)[1])
    print("Antibiotic resistance detection (TCUP) completed.")
    logger.debug("Finished antibiotic resistance detection: %s", ar_process_out)

    genome_mapping_out = genome_mapping_process.communicate()
    print("Genome mapping completed.")
    logger.debug("Finished genome mapping: %s", genome_mapping_out)

    taxcomp_process = run_taxonomic_composition(genome_mapping_output, options.taxref_db, options.annotation_db, taxcomp_output)
    taxcomp_process_out = taxcomp_process.communicate()
    print("Taxonomic composition estimation (TCUP) completed.")
    logger.debug("Finished taxonomic composition estimation: %s", taxcomp_process_out)


if __name__ == "__main__":
    main()

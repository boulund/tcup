#!/usr/bin/env python3.5
# encoding: utf-8
# Fredrik Boulund (c) 2015
# Proteotyping

from sys import argv, exit
from collections import defaultdict, OrderedDict, Counter
import shutil
import sqlite3
import gzip
import time
import os
import argparse
import logging
import re

try:
    from utils import find_files, grouper, existing_file
    from proteotyping_db import NCBITaxa_mod as NCBITaxa
except ImportError:
    from proteotyping.utils import find_files, grouper, existing_file
    from proteotyping.proteotyping_db import NCBITaxa_mod as NCBITaxa


def parse_commandline(argv):
    """Parse commandline arguments"""

    desc = """Proteotyping pipeline. (c) Fredrik Boulund 2015."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("FILE", nargs="+",
            type=existing_file,
            help="BLAT output file.")
    parser.add_argument("--proteodb", dest="proteodb", metavar="DB", type=str,
            default="proteodb.sql",
            help="Path to proteotyping sqlite3 db [%(default)s]")
    parser.add_argument("--sample-db", dest="sample_db", action="store_true",
            default=False,
            help="The supplied 'BLAT output file' is really a processed sample db with proteotyping results [%(default)s]")
    parser.add_argument("--taxonomic-rank", dest="taxonomic_rank", metavar="LVL", type=str,
            choices=["no rank", "subspecies", "species", "genus", "family", "order", "class", "phylum", "superkingdom"],
            default="family",
            help="Set the taxonomic level on which hits are grouped [%(default)s].")
    parser.add_argument("--min-matches", dest="min_matches", metavar="L", type=int,
            default=6,
            help="Minimum peptide matches (i.e. peptide length) [%(default)s].")
    parser.add_argument("-i", "--min-identity", dest="min_identity", metavar="I", type=float,
            default=90,
            help="Filter out hits with less than or equal to this percentage identity [%(default)s].")
    parser.add_argument("--min-coverage", dest="min_coverage", metavar="C", type=float,
            default=1,
            help="Proportion of fragment covered in alignment [%(default)s].")
    parser.add_argument("--max-pid-diff", dest="max_pid_diff", type=float, metavar="D",
            default=5.0,
            help="Maximum identity difference between highest and lowest hit for each peptide. Floating point between 0.0-100.0 [%(default)s].")
    parser.add_argument("--output", dest="output",
            default="",
            help="Write results to this filename [results/FILE.results].")


    devoptions = parser.add_argument_group("Developer options", "Voids warranty ;)")
    devoptions.add_argument("--loglevel", choices=["INFO", "DEBUG", "VERBOSE"],
            default="DEBUG",
            help="Set logging level [%(default)s].")
    devoptions.add_argument("--logfile", dest="logfile",
            default="proteotyping.log",
            help="Filename for log output [%(default)s].")
    devoptions.add_argument("--numCPUs", dest="numCPUs", type=int,
            default=16,
            help="Number of CPUs to utilize in parallel regions [%(default)s].")
    devoptions.add_argument("--leave-out", metavar="HEADERS", dest="leave_out",
            default="",
            help="Disregard any hits to sequence with HEADER when parsing and filtering blast8 output. Can be a list of comma separated FASTA headers (no spaces).")
    devoptions.add_argument("--blacklist", metavar="FILE", dest="blacklist",
            default=None,
            type=existing_file,
            help="File with sequence headers to blacklist (i.e. to ignore when parsing blast8 output).")

    if len(argv) < 2:
        parser.print_help()
        exit()

    options = parser.parse_args(argv[1:])
    logger = logging.getLogger()
    if options.loglevel == "VERBOSE":
        logger.setLevel(0)
    else:
        logger.setLevel(options.loglevel)
    fh = logging.FileHandler(options.logfile)
    ch = logging.StreamHandler()
    file_formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    console_formatter = logging.Formatter("%(levelname)s: %(message)s")
    fh.setFormatter(file_formatter)
    ch.setFormatter(console_formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)

    logging.info("----------========= LOGGING STARTED {} ========---------".format(time.strftime("%c")))

    return options


def prepare_blacklist(blacklist, additional=""):
    """ 
    Parses list of FASTA headers to blacklist.

    Assumes FASTA headers are split on first space character.
    
    :param filename:  filename of single column textfile with FASTA 
            headers to blacklist.
    :param additional:  string with comma separated headers to 
            blacklist.
    :return:  set of blacklisted FASTA headers.
    """
    blacklisted_seqs = set()
    if blacklist:
        with open(filename) as f:
            for line in f:
                blacklisted_seqs.add(line.strip())
    if additional:
        blacklisted_seqs.add(set(additional.split(",")))
    if blacklisted_seqs:
        logging.debug("Blacklisted sequences:\n%s", blacklisted_seqs)
    return blacklisted_seqs


def parse_blat_output(filename, min_identity, min_matches, 
            min_coverage, max_pid_diff, blacklisted_seqs):
    """
    Parses blat output.
    
    Filters out hits based on identity, fragment_length, 
    and skips blacklisted sequences.
    Columns in BLAT blast8 tabular format:
        0   Query id
        1   Subject id
        2   % identity
        3   alignment length
        4   mismatches
        5   gap openings
        6   q. start
        7   q. end
        8   s. start
        9   s. end
        10  e-value
        11  bit score
    """

    logging.info("Parsing and filtering hits from '{}'...".format(filename))

    hitlists = defaultdict(list)
    num_discarded = 0
    with open(filename) as blast8:
        for total_count, hit in enumerate(map(str.split, blast8), start=1):
            if hit[1] in blacklisted_seqs:
                logging.debug("Ignoring fragment %s hit to blacklisted %s", hit[0], hit[1])
                num_discarded += 1
                continue
            peptide, peptide_length = hit[0].rsplit("_", 1)
            alignment_length = int(hit[3])
            peptide_coverage = alignment_length / int(peptide_length)
            passing = (float(hit[2]) >= min_identity and 
                    alignment_length >= min_matches and 
                    peptide_coverage >= min_coverage)
            if passing:
                hitlists[hit[0]].append((hit[1], hit[2], hit[3], hit[8], hit[9]))
            else:
                num_discarded += 1
        logging.info("Parsed %s hits.", total_count)
        logging.info("Discarded %s hits based on primary criteria.", num_discarded)
        num_remain_pep = 0
        num_remain_hits = 0
        for peptide, hitlist in hitlists.items():
            max_pid = max(map(float, (h[1] for h in hitlist)))
            filtered = [h for h in hitlist if float(h[1]) >= max_pid - max_pid_diff]
            if len(filtered)>0:
                num_remain_pep += 1
                for hit in hitlist:
                    yield peptide, hit[0], float(hit[1]), int(hit[2]), int(hit[3]), int(hit[4])
                    num_remain_hits += 1
        logging.info("%s hits for %s peptides remain after relative filtering.", num_remain_hits, num_remain_pep)


class Proteotyping_DB_wrapper():
    """ 
    Wrapper for sqlite3 database.
    """

    def __init__(self, sample_db, ref_dbfile=None):
        self.dbfile = sample_db
        if ref_dbfile and os.path.isfile(self.dbfile):
            sleep_time = 5
            logging.warning("About to delete pre-existing sample DB: %s", self.dbfile)
            logging.warning("Press Ctrl+c to cancel within %i sec...", sleep_time)
            time.sleep(sleep_time)
            os.remove(self.dbfile)
            logging.debug("Copying reference DB (%s) to create new sample DB (%s)", ref_dbfile, self.dbfile)
            shutil.copy(ref_dbfile, self.dbfile)
        elif ref_dbfile:
            logging.error("Cannot connect to non-existent DB file: %s", ref_dbfile)
            exit(1)

        self.db = sqlite3.connect(self.dbfile)
        logging.info("Connected to sample DB %s", self.dbfile)
        try:
            logging.info("DB version: %s", self.db.execute("SELECT * FROM version").fetchall()[0])
        except sqlite3.OperationalError as msg:
            logging.error("Incorrect or missing DB %s: %s", self.dbfile, msg)
            exit(1)



    @staticmethod
    def lowest_common_ancestor(lineages, common_lineage=False):
        """
        Determine lowest common ancestors between a set of lineages.
        
        :param lineages: multiple lineages represented lists or tuples of ints. 
                The root node should be located at the end of each list. 
        :param common_lineage: boolean telling whether to return entire common
                lineage or only lowest common ancestor.
        :return: lowest common ancestor (list).
        """
       # if not (type(lineages[0]) is list or type(lineages[0]) is tuple):
       #     raise TypeError("lineages must be nested iterator"
       #                     "(e.g. list of lists or tuple of tuples)")
        
        occurrences = OrderedDict()
        for lineage in lineages:
            for taxid in lineage:
                if taxid in occurrences:
                    occurrences[taxid] += 1
                else:
                    occurrences[taxid] = 1
        common_ancestors = [taxid for taxid, occurrence in occurrences.items() 
                            if occurrence == len(lineages)]
        if common_lineage:
            return common_ancestors
        elif common_ancestors:  # Only return LCA if there is actually is one.
            return [common_ancestors[0]]
        else:
            return []


    def insert_blat_hits_into_db(self, blat_output, hits_per_chunk=50000):
        """
        Insert peptide hits from BLAT into db.
        """
        for values in grouper(hits_per_chunk, blat_output):
            self.db.executemany("INSERT INTO peptides VALUES (?, ?, ?, ?, ?, ?)", values)
        self.db.commit()

    
    def determine_discriminative_ranks(self):
        """
        Find the discriminative rank of each unique peptide in the db.
        """

        for peptide in self.db.execute("SELECT DISTINCT peptide FROM peptides"):
            cmd = """ SELECT track FROM species
              JOIN refseqs ON refseqs.taxid = species.taxid
              JOIN peptides ON peptides.target = refseqs.header AND peptides.peptide = ?
            """
            query = self.db.execute(cmd, peptide)
            tracks = [map(int, t[0].split(",")) for t in query.fetchall()]
            lca = self.lowest_common_ancestor(tracks)
            if lca[0] != 131567:  # taxid=131567 is "cellular organisms"
                # TODO: verbose
                #logging.debug("Peptide %s is discriminative at rank '%s' for %s", peptide[0], rank, spname)
                self.db.execute("INSERT INTO discriminative VALUES (?, ?)", (peptide[0], lca[0]))
        self.db.commit()


    def get_discriminative_at_rank(self, rank):
        """
        Retrieve all discriminative peptides at the specified rank.

        Stores a hard coded rank hierarchy based on the ranks seen 
        in NCBI Taxonomy. Note the problematic "no rank" assignment
        that can be encountered at different levels in the taxonomy, 
        here incorrectly positioned to always be at the lowest rank.
        """

        rank_hierarchy = """
        superkingdom
        kingdom
        subkingdom
        superphylum
        phylum
        subphylum
        class
        superclass
        subclass
        infraclass
        superorder
        order
        suborder
        infraorder
        parvorder
        superfamily
        family
        subfamily
        tribe
        subtribe
        genus
        subgenus
        species group
        species subgroup
        species
        subspecies
        varietas
        forma
        """.split()
        rank_hierarchy.append("no rank")
        ranks = {r: n for n, r in enumerate(rank_hierarchy)}
        rank_set = rank_hierarchy[ranks[rank]:]

        cmd = """SELECT peptide, rank, spname FROM discriminative
          JOIN species ON species.taxid = discriminative.taxid
          WHERE rank IN ({})
          ORDER BY rank""".format(",".join("?"*len(rank_set)))
        result = self.db.execute(cmd, rank_set).fetchall()
        return result

def print_discriminative_peptides(discriminative):
    """
    Print name and assignment of discriminative peptides.
    """
    print("-" * 60)
    print("Discriminative peptides")
    print("-" * 60)
    print("{:<20} {:<15} {:<30}".format("Peptide", "Rank", "Description"))
    for d in discriminative: 
        print("{:<20} {:<15} {:<30}".format(*d))


def print_peptides_per_spname(discriminative):
    """
    Print number of discriminative peptides per spname.
    """
    species_counts = Counter((species, rank) for pep, rank, species in discriminative)
    print("-" * 60)
    print("Discriminative peptides per spname")
    print("-" * 60)
    print("{:<6} {:<20} {:<40}".format("Count", "Rank", "Description"))
    for species, count in species_counts.most_common():
        print("{:<6} {:<20} {:<40}".format(count, species[1], species[0]))


def get_results_from_existing_db(options):
    """
    Retrieve results from an existing sample db.
    """
    for sample_db in options.FILE:
        refdb = Proteotyping_DB_wrapper(sample_db)
        disc = refdb.get_discriminative_at_rank(options.taxonomic_rank)

        print(refdb.dbfile)
        print_discriminative_peptides(disc)
        print_peptides_per_spname(disc)


def main(options):
    """
    Main function that runs the complete pipeline logic.
    """
    blacklisted_seqs = prepare_blacklist(options.blacklist, options.leave_out)

    for blat_file in options.FILE:
        refdb = Proteotyping_DB_wrapper(blat_file+".sqlite3", options.proteodb) 
        #refdb = Proteotyping_DB_wrapper(blat_file+".sqlite3") 
        blat_parser = parse_blat_output(blat_file, 
                options.min_identity, 
                options.min_matches, 
                options.min_coverage,
                options.max_pid_diff, 
                blacklisted_seqs)
        refdb.insert_blat_hits_into_db(blat_parser)
        refdb.determine_discriminative_ranks()
        disc = refdb.get_discriminative_at_rank(options.taxonomic_rank)

        print_discriminative_peptides(disc)
        print_peptides_per_spname(disc)


if __name__ == "__main__":

    options = parse_commandline(argv)

    if options.sample_db:
        get_results_from_existing_db(options)
    else:
        main(options)

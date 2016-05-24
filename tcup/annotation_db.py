#!/usr/bin/env python3.5
# encoding: utf-8
# Prepare an annotation database
# (c) Fredrik Boulund 2015

import sys
import os
import sqlite3
import argparse
import time
import logging

try: 
    from tcup.utils import read_fasta, find_files, grouper, existing_file
    from tcup.taxref_db import NCBITaxa_mod as Taxref_DB_wrapper
except ImportError:
    from utils import read_fasta, find_files, grouper, existing_file
    from taxref_db import NCBITaxa_mod as Taxref_DB_wrapper


def parse_args(argv):
    """
    Parse command line arguments.
    """

    desc = """Create an annotation database.
    (c) Fredrik Boulund 2015.
    """
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("taxref", type=existing_file, 
            help="Path to existing taxref DB.")
    parser.add_argument("annotation_dirs", nargs="+",
            help="Path to root dir(s) containing annotation files (*.gff).")
    parser.add_argument("--db-filename", dest="db_filename",
            default="annotation_db.sqlite3",
            help="Filename to save annotation DB to [%(default)s].")
    parser.add_argument("--glob-pattern-gff", dest="glob_pattern_gff",
            default="*.gff",
            help="Glob pattern for gff files [%(default)s].")
    parser.add_argument("--loglevel", choices=["INFO", "DEBUG"], 
            default="DEBUG", 
            help="Set logging level [%(default)s].")
    parser.add_argument("--logfile", 
            default=False,
            help="Log to file instead of STDOUT.")

    if len(argv) < 2:
        parser.print_help()
        exit()

    options = parser.parse_args()
    
    if options.logfile:
        logging.basicConfig(level=options.loglevel, filename=options.logfile)
    else:
        logging.basicConfig(level=options.loglevel)

    return options


def parse_gene_info(filename):
    """
    Parse tax_id, GeneID, Symbol, description from NCBI Gene (gene_info).

    Assumes that gene_info is structured as such:
    0  tax_id
    1  GeneID
    2  Symbol 
    3  LocusTag 
    4  Synonyms 
    5  dbXrefs 
    6  chromosome 
    7  map_location 
    8  description 
    9  type_of_gene 
    10 Symbol_from_nomenclature_authority 
    11 Full_name_from_nomenclature_authority 
    12 Nomenclature_status 
    13 Other_designations 
    14 Modification_date 
    (tab is used as a separator, pound sign - start of a comment)
    
    :param filename:  path to NCBI gene_info file.
    """
    with open(filename) as f:
        logging.info("Parsing %s", filename)
        f.readline() # Skip the header line
        for line in f:
            data = line.split("\t")
            taxid = int(data[0])
            geneid = int(data[1])
            symbol = data[2]
            description = data[8]
            yield (taxid, geneid, symbol, description)


def parse_gff_attributes(attributes_string):
    """
    Parse attributes field from GFF files.
    """

    attributes = attributes_string.split(";")
    return {k: v for k, v in (a.split("=") for a in attributes)}


def parse_gff(filename):
    """
    Parse gene annotations from GFF file.

    Parses GFF version 3 with the following columns:
    1   sequence
    2   source
    3   feature
    4   start
    5   end
    6   score
    7   strand
    8   phase
    9   attributes

    :param filename:  path to gff file.
    :return:  (sequence, start, end, attributes)
    """

    with open(filename) as f:
        logging.info("Parsing %s...", filename)
        line = f.readline()
        if not line.startswith("##gff-version 3"):
            logging.error("Parse Error: %s", filename)
            logging.error("Offending line: %s", line)
            raise Exception("Parse error: wrong gff version or not gff file: %s" % filename)
        while line.startswith("#"):
            line = f.readline()
        while line:
            try:
                sequence, source, feature, start, end, score, strand, phase, attributes_string = line.strip().split("\t")
                if feature == "CDS":
                    try:
                        product = attributes_string.split("product=")[1].split(";", 1)[0]
                    except IndexError:
                        logging.warning("Found no product info in attributes string: %s", attributes_string)
                        logging.warning("Entire line GFF line: %s", line)
                        product = "??"
                    yield (sequence, start, end, product, attributes_string)
            except ValueError:
                if line.startswith("#"):
                    # logging.debug("Skipping comment line: %s", line) # Too verbose!
                    pass
                else:
                    logging.error("Couldn't parse gff, the offending line was:\n%s", line)
            line = f.readline()


def parse_annotations(taxref_db, annotations_dir, pattern):
    """
    Recursively find gff files in refdir.

    """
    for gff_file in find_files(annotations_dir, pattern):
        for annotation_info in parse_gff(gff_file):
            try:
                header = taxref_db.find_refseq_header(annotation_info[0])
            except KeyError:
                continue
            yield (header, *annotation_info[1:])


class Annotation_DB_wrapper():
    """
    Wrapper for annotation DB.
    """

    def __init__(self, dbfile):
        if os.path.isfile(dbfile):
            logging.warning("About to overwrite existing DB: %s", dbfile)
            time.sleep(1)
            os.remove(dbfile)
        self.dbfile = dbfile
        self.db = sqlite3.connect(dbfile)

        self.db.execute("CREATE TABLE annotations (header TEXT, start INT, end INT, product TEXT, features TEXT)")

    def insert_annotations(self, annotation_data):
        """
        Insert annotation data into the DB.

        :param annotation_data:  Iterator of tuples with annotation info.
                (header, start, end, features)
        """

        self.db.executemany("INSERT INTO annotations VALUES (?,?,?,?,?)", annotation_data)
        records = self.db.execute("SELECT Count(*) FROM annotations").fetchone()[0]
        logging.info("Inserted %s annotation records.", records)
        self.db.commit()

    def create_indexes(self):
        """
        Create indexes for start, end, and a double index start_end.
        """
        self.db.execute("CREATE INDEX i_starts ON annotations(start)")
        self.db.execute("CREATE INDEX i_ends ON annotations(end)")
        self.db.execute("CREATE INDEX i_starts_ends ON annotations(start, end)")





def main():
    """
    Main function.
    """
    options = parse_args(sys.argv)

    taxref_db = Taxref_DB_wrapper(options.taxref)
    annotation_db = Annotation_DB_wrapper(options.db_filename)
    for annotations_dir in options.annotation_dirs:
        annotation_generator = parse_annotations(taxref_db, annotations_dir, options.glob_pattern_gff)
        annotation_db.insert_annotations(annotation_generator)
    annotation_db.create_indexes()


if __name__ == "__main__":
    main()

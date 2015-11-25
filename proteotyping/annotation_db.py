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
    from proteotyping.utils import read_fasta, find_files, grouper, existing_file
    from proteotyping.proteotyping_db import NCBITaxa_mod as Proteotyping_DB_wrapper
except ImportError:
    from utils import read_fasta, find_files, grouper, existing_file
    from proteotyping_db import NCBITaxa_mod as Proteotyping_DB_wrapper


def parse_args(argv):
    """
    Parse command line arguments.
    """

    desc = """Create an annotation database.
    (c) Fredrik Boulund 2015.
    """
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("gene_info",
            help="Path to NCBI Gene 'gene_info' file.")
    parser.add_argument("proteodb", type=existing_file,
            help="Path to existing proteotyping DB.")
    parser.add_argument("annotation_dirs", nargs="+",
            help="Path to root dir(s) containing annotation files (*.gff).")
    parser.add_argument("--db-filename", dest="db_filename",
            default="annotationdb.sql",
            help="Filename to save annotation DB to [%(default)s].")
    parser.add_argument("--glob-pattern-gff", dest="glob_pattern_gff",
            default="*.gff",
            help="Glob pattern for gff files [%(default)s].")
    parser.add_argument("--db-comment",
            help="A comment to be added to the database's version table")
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
            raise Exception("Parse error: wrong gff version or not gff file.")
        while line.startswith("#"):
            line = f.readline()
        while line:
            try:
                sequence, source, feature, start, end, score, strand, phase, attributes_string = line.split("\t")
                if feature == "gene" or feature == "protein":
                    yield (sequence, start, end, attributes_string)
            except ValueError:
                if line.startswith("#"):
                    logging.debug("Passing comment line: %s", line)
                    pass
                else:
                    logging.error("Couldn't parse gff, the offending line was:\n%s", line)
            line = f.readline()


def parse_annotations(proteodb, annotations_dir, pattern):
    """
    Recursively find gff files in refdir.

    """
    for gff_file in find_files(annotations_dir, pattern):
        for annotation_info in parse_gff(gff_file):
            header = proteodb.find_refseq_header(annotation_info[0])
            yield (header, *annotation_info[1:])


class Annotation_DB_wrapper():
    """
    Wrapper for annotation DB.
    """

    def __init__(self, dbfile, comment=None, create_db=False):
        if create_db:
            if os.path.isfile(dbfile):
                logging.warning("About to overwrite existing DB: %s", dbfile)
                time.sleep(1)
                os.remove(dbfile)
            self.dbfile = dbfile
            self.db = sqlite3.connect(dbfile)

            self.db.execute("DROP TABLE IF EXISTS annotations")
            self.db.execute("DROP TABLE IF EXISTS gene_info")
            self.db.execute("DROP TABLE IF EXISTS version")
            self.db.execute("CREATE TABLE annotations (header TEXT, start INT, end INT, attributes)")
            self.db.execute("CREATE TABLE gene_info (taxid INT, geneid INT, symbol TEXT, description TEXT)")
            self.db.execute("CREATE TABLE version (comment TEXT)")
            if comment:
                self.db.execute("INSERT INTO version VALUES (?)", (comment,))
                self.db.commit()
        else:
            if os.path.isfile(dbfile):
                self.dbfile = dbfile
                self.db = sqlite3.connect(dbfile)
                tables = self.db.execute("SELECT * FROM sqlite_master").fetchall()
                logging.debug("Connected to annotation db %s with tables: %s", dbfile, tables)
            else:
                logging.error("Found no annotation db at: %s", dbfile)
                exit(1)

    def insert_annotations(self, annotation_data):
        """
        Insert annotation data into the DB.

        :param annotation_data:  Iterator of tuples with annotation info.
                (header, start, end, attributes)
        """

        self.db.executemany("INSERT INTO annotations VALUES (?,?,?,?)", annotation_data)
        records = self.db.execute("SELECT Count(*) FROM annotations").fetchone()[0]
        logging.info("Table annotations now contains %s records.", records)
        # TODO: Create index on start and end columns
        self.db.commit()

    def insert_gene_info(self, gene_info):
        """
        Insert gene_info into the DB.
        
        :param gene_info:  Iterator with tuples of gene_info,
                (taxid, geneid, symbol, description)
        """
        self.db.executemany("INSERT INTO gene_info VALUES (?,?,?,?)", gene_info)
        records = self.db.execute("SELECT Count(*) FROM gene_info").fetchone()[0]
        logging.info("Table gene_info now contains %s records.", records)
        self.db.commit()

    def get_hits_to_annotated_regions(self, hits):
        """
        Return a list of annotated regions hit by peptides.
        
        :param hits:  Iterator with hits to reference sequences,
                (header, start, stop).
        :return:  Nested list of annotated regions hit by peptides.
        """

        for header, hit_start, hit_stop in hits:
            cmd = """
            WITH target_annotations AS 
            (SELECT * FROM annotations WHERE header = ?)
              SELECT attributes FROM target_annotations WHERE 
                start < ? AND end > ? OR
                start > ? AND end < ?
            """
            res = self.db.execute(cmd, (header, hit_start, hit_start, hit_stop, hit_stop))




if __name__ == "__main__":

    options = parse_args(sys.argv)

    proteodb = Proteotyping_DB_wrapper(options.proteodb)
    annotation_db = Annotation_DB_wrapper(options.db_filename, 
            options.db_comment, 
            create_db=True)
    for annotations_dir in options.annotation_dirs:
        annotation_generator = parse_annotations(proteodb, 
                annotations_dir, 
                options.glob_pattern_gff)
        annotation_db.insert_annotations(annotation_generator)

    annotation_db.insert_gene_info(parse_gene_info(options.gene_info))


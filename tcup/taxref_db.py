#!/usr/bin/env python3.5
# encoding: utf-8
#
#  ---------------------------------------------------------- 
#  This file is part of TCUP: http://tcup.readthedocs.org
#  ---------------------------------------------------------- 
#
#  Copyright (c) 2016, Fredrik Boulund <fredrik.boulund@chalmers.se>
#  
#  Permission to use, copy, modify, and/or distribute this software for any
#  purpose with or without fee is hereby granted, provided that the above
#  copyright notice and this permission notice appear in all copies.
#  
#  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
#  REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
#  FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
#  INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
#  LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
#  OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
#  PERFORMANCE OF THIS SOFTWARE.

import sys
from sys import argv, exit
from collections import OrderedDict
from functools import lru_cache
import requests
import sqlite3
import itertools
import time
import re
import logging
import argparse
import fnmatch
import os
import gzip

from ete3 import NCBITaxa

try: 
    from tcup.utils import read_fasta, find_files, grouper, existing_file
except ImportError:
    from utils import read_fasta, find_files, grouper, existing_file


def parse_commandline(argv):
    """
    Parse commandline arguments.
    """

    desc = """Prepare a TCUP "taxref" taxonomy reference database. 
    Uses header->taxid mappings to create a ready-to-use
    TCUP taxonomy reference database (taxref) to be filled in with
    sample data.  Fredrik Boulund (c) 2016."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("header_mappings", nargs="+",
            help="Path(s) to or filename(s) of tab-delimited text file with header->taxid mappings")
    parser.add_argument("--dbfile", type=str, dest="dbfile",
            default="taxref.sqlite3", 
            help="Filename to write the TCUP taxref database to [%(default)s].")
    parser.add_argument("--db-refseq-ver", dest="refseq_ver", type=str,
            default="",
            help="Specify RefSeq version, e.g. '2015-11-15'.")
    parser.add_argument("--db-comment", dest="comment", type=str,
            default="",
            help="A database creation comment added to the SQLite3 database.")

    update_db = parser.add_argument_group("update/extend existing taxref DB")
    update_db.add_argument("--add-sequences",
            action="store_true",
            default=False,
            help="""Read additional sequences from header_mappings to add to 
                existing taxref DB. 
                WARNING: Will overwrite existing taxref DB specified by --dbfile.""")

    other = parser.add_argument_group("other")
    other.add_argument("--loglevel", choices=["INFO", "DEBUG"], 
            default="INFO", 
            help="Set logging level [%(default)s].")
    other.add_argument("--logfile", 
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


class NCBITaxa_mod(NCBITaxa):
    """
    Extended/improved version of ete3.NCBITaxa.

    Improved functionality:
    - Added ability to select a custom dbfile location, instead of the default.
    - Made default dbfile location in user home dir cross-platform. 
    """

    def __init__(self, dbfile=False):
        if not dbfile:
            homedir = os.path.expanduser("~")
            self.dbfile = os.path.join(homedir, ".etetoolkit", "taxa.sqlite")
        elif dbfile and not os.path.exists(dbfile):
            logging.info("Downloading NCBI Taxonomy database")
            self.dbfile = dbfile
            self.update_taxonomy_database()
        else:
            self.dbfile = dbfile

        self.db = None
        self._connect()
 
    def expand_taxonomy_db(self, refseq_ver, taxonomy_ver, comment):
        """
        Prepare taxonomy DB for use with TCUP.

        Expands the ETE3-based taxonomy database with additional tables for
        storing what nodes are associated with reference genome sequence in the
        taxonomy. 

        :param refseq_ver:  String specifying what version of RefSeq was
                when creating the db.
        :param taxonomy_ver:  String specifying what version of NCBI
                Taxonomy was used when creating the db.
        :param comment:  String containing a comment describing the db.
        :return:  None.
        """
        
        creation_date = time.strftime("%Y-%m-%d")
        self.db.execute("DROP TABLE IF EXISTS version")
        self.db.execute("DROP TABLE IF EXISTS refseqs")
        self.db.execute("CREATE TABLE refseqs(header TEXT PRIMARY KEY, taxid INT)")
        self.db.execute("CREATE TABLE version(created TEXT, refseq TEXT, taxonomy TEXT, comment TEXT)")
        self.db.execute("INSERT INTO version VALUES (?, ?, ?, ?)", (creation_date, refseq_ver, taxonomy_ver, comment))
        self.db.commit()

    def insert_refseqs_into_db(self, refseqs):
        """
        Insert reference sequences into DB.

        :param refseqs:  Nested list/tuple with sequence header taxid pairs.
        """

        try:
            self.db.executemany("INSERT INTO refseqs VALUES (?, ?)", refseqs)
        except sqlite3.IntegrityError as e:
            logging.error("Error inserting header mappings. "
                          "Duplicate headers detected! "
                          "No mappings have been commited to DB, "
                          "remove duplicates and try again!")
            exit(2)
        self.db.commit()

    def create_refseq_indexes(self):
        """
        Create indexes for refseq_headers and refseq_taxids columns.
        """
        self.db.execute("CREATE INDEX i_refseq_headers ON refseqs (header)")
        self.db.execute("CREATE INDEX i_refseq_taxids ON refseqs (taxid)")
    
    @lru_cache(maxsize=2)
    def find_refseq_header(self, sequence_identifier):
        """
        Find a reference sequence header from the DB using a substring.
        """
        try:
            header = self.db.execute("SELECT header FROM refseqs WHERE header LIKE ?", ("%"+sequence_identifier+"%",)).fetchone()[0]
        except TypeError:
            logging.warning("Found no header for %s in database", sequence_identifier)
            raise KeyError("Found no header for %s in database" % sequence_identifier)
            # TODO: Raise custom exception?
        return header

    def extend_taxonomy_db(self, species_info_tuples):
        """
        Insert additional taxonomic nodes to taxonomy

        :param dbfile: path to dbfile.
        :param species_info_tuples: List/tuple of tuples containing 
                (taxid, parent, spname, common, rank, track), e.g.
                (12908, 1, "Viruses", "", "superkingdom", "10239,1").
        :return: None
        """

        print("NOT YET IMPLEMENTED") # TODO: Implement taxonomy db extension
        #con.executemany("INSERT INTO species VALUES (?, ?, ?, ?, ?, ?)", (taxid, parent, spname, common, rank, track))

    def dump_db(self, outputfile):
        """
        Dump entire DB in SQL text format.
        """

        if not outputfile.endswith(".gz"):
            outputfilename = outputfile+".gz"
        else:
            outputfilename = outputfile
        with gzip.GzipFile(outputfilename, "w") as out:
            logging.debug("Writing gzipped SQL DB dump to %s...", outputfilename)
            for line in self.db.iterdump():
                out.write(bytes("%s\n" % line, "utf-8"))
            logging.debug("Finished writing SQL DB dump.")


def parse_refseqs(filename):
    """
    Parse refseq:taxid mappings from tab or space delimited text file.
    """

    logging.debug("Opening %s to parse header mappings", filename)
    with open(filename) as f:
        for line in f:
            try:
                header, taxid = line.split()
            except ValueError:
                debug.error("Cannot parse line: %s", line)
                exit()
            yield header, int(taxid)


def prepare_db(dbfile, refseqs, refseq_ver, comment):
    """
    Prepare a reference DB based on the taxonomy from ETE3
    NCBITaxa, expanded with information on the nodes associated with reference
    genome sequences (e.g. from NCBI RefSeq).
    """
    
    if os.path.isfile(dbfile):
        logging.error("File '%s' already exists!", dbfile)
        exit()
    n = NCBITaxa_mod(dbfile)
    taxonomy_ver = time.strftime("%Y-%m-%d")
    n.expand_taxonomy_db(taxonomy_ver, refseq_ver, comment)
    for refseqs_file in refseqs:
        n.insert_refseqs_into_db(parse_refseqs(refseqs_file))
    n.create_refseq_indexes()


def add_sequences_to_existing_db(dbfile, header_mappings, refseq_ver, comment):
    """
    Add additional sequences to existing taxref DB file.
    """

    if not os.path.isfile(dbfile):
        logging.error("File '%s' not found!", dbfile)
        exit(1)
    n = NCBITaxa_mod(dbfile)
    for header_mappings_file in header_mappings:
        n.insert_refseqs_into_db(parse_refseqs(header_mappings_file))





def main():
    """
    Main function.
    """
    options = parse_commandline(argv)

    if options.add_sequences:
        add_sequences_to_existing_db(options.dbfile,
                options.header_mappings,
                options.refseq_ver,
                options.comment)
    else:
        prepare_db(options.dbfile, 
                options.header_mappings, 
                options.refseq_ver, 
                options.comment)

if __name__ == "__main__":
    main()

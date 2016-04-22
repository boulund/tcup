#!/usr/bin/env python3.5
# Fredrik Boulund
# (c) 2015-11-15

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


def gi_taxid_generator(gi_taxid_dmp):
    """
    Generate (gi, taxid) tuples from gi_taxid_dmp.
    """
    with open(gi_taxid_dmp) as f:
        for line in f:
            gi, taxid = line.split()
            yield int(gi), int(taxid)


def efetch_taxid(accno):
    payload = {"db": "nuccore", 
               "id": accno,
               "rettype": "fasta",
               "retmode": "xml"}
    logging.debug("Requesting taxid for accno %s via Entrez E-utils", accno)
    xml = requests.get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=payload)
    taxid = int(xml.text.split("taxid>")[1].split("<")[0])
    if taxid:
        return taxid
    else:
        return None


def create_header_taxid_mappings(refdirs, pattern, gi_taxid):
    """
    Generates a list of header:taxid pairs using gi numbers in FASTA headers.
    
    Reads gi numbers from FASTA headers (e.g. gi|52525|). The sequence
    header is split on the first space character. 
    """
    logging.getLogger("requests").setLevel(logging.WARNING)
    for refdir in refdirs:
        for fasta_file in find_files(refdir, pattern):
            for seqinfo in read_fasta(fasta_file):
                header = seqinfo[0].split()[0]
                gi = int(header.split("gi|")[1].split("|")[0])
                taxid = gi_taxid[gi]
                if taxid:
                    yield header, taxid
                else:
                    logging.debug("Found no taxid mapping for gi: %s in the taxdump db", header)
                    accno = header.split("ref|")[1].split("|")[0]
                    taxid = efetch_taxid(accno)
                    if taxid:
                        yield header, taxid
                    else:
                        logging.debug("Found no taxid for accno via NCBI E-utils: %s", accno)


class Taxdump_DB_wrapper():
    def __init__(self, sqlite3db_file, gi_taxid_dmp, rows_per_chunk=200000):
        if sqlite3db_file is None:
            self.sqlite3db_file = "gi_taxid.db"
        else:
            self.sqlite3db_file = sqlite3db_file
        if os.path.isfile(self.sqlite3db_file):
            logging.debug("Found previous DB file %s", self.sqlite3db_file)
            self.con = sqlite3.connect(self.sqlite3db_file)
        else:
            logging.debug("Found no previous DB file, creating new: %s", self.sqlite3db_file)
            if not gi_taxid_dmp: 
                raise Exception("Parameter 'gi_taxid_dmp' required to create new DB.")
            self.con = self.create_gi_taxid_db(self.sqlite3db_file, gi_taxid_dmp, rows_per_chunk)
    
    def create_gi_taxid_db(self, sqlite3db_file, gi_taxid_dmp, rows_per_chunk):
        """
        Create a (huge) sqlite3 DB with gi:taxid mappings.
        Consumes about 13 GiB storage and takes >1500 secs to make.

        Got some wild ideas on performance optimization from:
        http://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
        """

        logging.info("Creating header:taxid mappings...")

        con = sqlite3.connect(sqlite3db_file)
        con.execute("PRAGMA journal_mode = MEMORY")
        con.execute("PRAGMA synchronous = OFF")
        con.execute("CREATE TABLE gi_taxid(gi INT PRIMARY KEY, taxid INT)")
        tic = time.time()
        for pairs in grouper(rows_per_chunk, gi_taxid_generator(gi_taxid_dmp)):
            con.executemany("INSERT INTO gi_taxid VALUES (?, ?)", pairs)
        con.commit()
        toc = time.time()-tic
        num_mappings = con.execute("SELECT Count(*) FROM gi_taxid").fetchone()[0]
        logging.debug("Inserted %i mappings in %2.2f seconds.", num_mappings, toc)
        return con

    def __getitem__(self, key):
        taxid = self.con.execute("SELECT taxid FROM gi_taxid WHERE gi = ?", (key,)).fetchone()
        try:
            return taxid[0]
        except TypeError:
            return 

    def __len__(self):
        num_mappings = self.con.execute("SELECT Count(*) FROM gi_taxid").fetchone()[0]
        return num_mappings


def create_header_taxid_file(sqlite3db, refdirs, gi_taxid_dmp, outfilename, globpattern_fasta):
    """
    Create a two column text file containing FASTA header->taxids mappings.

    :param sqlite3db:  path to sqlite3 database storing gi->taxid mappings.
    :param refdirs:  list of paths to NCBI RefSeq base directory.
    :param gi_taxid_dmp:  path to file 'gi_taxid_nucl.dmp' from NCBI 
            taxonomy dump.
    :param outfile:  file to write header->taxid mappings to.
    :return:  None.
    """

    tic = time.time()
    gi_taxid_db = Taxdump_DB_wrapper(sqlite3db, gi_taxid_dmp)
    logging.debug("Database with %i entries created/loaded in %s seconds.", len(gi_taxid_db), time.time()-tic)

    tic = time.time()
    logging.debug("Parsing FASTA files under %s ...", refdirs)
    logging.debug("Writing results to %s...", outfilename)
    header_taxids = create_header_taxid_mappings(refdirs, globpattern_fasta, gi_taxid_db)
    with open(outfilename, "w") as outfile:
        for count, header_taxid in enumerate(header_taxids, start=1):
            outfile.write("{}\t{}\n".format(*header_taxid))
    logging.debug("Parsed and wrote %i header:taxid mappings in %s seconds.", count, time.time()-tic)


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

        self.db.executemany("INSERT INTO refseqs VALUES (?, ?)", refseqs)
        self.db.execute("CREATE INDEX i_refseq_headers ON refseqs (header)")
        self.db.exceute("CREATE INDEX i_refseq_taxids ON refseqs (taxid)")
        self.db.commit()
    
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

    with open(filename) as f:
        for line in f:
            header, taxid = line.split()
            yield header, int(taxid)


def prepare_db(dbfile, refseqs, taxonomy_ver, refseq_ver, comment):
    """
    Prepare a reference DB based on the taxonomy from ETE3
    NCBITaxa, expanded with information on the nodes associated with reference
    genome sequences (e.g. from NCBI RefSeq).
    """
    
    n = NCBITaxa_mod(dbfile)
    n.expand_taxonomy_db(taxonomy_ver, refseq_ver, comment)
    for refseqs_file in refseqs:
        n.insert_refseqs_into_db(parse_refseqs(refseqs_file))


def parse_commandline(argv):
    """
    Parse commandline arguments.
    """

    desc = """Prepare a TCUP reference database. 
    Step 1) Prepare a file with FASTA header->taxid mappings using sub-command
    "header_mappings".
    Step 2) Use the header->taxid mappings to create a ready-to-use
    TCUP taxonomy reference database (taxref) to be filled in with
    sample data.  Fredrik Boulund (c) 2015."""
    parser = argparse.ArgumentParser(description=desc)
    subparsers = parser.add_subparsers(dest="subcommand", help="Choose a sub-command.")

    parser_refseqs = subparsers.add_parser("header_mappings", 
            help="Prepare a list of 'sequence_header->taxid' mappings based on gi:taxid mappings from NCBI RefSeq sequences and NCBI Taxonomy taxdump.")
    parser_taxref_db = subparsers.add_parser("taxref_db",
            help="Prepare a TCUP taxonomy reference sequence (taxref) database based on NCBI Taxonomy.")
    
    parser_refseqs.add_argument("refdirs",  nargs="+",
            help="Path to NCBI RefSeq dir with sequences in FASTA format (*.fna). Walks subfolders.")
    parser_refseqs.add_argument("gi_taxid_dmp",
            help="Path to NCBI Taxonomy's 'gi_taxid_nucl.dmp'.")
    parser_refseqs.add_argument("--gi-taxid-db", dest="sqlite3db", type=existing_file,
            help="Specify a premade sqlite3 database with a gi_taxid(gi int, taxid int) table instad of creating a new one.")
    parser_refseqs.add_argument("--globpattern-fasta", dest="globpattern_fasta", type=str, metavar="'GLOB'",
            default="*.fna",
            help="Glob pattern for identifying FASTA files [%(default)s].")
    parser_refseqs.add_argument("-o", "--outfile", dest="outfile", metavar="FILE",
            default="header_taxid_mappings.tab",
            help="Output filename for header->taxid mappings [%(default)s].")
    parser_refseqs.add_argument("--loglevel", choices=["INFO", "DEBUG"], 
            default="DEBUG", 
            help="Set logging level [%(default)s].")
    parser_refseqs.add_argument("--logfile", 
            default=False,
            help="Log to file instead of STDOUT.")

    parser_taxref_db.add_argument("header_mappings", nargs="+",
            help="Path(s) to or filename(s) of two column tab-delimited text file with header->taxid mappings")
    parser_taxref_db.add_argument("--dbfile", type=str, dest="dbfile",
            default="taxref.sqlite3", 
            help="Filename to write the TCUP taxref database to [%(default)s].")
    parser_taxref_db.add_argument("--db-taxonomy-ver", dest="taxonomy_ver", type=str,
            default="",
            help="Specify Taxonomy version, e.g. '2015-11-15'.")
    parser_taxref_db.add_argument("--db-refseq-ver", dest="refseq_ver", type=str,
            default="",
            help="Specify RefSeq version, e.g. '2015-11-15'.")
    parser_taxref_db.add_argument("--db-comment", dest="comment", type=str,
            default="",
            help="A database creation comment added to the SQLite3 database.")
    parser_taxref_db.add_argument("--loglevel", choices=["INFO", "DEBUG"], 
            default="DEBUG", 
            help="Set logging level [%(default)s].")
    parser_taxref_db.add_argument("--logfile", 
            default=False,
            help="Log to file instead of STDOUT.")

    if len(argv) < 3:
        parser.print_help()
        exit()

    options = parser.parse_args()
    
    if options.logfile:
        logging.basicConfig(level=options.loglevel, filename=options.logfile)
    else:
        logging.basicConfig(level=options.loglevel)

    return options



def main():
    """
    Main function.
    """
    options = parse_commandline(argv)

    if options.subcommand == "header_mappings":
        create_header_taxid_file(options.sqlite3db, 
                options.refdir, 
                options.gi_taxid_dmp, 
                options.outfile,
                options.globpattern_fasta)
    elif options.subcommand == "taxref_db":
        prepare_db(options.dbfile, 
                options.header_mappings, 
                options.taxonomy_ver, 
                options.refseq_ver, 
                options.comment)

if __name__ == "__main__":
    main()

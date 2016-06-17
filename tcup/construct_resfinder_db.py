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

from sys import argv, exit
import difflib
import argparse
import logging
import sqlite3
import os
import re


def parse_commandline():
    """
    Parse commandline and set logging level.
    """

    desc = """Construct ResFinder sqlite3 DB. Fredrik Boulund 2016"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-s", "--sequences", dest="sequences",
            required=True,
            help="File with sequences in FASTA format.")
    parser.add_argument("-n", "--notes", dest="notes",
            required=True,
            help="File called 'notes.txt' from ResFinder.")
    parser.add_argument("-d", "--db", dest="dbfile",
            default="resfinder.sqlite3",
            help="Database filename for sqlite3 database [%(default)s].")
    parser.add_argument("--use-closest-match", dest="use_closest_match",
            action="store_true",
            default=False,
            help="Use closest match if unable to map gene symbol parsed from FASTA header to ResFinder's notes.txt [%(default)s].")
    parser.add_argument("--loglevel", 
        choices=["INFO", "DEBUG"],
        default="INFO",
        help="Set logging level [%(default)s].")

    if len(argv)<2:
        parser.print_help()
        exit()

    options = parser.parse_args()

    logging.basicConfig(level=options.loglevel)
    return options


def parse_headers_from_fasta(filename):
    """
    Parse headers from fasta file.

    Shows a warning if identical headers are seen.
    """
    already_seen = set()
    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                if header in already_seen:
                    logging.warning("Header %s seen (at least) twice.", header)
                yield header
                already_seen.add(header)


def parse_resfinder_notes(filename):
    """
    Parse database mappings from ResFinder notes.txt.
    """
    mappings = {}
    with open(filename) as f:
        for count, line in enumerate(f, start=1):
            if line.startswith("#"):
                continue
            try:
                symbol, resistance_type, extra = line.strip().split(":")
            except ValueError:
                logging.error("Line %i cannot be parsed: %s", count, line)
            mappings[symbol] = (resistance_type, extra)
    return mappings


def guess_family(h):
    """ 
    Guess the AR family of a sequence header in ResFinder. 

    The last if statement and final return statement can take care
    of most cases (such as blaXXX, QnrXXX, dfrXX), except e.g. 
    blaTEM members with capital characters on the end (e.g. blaTEM-1A).
    Adding capital characters to the rstrip at the bottom would remove 
    an 'A' too much from families such as blaFONA and several others. 
    Annoying...
    """
    if h.startswith("sul"):
        return h[0:4]  # e.g. sul1, sul2, sul3
    if h.startswith("blaTEM-"):
        return "blaTEM"
    if h.count("-") > 1:
        h = h.rsplit("-", 1)[0]
    return h.rstrip("0123456789-")


def merge_fasta_headers_and_notes(headers, notes, use_closest_match=False):
    """
    Merge FASTA headers and information from ResFinder notes.txt.

    Yields a stream of (header, symbol, family, class, extra) tuples.
    """
    missing_mappings = []

    for header in headers:
        symbol = header.split("_", 1)[0]
        family = guess_family(symbol)
        logging.debug("Guessed family %s for header %s", family, header)  # TODO: Verbose
        try:
            yield (header, symbol, family, notes[symbol][0], notes[symbol][1])
            continue
        except KeyError:
            if use_closest_match:
                closest_matching_symbol = difflib.get_close_matches(symbol, notes.keys(), n=1, cutoff=0.6)[0]
                if len(closest_matching_symbol) > 0:
                    logging.warning("Using '%s' as closest matching symbol for header '%s'", closest_matching_symbol, header)
                    infostring = " WARNING: AUTOMATICALLY ASSIGNED AS CLOSEST MATCH FOR HEADER {}".format(header)
                    yield (header, 
                           symbol, 
                           family, 
                           notes[closest_matching_symbol][0], 
                           notes[closest_matching_symbol][1]+infostring)
                    continue
            else:
                logging.warning("Symbol '%s' from header '%s' not found in %s", symbol, header, options.notes)
            missing_mappings.append(header)
    if len(missing_mappings) > 0:
        logging.warning("You should manually enter the %s missing mappings into db.", len(missing_mappings))
    return missing_mappings


def create_resfinder_sqlite3_db(dbfile, mappings):
    """
    Create and fill an sqlite3 DB with ResFinder mappings.

    Expects mappings to be a list of tuples:
      (header, symbol, family, class, extra)
    """

    logging.info("Creating sqlite3 db: %s ...", dbfile)
    if os.path.isfile(dbfile):
        logging.warning("Overwriting previously existing dbfile: %s")
        os.remove(dbfile)
        logging.debug("Removed pre-existing dbfile: %s")
    con = sqlite3.connect(dbfile)
    con.execute("CREATE TABLE resfinder(header TEXT PRIMARY KEY, symbol TEXT, family TEXT, class TEXT, extra TEXT)")
    con.executemany("INSERT INTO resfinder VALUES (?,?,?,?,?)", mappings)
    num_mappings = con.execute("SELECT Count(*) FROM resfinder").fetchone()[0]
    con.commit()
    logging.debug("Inserted %i mappings in to sqlite3 DB", num_mappings)
    return con


def main():
    """
    Main.
    """

    options = parse_commandline()
    headers = set(parse_headers_from_fasta(options.sequences))
    db_mappings = parse_resfinder_notes(options.notes)
    mapping_tuples = merge_fasta_headers_and_notes(headers, db_mappings, options.use_closest_match)
    create_resfinder_sqlite3_db(options.dbfile, mapping_tuples)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3.5
# Parse AR hits from pBLAT blast8
# Fredrik Boulund 2016

from sys import argv, exit, stdout
from collections import defaultdict, Counter
from itertools import groupby
import argparse
import logging
import sqlite3
from os import path, makedirs


def parse_commandline():
    """
    Parse commandline.
    """

    desc = """Parse AR hits from pBLAT blast8. Fredrik Boulund 2016"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FILE", type=str, nargs="+",
            help="pBLAT blast8 file to parse.")
    parser.add_argument("-r", "--resfinder", dest="resfinder",
            help="ResFinder sqlite3 database.")
    parser.add_argument("-m", "--min-identity", dest="min_identity",
            type=float,
            default=100.00,
            help="Minimum identity for pBLAT matches [%(default)s].")
    parser.add_argument("-o", "--output", dest="output", metavar="OUTFILE",
            help="Write output to OUTFILE.")
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


def parse_blat_output(filename, min_identity=100, max_pid_diff=0):
    """
    Filter out the best hit(s) for each query sequence.
    
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
    with open(filename) as blast8:
        for total_count, hit in enumerate(map(str.split, blast8), start=1):
            query = hit[0]
            target = hit[1]
            pid = float(hit[2])
            if pid < min_identity:
                continue
            hitlists[query].append((target, pid))

        try:
            logging.info("Parsed %s hits.", total_count)
        except UnboundLocalError:
            logging.error("Parsed no hits from %s; file is empty? Exiting...", filename)
            exit()
        num_remain_pep = 0
        num_remain_hits = 0
        for peptide, hitlist in hitlists.items():
            max_pid = max(h[1] for h in hitlist)
            filtered = [h for h in hitlist if h[1] >= (max_pid - max_pid_diff)]
            if len(filtered)>0:
                num_remain_pep += 1
                for hit in filtered:
                    yield peptide, hit[0], hit[1]
                    num_remain_hits += 1
        logging.info("%s hits for %s peptides remain after relative filtering.", num_remain_hits, num_remain_pep)


class ResFinderDB():
    """
    Convenience wrapper for ResFinder sqlite3 database.
    """

    def __init__(self, filename):
        self.dbfile = filename
        self.db = sqlite3.connect(self.dbfile)
        logging.debug("Connected to sqlite3 db: %s", self.dbfile)
        

    def __getitem__(self, key):
        """ Overload getter: retrieve row from db using ResFinderDB["header"] """
        header = self.db.execute("SELECT * FROM resfinder WHERE header = ?", (key,)).fetchone()
        return header


def best_matching_family_per_peptide(blast8file, min_identity, resfinder_db):
    """
    Yield the best matching family per peptide.
    """

    resfinder = ResFinderDB(resfinder_db)

    for query, hits in groupby(parse_blat_output(blast8file, min_identity, max_pid_diff=0), key=lambda x: x[0]):
        matching_targets = (target for _, target, _ in hits)
        matching_families = set(resfinder[target][2] for target in matching_targets)
        if len(matching_families) == 1:
            yield list(matching_families)[0]  # set is not indexable


def main(options):
    """
    Main.
    """

    if options.output:
        outfilehandle = open(options.output, 'w')
    else:
        outfilehandle = stdout

    with outfilehandle as outfile:
        for blast8file in options.FILE:
            family_counter = Counter(best_matching_family_per_peptide(blast8file, options.min_identity, options.resfinder))
            print("Results for {}".format(blast8file),  file=outfile)
            print("{:<10} {}".format("Count", "Family"), file=outfile)
            for family, count in family_counter.most_common():
                print("{:<10} {}".format(count, family), file=outfile)


if __name__ == "__main__":
    options = parse_commandline()
    main(options)

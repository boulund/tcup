#!/usr/bin/env python3.5
# Parse AR hits from pBLAT blast8
# Fredrik Boulund 2016

from sys import argv, exit
from collections import defaultdict
from itertools import groupby
import argparse
import logging
from os import path, makedirs


def parse_commandline():
    """
    Parse commandline.
    """

    desc = """Parse AR hits from pBLAT blast8. Fredrik Boulund 2016"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FILE", type=str, nargs="+",
        help="Filename of pBLAT blast8 file to parse.")
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


def parse_blat_output(filename, max_pid_diff=5):
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


def main(options):
    """
    Main.
    """
    for blast8file in options.FILE:
        for query, hits in groupby(parse_blat_output(blast8file, max_pid_diff=0), key=lambda x: x[0]):
            print(query)
            for q, target, pid in hits:
                print("  {}".format(target))
                pass


if __name__ == "__main__":
    options = parse_commandline()
    main(options)

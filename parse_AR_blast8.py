#!/usr/bin/env python3.5
# Parse AR hits from pBLAT blast8
# Fredrik Boulund 2016

from sys import argv, exit, stdout
from collections import defaultdict, Counter
from itertools import groupby, chain
import argparse
import logging
import sqlite3
from os import path, makedirs


def parse_commandline():
    """
    Parse commandline.
    """

    desc = """Parse AR hits from pBLAT blast8 using ResFinder. Fredrik Boulund 2016"""

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FILE", type=str, nargs="+",
            help="pBLAT blast8 file to parse.")
    parser.add_argument("-r", "--resfinder", dest="resfinder",
            required=True,
            help="ResFinder sqlite3 database.")
    parser.add_argument("-m", "--min-identity", dest="min_identity",
            type=float,
            default=100.00,
            help="Minimum identity for pBLAT matches [%(default)s].")
    parser.add_argument("-M", "--max-pid-difference", dest="max_pid_difference",
            type=float,
            default=0.0,
            help="Maximum percent identity difference for discriminative peptides [%(default)s].")
    parser.add_argument("-p", "--print-non-discriminative", dest="print_non_discriminative",
            action="store_true",
            default=False,
            help="Print all gene families matched by discriminative peptides [%(default)s].")
    parser.add_argument("-o", "--output", dest="output", metavar="OUTFILE",
            help="Write output to OUTFILE.")
    parser.add_argument("-k", "--keep-going", dest="keep_going",
            action="store_true",
            default=False,
            help="Keep going even if an blast8 input file is empty [%(default)s].")
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


def parse_blat_output(filename, min_identity=100, keep_going=False, max_pid_diff=0):
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
    global TOTAL_PEPTIDES

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
            TOTAL_PEPTIDES = total_count
        except UnboundLocalError:
            if not keep_going:
                logging.error("Parsed no hits from %s; file is empty? Exiting...", filename)
                exit()
            else:
                logging.warning("Parsed no hits from %s: file is empty? Continuing because --keep-going is set.", filename)
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
        """ 
        Overload getter: retrieve family from db using ResFinderDB["header"] 
        """
        get_family_cmd = "SELECT family FROM resfinder WHERE header = ?"
        family = self.db.execute(get_family_cmd, (key,)).fetchone()[0]
        return family
    

def best_matching_family_per_peptide(blast8file, min_identity, resfinder_db, keep_going, max_pid_difference):
    """
    Yield the best matching families per peptide.
    """

    resfinder = ResFinderDB(resfinder_db)

    for query, hits in groupby(parse_blat_output(blast8file, 
                                min_identity, 
                                keep_going, 
                                max_pid_difference), 
                                key=lambda x: x[0]):
        matching_targets = (target for _, target, _ in hits)
        matching_families = set(resfinder[target] for target in matching_targets)
        yield matching_families


def main(options):
    """
    Main.
    """
    global TOTAL_PEPTIDES

    if options.output:
        outfilehandle = open(options.output, 'w')
    else:
        outfilehandle = stdout

    with outfilehandle as outfile:
        for blast8file in options.FILE:
            families = list(best_matching_family_per_peptide(blast8file, 
                                            options.min_identity, 
                                            options.resfinder,
                                            options.keep_going,
                                            options.max_pid_difference))

            discriminative_families = [fam.pop() for fam in families if len(fam) == 1]
            discriminative_family_counter = Counter(discriminative_families)

            non_discriminative_families = list(chain(*(fam for fam in families if len(fam) > 1)))  # Flatten nested list
            non_discriminative_counter = Counter(non_discriminative_families)

            all_matched_families = set().union(discriminative_family_counter, non_discriminative_counter)
            total_discriminative_peptides = sum(discriminative_family_counter.values())
            hitcounts = {}
            for family in all_matched_families:
                try:
                    hitcounts[family] = discriminative_family_counter[family] + non_discriminative_counter[family]
                except IndexError:
                    hitcounts[family] = non_discriminative_counter[family]
            sorted_families = sorted(hitcounts, key=hitcounts.get, reverse=True)

            print("-"*70,  file=outfile)
            print("Results for {} discriminative peptides in {}".format(total_discriminative_peptides, blast8file), file=outfile)
            print("{:<6} {:<6} {:>6}  {}".format("Disc.", "Hits", "%", "Family"), file=outfile)

            if options.print_non_discriminative:
                for family in sorted_families:
                    print("{:<6} {:<6} {:>6.4f}  {}".format(discriminative_family_counter[family],
                                    hitcounts[family], 
                                    hitcounts[family]/TOTAL_PEPTIDES*100,
                                    family), file=outfile)
            else:
                for family, count in discriminative_family_counter.most_common():
                    print("{:<6} {:<6} {:>6.4f}  {}".format(count,
                                    hitcounts[family], 
                                    hitcounts[family]/TOTAL_PEPTIDES*100,
                                    family), file=outfile)
           


if __name__ == "__main__":
    # I dislike using a global value like this, but I can't
    # figure out a better way to pass this value out of the 
    # blast8-parser generator.
    global TOTAL_PEPTIDES

    options = parse_commandline()

    main(options)

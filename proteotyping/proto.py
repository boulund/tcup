#!/usr/bin/env python3.5
# encoding: utf-8
# Fredrik Boulund (c) 2015
# Proteotyping

from sys import argv, exit
from collections import defaultdict, OrderedDict, Counter
from itertools import groupby
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
    from annotation_db import Annotation_DB_wrapper
except ImportError:
    from proteotyping.utils import find_files, grouper, existing_file
    from proteotyping.annotation_db import Annotation_DB_wrapper


def parse_commandline(argv):
    """Parse commandline arguments"""

    desc = """Proteotyping pipeline. (c) Fredrik Boulund 2015."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("FILE", nargs="+",
            type=existing_file,
            help="BLAT output file.")
    parser.add_argument("--proteodb", dest="proteodb", metavar="DB", type=str,
            default="proteodb.sql",
            help="Path to proteotyping sqlite3 DB [%(default)s].")
    parser.add_argument("--annotation-db", dest="annotation_db", metavar="FILE",
            default="annotationdb.sql",
            help="Path to annoation sqlite3 DB [%(default)s].")
    parser.add_argument("--sample-db", dest="sample_db", action="store_true",
            default=False,
            help="The supplied 'BLAT output file' is really a processed sample db with proteotyping results [%(default)s]")
    parser.add_argument("--taxonomic-rank", dest="taxonomic_rank", metavar="LVL", type=str,
            choices=["no rank", "subspecies", "species", "genus", "family", "order", "class", "phylum", "superkingdom"],
            default="family",
            help="Set the taxonomic level on which hits are grouped [%(default)s].")
    parser.add_argument("--print-all-discriminative-peptides", 
            dest="print_all_discriminative_peptides", 
            action="store_true",
            default=False,
            help="Print all discriminative peptides [%(default)s].")
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
                    yield peptide, hit[0], int(hit[3]), int(hit[4]), float(hit[1]), int(hit[2])
                    num_remain_hits += 1
        logging.info("%s hits for %s peptides remain after relative filtering.", num_remain_hits, num_remain_pep)


class Proteotyping_DB_wrapper():
    """ 
    Wrapper for sqlite3 database.

    Stores a hard coded rank hierarchy based on the ranks seen 
    in NCBI Taxonomy. Note the problematic "no rank" assignment
    that can be encountered at different levels in the taxonomy, 
    here incorrectly positioned to always be at the lowest rank.
    """

    def __init__(self, sample_db, ref_dbfile=None):
        # Setup SQLite DB references
        self.dbfile = sample_db
        if ref_dbfile and os.path.isfile(self.dbfile):
            sleep_time = 5
            logging.warning("About to delete pre-existing sample DB: %s", self.dbfile)
            logging.warning("Press Ctrl+c to cancel within %i sec...", sleep_time)
            time.sleep(sleep_time)
            os.remove(self.dbfile)
            logging.debug("Copying reference DB (%s) to create new sample DB (%s)", ref_dbfile, self.dbfile)
            shutil.copy(ref_dbfile, self.dbfile)
        elif ref_dbfile and not os.path.isfile(ref_dbfile):
            logging.error("Cannot connect to non-existent DB file: %s", ref_dbfile)
            exit(1)
        elif ref_dbfile:
            logging.debug("Copying ref_dbfile %s to sample db location: %s", ref_dbfile, sample_db)
            shutil.copy(ref_dbfile, self.dbfile)

        self.db = sqlite3.connect(self.dbfile)
        logging.info("Connected to sample DB %s", self.dbfile)
        try:
            logging.info("DB version: %s", self.db.execute("SELECT * FROM version").fetchall()[0])
        except sqlite3.OperationalError as msg:
            logging.error("Incorrect or missing DB %s: %s", self.dbfile, msg)
            exit(1)

        # Define NCBI Taxonomy rank hierachy
        self.rank_hierarchy = ["superkingdom",
                "kingdom",
                "subkingdom",
                "superphylum",
                "phylum",
                "subphylum",
                "class",
                "superclass",
                "subclass",
                "infraclass",
                "superorder",
                "order",
                "suborder",
                "infraorder",
                "parvorder",
                "superfamily",
                "family",
                "subfamily",
                "tribe",
                "subtribe",
                "genus",
                "subgenus",
                "species group",
                "species subgroup",
                "species",
                "subspecies",
                "varietas",
                "forma",
                "no rank"]
        self.ranks = {r: n for n, r in enumerate(self.rank_hierarchy)}

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
            tracks = [list(map(int, t[0].split(","))) for t in query.fetchall()]
            lca = self.lowest_common_ancestor(tracks)
            try:
                if lca[0] != 131567:  # taxid=131567 is "cellular organisms"
                    # TODO: verbose
                    #logging.debug("Peptide %s is discriminative at rank '%s' for %s", peptide[0], rank, spname)
                    self.db.execute("INSERT INTO discriminative VALUES (?, ?)", (peptide[0], lca[0]))
            except IndexError:
                logging.warning("Found no LCA for %s with tracks %s", peptide, [list(track) for track in tracks])

            # Find the track lineage of the LCA to increment discriminative_count on all
            # lineage members.
            lca_lineage_query = self.db.execute("SELECT track FROM species WHERE taxid = (?)", lca).fetchone()
            lca_lineage = list(map(int, lca_lineage_query[0].split(",")))
            update_cmd = "UPDATE species SET discriminative_count = discriminative_count + 1 WHERE taxid IN ({})"
            self.db.execute(update_cmd.format(",".join("?"*len(lca_lineage))), lca_lineage)

        self.db.commit()



    def get_discriminative_at_rank(self, rank):
        """
        Retrieve all discriminative peptides at the specified rank.
        """

        rank_set = self.rank_hierarchy[self.ranks[rank]:]

        cmd = """SELECT peptide, rank, spname FROM discriminative
          JOIN species ON species.taxid = discriminative.taxid
          WHERE rank IN ({})
          ORDER BY rank""".format(",".join("?"*len(rank_set)))
        result = self.db.execute(cmd, rank_set).fetchall()
        return result

    def get_discriminative_counts(self):
        """
        Retrieve cumulative peptide counts across all ranks.
        """

        cmd = """SELECT discriminative_count, rank, spname FROM species 
            WHERE discriminative_count > 0"""
        result = self.db.execute(cmd).fetchall()
        return result
    
    def get_peptide_hits(self, rank):
        """
        Retrieve all hits at or below a given rank, sorted by target sequence.

        :param rank:  Rank at and below which to return hits.
        :return:  List of hits to reference sequences, 
                sorted by reference sequence headers.
        """
        
        pass


def print_cumulative_discriminative_counts(disc_peps_per_rank, ranks):
    """
    Print sorted lists of discriminative peptide counts across all ranks.
    """

    sorted_disc_peps_per_rank = sorted(disc_peps_per_rank, key=lambda entry: ranks[entry[1]])
    # Removed this printout because it is too messy and non-informative.
    #for rank, group in groupby(sorted_disc_peps_per_rank, key=lambda entry: ranks[entry[1]]):
    #    for count, rank, spname in sorted(group, reverse=True):
    #        if spname in ("root", "cellular organisms"):
    #            continue
    #        print("{:<6} {:<20} {:<40}".format(count, rank, spname))
    for count, rank, spname in sorted(disc_peps_per_rank, reverse=True):
        if spname in ("root", "cellular organisms"):
            continue
        print("{:<6} {:<20} {:<40}".format(count, rank, spname))



def print_discriminative_peptides(discriminative):
    """
    Print name and assignment of discriminative peptides.
    """
    print("Discriminative peptides".center(60, "-"))
    print("{:<20} {:<15} {:<30}".format("Peptide", "Rank", "Description"))
    for d in discriminative: 
        print("{:<20} {:<15} {:<30}".format(*d))


def print_peptides_per_spname(discriminative):
    """
    Print number of discriminative peptides per spname.
    """
    species_counts = Counter((species, rank) for pep, rank, species in discriminative)
    print("Discriminative peptides per spname".center(60, "-"))
    print("{:<6} {:<20} {:<40}".format("Count", "Rank", "Description"))
    for species, count in species_counts.most_common():
        print("{:<6} {:<20} {:<40}".format(count, species[1], species[0]))


def get_results_from_existing_db(sample_databases,
        annotation_db, 
        taxonomic_rank="family",
        print_all_discriminative_peptides=False,
        print_cumulative_counts=True):
    """
    Retrieve results from an existing sample db.
    """

    annotation_db = Annotation_DB_wrapper(annotation_db)

    for sample_db in sample_databases:
        if type(sample_db) == Proteotyping_DB_wrapper:
            refdb = sample_db
        else: 
            refdb = Proteotyping_DB_wrapper(sample_db)
        disc = refdb.get_discriminative_at_rank(taxonomic_rank)

        print(refdb.dbfile.center(60, "-"))
        if print_all_discriminative_peptides:
            print_discriminative_peptides(disc)
        print_peptides_per_spname(disc)

        print(refdb.dbfile.center(60, "-"))
        disc_peps_per_rank = refdb.get_discriminative_counts()
        if print_cumulative_counts:
            print_cumulative_discriminative_counts(disc_peps_per_rank, refdb.ranks)

        hits = [("gi|152977688|ref|NC_009655.1|", 2500, 2700)]
        annotation_db.get_hits_to_annotated_regions(hits)

        


def main(options):
    """
    Main function that runs the complete pipeline logic.
    """
    blacklisted_seqs = prepare_blacklist(options.blacklist, options.leave_out)

    for blat_file in options.FILE:
        refdb = Proteotyping_DB_wrapper(blat_file+".sqlite3", options.proteodb) 
        blat_parser = parse_blat_output(blat_file, 
                options.min_identity, 
                options.min_matches, 
                options.min_coverage,
                options.max_pid_diff, 
                blacklisted_seqs)
        refdb.insert_blat_hits_into_db(blat_parser)
        refdb.determine_discriminative_ranks()

        get_results_from_existing_db([refdb],
                options.annotation_db,
                options.taxonomic_rank,
                options.print_all_discriminative_peptides)


if __name__ == "__main__":

    options = parse_commandline(argv)

    if options.sample_db:
        get_results_from_existing_db(options.FILE,
                options.annotation_db,
                options.taxonomic_rank,
                options.print_all_discriminative_peptides)
    else:
        main(options)

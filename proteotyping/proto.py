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
import xlsxwriter

try:
    from utils import find_files, grouper, existing_file
except ImportError:
    from proteotyping.utils import find_files, grouper, existing_file


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
    parser.add_argument("--annotation-db", dest="annotation_db_file", metavar="FILE",
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
    parser.add_argument("--print-annotations", dest="print_annotations", action="store_true",
            default=False,
            help="Print all annotated regions hit by discriminative fragments [%(default)s].")
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
    parser.add_argument("--write-xlsx", action="store_true",
            default=False,
            help="Write results to Excel xslx file [%(default)s].")
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


class Sample_DB_wrapper():
    """ 
    Wrapper for an Sqlite3 database for storing proteotyping results.

    This class also stores a hard coded rank hierarchy based on the ranks seen
    in NCBI Taxonomy. Note the problematic "no rank" assignment that can be
    encountered at different levels in the taxonomy, here incorrectly
    positioned to always be at the lowest rank.

    The class has methods that allow attaching a proteotyping reference
    database and an annotation database.
    """

    def __init__(self, sample_db_filename, create_new=True):
        if os.path.isfile(sample_db_filename) and not create_new:
            self._connect_to_existing_sample_db(sample_db_filename)
        elif os.path.isfile(sample_db_filename) and create_new:
            sleep_time = 3
            logging.warning("About to delete pre-existing sample DB: %s", sample_db_filename)
            logging.warning("Press Ctrl+c to cancel within %i sec...", sleep_time)
            time.sleep(sleep_time)
            os.remove(sample_db_filename)
            self._create_new_sample_db(sample_db_filename)
        else:
            self._create_new_sample_db(sample_db_filename)

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

    def _connect_to_existing_sample_db(self, sample_db_filename):
        self.dbfile = sample_db_filename
        self.db = sqlite3.connect(self.dbfile)
        logging.debug("Connected to pre-existing sample DB %s", self.dbfile)
    
    def _create_new_sample_db(self, sample_db_filename):
        self.dbfile = sample_db_filename
        self.db = sqlite3.connect(self.dbfile)
        logging.info("Creating new sample DB %s", self.dbfile)
        self.db.execute("CREATE TABLE mappings(peptide TEXT, target TEXT, start INT, end INT, identity REAL, matches INT)")
        self.db.execute("CREATE TABLE peptides(peptide TEXT, discriminative_taxid INT)")
        self.db.execute("CREATE TABLE cumulative(taxid INT, count INT DEFAULT 0)")
        self.db.execute("CREATE TABLE rank_counts(count INT DEFAULT 0, rank TEXT)")


    def attach_proteotyping_ref_db(self, proteotyping_ref_db_file):
        self.db.execute("ATTACH ? as proteodb", (proteotyping_ref_db_file, ))
        logging.debug("Attached proteotyping taxref DB")
        self.db.execute("INSERT INTO cumulative (taxid) SELECT taxid from proteodb.species")
        self.db.commit()
    
    def attach_annotation_db(self, annotation_db_file):
        self.db.execute("ATTACH ? as annotationdb", (annotation_db_file, ))
        logging.debug("Attached annotation DB")

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
            self.db.executemany("INSERT INTO mappings VALUES (?, ?, ?, ?, ?, ?)", values)
        self.db.execute("INSERT INTO peptides (peptide) SELECT DISTINCT peptide FROM mappings")
        self.db.execute("CREATE INDEX i_peptides_disc ON peptides(discriminative_taxid)")
        self.db.execute("CREATE INDEX i_mappings_targets ON mappings(target)")
        self.db.commit()

    
    def determine_discriminative_ranks(self):
        """
        Find the discriminative rank of each unique peptide in the db.
        """

        for peptide in self.db.execute("SELECT DISTINCT peptide FROM peptides"):
            cmd = """ SELECT track FROM proteodb.species
              JOIN proteodb.refseqs ON refseqs.taxid = species.taxid
              JOIN mappings ON mappings.target = proteodb.refseqs.header AND mappings.peptide = ?
            """
            query = self.db.execute(cmd, peptide)
            tracks = [list(map(int, t[0].split(","))) for t in query.fetchall()]
            if not tracks:
                # Find track based on another taxid, perhaps has the target
                # taxid has been merged into another.
                logging.debug("Found no track for %s in proteodb.species, trying to see if taxid has been merged", peptide)
                taxid_old_cmd = """SELECT track FROM proteodb.species
                  JOIN proteodb.refseqs ON refseqs.taxid = proteodb.merged.taxid_old
                  JOIN proteodb.merged ON merged.taxid_new = species.taxid
                  JOIN mappings ON mappings.target = proteodb.refseqs.header
                  AND mappings.peptide = ?
                """
                query = self.db.execute(taxid_old_cmd, peptide)
                tracks = [list(map(int, t[0].split(","))) for t in query.fetchall()]
                if tracks:
                    logging.debug("Found track for %s, taxid was merged", peptide)
            lca = self.lowest_common_ancestor(tracks)
            try:
                if lca[0] != 131567:  # taxid=131567 is "cellular organisms"
                    # TODO: verbose
                    #logging.debug("Peptide %s is discriminative at rank '%s' for %s", peptide[0], rank, spname)
                    self.db.execute("UPDATE peptides SET discriminative_taxid = ? WHERE peptide = ?", (lca[0], peptide[0]))
            except IndexError:
                logging.warning("Found no LCA for %s", peptide)

            # Find the track lineage of the LCA to increment the count of
            # number of discriminative fragments beneath that node on all
            # lineage members.
            lca_lineage_query = self.db.execute("SELECT track FROM proteodb.species WHERE taxid = (?)", lca).fetchone()
            lca_lineage = list(map(int, lca_lineage_query[0].split(",")))
            update_cmd = "UPDATE cumulative SET count = count + 1 WHERE taxid IN ({})"
            self.db.execute(update_cmd.format(",".join("?"*len(lca_lineage))), lca_lineage)

        self.db.commit()

    def count_discriminative_per_rank(self):
        """
        Count the number of discriminative peptides per rank and store in
        sample DB.
        """

        logging.debug("Counting discriminative peptides per rank...")
        cmd = """INSERT INTO rank_counts
          SELECT count(rank), rank 
          FROM peptides
          JOIN proteodb.species ON peptides.discriminative_taxid = proteodb.species.taxid
          GROUP BY rank
          ORDER BY count(rank) DESC
          """
        self.db.execute(cmd)
        self.db.commit()

    def get_rank_counts(self):
        """
        Returns a dictionary with rank:counts mappings.
        """
        return {r: c for c, r in self.db.execute("SELECT * FROM rank_counts").fetchall()}

    def get_discriminative_at_rank(self, rank):
        """
        Retrieve all discriminative peptides at and below the specified rank.
        """

        rank_set = self.rank_hierarchy[self.ranks[rank]:]

        cmd = """SELECT peptide, rank, spname FROM peptides
          JOIN proteodb.species ON proteodb.species.taxid = peptides.discriminative_taxid
          WHERE rank IN ({})
          ORDER BY rank""".format(",".join("?"*len(rank_set)))
        logging.debug("Retrieving discriminative peptides at and below rank %s...", rank)
        result = self.db.execute(cmd, rank_set).fetchall()
        return result

    def get_discriminative_counts_at_rank(self, rank):
        """
        Retrieve cumulative peptide counts across all ranks.
        """

        rank_set = self.rank_hierarchy[self.ranks[rank]:]

        cmd = """SELECT cumulative.count, count(peptide), rank, spname FROM peptides
        JOIN proteodb.species, cumulative 
        ON proteodb.species.taxid = peptides.discriminative_taxid
        AND cumulative.taxid = proteodb.species.taxid
        WHERE proteodb.species.rank in ({})
        AND cumulative.count > 0
        GROUP BY rank, spname
        ORDER BY cumulative.count DESC
        """.format(",".join("?"*len(rank_set)))

        logging.debug("Retrieving cumulative counts at and below rank %s...", rank)
        result = self.db.execute(cmd, rank_set).fetchall()
        return result
    
    def get_peptide_hits(self, rank):
        """
        Retrieve all hits at or below a given rank, sorted by target sequence.

        :param rank:  Rank at and below which to return hits.
        :return:  List of hits to reference sequences, 
                sorted by reference sequence headers.
        """
        pass
    
    def get_hits_to_annotated_regions(self):
        """
        Retrieve all annotated regions matched by any peptide.
        """

        cmd = """SELECT target, spname, product, features FROM annotationdb.annotations 
          JOIN mappings ON mappings.target = annotationdb.annotations.header
            AND (mappings.start BETWEEN annotationdb.annotations.start 
              AND annotationdb.annotations.end
            OR mappings.end BETWEEN annotationdb.annotations.start 
              AND annotationdb.annotations.end
            OR (mappings.start <= annotationdb.annotations.start 
              AND mappings.end >= annotationdb.annotations.end
              )
            )
          JOIN proteodb.refseqs ON mappings.target = proteodb.refseqs.header
          JOIN proteodb.species ON proteodb.refseqs.taxid = proteodb.species.taxid
        """
        logging.debug("Retrieving all annotated regions matched by any discriminative peptide...")
        result = self.db.execute(cmd).fetchall()
        return result


def print_cumulative_discriminative_counts(disc_peps_per_rank, rank_counts):
    """
    Print sorted lists of discriminative peptide counts across all ranks.
    """

    print("Discriminative peptides per spname".center(60, "-"))
    print("{:<10} {:<6} {:>6} {:<20} {:<40}".format("Cumulative", "Count", "%", "Rank", "Description"))
    for cum_count, count, rank, spname in disc_peps_per_rank:
        if spname in ("root", "cellular organisms"):
            continue
        percentage = count/rank_counts[rank] * 100
        print("{:<10} {:<6} {:>6.2f} {:<20} {:<40}".format(cum_count, count, percentage, rank, spname))


def print_discriminative_peptides(discriminative):
    """
    Print name and assignment of discriminative peptides.
    """
    print("Discriminative peptides".center(60, "-"))
    print("{:<20} {:<15} {:<30}".format("Peptide", "Rank", "Description"))
    for d in discriminative: 
        print("{:<20} {:<15} {:<30}".format(*d))


def print_annotation_hits(hits):
    """
    Print hits to annotated genome regions.
    """
    print("Hits to annotated genome regions".center(60, "-"))
    print("{:<34}\t{:<40}\t{:<30}\t{:<}".format("Genome sequence", "Spname", "Product", "Features"))
    for seq, spname, product, features in hits:
        print("{:<34}\t{:<40}\t{:<30}\t{:<}".format(seq, spname, product, features))

def write_results_xlsx(disc_peps_per_rank, rank_counts, hits, results_filename):
    """
    Write results to an Excel xslx file.
    """

    workbook = xlsxwriter.Workbook(results_filename)
    worksheet_composition= workbook.add_worksheet("Taxonomic composition")

    percentage_format = workbook.add_format({"num_format": "0.00%"})

    worksheet_composition.write(0, 0, "Cumulative")
    worksheet_composition.write(0, 1, "Discriminative count")
    worksheet_composition.write(0, 2, "Percentage")
    worksheet_composition.write(0, 3, "Rank")
    worksheet_composition.write(0, 4, "Description")
    worksheet_composition.set_column(0, 0, 9.0)
    worksheet_composition.set_column(1, 1, 17.0)
    worksheet_composition.set_column(2, 2, 9.0)
    worksheet_composition.set_column(3, 3, 11.0)
    worksheet_composition.set_column(4, 4, 40.0)
    for row, data in enumerate(disc_peps_per_rank, start=1):
        cum_count, count, rank, spname = data
        if spname in ("root", "cellular organisms"):
            continue
        percentage = count/rank_counts[rank]
        worksheet_composition.write(row, 0, cum_count)
        worksheet_composition.write(row, 1, count)
        worksheet_composition.write(row, 2, percentage, percentage_format)
        worksheet_composition.write(row, 3, rank)
        worksheet_composition.write(row, 4, spname)

    worksheet_annotations = workbook.add_worksheet("Hits to annotated regions")
    worksheet_annotations.write(0, 0, "Genome sequence")
    worksheet_annotations.write(0, 1, "Species")
    worksheet_annotations.write(0, 2, "Product")
    worksheet_annotations.write(0, 3, "Features")
    worksheet_annotations.set_column(0, 0, 34.0)
    worksheet_annotations.set_column(1, 1, 40.0)
    worksheet_annotations.set_column(2, 2, 45.0)
    for row, data in enumerate(hits, start=1):
        seq, spname, product, features = data
        worksheet_annotations.write(row, 0, seq)
        worksheet_annotations.write(row, 1, spname)
        worksheet_annotations.write(row, 2, product)
        worksheet_annotations.write(row, 3, features)
    
    workbook.close()

def get_results_from_existing_db(sample_databases,
        annotation_db_file, 
        taxonomic_rank="family",
        print_all_discriminative_peptides=False,
        print_annotations=False,
        write_xlsx=True):
    """
    Retrieve results from existing sample database(s).
    """

    for sample_db in sample_databases:
        sample_db.attach_annotation_db(annotation_db_file)

        disc = sample_db.get_discriminative_at_rank(taxonomic_rank)
        disc_peps_per_rank = sample_db.get_discriminative_counts_at_rank(taxonomic_rank)
        rank_counts = sample_db.get_rank_counts()

        print(sample_db.dbfile.center(60, "-"))
        if print_all_discriminative_peptides:
            print_discriminative_peptides(disc)

        print_cumulative_discriminative_counts(disc_peps_per_rank, rank_counts)

        if print_annotations:
            hits = sample_db.get_hits_to_annotated_regions()
            print_annotation_hits(hits)

        if write_xlsx:
            hits = sample_db.get_hits_to_annotated_regions()
            xlsx_filename = sample_db.dbfile.split(".", 1)[0] + ".xlsx"
            write_results_xlsx(disc_peps_per_rank, rank_counts, hits, xlsx_filename)



        


def main(options):
    """
    Main function that runs the complete pipeline logic.
    """
    blacklisted_seqs = prepare_blacklist(options.blacklist, options.leave_out)

    for blat_file in options.FILE:
        if blat_file.endswith(".sqlite3"):
            sleep_time = 5
            logging.warning("{} ends with .sqlite3, did you mean to use --sample-db?".format(blat_file))
            logging.warning("Press Ctrl+c to cancel within %i sec...", sleep_time)
            time.sleep(sleep_time)
        sample_db_filename = blat_file + ".sqlite3"
        sample_db = Sample_DB_wrapper(sample_db_filename)
        blat_parser = parse_blat_output(blat_file, 
                options.min_identity, 
                options.min_matches, 
                options.min_coverage,
                options.max_pid_diff, 
                blacklisted_seqs)
        sample_db.insert_blat_hits_into_db(blat_parser)

        sample_db.attach_proteotyping_ref_db(options.proteodb)
        sample_db.determine_discriminative_ranks()
        sample_db.count_discriminative_per_rank()

        get_results_from_existing_db([sample_db],
                options.annotation_db_file,
                options.taxonomic_rank,
                options.print_all_discriminative_peptides,
                options.print_annotations,
                options.write_xlsx)


if __name__ == "__main__":

    options = parse_commandline(argv)

    if options.sample_db:
        sample_databases = [Sample_DB_wrapper(filename, create_new=False) for filename in options.FILE]
        for sample_db in sample_databases:
            sample_db.attach_proteotyping_ref_db(options.proteodb)
        get_results_from_existing_db(sample_databases,
                options.annotation_db_file,
                options.taxonomic_rank,
                options.print_all_discriminative_peptides,
                options.print_annotations,
                options.write_xlsx)
    else:
        main(options)

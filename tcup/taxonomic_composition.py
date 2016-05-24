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

from sys import argv, exit, stdout
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
    from tcup.utils import find_files, grouper, existing_file


def parse_commandline(argv):
    """Parse commandline arguments"""

    desc = """TCUP: Typing and Characterization using Proteomics. (c) Fredrik Boulund 2016."""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("FILE", nargs="+",
            type=existing_file,
            help="BLAT output file.")
    parser.add_argument("--taxref-db", dest="taxref_db", metavar="TAXREFDB", 
            type=existing_file,
            default="taxref_db.sqlite3",
            help="Path to taxonomy reference sqlite3 DB [%(default)s].")
    parser.add_argument("--annotation-db", dest="annotation_db_file", metavar="ANNOTATIONDB", 
            type=existing_file,
            default="annotation_db.sqlite3",
            help="Path to annotation sqlite3 DB [%(default)s].")
    parser.add_argument("--sample-db", dest="sample_db", metavar="SAMPLEDB",
            default=False,
            help="Filename of sqlite3 db created for the sample [<sample.blast8>.sqlite3].")
    parser.add_argument("--pre-existing", dest="pre_existing", action="store_true",
            default=False,
            help="The supplied 'BLAT output file' is really an sqlite3 sample db.")
    parser.add_argument("--taxonomic-rank", dest="taxonomic_rank", metavar="LVL", type=str,
            choices=["no rank", "subspecies", "species", "genus", "family", "order", "class", "phylum", "superkingdom"],
            default="family",
            help="Set the taxonomic level on which hits are grouped [%(default)s].")
    parser.add_argument("--blacklist", metavar="FILE", dest="blacklist",
            default=None,
            type=existing_file,
            help="File with sequence headers to blacklist (i.e. to ignore when parsing blast8 output).")
    parser.add_argument("--species-normalization", metavar="FILE", dest="species_normalization",
            default="",
            help="Tab separated file with species normalization factors. Using your own will overwrite built-in values.")
    parser.add_argument("--write-discriminative-peptides", metavar="FILE",
            dest="write_discriminative_peptides", 
            default=False,
            help="Write all discriminative peptides to FILE.")
    parser.add_argument("--print-annotations", dest="print_annotations", action="store_true",
            default=False,
            help="Print all annotated regions hit by discriminative fragments [%(default)s].")
    parser.add_argument("--min-matches", dest="min_matches", metavar="L", type=int,
            default=6,
            help="Minimum peptide matches (i.e. peptide length) [%(default)s].")
    parser.add_argument("-i", "--min-identity", dest="min_identity", metavar="I", type=float,
            default=90.0,
            help="Filter out hits with less than or equal to this percentage identity [%(default)s].")
    parser.add_argument("--min-coverage", dest="min_coverage", metavar="C", type=float,
            default=1.0,
            help="Proportion of fragment covered in alignment [%(default)s].")
    parser.add_argument("--max-pid-diff", dest="max_pid_diff", type=float, metavar="D",
            default=5.0,
            help="Maximum identity difference between highest and lowest hit for each peptide. Floating point between 0.0-100.0 [%(default)s].")
    parser.add_argument("--inclusion-threshold", dest="inclusion_threshold", metavar="c", type=int,
            default=0,
            help="Minimum number of discriminative peptides required to report taxa in output [%(default)s].")
    parser.add_argument("--write-xlsx", dest="write_xlsx", metavar="XLSX_FILE",
            default="",
            help="Write results to Excel file.")
    parser.add_argument("--output", dest="output", metavar="FILE",
            default=False,
            help="Write results to FILE instead of STDOUT.")


    devoptions = parser.add_argument_group("Developer options", "Voids warranty ;)")
    devoptions.add_argument("--loglevel", choices=["INFO", "DEBUG", "VERBOSE"],
            default="DEBUG",
            help="Set logging level [%(default)s].")
    devoptions.add_argument("--logfile", dest="logfile",
            default="tcup.log",
            help="Filename for log output [%(default)s].")
    devoptions.add_argument("--leave-out", metavar="HEADERS", dest="leave_out",
            default="",
            help="Disregard any hits to sequence with HEADER when parsing and filtering blast8 output. Can be a list of comma separated FASTA headers (no spaces).")

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
        with open(blacklist) as f:
            for line in f:
                blacklisted_seqs.add(line.strip())
    if additional:
        [blacklisted_seqs.add(header) for header in additional.split(",")]
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
        total_count = 0
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
        if not total_count:
            logging.error("Parsed no hits from %s, file empty? Exiting...", filename)
            exit(2)
        logging.info("Parsed %s hits.", total_count)
        logging.info("Discarded %s hits based on primary criteria.", num_discarded)
        num_remain_pep = 0
        num_remain_hits = 0
        for peptide, hitlist in hitlists.items():
            max_pid = max(map(float, (h[1] for h in hitlist)))
            filtered = [h for h in hitlist if float(h[1]) >= (max_pid - max_pid_diff)]
            if len(filtered)>0:
                num_remain_pep += 1
                for hit in filtered:
                    yield peptide, hit[0], int(hit[3]), int(hit[4]), float(hit[1]), int(hit[2])
                    num_remain_hits += 1
        logging.info("%s hits for %s peptides remain after relative filtering.", num_remain_hits, num_remain_pep)


class Sample_DB_wrapper():
    """ 
    Wrapper for an Sqlite3 database for storing TCUP results.

    This class also stores a hard coded rank hierarchy based on the ranks seen
    in NCBI Taxonomy. Note the problematic "no rank" assignment that can be
    encountered at different levels in the taxonomy, here incorrectly
    positioned to always be at the lowest rank.

    The class has methods that allow attaching a TCUP reference
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
        # Commented lines are not included in results printout.
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
                #"subfamily",
                #"tribe",
                #"subtribe",
                "genus",
                #"subgenus",
                #"species group",
                #"species subgroup",
                "species",
                #"subspecies",
                #"varietas",
                #"forma",
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
        self.db.execute("CREATE TABLE peptides(peptide TEXT PRIMARY KEY, discriminative_taxid INT)")
        self.db.execute("CREATE TABLE cumulative(taxid INT PRIMARY KEY, count INT DEFAULT 0)")
        self.db.execute("CREATE TABLE rank_counts(rank TEXT PRIMARY KEY, count INT DEFAULT 0)")


    def attach_taxref_db(self, taxref_db_file):
        self.db.execute("ATTACH ? as taxref", (taxref_db_file, ))
        logging.debug("Attached taxref DB")
        self.db.execute("INSERT OR IGNORE INTO cumulative (taxid) SELECT taxid from taxref.species")
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
            cmd = """ SELECT track FROM taxref.species
              JOIN taxref.refseqs ON refseqs.taxid = species.taxid
              JOIN mappings ON mappings.target = taxref.refseqs.header AND mappings.peptide = ?
            """
            query = self.db.execute(cmd, peptide)
            tracks = [list(map(int, t[0].split(","))) for t in query.fetchall()]
            if not tracks:
                # Find track based on another taxid, perhaps has the target
                # taxid has been merged into another.
                logging.debug("Found no track for %s in taxref.species, trying to see if taxid has been merged", peptide)
                taxid_old_cmd = """SELECT track FROM taxref.species
                  JOIN taxref.refseqs ON refseqs.taxid = taxref.merged.taxid_old
                  JOIN taxref.merged ON merged.taxid_new = species.taxid
                  JOIN mappings ON mappings.target = taxref.refseqs.header
                  AND mappings.peptide = ?
                """
                query = self.db.execute(taxid_old_cmd, peptide)
                tracks = [list(map(int, t[0].split(","))) for t in query.fetchall()]
                if tracks:
                    logging.debug("Found track for %s, taxid was merged", peptide)
            lca = self.lowest_common_ancestor(tracks)
            if len(lca) == 0:
                logging.warning("Found no LCA for %s", peptide)
                continue
            if lca[0] != 131567:  # taxid=131567 is "cellular organisms"
                # TODO: verbose
                #logging.debug("Peptide %s is discriminative at rank '%s' for %s", peptide[0], rank, spname)
                self.db.execute("UPDATE peptides SET discriminative_taxid = ? WHERE peptide = ?", (lca[0], peptide[0]))

            # Find the track lineage of the LCA to increment the count of
            # number of discriminative fragments beneath that node on all
            # lineage members.
            lca_lineage_query = self.db.execute("SELECT track FROM taxref.species WHERE taxid = (?)", lca).fetchone()
            lca_lineage = list(map(int, lca_lineage_query[0].split(",")))

            # TODO: Change this into an UPSERT-like thing to remove the need for
            # including ALL taxids when attaching the taxref the first time,
            # and remove the need to delete all count = 0 below.
            update_cmd = "UPDATE cumulative SET count = count + 1 WHERE taxid IN ({})"
            self.db.execute(update_cmd.format(",".join("?"*len(lca_lineage))), lca_lineage)

        self.db.execute("DELETE FROM cumulative WHERE count = 0")
        self.db.commit()

    def count_discriminative_per_rank(self):
        """
        Count the number of discriminative peptides per rank and store in
        sample DB.
        """

        logging.debug("Counting discriminative peptides per rank...")
        cmd = """INSERT INTO rank_counts
          SELECT rank, count(rank)
          FROM peptides
          JOIN taxref.species ON peptides.discriminative_taxid = taxref.species.taxid
          GROUP BY rank
          ORDER BY count(rank) DESC
          """
        self.db.execute(cmd)
        self.db.commit()

    def get_rank_counts(self):
        """
        Returns a dictionary with rank:count mappings.
        """
        return {r: c for c, r in self.db.execute("SELECT * FROM rank_counts").fetchall()}
    
    def get_cumulative_rank_counts(self):
        """
        Returns a dictionary with cumulative rank:count mappings.
        """
        logging.debug("Getting cumulative rank:count mappings.")
        cmd = """SELECT rank, sum(count) FROM cumulative
          JOIN taxref.species
          ON taxref.species.taxid = cumulative.taxid
          WHERE cumulative.taxid != 1 
            AND cumulative.taxid != 131567
          GROUP BY rank
        """
        result = self.db.execute(cmd).fetchall()
        return {r: c for r, c in result}

    def get_discriminative_peptides_from_rank(self, rank):
        """
        Retrieve all discriminative peptides at and below the specified rank.
        """
        logging.debug("Retrieving discriminative peptides at and below rank %s...", rank)

        rank_set = self.rank_hierarchy[self.ranks[rank]:]

        cmd = """SELECT peptide, rank, spname FROM peptides
          JOIN taxref.species ON taxref.species.taxid = peptides.discriminative_taxid
          WHERE rank IN ({})
          ORDER BY rank""".format(",".join("?"*len(rank_set)))
        result = self.db.execute(cmd, rank_set).fetchall()
        return result

    def get_discriminative_counts_from_rank(self, rank):
        """
        Retrieve cumulative peptide counts for all ranks at or below 'rank'.
        """
        logging.debug("Getting cumulative counts at and below rank: '%s'...", rank)

        rank_set = self.rank_hierarchy[self.ranks[rank]:]

        cmd = """
        SELECT cumulative.count, COALESCE(disc_pep_count, 0), rank, species.spname
        FROM cumulative
        JOIN species
        ON cumulative.taxid = taxref.species.taxid
        LEFT JOIN
          (SELECT COUNT(peptide) AS disc_pep_count, taxref.species.spname 
           FROM peptides 
           JOIN taxref.species 
           ON taxref.species.taxid = peptides.discriminative_taxid
           GROUP BY rank, taxref.species.spname
          ) AS disc_per_spname
        ON disc_per_spname.spname = taxref.species.spname
        WHERE taxref.species.rank in ({})
        AND cumulative.count > 0
        GROUP BY rank, taxref.species.spname
        ORDER BY cumulative.count DESC
        """.format(",".join("?"*len(rank_set)))

        result = self.db.execute(cmd, rank_set).fetchall()
        return result
    
    def get_discriminative_hits_to_annotated_regions_from_rank(self, rank):
        """
        Retrieve all annotation regions matched by discriminative peptides 
        from given rank.
        """

        rank_set = self.rank_hierarchy[self.ranks[rank]:]

        cmd = """
        SELECT spname, disc.rank, Count(product) as prod_count, product, features FROM 
          (SELECT peptide, rank FROM peptides
                  JOIN taxref.species ON taxref.species.taxid = peptides.discriminative_taxid
                  WHERE rank IN ({})) as disc
        JOIN mappings 
            ON mappings.peptide = disc.peptide
        JOIN annotationdb.annotations 
            ON annotationdb.annotations.header = mappings.target
                AND (mappings.start BETWEEN 
                     annotationdb.annotations.start AND annotationdb.annotations.end 
                     OR 
                     mappings.end BETWEEN 
                     annotationdb.annotations.start AND annotationdb.annotations.end
                     OR (mappings.start <= annotationdb.annotations.start 
                      AND mappings.end >= annotationdb.annotations.end
                     )
                )
        JOIN taxref.refseqs ON mappings.target = taxref.refseqs.header
        JOIN taxref.species ON taxref.refseqs.taxid = taxref.species.taxid
        GROUP BY spname, product
        ORDER BY prod_count DESC
        """.format(",".join("?"*len(rank_set)))

        logging.debug("Getting annotated regions matched by discriminative peptides at rank %s", rank)
        result = self.db.execute(cmd, rank_set).fetchall()
        return result
    
    def get_hits_to_annotated_regions(self):
        """
        Retrieve all annotated regions matched by any peptide.
        """

        cmd = """
        SELECT spname, count(product) as prod_count, product, features 
        FROM annotationdb.annotations 
        JOIN mappings 
        ON mappings.target = annotationdb.annotations.header
          AND (mappings.start 
            BETWEEN annotationdb.annotations.start AND annotationdb.annotations.end
            OR mappings.end 
            BETWEEN annotationdb.annotations.start AND annotationdb.annotations.end
            OR (mappings.start <= annotationdb.annotations.start 
              AND mappings.end >= annotationdb.annotations.end
            )
          )
        JOIN taxref.refseqs ON mappings.target = taxref.refseqs.header
        JOIN taxref.species ON taxref.refseqs.taxid = taxref.species.taxid
        GROUP BY spname, product
        ORDER BY prod_count DESC
        """
        logging.debug("Retrieving all annotated regions matched by any peptide...")
        result = self.db.execute(cmd).fetchall()
        return result


def parse_normalization_factors(filename):
    """
    Parse a tab separated file with per-species normalization factors.

    Returns a dictionary with per-species normalization factors.
    """
    normalization_factors = {}

    # TODO: The following if-statement is in place to enable
    # future development of TCUP where a file with default 
    # normalization factors can automatically be loaded if
    # the user does not specify one. 
    if not filename:
        logging.debug("Not parsing normalization factors.")
        return normalization_factors  # TODO: Remove this return statement in the future
        # Import and read the default normalization factor file 
        # included in the TCUP package distribution.
        import pkg_resources
        filename = pkg_resources.resource_filename("tcup", "species_normalization.tab")
    logging.debug("Parsing species normalization factors from %s", filename)
    with open(filename) as f:
        for species_line in f:
            try:
                species, factor = species_line.strip().split("\t")
                normalization_factors[species] = float(factor)
            except ValueError:
                logging.error("Could not parse %s. The offending line was:\n%s", filename, species_line)
                logging.warning("Continuing without normalization factors.")
    logging.debug("Available normalization factors: %s", normalization_factors)
    return normalization_factors


def compute_corrected_and_normalized_cumulative_percentages(disc_peps_per_rank, rank_counts, normalization_factors):
    """
    Compute corrected and normalized cumulative percentages for each species.

    Note that this only includes species with specified normalization factors!
    """

    percentages = {spname: (cum_count/rank_counts[rank], rank) for cum_count, count, rank, spname in disc_peps_per_rank}
    pre_normalized_corrected_percentages = {}
    for spname, percentage_rank in percentages.items():
        percentage, rank = percentage_rank
        if spname in normalization_factors:
            pre_normalized_corrected_percentages[spname] = percentage/normalization_factors[spname]
    pre_normalization_sum = sum(p for p in pre_normalized_corrected_percentages.values())
    normalized_percentages = {spname: p/pre_normalization_sum for spname, p in pre_normalized_corrected_percentages.items()}

    return normalized_percentages


def filter_cumulative_discriminative_counts(disc_peps_per_rank, rank_counts, min_count):
    """
    Filter out low-abundance taxa from the results listing.
    """

    filtered_disc_peps_per_rank = [taxon for taxon in disc_peps_per_rank if taxon[0] > min_count]
    filtered_rank_counts = defaultdict(int)
    for taxon in filtered_disc_peps_per_rank:
        filtered_rank_counts[taxon[2]] += taxon[0]
    num_filtered = len(disc_peps_per_rank) - len(filtered_disc_peps_per_rank)

    return filtered_disc_peps_per_rank, filtered_rank_counts, num_filtered


def print_cumulative_discriminative_counts(disc_peps_per_rank, rank_counts, normalization_factors, outfile):
    """
    Print sorted lists of discriminative peptide counts.
    """

    try:
        disc_peps_per_rank.sort(key=lambda c: c[0]/rank_counts[c[2]], reverse=True)
    except KeyError:
        # This happens when there are no discriminative peptides at/below the
        # requested level; normally only the taxonomic nodes "cellular organisms"
        # and "root" that both have official rank 'no rank', but should not be
        # included anyway.
        pass
    

    print("Discriminative peptides per spname".center(60, "-"), file=outfile)

    if normalization_factors:
        normalized_percentages = compute_corrected_and_normalized_cumulative_percentages(
                disc_peps_per_rank, rank_counts, normalization_factors)
        print("{:<10} {:<6} {:>6} {:>7} {:<20} {:<40}".format(
            "Cumulative", "Count", "%", "% corr.", "Rank", "Description"), file=outfile)
    else:
        print("{:<10} {:<6} {:>6} {:<20} {:<40}".format(
            "Cumulative", "Count", "%", "Rank", "Description"), file=outfile)

    for cum_count, count, rank, spname in disc_peps_per_rank:
        if spname in ("root", "cellular organisms"):
            continue
        percentage = cum_count/rank_counts[rank] * 100
        if normalization_factors:
            try:
                corrected_percentage = normalized_percentages[spname] * 100 
            except (KeyError, TypeError):
                corrected_percentage = float("NaN")
            print("{:<10} {:<6} {:>6.2f} {:>7.2f} {:<20} {:<40}".format(
                cum_count, count, percentage, corrected_percentage, rank, spname), file=outfile)
        else:
            print("{:<10} {:<6} {:>6.2f} {:<20} {:<40}".format(
                cum_count, count, percentage, rank, spname), file=outfile)



def write_discriminative_peptides(discriminative, outfilename):
    """
    Print name and assignment of discriminative peptides.
    """
    with open(outfilename, 'w') as outfile:
        print("Discriminative peptides".center(60, "-"), file=outfile)
        print("{:<20} {:<15} {:<30}".format("Peptide", "Rank", "Description"), file=outfile)
        for d in discriminative: 
            print("{:<20} {:<15} {:<30}".format(*d), file=outfile)


def print_annotation_hits(hits, outfile):
    """
    Print hits to annotated genome regions.
    """
    print("Hits to annotated genome regions".center(60, "-"), file=outfile)
    print("{:<40}\t{:<10}\t{:<7}\t{:<45}\t{:<}".format("Spname", "Disc. level", "Count", "Product", "Features"), file=outfile)
    for spname, disc_level, count, product, features in hits:
        print("{:<40}\t{:<10}\t{:<7}\t{:<45}\t{:<}".format(spname, disc_level, count, product, features), file=outfile)


def write_results_xlsx(disc_peps_per_rank, rank_counts, hits, results_filename, normalization_factors):
    """
    Write results to an Excel xlsx file.
    """
    workbook = xlsxwriter.Workbook(results_filename)
    worksheet_composition= workbook.add_worksheet("Taxonomic composition")

    percentage_format = workbook.add_format({"num_format": "0.00%"})

    if normalization_factors:
        normalized_percentages = compute_corrected_and_normalized_cumulative_percentages(
                disc_peps_per_rank, rank_counts, normalization_factors)
        worksheet_composition.write(0, 0, "Cumulative")
        worksheet_composition.write(0, 1, "Discriminative count")
        worksheet_composition.write(0, 2, "Percentage")
        worksheet_composition.write(0, 3, "% Normalized")
        worksheet_composition.write(0, 4, "Rank")
        worksheet_composition.write(0, 5, "Description")
        worksheet_composition.set_column(0, 0, 9.0)
        worksheet_composition.set_column(1, 1, 17.0)
        worksheet_composition.set_column(2, 2, 9.0)
        worksheet_composition.set_column(3, 3, 11.0)
        worksheet_composition.set_column(4, 4, 11.0)
        worksheet_composition.set_column(5, 5, 40.0)
    else:
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
    
    row_adjustment = 0
    for row, data in enumerate(disc_peps_per_rank, start=1):
        row = row - row_adjustment
        cum_count, count, rank, spname = data
        if spname == "root" or spname =="cellular organisms":
            row_adjustment += 1
            continue
        if normalization_factors:
            percentage = cum_count/rank_counts[rank]
            try:
                corrected_percentage = normalized_percentages[spname] * 100
            except (KeyError, TypeError):
                corrected_percentage = -0.0 
            worksheet_composition.write(row, 0, cum_count)
            worksheet_composition.write(row, 1, count)
            worksheet_composition.write(row, 2, percentage, percentage_format)
            worksheet_composition.write(row, 3, corrected_percentage, percentage_format)
            worksheet_composition.write(row, 4, rank)
            worksheet_composition.write(row, 5, spname)
        else:
            percentage = cum_count/rank_counts[rank]
            worksheet_composition.write(row, 0, cum_count)
            worksheet_composition.write(row, 1, count)
            worksheet_composition.write(row, 2, percentage, percentage_format)
            worksheet_composition.write(row, 3, rank)
            worksheet_composition.write(row, 4, spname)
    try:
        worksheet_composition.autofilter(0, 0, row, 4)
    except NameError:
        logging.debug("Not applying autofilter: No rows in xlsx file")

    worksheet_annotations = workbook.add_worksheet("Hits to annotated regions")
    worksheet_annotations.write(0, 0, "Species")
    worksheet_annotations.write(0, 1, "Disc. level")
    worksheet_annotations.write(0, 2, "Count")
    worksheet_annotations.write(0, 3, "Product")
    worksheet_annotations.write(0, 4, "Features")
    worksheet_annotations.set_column(0, 0, 53.0)
    worksheet_annotations.set_column(1, 1, 10.0)
    worksheet_annotations.set_column(2, 2, 8.0)
    worksheet_annotations.set_column(3, 3, 53.0)
    for row, data in enumerate(hits, start=1):
        spname, disc_level, product_count, product, features = data
        worksheet_annotations.write(row, 0, spname)
        worksheet_annotations.write(row, 1, disc_level)
        worksheet_annotations.write(row, 2, product_count)
        worksheet_annotations.write(row, 3, product)
        worksheet_annotations.write(row, 4, features)
    
    workbook.close()


def get_results_from_existing_db(sample_databases,
        annotation_db_file, 
        taxonomic_rank="family",
        discpeps_file=False,
        print_annotations=False,
        write_xlsx="",
        normalization_factors="",
        outfile=stdout,
        inclusion_threshold=0):
    """
    Retrieve and print results from existing sample database(s).
    """

    for sample_db in sample_databases:
        sample_db.attach_annotation_db(annotation_db_file)

        disc_peps_per_rank = sample_db.get_discriminative_counts_from_rank(taxonomic_rank)
        rank_counts = sample_db.get_cumulative_rank_counts()

        if inclusion_threshold:
            disc_peps_per_rank, rank_counts, num_filtered = filter_cumulative_discriminative_counts(disc_peps_per_rank, rank_counts, inclusion_threshold)
            logging.info("%s taxa were below inclusion threshold (%s) and thus discarded from the results.", num_filtered, inclusion_threshold)

        print(sample_db.dbfile.center(60, "-"), file=outfile)
        print_cumulative_discriminative_counts(disc_peps_per_rank, rank_counts, normalization_factors, outfile)

        if write_xlsx or print_annotations:
            hits = sample_db.get_discriminative_hits_to_annotated_regions_from_rank(taxonomic_rank)
        if write_xlsx:
            xlsx_filename = write_xlsx
            write_results_xlsx(disc_peps_per_rank, rank_counts, hits, xlsx_filename, normalization_factors)
        if print_annotations:
            print_annotation_hits(hits, outfile)

        if discpeps_file:
            disc = sample_db.get_discriminative_peptides_from_rank(taxonomic_rank)
            write_discriminative_peptides(disc, discpeps_file)


def run_complete_pipeline(options):
    """
    The complete pipeline logic.
    """
    blacklisted_seqs = prepare_blacklist(options.blacklist, options.leave_out)
    normalization_factors = parse_normalization_factors(options.species_normalization)

    for blat_file in options.FILE:
        if blat_file.endswith(".sqlite3"):
            sleep_time = 5
            logging.warning("{} ends with .sqlite3, did you mean to use --sample-db?".format(blat_file))
            logging.warning("Press Ctrl+c to cancel within %i sec...", sleep_time)
            time.sleep(sleep_time)

        if options.sample_db:
            sample_db_filename = options.sample_db
        else:
            sample_db_filename = blat_file + ".sqlite3"
        sample_db = Sample_DB_wrapper(sample_db_filename)

        blat_parser = parse_blat_output(blat_file, 
                options.min_identity, 
                options.min_matches, 
                options.min_coverage,
                options.max_pid_diff, 
                blacklisted_seqs)
        sample_db.insert_blat_hits_into_db(blat_parser)

        sample_db.attach_taxref_db(options.taxref_db)
        sample_db.determine_discriminative_ranks()
        sample_db.count_discriminative_per_rank()

        if options.output:
            output = open(options.output, 'w')
        else:
            output = stdout
        with output:
            get_results_from_existing_db([sample_db],
                    options.annotation_db_file,
                    options.taxonomic_rank,
                    options.write_discriminative_peptides,
                    options.print_annotations,
                    options.write_xlsx,
                    normalization_factors,
                    output,
                    options.inclusion_threshold)

def main():
    """
    Main function.
    """

    options = parse_commandline(argv)
    normalization_factors = parse_normalization_factors(options.species_normalization)

    if options.pre_existing:
        sample_databases = [Sample_DB_wrapper(filename, create_new=False) for filename in options.FILE]
        for sample_db in sample_databases:
            sample_db.attach_taxref_db(options.taxref_db)
        if options.output:
            output = open(options.output, 'w')
        else:
            output = stdout
        with output:
            get_results_from_existing_db(sample_databases,
                options.annotation_db_file,
                options.taxonomic_rank,
                options.write_discriminative_peptides,
                options.print_annotations,
                options.write_xlsx,
                normalization_factors,
                output,
                options.inclusion_threshold)
    else:
        run_complete_pipeline(options)

if __name__ == "__main__":
    main()

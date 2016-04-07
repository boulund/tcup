#!/usr/bin/env python3.5
# Fredrik Boulund
# 2015-11-15
# Utility functions

import sqlite3
import itertools
import time
import re
import logging
import argparse
import fnmatch
import os


def existing_file(f):
    """ 
    Check that a file exists.

    :param f:  filename or path to file.
    :return:  None if f is None, filename or path to file if it exists.
    :raises ArgumentError:  if file does not exist. 
    """
    if f is None:
        return None
    elif not os.path.exists(f) or not os.path.isfile(f):
        raise argparse.ArgumentError("{} does not exist.".format(f))
    return f


def read_fasta(filename, keep_formatting=True):
    """
    Generator that reads sequence records from FASTA file.

    NOTE: This is a generator, it yields after each sequence record. 
    Usage example::
    for header, seq in read_fasta(filename):
        print ">"+header
        print seq

    :param filename: string with path to FASTA file.
    :param keep_formatting: boolean determining whether to remove line
            breaks from headers and sequences. Default: True.
    :return: yields (header, sequence) tuples.
    """

    with open(filename) as fasta:
        line = fasta.readline().rstrip()
        if not line.startswith(">"):
            raise IOError("Not FASTA format? First line didn't start with '>'")
        if keep_formatting:
            sep = "\n"
        else:
            sep = ""
        first = True
        seq = []
        header = ""
        while fasta:
            if line == "": #EOF
                yield header, sep.join(seq)
                break
            elif line.startswith(">") and not first:
                yield header, sep.join(seq)
                header = line.rstrip()[1:]
                seq = []
            elif line.startswith(">") and first:
                header = line.rstrip()[1:]
                first = False
            else:
                seq.append(line.rstrip())
            line = fasta.readline()


def find_files(directory, pattern):
    """
    Recursively search a dir with a glob pattern yielding paths to matching filenames.

    :param directory:  Path to directory from which to start walking.
    :param pattern:  a glob pattern used select what files to find.
    :return:  yields path to each filename that matches the glob pattern.
    """
    for root, subfolders, files in os.walk(directory, followlinks=True):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


def grouper(n, iterable):
    """
    Groups an iterable into n-sized chunks.

    :param n:  number of items per chunk.
    :param iterable:  an iterable to divide into chunks of size n.
    :return:  n-sized chunks from the iterable.
    """
    it = iter(iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           return
       yield chunk


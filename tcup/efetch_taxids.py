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
from requests import get
from multiprocessing import Pool


def parse_accnos(filename):
    """
    Parse accnos from file.
    """

    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                if "ref|" in line:
                    accno = line.split("ref|")[1].split("|", 1)[0]
                else:
                    accno = line.split()[0][1:]
            else:
                accno = line.strip()
            yield accno


def efetch_taxid(accno, retry_attempts=2):
    """
    Use NCBI E-utils efetch to get taxid for accno.
    Returns an empty string if no match is found for accno.
    """
    payload = {"db": "nuccore", 
               "id": accno,
               "rettype": "fasta",
               "retmode": "xml"}

    taxid = ""
    for attempt in range(0,retry_attempts):
        xml = get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=payload)
        try:
            taxid = xml.text.split("taxid>")[1].split("<")[0]
            break
        except IndexError:
            continue
    if taxid:
        return accno, taxid
    else:
        return accno, ""


if __name__ == "__main__":
    if len(argv) < 2:
        print("usage: efetch_taxids.py FILE")
        print("Fredrik Boulund 2016")
        print("Fetch TAXIDs from Genbank using ACCNOs.")
        print("FILE should contain one ACCNO (or header) per line.")
        exit()

    accnos_file = argv[1]
    accnos = list(parse_accnos(accnos_file))

    p = Pool(20)
    taxids = p.map(efetch_taxid, accnos)

    for accno, taxid in taxids:
        print("{}\t{}".format(accno, taxid))

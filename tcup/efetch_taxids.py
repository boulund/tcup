#!/usr/bin/env python3.5
# Fredrik Boulund 2015
# Fetch taxids from GenBank via ACCNO
# Input is list of ACCNO 

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

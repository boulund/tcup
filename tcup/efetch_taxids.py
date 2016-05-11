#!/usr/bin/env python3.5
# Fetch taxids from GenBank by ACCNO
# Simple utility script without any error handling.
# Fredrik Boulund 2016

from sys import argv, exit
from requests import get

if len(argv) < 2:
    print("usage: efetch_taxids.py FILE")
    print("Fetch TAXIDs from Genbank using ACCNOs.")
    print("FILE should contain one ACCNO per line.")
    print("Fredrik Boulund 2016")
    exit()

with open(argv[1]) as f:
    for line in f:
        if line.startswith(">"):
            accno = line.split()[0][1:]
        else:
            accno = line.strip()

        payload = {"db": "nuccore", 
                   "id": accno,
                   "rettype": "fasta",
                   "retmode": "xml"}
        xml = get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=payload)
        try:
            taxid = xml.text.split("taxid>")[1].split("<")[0]
        except IndexError:
            print("{}\t{}".format(accno, ""))
        print("{}\t{}".format(accno, taxid))

# Proteotyping

Author: Fredrik Boulund
License: MIT

This project performs *proteotyping* using peptide data produced with tandem
mass spectrometry (MS/MS). 


# Installation 
I recommend using Anaconda Python 3.5. Download and install Anaconda Python 3.5
and create an anaconda environment using the following command:

    $ conda create -name proteotyping_env python=3.5 .... [more packages!]

## Dependencies
The proteotyping program depends on the following excellent Python packages,
most of which are easily installable via `conda` and `pip`.

 * ETE3

# How to run the program
In order to run the program, some databases need to be prepared.

## Building databases
Construct reference databases containing an expanded NCBI Taxonomy, and
annotations of all reference genome sequences.

## Use BLAT to map peptides to reference genome sequences
Run BLAT to map peptides in the sample to reference genomes sequences in FASTA
format. Make sure to store the output in `blast8` tabular format.

## Run proteotyping on the mapping results
The proteotyping program requires the following input:

 * [sample.blast8]: BLAT output file (blast8 format)
 * [proteodb.sql]: proteotyping database (expanded NCBI Taxonomy constructed in previous step)
 * [annotationdb.sql]: annotation database (constructed in previous step)

Optionally, a file containing a set of sequence headers to blacklist can be
specified using the `--blacklist` command line option.


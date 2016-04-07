Running Proteotyping
====================
Before running |name|, ensure that the requisite databases have been created
(see :doc:`preparing_databases`). The input data to |name| are mappings of
peptides from against the reference genomes database you used to create the
`taxref DB` (see :doc:`preparing_databases`). 

.. note::
    **IMPORTANT:** The FASTA headers of the peptides need to contain the
    peptide length. It should be encoded at the end of the header as a single 
    integer after an underscore character ("_") at the end of the first space
    separated part of the header. E.g. the following header::

        >peptidename_15 some other information 

    belongs to a peptide called `peptidename` that is 15 amino acids long.

The mappings of peptides to reference genome sequences must be in blast8
tabular format (without column headers or comments). We recommend using 
`BLAT`_ or `pBLAT`_. 

.. _BLAT: https://genome.ucsc.edu/FAQ/FAQblat.html
.. _pBLAT: http://icebert.github.io/pblat/


Taxonomic composition
*********************
A typical invocation might look something like this::

    taxonomic_composition.py \
        --taxref-db <TAXREF_DB> \
        --annotation-db <ANNOTATION_DB> \
        --write-xlsx <OUTPUT_XLSX_FILENAME> \
        --output <OUTPUT_TXT_FILENAME>
        <INPUT_BLAST8_FILENAME>

where ``<TAXREF_DB>`` is the path to the taxonomy reference database file,
``<ANNOTATION_DB>`` is the path to the annotation database file,
``<OUTPUT_XLSX_FILENAME>`` is desired filename for the output in Excel format,
``<OUTPUT_TXT_FILENAME>`` is the desired filename for the output in text
format, and ``<INPUT_BLAST8_FILENAME>`` is the filename of a file in BLAST8
format containing the mapping results of sample peptides against the reference
genomes.


Antibiotic resistance
*********************
A typical invocation might look something like this::

   parse_AR_blast8.py \
       --resfinder <RESFINDER_DB> \
       --output <OUTPUT_TXT_FILENAME> \
       <INPUT_BLAST8_FILENAME>

where ``<RESFINDER_DB>`` is the path to the ResFinder sqlite3 database file,
``<OUTPUT_TXT_FILENAME>`` is the desired filename for the text format output,
and ``<INPUT_BLAST8_FILENAME>`` is the filename of a file in BLAST8 format
containing the mapping results of peptides against the ResFinder database.

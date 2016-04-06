Running Proteotyping
====================
Before running |name|, ensure that the requisite databases have been
created (see :doc:`database_preparation`).


Taxonomic composition
*********************
A typical invocation might look something like this::

    taxonomic_composition.py \
        --taxref-db <TAXREF_DB> \
        --annotation-db <ANNOTATION_DB> \
        --write-xlsx <OUTPUT_XLSX_FILENAME> \
        --output <OUTPUT_TXT_FILENAME>
        <INPUT_BLAST8_FILENAME>

where `<TAXREF_DB>` is the path to the taxonomy reference database file,
`<ANNOTATION_DB>` is the path to the annotation database file,
`<OUTPUT_XLSX_FILENAME>` is desired filename for the output in Excel format, 
`<OUTPUT_TXT_FILENAME>` is the desired filename for the output in text format,
and `<INPUT_BLAST8_FILENAME>` is the filename of a file in BLAST8 format
containing the mapping results of peptides against the reference genomes.


Antibiotic resistance
*********************
A typical invocation might look something like this::

   parse_AR_blast8.py \
       --resfinder <RESFINDER_DB> \
       --output <OUTPUT_TXT_FILENAME> \
       <INPUT_BLAST8_FILENAME>

where `<RESFINDER_DB>` is the path to the ResFinder sqlite3 database file,
`<OUTPUT_TXT_FILENAME>` is the desired filename for the output in text format,
and `<INPUT_BLAST8_FILENAME>` is the filename of a file in BLAST8 format
containing the mapping results of peptides against the ResFinder database.

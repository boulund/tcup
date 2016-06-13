Running |name|
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
`BLAT`_ or `pBLAT`_. On Windows, `BLAST`_ can be used.

.. _BLAT: https://genome.ucsc.edu/FAQ/FAQblat.html
.. _pBLAT: http://icebert.github.io/pblat/
.. _BLAST: http://www.ncbi.nlm.nih.gov/books/NBK52637/

.. note::
    All examples of command invocations in this section are aimed towards
    Linux. Where notable differences exist, they will be pointed out.


BLAT/pBLAT
**********
The peptides in the sample need to mapped to the reference genomes. We
recommend the use of `BLAT`_ for this, as it performs very well with translated
mapping to large genome sequences. It is important to tell BLAT to do
protein-to-translated-dna mapping with the ``-t=dnax -q=prot`` arguments.  We
have successfully used the following command line options with tandem mass
spectrometry data::

    blat <ref_db.fasta> <sample.fasta> -t=dnax -q=prot -minScore=10 -stepSize=5 -tileSize=5 -minIdentity=90 -out=blast8 <outfilename.blast8>

For the antibiotic resistance detection, the high sensitivity afforded by the
above settings are not required, also the mapping is now protein-to-protein.
Instead of using the above command line for antibiotic resistance detection, we
recommend using the default settings, like so::

    blat <resfinder.fasta> <sample.fasta> -out=blast8 <outfilename.blast8>


BLAST
*****
We do not recommend using BLAST, but if no other option is available, it can be done
using the following settings for taxonomic composition estimation::

    tblastn.exe -query <query> -db <db> -out <outfile> -outfmt 6

For antibiotic resistance detection, we need to do protein-to-protein mapping using
the following command line::

    blastp.exe -query <query> -db <db> -out <outfile> -outfmt 6

.. note::

    |name| has not yet been optimized for use with BLAST, and the results might be
    be unreliable if BLAST is used for peptide to genome mapping. 


Taxonomic composition
*********************
A typical invocation might look something like this::

    taxonomic_composition \
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

Running ``taxonomic_composition --help`` will show a full list of all available
options and their default values.


Antibiotic resistance
*********************
A typical invocation might look something like this::

   antibiotic_resistance \
       --resfinder <RESFINDER_DB> \
       --output <OUTPUT_TXT_FILENAME> \
       <INPUT_BLAST8_FILENAME>

where ``<RESFINDER_DB>`` is the path to the ResFinder sqlite3 database file,
``<OUTPUT_TXT_FILENAME>`` is the desired filename for the text format output,
and ``<INPUT_BLAST8_FILENAME>`` is the filename of a file in BLAST8 format
containing the mapping results of peptides against the ResFinder database.

Running ``antibiotic_resistance --help`` will show a full list of all available
options and their default values.

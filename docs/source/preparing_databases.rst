Preparing databases for use with |name|
=======================================

|name| relies on several reference databases that need to be prepared in
order to analyze a sample.

For creating the databases required for taxonomic composition estimation, the
following reference data sources are required:

* Reference genome sequences (FASTA) 
* Reference genome annotations (GFF) 
* NCBI Taxonomy taxdump (tab separated)  *[optional]*

The information from these reference data sources is combined into two SQLite3
databases that the taxonomic composition estimation program uses when
estimating the taxonomic composition of a sample. These SQLite3 databases are
referred to as `taxref DB` and `annotation DB` throughout this documentation.
How to prepare these reference databases are described in `Taxref DB`_ and
`Annotation DB`_.

Additionally, for creating the database used for antibiotic resistance
detection, the following reference data sources are required:

* Antibiotic resistance gene database 

Note that the detection of expressed antibiotic resistance proteins and the
estimation of taxonomic composition are two independent programs. This means
that if you only require one of them, you only need to create the databases
associated with the program you want to use.
  

Taxref DB
*********
The `taxref DB` is central to determining the taxonomic composition of samples.
The database is built from `NCBI Taxonomy`_ and contains information linking
the reference genome sequences to nodes in the taxonomy. 

To create the `taxref DB`, a single user-created input file is required, called
`header_mappings`. This file is a tab separated text file with two columns::

    <FASTA_HEADER>  <TAXID> 

.. note::
    <FASTA_HEADER> should only contain the first part of the FASTA header up to
    the first space (also excluding the ">" angle bracket symbol).
    For example, for a FASTA header like this::

        >gi|158333233|ref|NC_009925.1| Acaryochloris marina MBIC11017 chromosome, complete genome

    the corresponding line in the `header_mappings` file should be::
        
        gi|158333233|ref|NC_009925.1| <TAB> 329726
 
    
Header mappings can be created for NCBI RefSeq sequences using the
`taxref_db.py header_mappings` subprogram described in `Header mappings`_
below. It is also possible to create this mapping yourself, however you want
to. The only important thing to note is to make sure that the first part of the
FASTA header (up to the first space) is mapped to a valid NCBI taxid. As the
accuracy of the taxonomic composition estimation is wholly dependent on correct
assignments in the `taxref DB`, it is imperative that these assignments are as
correct as possible.

The simplest invocation to create a `taxref DB` looks like this::
    
    taxref_db.py  taxref_db  <HEADER_MAPPINGS>

where ``<HEADER_MAPPINGS>`` is the path to a `header_mappings` file as
described above. The program allows multiple header mapping files to be
specified on the command line. The program will automatically download the most
recent version of NCBI Taxonomy to build the database from.


Header mappings
---------------
A `header_mappings` file is required to create a `taxref DB`.  
For NCBI RefSeq sequences, it is can be created like this::

    taxref_db.py  header_mappings  <REFSEQ_DIR>  <GI_TAXID_DUMP>

where ``<REFSEQ_DIR>`` is the path to a directory containing NCBI RefSeq
sequences in FASTA format (default having an ``*.fna`` extension), and
``<GI_TAXID_DUMP>`` is the path to a tab separated file with two columns
containing `gi` to `taxid` mappings.

.. note::
    Since early 2016, NCBI are phasing out the use of `gi` numbers. The
    ``taxref_db.py`` program will automatically search for taxid mappings for a
    given RefSeq sequence if the FASTA header contains a sequence accession
    number, e.g. ``ref|NC_009925.1|``. However, expect this to be quite slow 
    as it has to make queries for all accession numbers of the network.


.. _`NCBI Taxonomy`: http://www.ncbi.nlm.nih.gov/taxonomy


Annotation DB
*************
The `annotation DB` is used after estimating the taxonomic composition to
present what annotated regions in the reference genome sequences were matched
by discriminative peptides in the sample. The annotation matches are presented 
in the xlsx (Excel) output, and can optionally be printed in the text file output
using the ``--print-annotations`` flag to ``taxonomic_composition.py``. 

To create the `annotation DB`, the database creation program ``annotation_db.py``
parses GFF files (General Feature Format). The simplest invocation to create
an `annotation DB` looks like this::

    annotation_db.py  <TAXREF_DB>  <GFF_DIR>

where ``<TAXREF_DB>`` is the path to a pre-made `Taxref DB`_, and ``<GFF_DIR>``
is the path to a directory containing GFF files to parse. The program allows
multiple directories to be specified on the command line, and will recursively
search subdirectories for GFF files as well.

.. note::
    Creating an `annotation_DB` requires a pre-made `taxref DB`. Also note that
    the program won't know what to do if it encounters GFF files for reference
    genome sequences that were not included in the creation of the `taxref DB`. 
    


Antibiotic resistance DB
************************
A modified version of `ResFinder`_ suitable for use with |name| is available
for download `here`_. 

.. _ResFinder: https://cge.cbs.dtu.dk//services/ResFinder/
.. _here: https://bioinformatics.math.chalmers.se/proteotyping/resfinder_20160304.zip



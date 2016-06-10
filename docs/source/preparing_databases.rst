Preparing databases for use with |name|
=======================================

|name| relies on several reference databases that need to be prepared in
order to analyze a sample.

For creating the databases required for taxonomic composition estimation, the
following reference data sources are required:

* Reference genome sequences (FASTA) 
* Reference genome annotations (GFF) 

The headers of the genome sequences along with their taxonomic affiliations are
combined with NCBI Taxonomy to form an SQLite3 database used to map genome sequences
to the taxonomy, along with species names etc. The genome annotations are put into
a separate SQLite3 database, used to retrieve annotation information if required.
These SQLite3 databases are referred to as `taxref DB` and `annotation DB`
throughout this documentation.  How to prepare these reference databases are
described in `Taxref DB`_ and `Annotation DB`_.

Additionally, for creating the database used for antibiotic resistance
detection, some kind of reference data sources for antibiotic resistance
genes is required, we recommend something like the publicly available `ResFinder`_
database. A version of the `ResFinder`_ database suitable for use with |name| is
available from our `website`_.


.. _ResFinder: https://cge.cbs.dtu.dk//services/ResFinder/
.. _website: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/resfinder.zip


Note that the detection of expressed antibiotic resistance proteins and the
estimation of taxonomic composition are actually two independent programs in
the |name| package. This means that if you only require one of them, you only
need to create the databases associated with the program you want to use.


.. note::
    All examples of commands below are for Linux environments. Slight changes
    to the displayed commands might be required if working in a Windows environment.
  

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

    the corresponding line in the `header_mappings` file should be (without the
    extra spaces around the TAB character)::
        
        gi|158333233|ref|NC_009925.1| <TAB> 329726
 
    
Header mappings can be created for NCBI RefSeq sequences using the
`efetch_taxids.py` subprogram described in `Header mappings`_
below. It is also possible to create this mapping yourself, however you want
to. The only important thing to note is to make sure that the first part of the
FASTA header (up to the first space) is mapped to a valid NCBI taxid. As the
accuracy of the taxonomic composition estimation is wholly dependent on correct
assignments in the `taxref DB`, it is imperative that these assignments are as
correct as possible.

The simplest invocation to create a `taxref DB` looks like this::
   
    taxref_db  <HEADER_MAPPINGS>

where ``<HEADER_MAPPINGS>`` is the path to a `header_mappings` file as
described above. The program allows multiple header mapping files to be
specified on the command line. The program will automatically download the most
recent version of NCBI Taxonomy to build the database from. Run ``taxref_db
--help`` for a listing of all the available options. 


Header mappings
---------------
A `header_mappings` file is required to create a `taxref DB`.  
For NCBI RefSeq sequences, it can be created using ``efetch_taxids.py`` to look
up taxids via NCBI E-utils. This example works on Linux only::

    $ grep ">" path/to/reference_genomes.fasta > reference_genome_headers.txt
    $ efetch_taxids reference_genome_headers.txt > header_mappings.tab

It is also perfectly possible to manually create this file, or use some other
method for pairing up the sequence headers to their respective taxid. Learn more
about taxids on `NCBI Taxonomy`_. 

.. _`NCBI Taxonomy`: http://www.ncbi.nlm.nih.gov/taxonomy


Annotation DB
*************
The `annotation DB` is used after estimating the taxonomic composition to
present what annotated regions in the reference genome sequences were matched
by discriminative peptides in the sample. The annotation matches are presented 
in the xlsx (Excel) output, and can optionally be printed in the text file output
using the ``--print-annotations`` flag to ``taxonomic_composition``. 

To create the `annotation DB`, the database creation program ``annotation_db``
parses GFF files (General Feature Format). The simplest invocation to create
an `annotation DB` looks like this::

    annotation_db  <TAXREF_DB>  <GFF_DIR>

where ``<TAXREF_DB>`` is the path to a pre-made `Taxref DB`_, and ``<GFF_DIR>``
is the path to a directory containing GFF files to parse. The program allows
multiple directories to be specified on the command line, and will recursively
search subdirectories for GFF files as well.

.. note::
    Creating an `annotation DB` requires that you have already completed a
    `taxref DB`. Also note that the program won't do anything with GFF files
    for reference genome sequences that were not included in the creation of
    the `taxref DB`. 
    


Antibiotic resistance DB
************************
To create your own antibiotic resistance gene database for use with |name|. Use
the ``construct_resfinder_db.py`` program supplied with |name|.  A typical
invocation might look like this::

    $ construct_resfinder_db.py --sequences <SEQUENCES> --notes <NOTES>

where ``<SEQUENCES>`` is the path to a FASTA file with antibiotic resistance
gene sequences. ``<NOTES>`` is the path to the ``notes.txt`` file included with
the `ResFinder`_ distribution.  

.. note:: 
    A manually curated version of `ResFinder`_ suitable for use with |name| is
    available for download from our `website`_. 

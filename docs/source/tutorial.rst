Tutorial
========
In this tutorial we will run |name| on a sample containing an antibiotic
resistant *E. coli*. The sample input file (peptides from bottom-up tandem mass
spectrometry) is available for download from here:
http://bioinformatics.math.chalmers.se/tcup/tutorial/tcup_tutorial_sample.fasta.zip


Preparations
************
First, ensure you have installed |name| according to the instructions in the
:doc:`install` section. Then, to run |name|, you need the following databases:

* `Reference genome database`_. BLAST or BLAT reference databases [FASTA or
  BLASTDB]. About 21 GB uncompressed.
* `Taxonomy/Reference`_ ("taxref") database. Linking reference genome sequences
  to taxonomy [sqlite3]. About 340 MB uncompressed.
* `Annotation database`_. Contains annotations for all genome sequences
  included in the reference genome database [sqlite3]. About 5 GB uncompressed.
* `Resistance gene database`_. Based on `ResFinder`_ [FASTA or BLASTDB and
  sqlite3]. About 1 MB uncompressed.

Example versions of these databases can be downloaded from the links in the
list above.  Note that the databases are quite large. The total download size
is approximately 10.3 GB.  There are instructions on how to create your own
databases in the :doc:`preparing_databases` section. Do not forget to download
the sample (linked in the section above).

.. Download sizes:
   9.1GB reference_genomes.zip
   1.2GB annotation_db.zip
   0.2MB resfinder.zip
   83 MB taxref.sqlite3.zip
   == 10.3GB


.. _Reference genome database: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/reference_genomes.zip
.. _Taxonomy/Reference: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/taxref.sqlite3.zip 
.. _Annotation database: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/annotation_db.sqlite3.zip
.. _Resistance gene database: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/resfinder.zip

.. _ResFinder: https://cge.cbs.dtu.dk//services/ResFinder/


Running |name|
**************
A simple way to run |name| is to use the included ``run_tcup.py`` script that
is installed as part of the |name| package. Running ``run_tcup.py`` without any
arguments (or with ``-h``/``--help``) produces this helpful output::

    $ run_tcup.py
    usage: run_tcup.py [-h] -t TAXREF_DB -a ANNOTATION_DB -r RESISTANCE_DB
                       FASTA GENOME_DB RESISTANCE_DB
    
    TCUP wrapper; align peptides to reference databases in parallel and run TCUP
    on alignment results. Fredrik Boulund 2016
    
    positional arguments:
      FASTA                 FASTA file with peptides from tandem MS.
      GENOME_DB             Path to reference bacterial genome db (FASTA or
                            blastdb format depending on OS).
      RESISTANCE_DB         Path to antibiotic resistance gene db (FASTA or
                            blastdb format depending on OS.
    
    optional arguments:
      -h, --help            show this help message and exit
    
    Taxonomic composition:
      -t TAXREF_DB, --taxref-db TAXREF_DB
                            Path to taxref db (sqlite3).
      -a ANNOTATION_DB, --annotation-db ANNOTATION_DB
                            Path to annotation db (sqlite3).
    
    Antibiotic resistance:
      -r RESISTANCE_DB, --resistance-db RESISTANCE_DB
                            Path to resistance db (sqlite3).

As indicated by the help text, we need to supply a number of arguments to
``run_tcup.py`` in order to run |name|.  The program requires the ``FASTA``
file containing the sample peptides as the first argument. The second argument
is the path to the ``GENOME_DB``. The third argument is the path to the
``RESISTANCE_DB``. In addition to these positional arguments, three additional
arguments are required: ``-t`` specifies the path to the ``taxref.sqlite3``
file, ``-a`` is the path to the ``annotation_db.sqlite3`` file, and ``-r`` is
the path to the ``resfinder.sqlite3`` file.

Assuming you have downloaded and extracted all the databases listed above, and
downloaded the tutorial sample FASTA file, into a folder like this::

   tutorial/
       tcup_tutorial_sample.fasta
       databases/
           annotation_db.sqlite3
           reference_genomes.fasta
           reference_genomes.fasta.00.nhr
           reference_genomes.fasta.00.nin
           reference_genomes.fasta.00.nsq
           reference_genomes.fasta.01.nhr
           reference_genomes.fasta.01.nin
           reference_genomes.fasta.01.nsq
           reference_genomes.fasta.02.nhr
           reference_genomes.fasta.02.nin
           reference_genomes.fasta.02.nsq
           reference_genomes.fasta.03.nhr
           reference_genomes.fasta.03.nhr
           reference_genomes.fasta.03.nin
           reference_genomes.fasta.04.nsq
           reference_genomes.fasta.04.nin
           reference_genomes.fasta.04.nsq
           reference_genomes.fasta.nal
           resfinder.fasta
           resfinder.fasta.00.nhr
           resfinder.fasta.00.nin
           resfinder.fasta.00.nsq
           resfinder.fasta.nal
           resfinder.sqlite3
           taxref.sqlite3


Running |name| is now as easy as::

   [tutorial]$ run_tcup.py tcup_tutorial_sample.fasta databases/reference_genomes.fasta databases/resfinder.fasta -t databases/taxref.sqlite3 -a databases/annotation_db.sqlite3 -r databases/resfinder.sqlite3

The above command will produce the output files::

    tcup_tutorial_sample.fasta.genomes.blast8
    tcup_tutorial_sample.fasta.ar.blast8
    tcup_tutorial_sample.fasta.antibiotic_resistance.txt
    tcup_tutorial_sample.fasta.taxonomic_composition.txt
    tcup_tutorial_sample.fasta.taxonomic_composition.xslx

In the next section we will look at the output from the analysis.

Analysis of the results
***********************



Congratulations, you have now completed the tutorial. The is more detailed
information on how to use |name| in the :doc:`running` section.



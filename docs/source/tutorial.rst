Tutorial
========
In this tutorial we will run |name| on a sample containing an antibiotic
resistant *E. coli*. The sample input file (peptides from bottom-up tandem mass
spectrometry) is available for download from here:
http://bioinformatics.math.chalmers.se/tcup/tutorial/tcup_tutorial_sample.fasta.zip


Preparations
************
Follow the installation instructions in :doc:`installation` to install |name|. 
When you have completed the installation, download the following databases that
are required to run |name|:

* `Reference genome database`_. BLAST or BLAT reference databases [FASTA or
  BLASTDB]. About 21 GB uncompressed.
* `Taxonomy`_ ("taxref") database. Links reference genome sequences
  to taxonomy structure [sqlite3]. About 340 MB uncompressed.
* `Annotation database`_. Contains annotations for all genome sequences
  included in the reference genome database [sqlite3]. About 5 GB uncompressed.
* `Resistance gene database`_. Based on `ResFinder`_ [FASTA or BLASTDB and
  sqlite3]. About 1 MB uncompressed.

Example versions of these databases can be downloaded from the links in the
list above.  Note that the databases are quite large. The total download size
is approximately 10.3 GB.  There are instructions on how to create your own
databases in :doc:`preparing_databases`. Do not forget to download the sample
(linked in the section above).

.. Download sizes:
   9.1GB reference_genomes.zip
   1.2GB annotation_db.zip
   0.3MB resfinder.zip
   83 MB taxref.zip
   == 10.3GB

.. _Reference genome database: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/reference_genomes.zip
.. _Taxonomy: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/taxref.zip 
.. _Annotation database: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/annotation_db.zip
.. _Resistance gene database: http://bioinformatics.math.chalmers.se/tcup/tutorial/databases/resfinder.zip

.. _ResFinder: https://cge.cbs.dtu.dk//services/ResFinder/

You also need to install either `BLAT`_ (Linux) or `BLAST`_ (Windows),
depending on your OS. Please refer to their respective installation
instructions for installing them.

.. _BLAT: https://genome.ucsc.edu/FAQ/FAQblat.html
.. _BLAST: http://www.ncbi.nlm.nih.gov/books/NBK52637/




Running |name|
**************
A simple way to run |name| is to use the included ``run_tcup.py`` script that
is installed as part of the |name| package. Running ``run_tcup`` without any
arguments (or with ``-h``/``--help``) produces this helpful output::

    $ run_tcup.py
    usage: run_tcup.py [-h] -t TAXREF_DB -a ANNOTATION_DB -r RESISTANCE_DB        
                       SAMPLE GENOME_DB RESISTANCE_DB                             
                                                                                  
    TCUP wrapper; align peptides to reference databases in parallel and run TCUP  
    on alignment results. Fredrik Boulund 2016                                    
                                                                                  
    positional arguments:                                                         
      SAMPLE                FASTA file with peptides from tandem MS.              
      GENOME_DB             Reference bacterial genome db (FASTA or blastdb format
                            depending on OS).                                     
      RESISTANCE_DB         Antibiotic resistance gene db (FASTA or blastdb format
                            depending on OS).                                     
                                                                                  
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
``run_tcup`` in order to run |name|.  The program requires the ``SAMPLE``
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
           reference_genomes.00.nhr
           reference_genomes.00.nin
           reference_genomes.00.nsq
           reference_genomes.01.nhr
           reference_genomes.01.nin
           reference_genomes.01.nsq
           reference_genomes.01.fasta
           reference_genomes.02.nhr
           reference_genomes.02.nin
           reference_genomes.02.nsq
           reference_genomes.02.fasta
           reference_genomes.03.nhr
           reference_genomes.03.nhr
           reference_genomes.03.nin
           reference_genomes.04.nsq
           reference_genomes.04.nin
           reference_genomes.04.nsq
           reference_genomes.nal
           resfinder.fasta
           resfinder.phr
           resfinder.pin
           resfinder.psq
           resfinder.sqlite3
           taxref.sqlite3


To run |name| on Windows, type the following command line (without linebreaks)::

   > run_tcup.exe
        -t databases\taxref.sqlite3 
        -a databases\annotation_db.sqlite3 
        -r databases\resfinder.sqlite3 
        tcup_tutorial_sample.fasta 
        databases\reference_genomes 
        databases\resfinder

To run |name| on Linux, type the following command line (without linebreaks)::

   $ run_tcup 
        -t databases/taxref.sqlite3 
        -a databases/annotation_db.sqlite3 
        -r databases/resfinder.sqlite3 
        tcup_tutorial_sample.fasta 
        databases/reference_genomes
        databases/resfinder.fasta 

After completing, |name| will produce the following output files::

    tcup_tutorial_sample.fasta.genomes.blast8
    tcup_tutorial_sample.fasta.ar.blast8
    tcup_tutorial_sample.fasta.antibiotic_resistance.txt
    tcup_tutorial_sample.fasta.taxonomic_composition.txt
    tcup_tutorial_sample.fasta.taxonomic_composition.xslx


.. note::

    TCUP is actually not intended to be run via the 'run_tcup' script as
    described in this section. The script is provided as a convenience to easily
    try out TCUP to see how it works, but for real world use of TCUP, please refer
    to :doc:`running`.

In the next section we will analyze the output from |name|.

Analysis of the results
***********************
.. note::
    
    NCBI BLAST produces more false positives than BLAT, and TCUP has only been
    optimized for use with BLAT at this time. The use of BLAST together with
    TCUP to determine taxonomic composition or expressed antibiotic resistance
    peptides is currently not recommended. Thus, if you are running TCUP on
    Windows, keep in mind that the results likely will contain a high number of
    false positive assignments, both for taxonomic affiliation and antibiotic
    resistance detection.

Taxonomic composition
---------------------
First off, let's have a look at the taxonomic composition of the sample. The
taxonomic composition estimation is presented in two formats: plain text and 
as an Excel spreadsheet. They both contain the same information regarding the 
taxonomic composition estimation of the sample, but the Excel file also includes 
a sheet with information on hits to annotated regions of the reference sequences.

The table in the first sheet of ``tcup_tutorial_sample.fasta.taxonomic_composition.xslx``
shows columns containing::

    Cumulative  Count   Percentage  Rank    Spname

The leftmost column, ``Cumulative``, shows the number of peptides that are
discriminative at the taxonomic rank specified in the ``Rank`` column. This
forms a cumulative sum as you look at ranks higher up in the taxonomic
hierarchy. If e.g. the rank of superkingdom was included in the results, it
would contain the total cumulative sum of the number of discriminative peptides
at all taxa in the bacterial tree. 

The ``Percentage`` column shows the relative proportion of peptides classified
to the species given in the ``Spname`` column. This number is relative to all
other entries of the same taxonomic rank, e.g. the sum of all the percentages
across all *species* would sum to 100%.

The Excel format makes it easy to use the filtering functions in Excel to look
at the most interesting parts of the results, e.g. to filter out only matches
to the *genus* or *species* levels. 

The second sheet in the Excel file contains a listing of all hits to regions in
the reference genome sequences that were matched by any discriminative peptide.



Antibiotic resistance
---------------------
Second, let's have a look at the antibiotic resistance results. These are presented
in a text file. The output contains four columns::

    Disc.  Hits   %    Family

The first column, ``Disc.``, shows the number of discriminative peptides that matched
to the resistance gene family listed in the ``Family`` column. The ``Hits`` column shows
how many separate matches the discriminative peptides produced to the family in question.
The ``%`` column shows the proportion of peptides that matched to each family.


Congratulations, you have now completed the tutorial. There is more detailed
information on how to use |name| in the :doc:`running` section.



Introduction and overview
=========================
|name| can estimate the taxonomic composition of a sample of microorganisms
using peptides generated from mass spectrometry.  The programs contained in the
|name| package are meant to be used in a larger proteotyping workflow like
this:

.. image:: img/overview.*
    :alt: Overview of how |name| fits into the greater picture of the proteotyping workflow.
    :align: center

In the above picture, the components of |name| are represented as the two boxes
inside the dashed light green area in the lower center of the picture. Before
|name| can be used to estimate the taxonomic composition or detect expressed
antibiotic resistance proteins in a sample, the following steps must be
performed after the mass spectrometry:

1. Match tandem mass spectra with peptide sequences (i.e. protein identification).
2. Map the identified peptides to a database of reference genome sequences.

Step 1 can be performed using any mass spectrometry search engine, such as
`X!Tandem`_, `Mascot`_, `MyriMatch`_, `Comet`_, `Tide`_, etc. It is recommended
to use the largest possible database of non-redundant protein sequences in this
step, to allow for maximal similarity between the identified peptide sequences
and the actual peptides in the sample.  
Step 2 can be performed using any sequence aligner/mapper capable of mapping
protein sequences to nucleotide references, producing results in blast8 tabular
format. We recommend using `BLAT`_ (or `pBLAT`_).

|name| uses several reference data sources (the two main ones depicted in the
picture above) to create its estimate of a sample's taxonomic composition and
detecting expressed antibiotic resistance proteins. More information on how to
prepare these reference databases is available in :doc:`preparing_databases`.


.. _X!Tandem: http://www.thegpm.org/tandem/
.. _Mascot: http://www.matrixscience.com/
.. _MyriMatch: https://medschool.vanderbilt.edu/msrc-bioinformatics/software
.. _Comet: http://comet-ms.sourceforge.net/
.. _Tide: https://noble.gs.washington.edu/proj/tide/
.. _BLAT: https://genome.ucsc.edu/FAQ/FAQblat.html
.. _pBLAT: http://icebert.github.io/pblat/ 

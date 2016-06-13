Example output
==============
This section contains some examples of what the output of |name| looks
like. |name| can output results in several formats. 


Taxonomic composition
*********************
The main output format intended for human consumption for taxonomic composition
estimation is ``xlsx``, i.e.  Microsoft Excel format. The Excel file is not
produced by default, and has to be requested on the command line using
``--write-xlsx XLSX_FILE``. The Excel output contains two sheets: *Taxonomic
composition* and *Hits to annotated regions*. These two outputs are described
in detail below.

Taxonomic composition
---------------------
A summary table of all discriminative peptide assignments are listed in this
sheet. The output can look like this (shortened for brevity)::

    Cumulative  Discriminative count    Percentage  Rank    Description
    921         738                     98.29%      genus   Streptococcus
    953         32                      97.94%      family  Streptococcaceae
    171         157                     94.34%      species Streptococcus pneumoniae
    7           7                       25.93%      no rank Streptococcus pneumoniae 70585
    5           5                       18.52%      no rank Streptococcus pneumoniae CGSP14
    2           2                       7.41%       no rank Streptococcus oralis Uo5
    1           1                       3.70%       no rank [Eubacterium] eligens ATCC 27750
    1           1                       3.70%       no rank Campylobacter concisus 13826
    1           0                       3.70%       no rank Clostridiales incertae sedis
    1           1                       3.70%       no rank Mycobacterium rhodesiae NBB3
    [...]

The columns, from left to right, are:

 * ``Cumulative``: The number of discriminative fragments at the reported Rank.
   These figures also include all peptides below the indicated Rank (hence
   *Cumulative*).
 * ``Discriminative count``: The number of discriminative peptides only at the
   indicated Rank. This figure does **not** include peptides on any other Rank.
 * ``Percentage``: The proportion of all discriminative peptides at the
   indicated Rank that belong to this identification. Note that these numbers
   are only comparable within a given Rank; i.e. the sum to 100% within each
   Rank and have no meaning when comparing two different ranks.
 * ``Rank``: The rank at which peptides were discriminative.
 * ``Description``: The strain, species, genus, family, etc. 

Looking at the results displayed above, we can see that this sample most likely
contains *Streptococcus pneumoniae*, as indicated by the strong signal at
species level (94.34%). This is further strengthened by the strong indicators
at both genus and family level. In this shortened output we see some peptides
that were discriminative to different strains of *S. pneumoniae*, but also two
peptides discriminative to *S. oralis*. Towards the end of the output, several
spurious matches to organisms that were probably not in the sample are also
visible.

.. note::

    The taxonomic composition output can also be written to a text file,
    using the ``--output`` command line option. If this option is not used,
    this output is printed to STDOUT.


Hits to annotated regions
-------------------------
Results in the *Hits to annotated regions* sheet are somewhat harder to 
analyze than the results in the *Taxonomic composition* sheet. In this 
sheet, the output might look like this (shortened for brevity)::

    Species                                 Disc. level Count   Product                             Features
    Streptococcus pneumoniae AP200          family      151     hypothetical protein                ID=cds741;Name=YP_003876344.1;Parent=gene747;Note=Operon 340 Gene 1 protein supported gi%7C148997148%7Cref%7CZP_01824802.1%7C ribosomal protein S16;Dbxref=Genbank:YP_003876344.1,GeneID:9726765;gbkey=CDS;product=hypothetical protein;protein_id=YP_003876344.1;transl_table=11
    Streptococcus pneumoniae R6             genus       96      hypothetical protein                ID=cds26;Name=NP_357627.1;Parent=gene36;Dbxref=Genbank:NP_357627.1,GeneID:933816;gbkey=CDS;product=hypothetical protein;protein_id=NP_357627.1;transl_table=11
    Streptococcus pneumoniae CGSP14         genus       81      hypothetical protein                ID=cds32;Name=YP_001834750.1;Parent=gene38;Dbxref=Genbank:YP_001834750.1,GeneID:6216513;gbkey=CDS;product=hypothetical protein;protein_id=YP_001834750.1;transl_table=11
    Streptococcus pneumoniae ST556          genus       78      hypothetical protein                ID=cds106;Name=YP_006252305.1;Parent=gene112;Dbxref=Genbank:YP_006252305.1,GeneID:12900645;gbkey=CDS;product=hypothetical protein;protein_id=YP_006252305.1;transl_table=11
    Streptococcus pneumoniae TIGR4          genus       78      hypothetical protein                ID=cds29;Name=NP_344583.1;Parent=gene38;Dbxref=Genbank:NP_344583.1,GeneID:929780;gbkey=CDS;product=hypothetical protein;protein_id=NP_344583.1;transl_table=11
    Streptococcus mitis B6                  genus       35      Translation elongation factor TU    ID=cds854;Name=YP_003446034.1;Parent=gene929;Dbxref=Genbank:YP_003446034.1,GeneID:8797981;gbkey=CDS;product=Translation elongation factor TU;protein_id=YP_003446034.1;transl_table=11
    Streptococcus pneumoniae 670-6B         genus       35      translation elongation factor Tu    ID=cds793;Name=YP_003878976.1;Parent=gene805;Dbxref=Genbank:YP_003878976.1,GeneID:9729530;gbkey=CDS;product=translation elongation factor Tu;protein_id=YP_003878976.1;transl_table=11
    Streptococcus pneumoniae 70585          genus       35      elongation factor Tu                ID=cds1437;Name=YP_002740757.1;Parent=gene1488;Note=EF-Tu%3B promotes GTP-dependent binding of aminoacyl-tRNA to the A-site of ribosomes during protein biosynthesis%3B when the tRNA anticodon matches the mRNA codon%2C GTP hydrolysis results%3B the inactive EF-Tu-GDP leaves the ribosome and release of GDP is promoted by elongation factor Ts%3B many prokaryotes have two copies of the gene encoding EF-Tu;Dbxref=Genbank:YP_002740757.1,GeneID:7684439;gbkey=CDS;product=elongation factor Tu;protein_id=YP_002740757.1;transl_table=11
    Streptococcus pneumoniae A026           genus       35      elongation factor Tu                ID=cds1209;Name=YP_008730608.1;Parent=gene1336;Note=EF-Tu%3B promotes GTP-dependent binding of aminoacyl-tRNA to the A-site of ribosomes during protein biosynthesis%3B when the tRNA anticodon matches the mRNA codon%2C GTP hydrolysis results%3B the inactive EF-Tu-GDP leaves the ribosome and release of GDP is promoted by elongation factor Ts%3B many prokaryotes have two copies of the gene encoding EF-Tu;Dbxref=Genbank:YP_008730608.1,GeneID:17439784;gbkey=CDS;gene=tuf;product=elongation factor Tu;protein_id=YP_008730608.1;transl_table=11
    Streptococcus pneumoniae ATCC 700669    genus       35      elongation factor Tu                ID=cds1258;Name=YP_002511360.1;Parent=gene1394;Note=EF-Tu%3B promotes GTP-dependent binding of aminoacyl-tRNA to the A-site of ribosomes during protein biosynthesis%3B when the tRNA anticodon matches the mRNA codon%2C GTP hydrolysis results%3B the inactive EF-Tu-GDP leaves the ribosome and release of GDP is promoted by elongation factor Ts%3B many prokaryotes have two copies of the gene encoding EF-Tu;Dbxref=Genbank:YP_002511360.1,GeneID:7329336;gbkey=CDS;product=elongation factor Tu;protein_id=YP_002511360.1;transl_table=11

    [...]

The columns, left to right, are:

 * ``Species``: The species/strain name of the genome to which the matches
   occurred.
 * ``Disc. level``: The taxonomic rank at which the peptide that matched this
   genome was discriminative.
 * ``Count``: The number of hits/matches accrued by this annotated region
   (total across the entire sample). Note that these counts are not quantitative, because 
   each discriminative peptide can produce multiple matches to several
   different annotated regions in several different reference genome sequences.
 * ``Product``: Annotation for the protein product of this annotated region.
 * ``Features``: The raw feature string for this annotated region.

As we can see in this output, there are several hits to regions in different
*S. pneumoniae* genomes, to regions annotated with *hypothetical protein* and
*elongation factor Tu*. This sheet provides information that can be useful to
identify potential biomarkers or help understand what regions of the studied
bacteria are discriminative, and at what taxonomic level.

.. note::

    The *Hits to annotated regions* output can also be printed to STDOUT
    using the --print-annotations command line options. However, this is
    **not** recommended as it is very verbose.


Discriminative peptides
-----------------------
The names of all discriminative peptides along with the rank at which they were
discriminative can be output to a text file using the
``--write-discriminative-peptides`` command line option.

Antibiotic resistance
*********************
Detected expressed antibiotic resistance peptides are summarized in a
text file. The output can look like this::

    ----------------------------------------------------------------------
    Results for 124 discriminative peptides in tcup_tutorial_sample.fasta.ar.blast8
    Disc.  Hits        %  Family
    64     84     0.5114  blaCTX-M
    45     45     0.2739  blaTEM
    10     10     0.0609  aac(3)-II
    5      5      0.0304  mph(A)

Here, we see that the sample likely contains four different families of
antibiotic resistance mechanisms: blaCTX-M, blaTEM, aac(3)-II, and mph(A). 
The columns, from left to right, are:

 * ``Disc.``: The number of peptides that were discriminative to this family,
   i.e. peptides that did not match any other family.
 * ``Hits``: The number of aligned regions to this family that the
   discriminative peptides produced. This number of sometimes higher than the number
   of discriminative peptides, as each discriminative peptide can sometimes match to
   several variants within a family, or sometimes even several positions in the same
   protein.
 * ``%``: The proportion of peptides in the sample that were discriminative to 
   this family.
 * ``Family``: The antibotic resistance gene family matched by discriminative
   peptides.



Change log
==========
All notable changes to this project will be documented in this file.
The project loosely follows `semantic versioning <www.semver.org>`_.
The following subheaders are used in the change log: 

* `Added` for new features.
* `Changed` for changes in existing functionality.
* `Deprecated` for once-stable features removed in upcoming releases.
* `Removed` for deprecated features removed in this release.
* `Fixed` for any bug fixes.
* `Security` to invite users to upgrade in case of vulnerabilities.

The topmost subheader is always `[wip]`, and lists what is currently being
worked on. Below the work in progress, a subheader for each notable release 
lists the main changes introduced by that release.

[wip]
*****

Added
-----
*

Fixed
-----
*


[Second release] 2017-04-21 version 1.1.0
****************************************
Added
-----
* Peptide sequence is now included in `--write-discriminative-peptides` output
* `--reassign-taxids` now allows to modify taxid assignments for refseqs in the
  taxref DB.
* `--peptide-fasta` to include peptide sequences in `taxonomic composition`
  output file.
* `--add-annotations` to add annotations to annotation DB. 
* New helper script: `get_peptide_sequences.py` to print discriminative
  peptides and their taxonomic assignments.

Fixed
-----
* Removed `Terrabacteria group` with `no rank` from output listing.
* No longer try to insert taxids into taxref DB when pre existing sample DB is used.


[First release] 2016-06-13 version 1.0.0
****************************************
Added
-----
* Taxonomic composition estimation with xslx output.
* Antibiotic resistance detection.

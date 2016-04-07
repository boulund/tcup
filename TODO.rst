TODO
====

* Clean up and improve documentation.
* Modify SQLite3 code around lines 355-370 in taxonomic_composition.py to
  something UPSERT-like, to remove the need for including ALL taxids in the 
  'cumulative' table. This is one of the slowest parts of the entire code now 
  because it isn't indexed.
* Flesh out the CONTRIBUTING document with more details.
* Pick a license.

# TODO

 * Fail early for missing database files (e.g. proteodb.sql or annotationdb.sql).
 * Finish installation instructions in README.
 * Modify SQLite3 code around lines 355-370 in taxonomic_composition.py to
   something UPSERT-like, to remove the need for including ALL taxids in the 
   'cumulative' table. This is one of the slowest parts of the entire code now 
   because it isn't indexed.

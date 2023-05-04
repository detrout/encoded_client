Changelog
=========

0.3.0
-----

Start adding IGVF portal compatibility

  - New server names!
  - New endpoints! (There are more types of files now)

This lead to a new process_endpoint_file to replace the old
process_file (now Deprecated) to specify where to upload files to.

The submission main() function now also has an --end-point argument to
specify the endpoint being posted to.

There's new updates to support computing metadata needed to submit
splitseq subpools to ENCODE.

Fixed a bug where the list of files on an experiment object kept
getting added to each type the object was accessed.

Also use default_analysis to figure out which set of analysis files
encoded.EncodeExperiment should be returning instead of just
assuming it was the first set.

0.2.0
-----

- Fix RDFa parsing when used with more recent rdflib.
- Support posting documents with mime types other than text/pdf
- Verify that aws command is available before posting metadata.
- Downgrade not a submitter to a warning, instead of an error.
  (Because of validating using records submitted by others)
- Remove xlrd dependency in favor of openpyxl to read current xlsx files.
- When submitting documents quote any path seperators since the DCC requires
  document names to be in a flat hierarchy.
- Add a server paremeter to Document class for validation.
- Implement qc metric parsers for star solo and scRNAseq quality
  control metrics

0.1.1
-----

Fixes to setup.cfg to improve pypi landing page.

0.1.0
-----

There was a great deal of work on improving organizing and collecting
metadata for submitting ENCODE scRNA-seq experiments.

encoded_client now supports the DCC_API_KEY and DCC_SECRET_KEY
environment variables used by encode_utils in addition to reading from
the .netrc file.

ENCODED.get_response has now been made more general, previously it was
what was setting the accept: application/json header but that ment I
could use it for downloading files. now get_json sets the accept
header.

There were many changes to make skipping less common modules if their
dependencies weren't installed

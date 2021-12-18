Changelog
=========

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

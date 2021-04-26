Introduction
============

This is the the set of classes I have been using to handle submitting
data to the ENCODE DCC. I would only expect the most current ENCODE 4
submission tools to be useful, though some of the code for the earlier
versions of the project are still present.

I believe there is another encoded client, but I believe this version
is the only to attempt to use the DCC's provided jsonschema to validate
objects before attempting to post them.

The rough workflow I've been following is to generate a metadata
workbook that contains the information for the main 5 data types
relevant for submission into separate sheets named after the
submission object types.

I populate that with information available in our LIMS and then add in
the additional information the DCC has started to request of us that
hasn't made it into our LIMS.

I then have a jupyter notebook that reads each of the sheets, and in
dry-run mode attempts to validate the sheet. Identifiers are cached
between sheets so later sheets can use IDs defined in earlier sheets.

This ends up making submissions easier as instead of having to copy
the returned IDs from one object type to the next object type, I can
define aliases ahead of time and use those.

The call ENCODED.post_sheet is expecting a pandas DataFrame. Pandas
can natively read Excel sheets, and I contributed the ability to read
OpenDocument .ods sheets which was merged in 0.25.

I also ported a project I found called "gcat" that can parse Google
sheets to Python 3 and have occasionally made fixes to my
`fork <https://github.com/detrout/gcat>`_.

I occasionally wonder if I should fix the current warnings, add some
tests, and release it to pypi.


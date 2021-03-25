Usage
=====

The client is actually useful for both query and submission, and
querying is safer I'll start illustrating use there.

Retrieving an experiment
------------------------

A basic use case to find information about an experiment.

.. testcode::

   from encoded_client.encoded import ENCODED

   server = ENCODED("www.encodeproject.org")
   experiment = server.get_json("ENCSR000AEG")

   print(experiment["description"])

.. testoutput::

   PolyA RNA-Seq from oligo-dT primed Total RNA on the GM12878 cell line


Searching
---------

Sometimes, however, you may not have a specific accession list of access
ids to search for and you need to look for objects.

.. testcode::

   from encoded_client.encoded import ENCODED

   server = ENCODED("www.encodeproject.org")
   query = server.search_jsonld(searchTerm="C2C12", limit=5)

   print("Results:", len(query["@graph"]))
   for row in query["@graph"]:
       print(row["@id"], row.get("description", "n/a"))

.. testoutput::

    Results: 5
    /annotations/ENCSR777WDG/ candidate Cis-Regulatory Elements in C3H C2C12 for GRCh38
    /annotations/ENCSR991QYA/ Functional validation of enhancers in C2C12 myoblast cells
    /annotations/ENCSR953JBP/ Functional validation of enhancers in C2C12 myocyte cells
    /annotations/ENCSR309JHZ/ candidate Cis-Regulatory Elements in C3H myocyte originated from C2C12 for GRCh38
    /documents/a5f5c35a-cdda-4a45-9742-22e69ff50c9c/ C2C12 cell culture, differentiation treatment, and cross-linking protocol


Using the "searchTerm" argument generates a URL exactly like entering
the query into the search box on the https://www.encodeproject.org
website. The limit argument is optional and the default is to return
25 records. I shortened it for this example.

Additionally when doing programatic access, you are likely to want to
use limit="all" instead.


Submitting
----------

Here's a block from one of my submission notebooks, preparing to
submit a Biosamples page.


.. code-block:: python

    import pandas
    from encoded_client.encoded import ENCODED, DCCValidator

    server = ENCODED("test.encodedcc.org")
    server.load_netrc()
    validator = DCCValidator(server)
    
    biosample = pandas.read_excel(
                   spreadsheet_name,
                   sheet_name='Biosamples',
                   header=0,
                   engine='odf')

    created = server.post_sheet('/biosamples/',
                                verbose=True,
                                dry_run=True,
                                validator=validator)

    if created:
        biosample.to_excel("biosamples-created.xlsx", index=False)

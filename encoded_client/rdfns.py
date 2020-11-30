"""Namespace definitions

All in one place to make import rdfns.* work safely
"""
from rdflib import Namespace
from rdflib.namespace import DC, RDF, RDFS, OWL, XMLNS, XSD

# standard ontology namespaces
VS = Namespace("http://www.w3.org/2003/06/sw-vocab-status/ns#")
WOT = Namespace("http://xmlns.com/wot/0.1/")

# internal ontologies
submissionOntology = Namespace(
    "http://jumpgate.caltech.edu/wiki/UcscSubmissionOntology#"
)
dafTermOntology = Namespace("http://jumpgate.caltech.edu/wiki/UcscDaf#")
libraryOntology = Namespace("http://jumpgate.caltech.edu/wiki/LibraryOntology#")
inventoryOntology = Namespace("http://jumpgate.caltech.edu/wiki/InventoryOntology#")
submissionLog = Namespace("http://jumpgate.caltech.edu/wiki/SubmissionsLog/")
geoSoftNS = Namespace("http://www.ncbi.nlm.nih.gov/geo/info/soft2.html#")
encode3NS = Namespace("http://jumpgate.caltech.edu/wiki/Encode3#")

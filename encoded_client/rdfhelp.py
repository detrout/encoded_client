"""Helper features for working with librdf
"""
from __future__ import print_function, absolute_import

import collections
from glob import glob
import six
from six.moves import urllib
import logging
import os
import sys
import types
from pkg_resources import resource_listdir, resource_string

from rdflib import ConjunctiveGraph, Graph, Literal, BNode, URIRef, Namespace
from rdflib.namespace import ClosedNamespace, DC, RDF, RDFS, OWL, XMLNS, XSD
from .rdfns import (
    encode3NS,
    dafTermOntology,
    geoSoftNS,
    libraryOntology,
    submissionOntology,
    VS,
    WOT,
)

import lxml.html
import lxml.html.clean

logger = logging.getLogger(__name__)


SCHEMAS_URL = "http://jumpgate.caltech.edu/phony/schemas"
INFERENCE_URL = "http://jumpgate.caltech.edu/phony/inference"

ISOFORMAT_MS = "%Y-%m-%dT%H:%M:%S.%f"
ISOFORMAT_SHORT = "%Y-%m-%dT%H:%M:%S"


def sparql_query(model, query_filename, output_format="text"):
    """Execute sparql query from file"""
    logger.info("Opening: %s" % (query_filename,))
    query_body = open(query_filename, "r").read()
    query = RDF.SPARQLQuery(query_body)
    results = query.execute(model)
    if output_format == "html":
        html_query_results(results)
    else:
        display_query_results(results)


def display_query_results(results):
    """A very simple display of sparql query results showing name value pairs"""
    for row in results:
        for k, v in row.items()[::-1]:
            print("{0}: {1}".format(k, v))
        print()


def html_query_results(result_stream):
    from django.conf import settings
    from django.template import Context, loader

    # I did this because I couldn't figure out how to
    # get simplify_rdf into the django template as a filter
    class Simplified(object):
        def __init__(self, value):
            self.simple = simplify_rdf(value)
            if value.is_resource():
                self.url = value
            else:
                self.url = None

    template = loader.get_template("rdf_report.html")
    results = []
    for row in result_stream:
        new_row = collections.OrderedDict()
        row_urls = []
        for k, v in row.items():
            new_row[k] = Simplified(v)
        results.append(new_row)
    context = Context(
        {
            "results": results,
        }
    )
    print(template.render(context))


def get_node_type(node):
    """Return just the base name of a XSD datatype:
    e.g. http://www.w3.org/2001/XMLSchema#integer -> integer
    """
    # chop off xml schema declaration
    value_type = node.datatype
    if value_type is None:
        return "string"
    else:
        return value_type.replace(str(XSD), "").lower()


def simplify_rdf(value):
    """Return a short name for a RDF object
    e.g. The last part of a URI or an untyped string.
    """
    if isinstance(value, Literal):
        name = value.value
    elif isinstance(value, BNode):
        name = "<BLANK>"
    elif isinstance(value, URIRef):
        name = split_uri(str(value))
    else:
        name = value
    return str(name)


def simplify_uri(uri):
    """Split off the end of a uri

    >>> simplify_uri('http://asdf.org/foo/bar')
    'bar'
    >>> simplify_uri('http://asdf.org/foo/bar#bleem')
    'bleem'
    >>> simplify_uri('http://asdf.org/foo/bar/')
    'bar'
    >>> simplify_uri('http://asdf.org/foo/bar?was=foo')
    'was=foo'
    """
    if isinstance(uri, Literal) and uri.datatype not in (XSD.anyURI,):
        raise ValueError("Literal terms must be of URI type")

    uri = str(uri)

    parsed = urllib.parse.urlparse(uri)
    if len(parsed.query) > 0:
        return parsed.query
    elif len(parsed.fragment) > 0:
        return parsed.fragment
    elif len(parsed.path) > 0:
        for element in reversed(parsed.path.split("/")):
            if len(element) > 0:
                return element
    raise ValueError("Unable to simplify %s" % (uri,))


def strip_namespace(namespace, term):
    """Remove the namespace portion of a term

    returns None if they aren't in common
    """
    if not isinstance(namespace, (URIRef, Namespace, ClosedNamespace)):
        raise ValueError("Requires a URIRef namespace")

    if isinstance(term, Literal) and term.datatype not in (XSD.anyURI,):
        raise ValueError("Term literals must be a URI type")
    elif not isinstance(term, URIRef):
        raise ValueError("Term must be a URI type")

    term_s = str(term)
    if not term_s.startswith(str(namespace)):
        return None
    return term_s.replace(str(namespace), "")


def load_into_model(model, parser_name, path, ns=None):
    if isinstance(ns, six.string_types):
        ns = URIRef(ns)

    if isinstance(path, URIRef):
        path = str(path)

    url_parts = list(urllib.parse.urlparse(path))
    if len(url_parts[0]) == 0 or url_parts[0] == "file":
        url_parts[0] = "file"
        url_parts[2] = os.path.abspath(url_parts[2])
    if parser_name is None or parser_name == "guess":
        parser_name = guess_parser_by_extension(path)
    url = urllib.parse.urlunparse(url_parts)
    logger.info("Opening {0} with parser {1}".format(url, parser_name))

    kwargs = {}
    if parser_name == "rdfa":
        kwargs["media_type"] = "text/html"

    model.parse(url, format=parser_name, publicID=ns, **kwargs)


def load_string_into_model(model, parser_name, data, ns=None):
    ns = fixup_namespace(ns)
    logger.debug(
        "load_string_into_model parser={0}, len={1}".format(parser_name, len(data))
    )

    model.parse(data=data, format=parser_name, publicID=ns)
    add_imports(model, ns)


def fixup_namespace(ns):
    if ns is None:
        ns = URIRef("http://localhost/")
    elif isinstance(ns, six.string_types):
        ns = URIRef(ns)
    elif not (isinstance(ns, URIRef)):
        errmsg = "Namespace should be string or uri not {0}"
        raise ValueError(errmsg.format(str(type(ns))))
    return ns


def add_imports(model, ns):
    for s, p, o in model.triples((None, OWL.imports, None)):
        if p == OWL.imports:
            model.remove((s, p, o))
            load_into_model(model, None, o, ns)


def add_default_schemas(model, schema_path=None):
    """Add default schemas to a model
    Looks for turtle files in either htsworkflow/util/schemas
    or in the list of directories provided in schema_path
    """
    schemas = resource_listdir(__name__, "schemas")
    for s in schemas:
        schema = resource_string(__name__, "schemas/" + s)
        if six.PY3:
            # files must be encoded utf-8
            schema = schema.decode("utf-8")
        namespace = "file://localhost/htsworkflow/schemas/" + s
        add_schema(model, schema, namespace)

    if schema_path:
        if type(schema_path) in types.StringTypes:
            schema_path = [schema_path]

        for path in schema_path:
            for pathname in glob(os.path.join(path, "*.turtle")):
                url = "file://" + os.path.splitext(pathname)[0]
                stream = open(pathname, "rt")
                add_schema(model, stream, url)
                stream.close()


def add_schema(model, schema, url):
    """Add a schema to a model.

    Main difference from 'load_into_model' is it tags it with
    a RDFlib context so I can remove them later.
    """
    if not isinstance(model, ConjunctiveGraph):
        raise ValueError("Schemas requires a graph that supports quads")

    context = URIRef(SCHEMAS_URL)
    tmpmodel = Graph()
    tmpmodel.parse(data=schema, format="turtle", publicID=url)
    for s, p, o in tmpmodel:
        model.add((s, p, o, context))


def remove_schemas(model):
    """Remove statements labeled with our schema context"""
    context = URIRef(SCHEMAS_URL)
    for quad in model.triples((None, None, None, context)):
        model.remove(quad)
        # model.remove_context(context)


def sanitize_literal(node):
    """Clean up a literal string"""
    if not isinstance(node, Literal):
        raise ValueError("sanitize_literal only works on Literals")

    s = node.value
    if len(s) > 0:
        element = lxml.html.fromstring(s)
        cleaner = lxml.html.clean.Cleaner(page_structure=False)
        element = cleaner.clean_html(element)
        if six.PY3:
            text = lxml.html.tostring(element, encoding=str)
        else:
            text = lxml.html.tostring(element)
        p_len = 3
        slash_p_len = 4

        value = text[p_len:-slash_p_len]
    else:
        value = ""
    args = {}
    if node.datatype is not None:
        args["datatype"] = node.datatype
    if node.language is not None:
        args["lang"] = node.language
    return Literal(value, **args)


def guess_parser(content_type, pathname):
    if content_type in ("application/rdf+xml",):
        return "rdfxml"
    elif content_type in ("application/x-turtle",):
        return "turtle"
    elif content_type in ("text/html",):
        return "rdfa"
    elif content_type is None or content_type in ("text/plain",):
        return guess_parser_by_extension(pathname)


def guess_parser_by_extension(pathname):
    _, ext = os.path.splitext(pathname)
    if ext in (".xml", ".rdf"):
        return "rdfxml"
    elif ext in (".html",):
        return "rdfa"
    elif ext in (".turtle",):
        return "turtle"
    return "guess"


def add_default_namespaces(model):
    """Return a serializer with our standard prefixes loaded"""
    model.bind("rdf", RDF)
    model.bind("rdfs", RDFS)
    model.bind("owl", OWL)
    model.bind("dc", DC)
    model.bind("xml", XMLNS)
    model.bind("xsd", XSD)
    model.bind("vs", VS)
    model.bind("wot", WOT)

    # should these be here, kind of specific to an application
    model.bind("htswlib", libraryOntology)
    model.bind("ucscSubmission", submissionOntology)
    model.bind("ucscDaf", dafTermOntology)
    model.bind("geoSoft", geoSoftNS)
    model.bind("encode3", encode3NS)
    return model


def get_turtle_header():
    """Return a turtle header with our typical namespaces"""
    empty = ConjunctiveGraph()
    add_default_namespaces(empty)
    turtle_header = []
    for term, urlterm in empty.namespaces():
        turtle_header.append("@prefix {}: <{}>.".format(term, urlterm))
    return "\n".join(turtle_header)


def dump_model(model, destination=None):
    if destination is None:
        destination = sys.stdout
    add_default_namespaces(model)
    model.serialize(destination, format="turtle")

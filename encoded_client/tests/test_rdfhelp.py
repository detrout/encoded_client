from __future__ import print_function

import os
from unittest import TestCase
from datetime import datetime
import six
import pytest

pytest.importorskip("rdflib")

from rdflib import ConjunctiveGraph, Literal, Namespace, URIRef, Graph

from ..rdfhelp import (
    add_default_schemas,
    DC,
    guess_parser,
    guess_parser_by_extension,
    load_string_into_model,
    OWL,
    remove_schemas,
    RDF,
    RDFS,
    strip_namespace,
    simplify_uri,
    sanitize_literal,
    XSD,
)


class TestRDFHelp(TestCase):
    def test_typed_node_boolean(self):
        node = Literal(True)
        self.assertTrue(node.value)
        self.assertEqual(str(node.datatype), "http://www.w3.org/2001/XMLSchema#boolean")

    def test_typed_node_string_node_attributes(self):
        node = Literal("hello")
        self.assertEqual(node.value, "hello")
        self.assertTrue(node.datatype is None)

    def test_typed_real_like(self):
        num = 3.14
        node = Literal(num)
        self.assertEqual(node.toPython(), num)
        self.assertEqual(type(node.toPython()), float)

    def test_typed_integer(self):
        num = 3
        node = Literal(num)
        self.assertEqual(node.toPython(), num)
        self.assertEqual(type(node.toPython()), type(num))

    def test_typed_node_string(self):
        s = "Argh matey"
        node = Literal(s)
        self.assertEqual(node.toPython(), s)
        self.assertTrue(isinstance(node.toPython(), six.text_type))

    def test_unicode_node_roundtrip(self):
        literal = u"\u5927"
        roundtrip = Literal(literal).toPython()
        self.assertTrue(isinstance(roundtrip, six.text_type))

    def test_datetime_no_microsecond(self):
        dateTimeType = XSD.dateTime
        short_isostamp = "2011-12-20T11:44:25"
        short_node = Literal(short_isostamp, datatype=dateTimeType)
        short_datetime = datetime(2011, 12, 20, 11, 44, 25)

        self.assertEqual(short_node.toPython(), short_datetime)
        self.assertEqual(Literal(short_datetime), short_node)
        self.assertEqual(Literal(short_datetime).toPython(), short_datetime)

    def test_datetime_with_microsecond(self):
        dateTimeType = XSD.dateTime
        long_isostamp = "2011-12-20T11:44:25.081776"
        long_node = Literal(long_isostamp, datatype=dateTimeType)
        long_datetime = datetime(2011, 12, 20, 11, 44, 25, 81776)

        self.assertEqual(long_node.toPython(), long_datetime)
        self.assertEqual(Literal(long_datetime), long_node)
        self.assertEqual(Literal(long_datetime).toPython(), long_datetime)

    def test_strip_namespace_uri(self):
        nsOrg = Namespace("example.org/example#")
        nsCom = Namespace("example.com/example#")

        term = "foo"
        node = nsOrg[term]
        self.assertEqual(strip_namespace(nsOrg, node), term)
        self.assertEqual(strip_namespace(nsCom, node), None)

    def test_strip_namespace_exceptions(self):
        nsOrg = Namespace("example.org/example#")
        nsCom = Namespace("example.com/example#")

        node = Literal("bad")
        self.assertRaises(ValueError, strip_namespace, nsOrg, node)
        self.assertRaises(ValueError, strip_namespace, nsOrg, nsOrg)
        self.assertRaises(ValueError, strip_namespace, nsOrg, str(node))

    def test_simplify_uri(self):
        DATA = [
            ("http://asdf.org/foo/bar", "bar"),
            ("http://asdf.org/foo/bar#bleem", "bleem"),
            ("http://asdf.org/foo/bar/", "bar"),
            ("http://asdf.org/foo/bar?was=foo", "was=foo"),
        ]

        for uri, expected in DATA:
            self.assertEqual(simplify_uri(uri), expected)

        for uri, expected in DATA:
            n = URIRef(uri)
            self.assertEqual(simplify_uri(n), expected)

        for uri, expected in DATA:
            n = Literal(URIRef(uri), datatype=XSD.anyURI)
            self.assertEqual(simplify_uri(n), expected)

        # decoding literals is questionable
        n = Literal("http://foo/bar")
        self.assertRaises(ValueError, simplify_uri, n)

    def test_owl_import(self):
        path, name = os.path.split(__file__)
        # loc = 'file://'+os.path.abspath(path)+'/'
        loc = os.path.abspath(path) + "/"
        model = Graph()
        fragment = """
@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix owl: <http://www.w3.org/2002/07/owl#> .

_:a owl:imports "{loc}extra.turtle" .
""".format(
            loc=loc
        )
        load_string_into_model(model, "turtle", fragment, loc)
        tc = URIRef("http://jumpgate.caltech.edu/wiki/TestCase")
        result = list(model.triples((tc, RDFS.label, None)))
        self.assertEqual(len(result), 1)
        self.assertEqual(str(result[0][2]), "TestCase")

    def test_sanitize_literal_text(self):
        self.assertRaises(ValueError, sanitize_literal, "hi")
        hello_text = "hello"
        hello_none = Literal(hello_text)
        self.assertEqual(str(sanitize_literal(hello_none)), hello_text)
        hello_str = Literal(hello_text, datatype=XSD["string"])
        hello_clean = sanitize_literal(hello_str)
        self.assertEqual(hello_clean.value, hello_text)

    def test_sanitize_literal_empty_string(self):
        value = ""
        value_node = Literal(value)
        self.assertEqual(str(sanitize_literal(value_node)), value)

    def test_sanitize_literal_html(self):
        hello = "hello <a onload='javascript:alert(\"foo\");' href='http://google.com'>google.com</a>, whats up?"
        hello_clean = 'hello <a href="http://google.com">google.com</a>, whats up?'
        hello_node = Literal(hello, datatype=XSD["string"])
        hello_sanitized = sanitize_literal(hello_node)
        self.assertEqual(hello_sanitized.value, hello_clean)

        hostile = "hi <b>there</b><script type='text/javascript>alert('boo');</script><a href='javascript:alert('poke')>evil</a> scammer"
        hostile_node = Literal(hostile)
        hostile_sanitized = sanitize_literal(hostile_node)
        # so it drops the stuff after the javascript link.
        # I suppose it could be worse
        hostile_result = """hi <b>there</b>"""
        self.assertEqual(str(hostile_sanitized), hostile_result)

    def test_guess_parser_from_file(self):
        DATA = [
            ("/a/b/c.rdf", "rdfxml"),
            ("/a/b/c.xml", "rdfxml"),
            ("/a/b/c.html", "rdfa"),
            ("/a/b/c.turtle", "turtle"),
            ("http://foo.bar/bleem.turtle", "turtle"),
        ]
        for path, parser in DATA:
            self.assertEqual(guess_parser_by_extension(path), parser)
            self.assertEqual(guess_parser(None, path), parser)

        DATA = [
            ("application/rdf+xml", "http://a.org/b/c", "rdfxml"),
            ("application/x-turtle", "http://a.org/b/c", "turtle"),
            ("text/html", "http://a.org/b/c", "rdfa"),
            ("text/html", "http://a.org/b/c.html", "rdfa"),
            ("text/plain", "http://a.org/b/c.turtle", "turtle"),
            ("text/plain", "http://a.org/b/c", "guess"),
        ]
        for contenttype, url, parser in DATA:
            self.assertEqual(guess_parser(contenttype, url), parser)


class TestRDFSchemas(TestCase):
    def test_rdf_schema(self):
        """Does it basically work?"""
        model = ConjunctiveGraph()
        self.assertEqual(len(model), 0)
        add_default_schemas(model)
        self.assertTrue(len(model) > 0)
        remove_schemas(model)
        self.assertEqual(len(model), 0)

    def test_included_schemas(self):
        model = ConjunctiveGraph()
        add_default_schemas(model)

        # rdf test
        s = [RDF, DC["title"], None]
        title = model.objects(RDF, DC["title"])
        self.assertTrue(title is not None)

        s = [RDF["Property"], RDF["type"], RDFS["Class"]]
        self.assertIn(s, model)

        # rdfs test
        s = [RDFS["Class"], RDF["type"], RDFS["Class"]]
        self.assertIn(s, model)

        s = [OWL["inverseOf"], RDF["type"], RDF["Property"]]
        self.assertIn(s, model)


def suite():
    from unittest import TestSuite, defaultTestLoader

    suite = TestSuite()
    suite.addTests(defaultTestLoader.loadTestsFromTestCase(TestRDFHelp))
    suite.addTests(defaultTestLoader.loadTestsFromTestCase(TestRDFSchemas))
    return suite


if __name__ == "__main__":
    from unittest import main

    main(defaultTest="suite")

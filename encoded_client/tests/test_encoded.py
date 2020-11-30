from __future__ import absolute_import, print_function

import datetime
import json
import jsonschema
import logging
import os
from unittest import TestCase

from ..encoded import (
    ENCODED,
    ENCODED_CONTEXT,
    ENCODED_NAMESPACES,
    DCCValidator,
    typed_column_parser,
)


class TestEncoded(TestCase):
    def setUp(self):
        self.encode = ENCODED("www.encodeproject.org")
        self.encode._user = {
            "@context": "/terms/",
            "@id": "/users/bc5b62f7-ce28-4a1e-b6b3-81c9c5a86d7a/",
            "@type": ["User", "Item"],
            "first_name": "Diane",
            "groups": [],
            "job_title": "Submitter",
            "lab": {
                "@id": "/labs/barbara-wold/",
                "@type": ["Lab", "Item"],
                "country": "USA",
                "institute_label": "Caltech",
                "institute_name": "California Institute of Technology",
                "pi": "/users/0598c868-0b4a-4c5b-9112-8f85c5de5374/",
                "schema_version": "4",
                "title": "Barbara Wold, Caltech",
                "uuid": "72d5666a-a361-4f7b-ab66-a88e11280937",
            },
            "last_name": "Trout",
            "schema_version": "5",
            "submits_for": [
                "/labs/barbara-wold/",
                "/labs/richard-myers/",
                "/labs/ali-mortazavi/",
            ],
            "uuid": "bc5b62f7-ce28-4a1e-b6b3-81c9c5a86d7a",
        }

        logging.disable(logging.WARNING)
        self.validator = DCCValidator(self.encode)
        for schema, filename in [
            ("library", "library.json"),
            ("biosample", "biosample.json"),
        ]:
            schema_file = os.path.join(os.path.dirname(__file__), filename)
            with open(schema_file, "rt") as instream:
                self.validator._schemas[schema] = json.load(instream)

    def tearDown(self):
        logging.disable(logging.NOTSET)

    def test_prepare_url(self):
        tests = [
            ("/experiments", "https://www.encodeproject.org/experiments"),
            (
                "/experiments/ENCLB045ZZZ",
                "https://www.encodeproject.org/experiments/ENCLB045ZZZ",
            ),
            (
                "https://www.encodeproject.org/experiments/ENCLB045ZZZ",
                "https://www.encodeproject.org/experiments/ENCLB045ZZZ",
            ),
            ("barbara-wold:11111", "https://www.encodeproject.org/barbara-wold:11111"),
        ]
        for url, result in tests:
            self.assertEqual(self.encode.prepare_url(url), result)

    def test_validate_library(self):
        """Test validation of a Library object"""
        obj = {
            u"aliases": [],
            u"alternate_accessions": [],
            u"award": u"/awards/U54HG006998/",
            u"biosample": u"/biosamples/ENCBS089RNA/",
            u"date_created": u"2014-01-14T19:44:51.061770+00:00",
            u"documents": [],
            u"extraction_method": u"Ambion mirVana",
            u"fragmentation_method": u"chemical (Nextera tagmentation)",
            u"lab": u"/labs/barbara-wold/",
            u"library_size_selection_method": u"SPRI beads",
            u"lysis_method": u"Ambion mirVana",
            u"nucleic_acid_term_name": u"polyadenylated mRNA",
            u"size_range": u">200",
            u"status": u"released",
            u"strand_specificity": False,
            u"treatments": [],
        }
        self.validator.validate(obj, "library")

        # test requestMethod
        obj["schema_version"] = u"2"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, obj, "library"
        )
        del obj["schema_version"]

        # test calculatedProperty
        obj["nucleic_acid_term_name"] = u"SO:0000871"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, obj, "library"
        )
        del obj["nucleic_acid_term_name"]

        # test permssionValidator
        obj["uuid"] = u"42c46028-708f-4347-a3df-2c82dfb021c4"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, obj, "library"
        )
        del obj["uuid"]

    def test_validate_biosample(self):
        bio = {
            "aliases": ["barbara-wold:c1_e12.5_mouse_limb_donor"],
            "award": "U54HG006998",
            "biosample_term_id": "UBERON:0002101",
            "biosample_term_name": "limb",
            "biosample_type": "tissue",
            "date_obtained": "2017-02-01",
            "description": "C57Bl6 wild-type embryonic mouse",
            "donor": "/mouse-donors/ENCDO956IXV/",
            "lab": "/labs/barbara-wold",
            "model_organism_age": "12.5",
            "model_organism_age_units": "day",
            "mouse_life_stage": "embryonic",
            "organism": "3413218c-3d86-498b-a0a2-9a406638e786",
            "source": "/sources/gems-caltech/",
            "starting_amount": 1,
            "starting_amount_units": "items",
        }

        # test aliases
        bio["aliases"] = ["barbara-wold:c1-donor"]
        self.validator.validate(bio, "biosample")
        del bio["aliases"]

        # test part_of
        bio["part_of"] = "barbara-wold:c1-donor"
        self.validator.validate(bio, "biosample")
        del bio["part_of"]

        # tests linkTo
        self.validator.validate(bio, "biosample")
        bio["organism"] = "/organisms/mouse/"

        bio["lab"] = "/labs/alkes-price/"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, bio, "biosample"
        )
        bio["lab"] = "/labs/barbara-wold"

        bio["organism"] = "7745b647-ff15-4ff3-9ced-b897d4e2983c"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, bio, "biosample"
        )
        bio["organism"] = "/organisms/human"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, bio, "biosample"
        )
        bio["organism"] = "/organisms/mouse/"

    def test_aliases(self):
        """Test that objects being validated can access previous ones by alias

        Some properties can require that an object exist, since we validate
        before submitting, we need to be able cache object by its alias and use
        that to retrieve it.
        """
        donor = {
            "aliases": ["barbara-wold:c1_e12.5_mouse_limb_donor"],
            "award": "U54HG006998",
            "biosample_term_id": "UBERON:0002101",
            "biosample_term_name": "limb",
            "biosample_type": "tissue",
            "date_obtained": "2017-02-01",
            "description": "C57Bl6 wild-type embryonic mouse",
            "donor": "/mouse-donors/ENCDO956IXV/",
            "lab": "/labs/barbara-wold",
            "model_organism_age": "12.5",
            "model_organism_age_units": "day",
            "mouse_life_stage": "embryonic",
            "organism": "3413218c-3d86-498b-a0a2-9a406638e786",
            "source": "/sources/gems-caltech/",
            "starting_amount": 1,
            "starting_amount_units": "items",
        }
        part = donor.copy()
        part["aliases"] = ["barbara-wold:A7"]
        part["part_of"] = "barbara-wold:c1_e12.5_mouse_limb_donor"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, part, "biosample"
        )

        self.validator.validate(donor, "biosample")
        self.validator.validate(part, "biosample")

    def test_create_context(self):
        linked_id = {"@type": "@id"}
        library = {"@id": "/libraries/1234", "@type": ["Library", "Item"]}

        url = self.encode.prepare_url(library["@id"])
        context = self.encode.create_jsonld_context(library, url)
        self.assertEqual(
            context["@vocab"], "https://www.encodeproject.org/profiles/library.json#"
        )
        self.assertEqual(context["award"], linked_id)
        self._verify_context(context, "Library")
        # namespaces not added yet.
        self.assertRaises(AssertionError, self._verify_namespaces, context)
        self.encode.add_jsonld_namespaces(context)
        self._verify_namespaces(context)

    def test_add_context(self):
        """Checking to make sure nested @base and @vocab urls are set correctly"""
        obj = {
            "nucleic_acid_term_name": "RNA",
            "accession": "ENCLB044ZZZ",
            "@id": "/libraries/ENCLB044ZZZ/",
            "schema_version": "1",
            "@type": ["Library", "Item"],
            "lysis_method": "Ambion mirVana",
            "nucleic_acid_term_id": "SO:0000356",
            "biosample": {
                "biosample_term_name": "GM12878",
                "description": "B-lymphocyte, lymphoblastoid, International HapMap Project - CEPH/Utah - European Caucasion, Epstein-Barr Virus",
                "accession": "ENCBS090RNA",
                "date_created": "2013-10-29T21:15:29.144260+00:00",
                "@id": "/biosamples/ENCBS090RNA/",
                "aliases": ["brenton-graveley:GM12878-2", "thomas-gingeras:191WC"],
                "organism": "/organisms/human/",
                "@type": ["Biosample", "Item"],
            },
        }

        bio_base = self.encode.prepare_url(obj["biosample"]["@id"])

        url = self.encode.prepare_url("/libraries/ENCLB044ZZZ/?format=json&embed=False")
        obj_type = self.encode.get_object_type(obj)
        schema_url = self.encode.get_schema_url(obj_type)
        self.encode.add_jsonld_context(obj, url)

        self.assertEqual(obj["biosample"]["@context"]["@base"], bio_base)
        self.assertEqual(obj["@context"]["@vocab"], schema_url)
        self._verify_context(obj["@context"], "Library")
        self._verify_namespaces(obj["@context"])
        self._verify_context(obj["biosample"]["@context"], "Biosample")
        self.assertEqual(
            obj["@context"]["rdf"], "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        )
        self.assertEqual(obj["@context"]["OBO"], "http://purl.obolibrary.org/obo/")

    def test_convert_search_to_jsonld(self):
        example = {
            "count": {"biosamples": 2},
            "portal_title": "ENCODE",
            "title": "Search",
            "notification": "Success",
            "filters": [],
            "@id": "/search/?searchTerm=wold",
            "@type": ["search"],
            "facets": [],
            "@graph": [
                {
                    u"@id": u"/biosamples/ENCBS125ENC/",
                    u"@type": [u"Biosample", u"Item"],
                    u"accession": u"ENCBS125ENC",
                    u"award.rfa": u"ENCODE2-Mouse",
                    u"biosample_term_name": u"myocyte",
                    u"biosample_type": u"in vitro differentiated cells",
                    u"characterizations.length": [],
                    u"constructs.length": [],
                    u"lab.title": u"Barbara Wold, Caltech",
                    u"life_stage": u"unknown",
                    u"organism.name": u"mouse",
                    u"source.title": u"Barbara Wold",
                    u"status": u"CURRENT",
                    u"treatments.length": [],
                },
                {
                    u"@id": u"/biosamples/ENCBS126ENC/",
                    u"@type": [u"Biosample", u"Item"],
                    u"accession": u"ENCBS126ENC",
                    u"award.rfa": u"ENCODE2-Mouse",
                    u"biosample_term_name": u"myocyte",
                    u"biosample_type": u"in vitro differentiated cells",
                    u"characterizations.length": [],
                    u"constructs.length": [],
                    u"lab.title": u"Barbara Wold, Caltech",
                    u"life_stage": u"unknown",
                    u"organism.name": u"mouse",
                    u"source.title": u"Barbara Wold",
                    u"status": u"CURRENT",
                    u"treatments.length": [],
                },
            ],
        }

        result = self.encode.convert_search_to_jsonld(example)
        for obj in result["@graph"]:
            self.assertNotIn("award.rfa", obj)

    def _verify_context(self, context, obj_type):
        for context_key in [None, obj_type]:
            for k in ENCODED_CONTEXT[context_key]:
                self.assertIn(k, context)
                self.assertEqual(ENCODED_CONTEXT[context_key][k], context[k])

    def _verify_namespaces(self, context):
        for k in ENCODED_NAMESPACES:
            self.assertIn(k, context)
            self.assertEqual(ENCODED_NAMESPACES[k], context[k])


class TestTypedColumnParser(TestCase):
    def test_skip(self):
        self.assertEqual(typed_column_parser("foo:skip", "value"), (None, None))

    def test_default(self):
        name, value = typed_column_parser("foo", "3")
        self.assertEqual(name, "foo")
        self.assertEqual(value, "3")
        self.assertIsInstance(value, str)

    def test_array(self):
        name, value = typed_column_parser("foo:array", "a, b, c")
        self.assertEqual(name, "foo")
        self.assertEqual(len(value), 3)
        self.assertEqual(value, ["a", "b", "c"])

        name, value = typed_column_parser("foo:array", "a,b,c")
        self.assertEqual(name, "foo")
        self.assertEqual(len(value), 3)
        self.assertEqual(value, ["a", "b", "c"])

        name, value = typed_column_parser("foo:array", "a, b,\nc")
        self.assertEqual(name, "foo")
        self.assertEqual(len(value), 3)
        self.assertEqual(value, ["a", "b", "c"])

    def test_boolean(self):
        for b in ["True", "False", "1", "0"]:
            name, value = typed_column_parser("foo:boolean", b)
            self.assertEqual(name, "foo")
            self.assertEqual(value, bool(b))

    def test_integer(self):
        name, value = typed_column_parser("foo:integer", "3")
        self.assertEqual(name, "foo")
        self.assertEqual(value, 3)

        self.assertRaises(ValueError, typed_column_parser, "foo:integer", "a")
        self.assertRaises(ValueError, typed_column_parser, "foo:integer", "3.14")

    def test_number(self):
        name, value = typed_column_parser("foo:number", "3")
        self.assertEqual(name, "foo")
        self.assertEqual(value, 3)
        self.assertIsInstance(value, int)

        name, value = typed_column_parser("foo:number", "3.14")
        self.assertEqual(name, "foo")
        self.assertEqual(value, 3.14)
        self.assertIsInstance(value, float)

        self.assertRaises(ValueError, typed_column_parser, "foo:number", "a")

    def test_timestamp(self):
        name, value = typed_column_parser("foo:date", "1999-1-1")
        self.assertEqual(name, "foo")
        self.assertEqual(value, "1999-1-1")

        name, value = typed_column_parser("foo:date", datetime.date(1999, 1, 1))
        self.assertEqual(name, "foo")
        self.assertEqual(value, "1999-01-01")

    def test_json(self):
        name, value = typed_column_parser("foo:json", '{"a": 3}')
        self.assertEqual(value, {"a": 3})

        name, value = typed_column_parser("foo:json", '[{"a": 3}]')
        self.assertEqual(value, [{"a": 3}])


if __name__ == "__main__":
    from unittest import main

    main()

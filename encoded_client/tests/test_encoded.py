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
    DuplicateAliasError,
    get_object_type,
    typed_column_parser,
    LOGGER,
)


def load_dcc_schemas(validator):
    for schema, filename in [
        ("file", "file.json"),
        ("library", "library.json"),
        ("biosample", "biosample.json"),
        ("star_solo_quality_metric", "star_solo_quality_metric.json"),
    ]:
        schema_file = os.path.join(os.path.dirname(__file__), filename)
        with open(schema_file, "rt") as instream:
            validator._schemas[schema] = json.load(instream)


def get_test_dcc_user():
    return {
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


class TestEncoded(TestCase):
    def setUp(self):
        self.encode = ENCODED("www.encodeproject.org")
        self.encode._user = get_test_dcc_user()

        logging.disable(logging.WARNING)
        self.validator = DCCValidator(self.encode)
        load_dcc_schemas(self.validator)

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
            "aliases": [],
            "alternate_accessions": [],
            "award": "/awards/U54HG006998/",
            "biosample": "/biosamples/ENCBS089RNA/",
            "date_created": "2014-01-14T19:44:51.061770+00:00",
            "documents": [],
            "extraction_method": "Ambion mirVana",
            "fragmentation_methods": ["chemical (Nextera tagmentation)"],
            "lab": "/labs/barbara-wold/",
            "library_size_selection_method": "SPRI beads",
            "lysis_method": "Clontech UltraLow for Illumina sequencing",
            "nucleic_acid_term_name": "polyadenylated mRNA",
            "size_range": ">200",
            "status": "released",
            "strand_specificity": "unstranded",
            "treatments": [],
        }
        self.validator.validate(obj, "library")

        # test requestMethod
        obj["schema_version"] = "2"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, obj, "library"
        )
        del obj["schema_version"]

        # test calculatedProperty
        obj["nucleic_acid_term_name"] = "SO:0000871"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, obj, "library"
        )
        del obj["nucleic_acid_term_name"]

        # test permssionValidator
        obj["uuid"] = "42c46028-708f-4347-a3df-2c82dfb021c4"
        self.assertRaises(
            jsonschema.ValidationError, self.validator.validate, obj, "library"
        )
        del obj["uuid"]

    def test_validate_biosample(self):
        bio = {
            "aliases": ["barbara-wold:c1_e12.5_mouse_limb_donor"],
            "award": "U54HG006998",
            "biosample_ontology": "/biosample-types/tissue_UBERON_0002101/",
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
        with self.assertLogs(LOGGER) as cm:
            self.validator.validate(bio, "biosample")

        self.assertEqual(
            cm.output,
            [
                "ERROR:encoded_client.encoded:/labs/alkes-price/ is not a lab in user "
                "Diane Trout can submits_for"
            ],
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
            "biosample_ontology": "/biosample-types/tissue_UBERON_0002101/",
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

    def test_alias_collision(self):
        # Make sure an alias can only be defined once
        donor = {
            "aliases": ["barbara-wold:c1_e12.5_mouse_limb_donor"],
            "award": "U54HG006998",
            "biosample_ontology": "/biosample-types/tissue_UBERON_0002101/",
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
        part["description"] = "C57Bl6 wild-type embryonic mouse collision"

        self.validator.validate(donor, "biosample")
        self.assertRaises(
            DuplicateAliasError, self.validator.validate, part, "biosample"
        )

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
        obj_type = get_object_type(obj)
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
                    "@id": "/biosamples/ENCBS125ENC/",
                    "@type": ["Biosample", "Item"],
                    "accession": "ENCBS125ENC",
                    "award.rfa": "ENCODE2-Mouse",
                    "biosample_term_name": "myocyte",
                    "biosample_type": "in vitro differentiated cells",
                    "characterizations.length": [],
                    "constructs.length": [],
                    "lab.title": "Barbara Wold, Caltech",
                    "life_stage": "unknown",
                    "organism.name": "mouse",
                    "source.title": "Barbara Wold",
                    "status": "CURRENT",
                    "treatments.length": [],
                },
                {
                    "@id": "/biosamples/ENCBS126ENC/",
                    "@type": ["Biosample", "Item"],
                    "accession": "ENCBS126ENC",
                    "award.rfa": "ENCODE2-Mouse",
                    "biosample_term_name": "myocyte",
                    "biosample_type": "in vitro differentiated cells",
                    "characterizations.length": [],
                    "constructs.length": [],
                    "lab.title": "Barbara Wold, Caltech",
                    "life_stage": "unknown",
                    "organism.name": "mouse",
                    "source.title": "Barbara Wold",
                    "status": "CURRENT",
                    "treatments.length": [],
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
        self.assertEqual(value, "1999-01-01")

        name, value = typed_column_parser("foo:date", "1/1/1999")
        self.assertEqual(name, "foo")
        self.assertEqual(value, "1999-01-01")

        name, value = typed_column_parser("foo:date", datetime.date(1999, 1, 1))
        self.assertEqual(name, "foo")
        self.assertEqual(value, "1999-01-01")

        self.assertRaises(ValueError, typed_column_parser, "foo:date", "asdf")

    def test_json(self):
        name, value = typed_column_parser("foo:json", '{"a": 3}')
        self.assertEqual(value, {"a": 3})

        name, value = typed_column_parser("foo:json", '[{"a": 3}]')
        self.assertEqual(value, [{"a": 3}])


# should have a guard for "requires network access".
class TestEncodeExperiment(TestCase):
    def test_experiment(self):
        server = ENCODED("www.encodeproject.org")

        e = server.get_experiment("/experiments/ENCSR000AEH/")
        self.assertIsNot(e._json, None)
        self.assertEqual(e.assay_term_name, "polyA plus RNA-seq")
        self.assertEqual(
            e.description, "RNA Evaluation Gm12878 Long Poly-A+ RNA-seq from Wold"
        )

        replicates = list(e.replicates)
        self.assertEqual(len(replicates), 2)

        replicate_files = {}
        file_formats_expected = set(["tsv", "fastq", "bam", "bigWig"])
        output_types_expected = set(
            [
                "reads",
                "signal of all reads",
                "signal of unique reads",
                "transcriptome alignments",
                "transcript quantifications",
                "gene quantifications",
                "alignments",
                "signal",
            ]
        )
        library_aliases_expected = set(["barbara-wold:13713", "barbara-wold:13714"])
        library_aliases_seen = set()
        for i, r in enumerate(replicates):
            for alias in r["library"]["aliases"]:
                library_aliases_seen.add(alias)
            files = list(r.files)
            replicate_files[i] = [f["@id"] for f in files]
            self.assertGreaterEqual(len(files), 12)
            file_formats_seen = set()
            output_types_seen = set()
            qc_seen = []
            qc_types_seen = set()
            for f in files:
                file_formats_seen.add(f["file_format"])
                output_types_seen.add(f["output_type"])
                for qc in f.quality_metrics:
                    qc_seen.append(qc)
                    qc_types_seen.add(get_object_type(qc_seen[-1]))

            self.assertGreaterEqual(len(qc_seen), 4)
            expected = {
                "GeneQuantificationQualityMetric",
                "SamtoolsFlagstatsQualityMetric",
                "MadQualityMetric",
                "GeneTypeQuantificationQualityMetric",
                "StarQualityMetric",
            }
            for e in expected:
                self.assertIn(e, qc_types_seen)

        self.assertEqual(library_aliases_seen, library_aliases_expected)
        self.assertEqual(file_formats_expected, file_formats_seen)
        self.assertEqual(output_types_expected, output_types_seen)

        self.assertEqual(
            len(set(replicate_files[0]).intersection(replicate_files[1])), 0
        )


if __name__ == "__main__":
    from unittest import main

    main()

import datetime
import hashlib
import logging
import os
from unittest import TestCase
from unittest.mock import patch

from .test_encoded import get_test_dcc_user, load_dcc_schemas
from ..encoded import (
    ENCODED,
    DCCValidator,
)
from ..metadata import (
    compute_alignment_derived_from,
    compute_count_matrix_derived_from,
    compute_alignment_alias,
    generate_star_solo_processed_metadata,
    generate_star_solo_processed_sheet,
)


def mock_stat():
    class MockStatResult:
        def stat(self):
            return os.stat_result(
                (33188, 300, 50, 1, 1000, 1000, 1000, 1635793563, 1635793563, 1635793563)
            )

    return MockStatResult()


def md5sum_string(s):
    return hashlib.md5(s.encode("utf-8")).hexdigest()


def get_sample_config1():
    return {
        "read1": ["ENCFF150FBF", "ENCFF385IAW"],
        "read2": ["ENCFF351VBS", "ENCFF503CCI"],
        "include_intron": True,
        "stranded": "Forward",
        "umi_length": 12,
        "experiment_accession": "ENCSR724KET",
        "description": "snRNA on human adrenal gland",
        "library_accession": "ENCLB002DZK",
        "genome_accession": "ENCFF795ZFF",
        "genome_assembly": "GRCh38",
        "genome_annotation": "V29",
        "lab": "/labs/barbara-wold/",
        "alias_prefix": "barbara-wold",
        "award": "UM1HG009443",
        "alignment_step_run": "barbara-wold:starsolo-alignment-step-run",
        "quantification_step_run": "barbara-wold:starsolo-quantification-step-run",
    }


def get_processed_files1():
    return {
        "alignments": "Aligned.sortedByCoord.out.bam",
        "sparse gene count matrix of unique reads": "GeneFull_Ex50pAS_Unique_filtered.tar.gz",
        "sparse gene count matrix of all reads": "GeneFull_Ex50pAS_EM_filtered.tar.gz",
        "unfiltered sparse gene count matrix of unique reads": "GeneFull_Ex50pAS_Unique_raw.tar.gz",
        "unfiltered sparse gene count matrix of all reads": "GeneFull_Ex50pAS_EM_raw.tar.gz",
        "unfiltered sparse splice junction count matrix of unique reads": "SJ_Unique_raw.tar.gz",
    }

class test_metadata(TestCase):
    def setUp(self):
        sctest = "encd-6248-78c147870-qing.demo.encodedcc.org"
        production = "www.encodeproject.org"
        self.encode = ENCODED(sctest)
        self.encode._user = get_test_dcc_user()

        logging.disable(logging.WARNING)
        self.validator = DCCValidator(self.encode)
        load_dcc_schemas(self.validator)

    def tearDown(self):
        logging.disable(logging.NOTSET)

    def test_compute_alignment_derived_from(self):
        read1 = ["ENCFF150FBF", "ENCFF385IAW"]
        read2 = ["ENCFF351VBS", "ENCFF503CCI"]
        # This is actually the bulk index not the single cell index
        star_index = [" ENCFF795ZFF"]

        expected_single = [
            "/files/{}/".format(star_index[0]),
            "/files/{}/".format(read1[0]),
            "/files/{}/".format(read1[1]),
        ]

        expected_paired = [
            "/files/{}/".format(star_index[0]),
            "/files/{}/".format(read1[0]),
            "/files/{}/".format(read2[0]),
            "/files/{}/".format(read1[1]),
            "/files/{}/".format(read2[1]),
        ]

        # single end
        values = compute_alignment_derived_from(star_index, read1)
        self.assertEqual(len(values), 3)
        self.assertEqual(values, expected_single)

        # paired end
        values = compute_alignment_derived_from(star_index, read1, read2)
        self.assertEqual(len(values), 5)
        self.assertEqual(values, expected_paired)

    def test_compute_count_matrix_derived_from(self):
        datestamp = "2021-11-22"
        alias = compute_alignment_alias("barbara-wold", "ENCLB002DZK", datestamp)
        self.assertEqual(alias, "barbara-wold:ENCLB002DZK_alignment_2021-11-22")

        derived_from = compute_count_matrix_derived_from(alias)
        self.assertEqual(len(derived_from), 1)
        self.assertEqual(derived_from, [alias])

    @patch("encoded_client.metadata.Path")
    @patch("encoded_client.metadata.make_md5sum", wraps=md5sum_string)
    def test_generate_star_solo_processed_metadata(self, mock_md5sum, mock_path):
        mock_path.return_value = mock_stat()

        config = get_sample_config1()
        alignment_step_run = config["alignment_step_run"]
        quantification_step_run = config["quantification_step_run"]
        processed_files = get_processed_files1()

        datestamp = datetime.datetime.now().strftime("%Y-%m-%d")
        alias = compute_alignment_alias("barbara-wold", "ENCLB002DZK", datestamp)
        expected_bam_derived_from = [
            "/files/ENCFF795ZFF/",
            "/files/ENCFF150FBF/",
            "/files/ENCFF351VBS/",
            "/files/ENCFF385IAW/",
            "/files/ENCFF503CCI/",
        ]

        metadata = generate_star_solo_processed_metadata(
            config, processed_files
        )
        self.assertEqual(len(metadata), len(processed_files))

        for record, term in zip(metadata, processed_files):
            self.assertIsNone(record["accession"])
            self.assertIsNone(record["uuid"])
            self.assertEqual(record["dataset"], "ENCSR724KET")
            self.assertEqual(record["award"], "UM1HG009443")
            self.assertEqual(record["lab"], "/labs/barbara-wold/")
            self.assertEqual(
                record["submitted_file_name"], processed_files[term]
            )
            self.assertEqual(record["output_type"], term)
            if record["file_format"] == "bam":
                self.assertEqual(
                    record["derived_from"], expected_bam_derived_from
                )
                self.assertEqual(record["step_run"], alignment_step_run)
            else:
                self.assertEqual(record["derived_from"], [alias])
                self.assertEqual(record["step_run"], quantification_step_run)

            hidden = record.copy()
            del hidden["accession"]
            self.validator.validate(hidden, "file")

    @patch("encoded_client.metadata.Path")
    @patch("encoded_client.metadata.make_md5sum", wraps=md5sum_string)
    def test_generate_star_solo_processed_sheet(self, mock_md5sum, mock_path):
        mock_path.return_value = mock_stat()

        config = get_sample_config1()
        processed_files = get_processed_files1()

        metadata = generate_star_solo_processed_sheet(
            config, processed_files
        )

        to_create = self.encode.post_sheet('/files/', metadata, dry_run=True, validator=self.validator)
        self.assertEqual(len(to_create), 6)
        for row in to_create:
            self.assertEqual(row["accession"], "would create")

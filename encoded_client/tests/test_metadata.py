import datetime
import hashlib
import os
from pprint import pprint
from unittest import TestCase
from unittest.mock import patch

from ..metadata import (
    compute_alignment_derived_from,
    compute_count_matrix_derived_from,
    compute_alignment_alias,
    generate_star_solo_processed_metadata,
)


def mock_stat():
    return os.stat_result(
        (33188, 300, 50, 1, 1000, 1000, 1000, 1635793563, 1635793563, 1635793563)
    )


def md5sum_string(s):
    return hashlib.md5(s.encode("utf-8")).hexdigest()


class test_metadata(TestCase):
    def test_compute_alignment_derived_from(self):
        read1 = ["ENCFF150FBF", "ENCFF385IAW"]
        read2 = ["ENCFF351VBS", "ENCFF503CCI"]
        star_index = ["TST123456"]

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
        mock_path.stat.return_value = mock_stat()

        config = {
            "read1": ["ENCFF150FBF", "ENCFF385IAW"],
            "read2": ["ENCFF351VBS", "ENCFF503CCI"],
            "include_intron": True,
            "stranded": "Forward",
            "umi_length": 12,
            "experiment_accession": "ENCSR724KET",
            "description": "snRNA on human adrenal gland",
            "library_accession": "ENCLB002DZK",
            "genome_accession": "ENCFF123456",
            "genome_assembly": "GRCh38",
            "genome_annotation": "V29",
            "lab": "/labs/my-lab/",
            "alias_prefix": "my-lab",
            "award": "award",
        }

        processed_files = {
            "alignments": "Aligned.sortedByCoord.out.bam",
            "sparse gene count matrix of unique reads": "GeneFull_Ex50pAS_Unique_filtered.tar.gz",
            "sparse gene count matrix of all reads": "GeneFull_Ex50pAS_EM_filtered.tar.gz",
            "unfiltered sparse gene count matrix of unique reads": "GeneFull_Ex50pAS_Unique_raw.tar.gz",
            "unfiltered sparse gene count matrix of all reads": "GeneFull_Ex50pAS_EM_raw.tar.gz",
            "unfiltered sparse splice junction count matrix of unique reads": "SJ_Unique_raw.tar.gz",
        }
        datestamp = datetime.datetime.now().strftime("%Y-%m-%d")
        alias = compute_alignment_alias("my-lab", "ENCLB002DZK", datestamp)
        expected_bam_derived_from = ",".join(
            [
                "/files/ENCFF123456/",
                "/files/ENCFF150FBF/",
                "/files/ENCFF351VBS/",
                "/files/ENCFF385IAW/",
                "/files/ENCFF503CCI/",
            ]
        )

        metadata = generate_star_solo_processed_metadata(
            config, processed_files
        )
        self.assertEqual(len(metadata), len(processed_files))

        pprint(metadata)
        for record, term in zip(metadata, processed_files):
            self.assertIsNone(record["accession"])
            self.assertIsNone(record["uuid"])
            self.assertEqual(record["dataset"], "ENCSR724KET")
            self.assertEqual(record["award"], "award")
            self.assertEqual(record["lab"], "/labs/my-lab/")
            self.assertEqual(
                record["submitted_file_name"], processed_files[term]
            )
            self.assertEqual(record["output_type"], term)
            if record["file_format"] == "bam":
                self.assertEqual(
                    record["derived_from:array"], expected_bam_derived_from
                )
            else:
                self.assertEqual(record["derived_from:array"], alias)

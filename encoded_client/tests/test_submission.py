import pandas
import tempfile
from unittest import TestCase
import pytest

from ..encoded import ENCODED
from ..submission import process_files, process_endpoint_files, main


class TestSubmission(TestCase):
    def setUp(self):
        self.example_files = pandas.DataFrame(
            {
                "uuid": [None],
                "accession": [None],
                "dataset": ["ENCSR411MUF"],
                "submitted_file_name": [
                    __file__
                ],
                "md5sum": ["d826c25c5f5bd3263f40df6556c87180"],
                "flowcell_details:json": [
                    '[{"machine": "http://jumpgate.caltech.edu/sequencer/8", "flowcell": "HGLTKBCX3", "lane": "1", "barcode": "TTACCGAC"}]'
                ],
                "read_length:integer": ["100"],
                "file_format": ["fastq"],
                "output_type": ["reads"],
                "run_type": ["single-ended"],
                "platform": ["/platforms/OBI%3A0002002/"],
                "replicate": ["barbara-wold:22588_b1_t1"],
                "lab": ["barbara-wold"],
                "award": ["UM1HG009443"],
            }
        )

    def test_process_files(self):
        server = ENCODED("test.encodedcc.org")
        self.assertIsNone(self.example_files.iloc[0]["uuid"])
        self.assertIsNone(self.example_files.iloc[0]["accession"])

        self.assertRaises(
            DeprecationWarning, process_files, server, self.example_files, dry_run=True)
        process_endpoint_files(server, "/files/", self.example_files, dry_run=True)

        self.assertEqual(self.example_files.iloc[0]["accession"], "would create")

    def test_sumission_main(self):
        pytest.importorskip("openpyxl")
        with tempfile.NamedTemporaryFile(suffix="_encoded.xlsx") as out:
            self.example_files.to_excel(out, sheet_name="File", index=False)

            main(["-s", "test.encodedcc.org", "-f", out.name, "--dry-run"])

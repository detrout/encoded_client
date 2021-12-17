from unittest import TestCase
import pandas
import tempfile
import pytest

from ..sheet import (
    open_book,
    save_book,
)


class TestSheet(TestCase):
    def test_spreadsheet_wrappers(self):
        pytest.importorskip("openpyxl")
        sheet = pandas.DataFrame({
            "uuid": ["fee130b0-aed1-4889-8703-4a520b70fc79", None, None],
            "accession": ["ENCSR000AEC", None, None],
            "submitted_filename": [
                "ENCSR000AEC/f1.fastq.gz",
                "ENCSR000AED/f1.fastq.gz",
                "ENCSR000AEF/f1.fastq.gz",
            ]
        })
        sheet_name = "File"

        with tempfile.NamedTemporaryFile(suffix="_sheet.xlsx") as out:
            sheet.to_excel(out.name, sheet_name=sheet_name, index=False)

            book1 = open_book(out.name)

        with tempfile.NamedTemporaryFile(suffix="_book.xlsx") as out:
            save_book(out.name, book1)

            book2 = open_book(out.name)

        sheet2 = book2.parse(sheet_name)
        for (row1_index, row1), (row2_index, row2) in zip(sheet.iterrows(), sheet2.iterrows()):
            for val1, val2 in zip(row1, row2):
                if pandas.isnull(val1) and pandas.isnull(val2):
                    pass
                else:
                    self.assertEqual(val1, val2)

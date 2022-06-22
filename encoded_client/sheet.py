"""Utilities to help streamline spreadsheet IO
"""

from pathlib import Path
import pandas


def open_book(spreadsheet_file, force_extension=None):
    """Open saved spreadsheet
    """
    spreadsheet_file = Path(spreadsheet_file)

    extension = spreadsheet_file.suffix if force_extension is None else force_extension

    engine = {
        ".ods": "odf",
        ".xls": "xlrd",
        ".xlsx": "openpyxl",
    }[extension]

    book = pandas.ExcelFile(spreadsheet_file, engine=engine)
    return book


def save_book(filename, book, replacements={}):
    with pandas.ExcelWriter(filename) as new_submission:
        for sheet_name in book.sheet_names:
            sheet = replacements.get(
                sheet_name,
                book.parse(sheet_name, header=0)
            )
            sheet.to_excel(
                new_submission,
                sheet_name=sheet_name,
                index=False
            )

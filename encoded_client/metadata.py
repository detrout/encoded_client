import datetime
import logging
import pandas
from pathlib import Path
from urllib.parse import urlparse
from encoded_client.hashfile import make_md5sum

logger = logging.getLogger(__name__)


def format_accession(container, accession):
    # ENCODE aliases have a : in them, if we have an alias, return it,
    # otherwise format the accession for our target container.
    if ":" in accession:
        return accession
    else:
        return "/{}/{}/".format(container, accession)


def compute_alignment_derived_from(index_accessions, read1, read2=None):
    if isinstance(index_accessions, str):
        index_accessions = [index_accessions]

    derived_from = []
    for accession in index_accessions:
        if ":" in accession or accession.startswith("/files"):
            derived_from.append(accession)
        else:
            derived_from.append(format_accession("files", accession))

    if read2 is not None:
        if len(read2) != len(read1):
            raise ValueError("Lengths of read lists must match")
        for pair in zip(read1, read2):
            for accession in pair:
                derived_from.append(format_accession("files", accession))
    else:
        for accession in read1:
            derived_from.append(format_accession("files", accession))

    return derived_from


def compute_dcc_file_accession_from_url(url):
    """The ENCODE DCC accession ID is part of the download url.

    This extracts the ID from the url and reformats it to match
    the ENCODE object id.
    """
    parsed = urlparse(url)
    # if we already have an ID the http/https will already be removed
    if parsed.scheme is None or len(parsed.scheme) == 0:
        return url
    path = parsed.path
    name = path[path.rindex("/") + 1 :]
    for suffix in [".txt.gz", ".tar.gz"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    dcc_id = format_accession("files", name)
    return dcc_id


def compute_count_matrix_derived_from(config, alignment):
    if "inclusion_list_url" in config:
        config.setdefault(
            "inclusion_accession",
            compute_dcc_file_accession_from_url(config["inclusion_list_url"]),
        )

    derived_from = [config["inclusion_accession"]]
    derived_from.append(format_accession("files", alignment))
    return derived_from


def compute_alignment_alias(alias_prefix, library, datestamp):
    return "{prefix}:{library}_alignment_{datestamp}".format(
        prefix=alias_prefix, library=library, datestamp=datestamp
    )


def generate_star_solo_processed_metadata(config, records):
    datestamp = datetime.datetime.now().strftime("%Y-%m-%d")
    output_type_to_file_type = {
        "alignments": "bam",
        "sparse gene count matrix of unique reads": "tar",
        "sparse gene count matrix of all reads": "tar",
        "unfiltered sparse gene count matrix of unique reads": "tar",
        "unfiltered sparse gene count matrix of all reads": "tar",
        "unfiltered sparse splice junction count matrix of unique reads": "tar",
    }

    alignment_alias = compute_alignment_alias(
        config["alias_prefix"], config["library_accession"], datestamp
    )
    rows = []
    for output_type in records:
        filename = records[output_type]
        file_type = output_type_to_file_type[output_type]

        if file_type == "bam":
            derived_from = compute_alignment_derived_from(
                config["genome_accession"], config["read1"], config["read2"]
            )
        elif file_type == "tar":
            derived_from = compute_count_matrix_derived_from(config, alignment_alias)

        obj = {
            "uuid": None,
            "accession": None,
            "dataset": config["experiment_accession"],
            "file_format": file_type,
            "output_type": output_type,
            "assembly": config["assembly"],
            "genome_annotation": config["genome_annotation"],
            "derived_from": derived_from,
            "md5sum": make_md5sum(filename),
            "file_size": Path(filename).stat().st_size,
            "submitted_file_name": str(filename),
            "award": config["award"],
            "lab": config["lab"],
        }
        if file_type == "bam":
            obj["step_run"] = config["alignment_step_run"]
            obj["aliases"] = [alignment_alias]
        elif file_type == "tar":
            obj["step_run"] = config["quantification_step_run"]
        else:
            logger.error("Unknown file type {}".format(file_type))
        rows.append(obj)

    return rows


def generate_star_solo_processed_sheet(config, records):
    rows = generate_star_solo_processed_metadata(config, records)

    sheet = pandas.DataFrame(rows)

    for array in ["aliases", "derived_from"]:
        if array in sheet.columns:
            sheet[array] = sheet[array].apply(to_array_sheet)

    sheet = sheet.rename(
        {
            "aliases": "aliases:array",
            "derived_from": "derived_from:array",
            "file_size": "file_size:integer",
        },
        axis=1,
    )

    return sheet


def to_array_sheet(value):
    nulls = pandas.isnull(value)
    if nulls is True:
        return None

    return ",".join(value)


####
# Metadata for encode subpools
#
def compute_subpool_matrix_alias(config, library_id, output_type, datestamp):
    """Compute an alias for mex archive files generated by subpools"""
    output_type_to_alias = {
        "sparse gene count matrix of unique reads": "filtered_unique",
        "sparse gene count matrix of all reads": "filtered_em",
        "unfiltered sparse gene count matrix of unique reads": "raw_unique",
        "unfiltered sparse gene count matrix of all reads": "raw_em",
        "unfiltered sparse splice junction count matrix of unique reads": "raw_splice",
    }

    return "{prefix}:{library_id}_{term}_{datestamp}".format(
        prefix=config["alias_prefix"],
        library_id=library_id,
        term=output_type_to_alias[output_type],
        datestamp=datestamp,
    )


def generate_star_solo_subpool_metadata(config, records, library_id=None):
    """Generate star solo metadata for a subpool run"""
    datestamp = datetime.datetime.now().strftime("%Y-%m-%d")
    output_type_to_file_type = {
        "alignments": "bam",
        "sparse gene count matrix of unique reads": "tar",
        "sparse gene count matrix of all reads": "tar",
        "unfiltered sparse gene count matrix of unique reads": "tar",
        "unfiltered sparse gene count matrix of all reads": "tar",
        "unfiltered sparse splice junction count matrix of unique reads": "tar",
    }

    alignment_alias = compute_alignment_alias(
        config["alias_prefix"], library_id, datestamp
    )
    rows = []
    for output_type in records:
        filename = records[output_type]
        file_type = output_type_to_file_type[output_type]

        if "library" in config:
            library = config["library"][library_id]
        else:
            library = config

        if file_type == "bam":
            derived_from = compute_alignment_derived_from(
                config["genome_accession"], library["read1"], library["read2"]
            )
            matrix_alias = None
        elif file_type == "tar":
            derived_from = compute_count_matrix_derived_from(config, alignment_alias)
            matrix_alias = compute_subpool_matrix_alias(
                config, library_id, output_type, datestamp
            )

        obj = {
            "uuid": None,
            "accession": None,
            "dataset": config["experiment_accession"],
            "file_format": file_type,
            "output_type": output_type,
            "assembly": config["assembly"],
            "genome_annotation": config["genome_annotation"],
            "derived_from": derived_from,
            "md5sum": make_md5sum(filename),
            "file_size": Path(filename).stat().st_size,
            "submitted_file_name": str(filename),
            "award": config["award"],
            "lab": config["lab"],
        }
        if file_type == "bam":
            obj["step_run"] = config["alignment_step_run"]
            obj["aliases"] = [alignment_alias]
        elif file_type == "tar":
            obj["step_run"] = config["quantification_step_run"]
            obj["aliases"] = [matrix_alias]
        else:
            logger.error("Unknown file type {}".format(file_type))
        rows.append(obj)

    return rows


def generate_star_solo_subpool_sheet(config, records, library_id=None):
    rows = generate_star_solo_subpool_metadata(config, records, library_id)

    sheet = pandas.DataFrame(rows)

    for array in ["aliases", "derived_from"]:
        if array in sheet.columns:
            sheet[array] = sheet[array].apply(to_array_sheet)

    sheet = sheet.rename(
        {
            "aliases": "aliases:array",
            "derived_from": "derived_from:array",
            "file_size": "file_size:integer",
        },
        axis=1,
    )

    return sheet

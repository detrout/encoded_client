import datetime
import logging
import pandas
from pathlib import Path
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


def compute_inclusion_id(inclusion_list_url):
    name = Path(inclusion_list_url).name
    if name.endswith(".txt.gz"):
        name = name[:-len(".txt.gz")]
    dcc_id = format_accession("files", name)
    return dcc_id


def compute_count_matrix_derived_from(inclusion_url, alignment):
    derived_from = [compute_inclusion_id(inclusion_url)]
    derived_from.append(format_accession("files", alignment))
    return derived_from


def compute_alignment_alias(alias_prefix, library, datestamp):
    return "{prefix}:{library}_alignment_{datestamp}".format(
        prefix=alias_prefix,
        library=library,
        datestamp=datestamp)


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

    alignment_alias = compute_alignment_alias(config['alias_prefix'], config['library_accession'], datestamp)
    rows = []
    for output_type in records:
        filename = records[output_type]
        file_type = output_type_to_file_type[output_type]

        if file_type == "bam":
            derived_from = compute_alignment_derived_from(config["genome_accession"], config['read1'], config['read2'])
        elif file_type == "tar":
            inclusion_id = compute_inclusion_id(config["inclusion_list_url"])
            derived_from = compute_count_matrix_derived_from(inclusion_id, alignment_alias)

        obj = {
            'uuid': None,
            'accession': None,
            'dataset': config["experiment_accession"],
            'file_format': file_type,
            'output_type': output_type,
            'assembly': config["assembly"],
            'genome_annotation': config["genome_annotation"],
            'derived_from': derived_from,
            'md5sum': make_md5sum(filename),
            'file_size': Path(filename).stat().st_size,
            'submitted_file_name': str(filename),
            'award': config["award"],
            'lab': config["lab"],
        }
        if file_type == "bam":
            obj['step_run'] = config["alignment_step_run"]
            obj['aliases'] = [alignment_alias]
        elif file_type == "tar":
            obj['step_run'] = config["quantification_step_run"]
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

    sheet = sheet.rename({
        "aliases": "aliases:array",
        "derived_from": "derived_from:array",
        "file_size": "file_size:integer",
    }, axis=1)

    return sheet


def to_array_sheet(value):
    nulls = pandas.isnull(value)
    if nulls is True:
        return None

    return ",".join(value)

import datetime
from pathlib import Path
from encoded_client.hashfile import make_md5sum


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


def compute_count_matrix_derived_from(alignment):
    return [format_accession("files", alignment)]


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
            derived_from = compute_count_matrix_derived_from(alignment_alias)

        obj = {
            'uuid': None,
            'accession': None,
            'dataset': config["experiment_accession"],
            'file_format': file_type,
            'output_type': output_type,
            'assembly': config["genome_assembly"],
            'genome_annotation': config["genome_annotation"],
            'step_run': config["step_run"],
            'derived_from:array': ",".join(derived_from),
            'md5sum': make_md5sum(filename),
            'file_size:integer': Path(filename).stat().st_size,
            'submitted_file_name': str(filename),
            #'quality_metrics:json':
            'award': config["award"],
            'lab': config["lab"],
        }
        if file_type == "bam":
            obj['aliases'] = [alignment_alias]
        rows.append(obj)

    return rows
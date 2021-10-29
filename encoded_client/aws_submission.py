"""Partially through ENCODE3 the DCC switched to needing to upload via AWS
"""

import logging
import json
import os
import pandas
import subprocess
import time

from .encoded import DCCValidator

logger = logging.getLogger(__name__)


def run_aws_cp(pathname, creds):
    env = os.environ.copy()
    env.update(
        {
            "AWS_ACCESS_KEY_ID": creds["access_key"],
            "AWS_SECRET_ACCESS_KEY": creds["secret_key"],
            "AWS_SECURITY_TOKEN": creds["session_token"],
        }
    )
    start = time.time()
    try:
        subprocess.check_call(
            ["aws", "s3", "cp", pathname, creds["upload_url"]], env=env
        )
    except subprocess.CalledProcessError as e:
        logger.error("Upload of %s failed with exit code %d", pathname, e.returncode)
        return
    else:
        end = time.time()
        logger.info("Upload of %s finished in %.2f seconds", pathname, end - start)


def process_files(server, files, dry_run):
    """Validate sheet and then upload files
    """
    logger.info('Validating metadata')
    validator = DCCValidator(server=server)

    logger.info('Uploading files')
    upload(server, validator, files, dry_run=dry_run)


def upload(server, validator, files, dry_run=True, retry=False):
    created = []
    print(files.columns)
    to_create = server.prepare_objects_from_sheet('/files/', files, validator=validator)
    for i, new_object in to_create:
        if new_object is not None and pandas.isnull(new_object.get('accession')):
            print('Would upload {}'.format(new_object['submitted_file_name']))
            posted_object = upload_file(server, validator, new_object, dry_run, retry)
            created.append(posted_object)

            if posted_object:
                accession = posted_object.get('accession')
                uuid = posted_object.get('uuid')
                if 'accession' in files.columns:
                    files.loc[i, 'accession'] = accession
                if 'uuid' in files.columns:
                    files.loc[i, 'uuid'] = uuid

    logger.info('created {} objects'.formt(len(created)))
    return created


def upload_file(encode, validator, metadata, dry_run=True, retry=False):
    """Upload a file to the DCC

    :Parameters:
      - encode: ENCODED instance pointing to server to upload to
      - validator: DCCValidator instance
      - dry_run: bool indicating if this is for real
      - retry: try uploading again.
    """
    if not isinstance(validator, DCCValidator):
        raise RuntimeError("arguments to upload_file changed")

    validator.validate(metadata, "file")

    file_name_fields = ["submitted_file_name", "pathname:skip", "pathname"]
    file_name_field = None
    for field in file_name_fields:
        if field in metadata and os.path.exists(metadata[field]):
            file_name_field = field

    if file_name_field is None:
        logger.error("Couldn't find file name to upload in metadata")
        logger.error(json.dumps(metadata, indent=4, sort_keys=True))
        return

    upload = make_upload_filename(metadata, encode)
    if retry or not os.path.exists(upload):
        logger.debug(json.dumps(metadata, indent=4, sort_keys=True))
        if not dry_run:
            item = post_file_metadata(encode, metadata, upload, retry)
            creds = item["upload_credentials"]
            run_aws_cp(metadata[file_name_field], creds)
            return item
        else:
            logger.info("Would upload %s", metadata[file_name_field])
            metadata["accession"] = "would create"
            return metadata
    else:
        logger.info("%s already uploaded", metadata[file_name_field])


def post_file_metadata(encode, metadata, upload, retry=False):
    """Post file metadata to ENCODE server

    :Paramters:
      - encode: ENCODED instance pointing to server to upload to
      - upload: name to store upload metadata cache
      - retry: bool if try, return saved metadata object,
               instead of posting a new object.
    """
    if not retry:
        response = encode.post_json("/file", metadata)
        logger.info(json.dumps(response, indent=4, sort_keys=True))
        with open(upload, "w") as outstream:
            json.dump(response, outstream, indent=4, sort_keys=True)

        item = response["@graph"][0]
    else:
        # Retry
        with open(upload, "r") as instream:
            item = json.load(instream)["@graph"][0]
    return item


def make_upload_filename(metadata, server=None):
    filename = metadata["submitted_file_name"].replace(os.path.sep, "_")
    if server is not None:
        extension = ".{}.upload".format(server.server)
    else:
        extension = ".upload"
    return filename + extension

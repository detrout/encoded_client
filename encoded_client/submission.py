"""Partially through ENCODE3 the DCC switched to needing to upload via AWS
"""
from argparse import ArgumentParser
import logging
import json
import os
import pandas
import subprocess
import time

from .encoded import ENCODED, DCCValidator, HTTPError
from .sheet import open_book, save_book

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
    return upload(server, validator, files, dry_run=dry_run)


def upload(server, validator, files, dry_run=True, retry=False):
    created = []
    to_create = server.prepare_objects_from_sheet('/files/', files, validator=validator)
    for i, new_object in to_create:
        if new_object is not None and pandas.isnull(new_object.get('accession')):
            logger.debug('Would upload {}'.format(new_object['submitted_file_name']))
            posted_object = upload_file(server, validator, new_object, dry_run, retry)
            created.append(posted_object)

            if posted_object:
                accession = posted_object.get('accession')
                uuid = posted_object.get('uuid')
                if 'accession' in files.columns:
                    files.loc[i, 'accession'] = accession
                if 'uuid' in files.columns:
                    files.loc[i, 'uuid'] = uuid

    logger.info('created {} objects'.format(len(created)))
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
        # was it already submitted?
        md5 = metadata["md5sum"]
        submitted = None
        try:
            submitted = encode.get_json("md5:{}".format(md5))
            logger.info("Previous object found: {}".format(submitted["accession"]))
            if submitted["status"] == "uploading":
                creds = encode.get_json("{}/upload/".format(submitted["@id"]))
                logger.info("Previous object has upload credentials: {}".format(submitted["accession"]))
                if "@graph" in creds:
                    creds = creds["@graph"]
                if isinstance(creds, list) and len(creds) == 1:
                    creds = creds[0]
                    if creds["@id"] == submitted["@id"]:
                        if "upload_credentials" in creds:
                            submitted["upload_credentials"] = creds["upload_credentials"]
                            with open(upload, "a") as outstream:
                                json.dump(submitted, outstream, indent=4, sort_keys=True)
                        return submitted
                    else:
                        logger.error("Credential id {} does not match object id {}".format(creds["@id"], submitted["@id"]))
        except HTTPError as e:
            if e.response.status_code == 404:
                pass
            else:
                print("{}".format(e))


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


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument(
        '-s',
        '--server',
        required=True,
        choices=[
            'www.encodeproject.org',
            'test.encodedcc.org',
        ],
        help='DCC Server to upload to'
    )
    parser.add_argument(
        '-f',
        '--spreadsheet-file',
        required=True,
        help='Metadata spreadsheet filename'
    )
    parser.add_argument(
        '-t',
        '--sheet-name',
        default='File',
        help='Override default file sheet name'
    )
    parser.add_argument('-o', '--output-file', help='Write sheet progress')
    parser.add_argument('-n', '--dry-run', action='store_true', default=False)
    args = parser.parse_args(cmdline)

    logging.basicConfig(level=logging.INFO)

    logging.info('Server: %s', args.server)
    logging.info('Spreadsheet: %s', args.spreadsheet_file)
    server = ENCODED(args.server)
    server.load_netrc()

    book = open_book(args.spreadsheet_file)
    files = book.parse(args.sheet_name, header=0)
    try:
        process_files(server, files, args.dry_run)
    except Exception as e:
        logger.error(e)

    if args.output_file:
        save_book(args.output_file, book, {'File': files})


if __name__ == "__main__":
    main()

"""Partially through ENCODE3 the DCC switched to needing to upload via AWS
"""

import logging
import json
import os
from pprint import pformat
import subprocess
import time

from requests.exceptions import HTTPError

from .submission import Submission
from .encoded import ENCODED, DCCValidator

LOGGER = logging.getLogger(__name__)


class AWSSubmission(Submission):
    def __init__(self, name, model, encode_host, lims_host):
        """Create a AWS based submission

        :Parameters:
          - `name`: Name of submission
          - `model`: librdf model reference
          - `host`: hostname for library pages.
        """
        super(AWSSubmission, self).__init__(name, model, lims_host)
        self.encode = ENCODED(encode_host)
        self.validator = DCCValidator(self.encode)
        self.encode.load_netrc()

        self._replicates = {}
        self._files = {}

    def check_upload(self, results_map):
        tocheck = []
        # phase one download data
        for an_analysis in self.analysis_nodes(results_map):
            for metadata in self.get_metadata(an_analysis):
                filename = make_upload_filename(metadata)
                if os.path.exists(filename):
                    with open(filename, "rt") as instream:
                        uploaded = json.load(instream)

                    if "submitted_file_name" not in uploaded:
                        LOGGER.warning("No submitted_file_name: %s", pformat(uploaded))
                    else:
                        tocheck.append(
                            {
                                "submitted_file_name": uploaded["submitted_file_name"],
                                "md5sum": uploaded["md5sum"],
                            }
                        )
                        self.update_replicates(uploaded)

        # phase 2 make sure submitted file is there
        md5sums = set((self._files[f]["md5sum"] for f in self._files))
        submitted_file_names = set(
            (self._files[f]["submitted_file_name"] for f in self._files)
        )
        errors_detected = 0
        for row in tocheck:
            error = []
            if row["submitted_file_name"] not in submitted_file_names:
                error.append("!file_name")
            if row["md5sum"] not in md5sums:
                error.append("!md5sum")
            if error:
                print(
                    "{} failed {} checks".format(
                        row["submitted_file_name"], ", ".join(error)
                    )
                )
                errors_detected += 1

        if not errors_detected:
            print("No errors detected")

    def update_replicates(self, metadata):
        replicate_id = metadata["replicate"]
        if replicate_id not in self._replicates:
            LOGGER.debug("Downloading %s", replicate_id)
            try:
                rep = self.encode.get_json(replicate_id)

                self._replicates[replicate_id] = rep
                for file_id in rep["experiment"]["files"]:
                    self.update_files(file_id)
            except HTTPError as err:
                print(err.response, dir(err.response))
                if err.response.status_code == 404:
                    print(
                        "Unable to find {} {}".format(
                            replicate_id, metadata["submitted_file_name"]
                        )
                    )
                else:
                    raise err

    def update_files(self, file_id):
        if file_id not in self._files:
            LOGGER.debug("Downloading %s", file_id)
            try:
                file_object = self.encode.get_json(file_id)
                self._files[file_id] = file_object
            except HTTPError as err:
                if err.response.status_code == 404:
                    print("unable to find {}".format(file_id))
                else:
                    raise err

    def upload(self, results_map, dry_run=False, retry=False):
        for an_analysis in self.analysis_nodes(results_map):
            for metadata in self.get_metadata(an_analysis):
                upload_file(self.encode, self.validator, metadata, dry_run, retry)

    def get_metadata(self, analysis_node):
        from django.template import loader

        # convert our model names to encode project aliases
        platform_alias = {"Illumina HiSeq 2500": "ENCODE:HiSeq2500"}
        run_type_alias = {
            "Single": "single-ended",
            "Paired": "paired-ended",
        }
        query_template = loader.get_template("aws_metadata.sparql")

        context = {
            "submission": str(analysis_node),
            "submissionSet": str(self.submissionSetNS[""]),
        }
        results = self.execute_query(query_template, context)
        LOGGER.info("scanned %s for results found %s", str(analysis_node), len(results))

        # need to adjust the results of the query slightly.
        for row in results:
            if "platform" in row:
                row["platform"] = platform_alias[row["platform"]]
            if "read_length" in row:
                row["read_length"] = int(row["read_length"])
            if "run_type" in row:
                row["run_type"] = run_type_alias[str(row["run_type"])]
            flowcell_details = {}
            for term in ["machine", "flowcell", "lane", "barcode"]:
                if term in row:
                    value = str(row[term])
                    flowcell_details[term] = value
                    del row[term]
            if len(flowcell_details) > 0:
                row["flowcell_details"] = [flowcell_details]

        return results


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
        LOGGER.error("Upload of %s failed with exit code %d", pathname, e.returncode)
        return
    else:
        end = time.time()
        LOGGER.info("Upload of %s finished in %.2f seconds", pathname, end - start)


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
        LOGGER.error("Couldn't find file name to upload in metadata")
        LOGGER.error(json.dumps(metadata, indent=4, sort_keys=True))
        return

    upload = make_upload_filename(metadata, encode)
    if retry or not os.path.exists(upload):
        LOGGER.debug(json.dumps(metadata, indent=4, sort_keys=True))
        if not dry_run:
            item = post_file_metadata(encode, metadata, upload, retry)
            creds = item["upload_credentials"]
            run_aws_cp(metadata[file_name_field], creds)
            return item
        else:
            LOGGER.info("Would upload %s", metadata[file_name_field])
            metadata["accession"] = "would create"
            return metadata
    else:
        LOGGER.info("%s already uploaded", metadata[file_name_field])


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
        LOGGER.info(json.dumps(response, indent=4, sort_keys=True))
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

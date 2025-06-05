"""Interface for the encoded software for IGVF and ENCODE object database

The classes here are intended to help manage the HTTP requsts and authentication
necessary to return the JSON-LD objects.
"""
from __future__ import print_function
from datetime import datetime
import pandas
import base64
from collections.abc import Sequence, Iterable, Mapping
import hashlib
import logging
from pathlib import Path
import json
import jsonschema
import os
from pprint import pformat
import re
import requests
from requests import HTTPError
from uuid import UUID

import six
from urllib.parse import (
    quote_plus,
    urljoin,
    urlparse,
    urlunparse
)

LOGGER = logging.getLogger(__name__)

ENCODED_CONTEXT = {
    # The None context will get added to the root of the tree and will
    # provide common defaults.
    None: {
        # terms in multiple encoded objects
        "award": {"@type": "@id"},
        "dataset": {"@type": "@id"},
        "description": "rdf:description",
        "documents": {"@type": "@id"},
        "experiment": {"@type": "@id"},
        "href": {"@type": "@id"},
        "lab": {"@type": "@id"},
        "library": {"@type": "@id"},
        "pi": {"@type": "@id"},
        "platform": {"@type": "@id"},
        "replicates": {"@type": "@id"},
        "submitted_by": {"@type": "@id"},
        "url": {"@type": "@id"},
    },
    # Identify and markup contained classes.
    # e.g. in the tree there was a sub-dictionary named 'biosample'
    # That dictionary had a term 'biosample_term_id, which is the
    # term that should be used as the @id.
    "Biosample": {
        "biosample_term_id": {"@type": "@id"},
    },
    "Experiment": {
        "assay_term_id": {"@type": "@id"},
        "files": {"@type": "@id"},
        "original_files": {"@type": "@id"},
    },
    # I tried to use the JSON-LD mapping capabilities to convert the lab
    # contact information into a vcard record, but the encoded model
    # didn't lend itself well to the vcard schema
    # 'lab': {
    #    "address1": "vcard:street-address",
    #    "address2": "vcard:street-address",
    #    "city": "vcard:locality",
    #    "state": "vcard:region",
    #    "country": "vcard:country"
    # },
    "Library": {"nucleic_acid_term_id": {"@type": "@id"}},
}

# FIXME: this needs to be initialized from rdfns
ENCODED_NAMESPACES = {
    # JSON-LD lets you define namespaces so you can used the shorted url syntax.
    # (instead of http://www.w3.org/2000/01/rdf-schema#label you can do
    # rdfs:label)
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
    "owl": "http://www.w3.org/2002/07/owl#",
    "dc": "htp://purl.org/dc/elements/1.1/",
    "xsd": "http://www.w3.org/2001/XMLSchema#",
    "vcard": "http://www.w3.org/2006/vcard/ns#",
    # for some namespaces I made a best guess for the ontology root.
    "EFO": "http://www.ebi.ac.uk/efo/",  # EFO ontology
    "OBO": "http://purl.obolibrary.org/obo/",  # OBO ontology
    "OBI": "http://purl.obolibrary.org/obo/OBI_",  # Ontology for Biomedical Investigations
    # OBI: available from http://svn.code.sf.net/p/obi/code/releases/2012-07-01/merged/merged-obi-comments.owl
    "SO": "http://purl.obolibrary.org/obo/SO_",  # Sequence ontology
    # SO: available from http://www.berkeleybop.org/ontologies/so.owl
    # NTR: New Term Request space for DCC to implement new ontology terms
}

ENCODED_SCHEMA_ROOT = "/profiles/"
SCHEMA_TYPE_OVERRIDES = {
    "aggregateseries": "aggregate_series",
    "analysissteprun": "analysis_step_run",
    "mousedonor": "mouse_donor",
    "publicationdata": "publication_data",
    "samtoolsflagstatsqualitymetric": "samtools_flagstats_quality_metric",
}

COLLECTION_TO_TYPE = {
    "/aggregate-series/": "AggregateSeries",
    "/analysis-step-runs/": "AnalysisStepRun",
    "/annotations/": "Annotation",
    "/award/": "Award",
    "/biosamples/": "Biosample",
    "/datasets/": "Dataset",
    "/documents/": "Document",
    "/experiments/": "Experiment",
    "/labs/": "Lab",
    "/libraries/": "Library",
    "/mouse-donors/": "MouseDonor",
    "/organisms/": "Organism",
    "/replicates/": "Replicate",
    "/star_quality_metric": "StarQualityMetric",
    "/files/": "File",
}

TYPE_TO_COLLECTION = {COLLECTION_TO_TYPE[k]: k for k in COLLECTION_TO_TYPE}


document_mime_type_default = {
    ".pdf": "application/pdf",
    ".tar": "application/x-tar",
    ".json": "application/json",
    ".gif": "image/gif",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".png": "image/png",
    ".tif": "image/tiff",
    ".tiff": "image/tiff",
    ".html": "text/html",
    ".txt": "text/plain",
    ".tsv": "text/tab-separated-values",
    ".yaml": "text/plain",  # not actually true, but it's what's supported
}


class DuplicateAliasError(jsonschema.ValidationError):
    """Exception to indicate that an alias was reused in this submission

    aliases on the portal must be unique
    """
    pass


def get_object_type(obj):
    """Return type for a encoded object

    Types in the encoded database are arrainged most specific to least
    specific. For instance an IGVF SequenceFile's type list is:

    ['SequenceFile', 'File', 'Item']

    And this will return the most specific type SequenceFile.

    :param obj: A dictionary containing a JSON-LD object
    :rtype str: The first type

    :raises ValueError: if object lacks an indication that this looks like a
                        JSON-LD object this will raise a ValueError.
    """
    obj_type = obj.get("@type")
    if not obj_type:
        raise ValueError("None type")
    if isinstance(obj_type, six.string_types):
        raise ValueError("@type should be a list, not a string")
    if not isinstance(obj_type, Sequence):
        raise ValueError("@type is not a sequence")
    return obj_type[0]


class ENCODED:
    """Programatic object orientated access to encoded

    encoded is the software powering the IGVF and ENCODE object database.
    It's also been used in some other projects as well.

    :param server: Domain name we are expecting to connect to.
    :type server: str
    :param context: A JSON-LD context when None uses ENCODED_CONTEXT as a default
    :type context: dict, optional
    :param namespaces: A dictionary of short CURIE names to urls, if None defaults to
                       ENCODED_NAMESPACES
    :type namespaces: dict, optional
    """

    def __init__(self, server, contexts=None, namespaces=None):
        self.server = server
        self.scheme = "https"
        self._session = requests.session()
        self.username = None
        self.password = None
        self.load_auth()
        self._user = None
        self.contexts = contexts if contexts else ENCODED_CONTEXT
        self.namespaces = namespaces if namespaces else ENCODED_NAMESPACES
        self.json_headers = {
            "content-type": "application/json",
            "accept": "application/json",
        }

    @property
    def auth(self):
        """A basic HTTP authentication tuple or None"""
        if self.username is None and self.password is None:
            return None

        return (self.username, self.password)

    def load_auth(self):
        """Abstract finding the access tokens

        This will first try looking for the DCC_API_KEY/DCC_SECRET_KEY
        environment variables, and fhat fails, then will look in a
        .netrc file.

        If the server and username are provided but no password, it
        will also attempt to use the python module keyring to retrieve
        the password.
        https://pypi.org/project/keyring/

        If self.username and self.password are already set this
        function does nothing.
        """
        if self.username is None:
            self.username = os.environ.get("DCC_API_KEY")
            self.password = os.environ.get("DCC_SECRET_KEY")

        if self.username is None:
            self.load_netrc()

        if self.username is not None and self.password is None:
            try:
                import keyring
                self.password = keyring.get_password(self.server, self.username)
            except ImportError:
                LOGGER.info("Unable to load password via keyring")
                pass

    def load_netrc(self):
        """Load authentication tokens from the .netrc file

        Uses the server name as the machine name for the look up for
        the username and password
        """
        netrc_path = Path("~/.netrc").expanduser()
        if netrc_path.exists():
            import netrc

            session = netrc.netrc()
            authenticators = session.authenticators(self.server)
            if authenticators is not None:
                self.username = authenticators[0]
                self.password = authenticators[2]

    def add_jsonld_context(self, tree, default_base):
        """Add contexts to various objects in the tree.

        This adds in some JSON-LD @context to the root and child
        objects if they were missing

        :param tree: is a nested tree of JSON-LD objects returned from
                     the object database.
        :type tree: dict
        :param default_base: if supplied allows setting the base url that relative
                             urls will be resolved against.
        :type default_base: str

        :return: Nothing, though as a side effect it modifies the
                 passed in objects in place
        """
        self.add_jsonld_child_context(tree, default_base)
        self.add_jsonld_namespaces(tree["@context"])

    def add_jsonld_child_context(self, obj, default_base):
        """Add JSON-LD context to the encoded JSON.

        This is recursive because some of the IDs were relative URLs
        and I needed a way to properly compute a the correct base URL.

        :param obj: An object retrieved from the portal
        :type obj: dict
        :param default_base: the default base url to use for this object
        :type default_base: str

        See https://www.w3.org/TR/json-ld/#the-context for information
        about the JSON-LD context

        """
        # pretend strings aren't iterable
        if isinstance(obj, six.string_types):
            return

        # recurse on container types
        if isinstance(obj, Sequence):
            # how should I update lists?
            for v in obj:
                self.add_jsonld_child_context(v, default_base)
            return

        if isinstance(obj, Mapping):
            for v in obj.values():
                self.add_jsonld_child_context(v, default_base)

        # we have an object. attach a context to it.
        if self._is_encoded_object(obj):
            context = self.create_jsonld_context(obj, default_base)
            if len(context) > 0:
                # this is a total hack for relese 33 of
                # encoded. They changed their model and
                # i'm not sure what to do about it.
                if obj.get("@context") == "/terms/":
                    del obj["@context"]
                obj.setdefault("@context", {}).update(context)

    def add_jsonld_namespaces(self, context):
        """Add shortcut namespaces to a context

        Only needs to be run on the top-most context to set the
        default url mappings.

        :param context: a dictionary that represent a JSON-LD context.
        :type context: dict

        """
        context.update(self.namespaces)

    def create_jsonld_context(self, obj, default_base):
        """Synthesize the context for a encoded type

        self.contexts[None] = default context attributes added to any type
        self.contexts[type] = context attributes for this type.

        :param obj: An object from that portal that might need to have more
                    details added to JSON-LD context.
        :type obj: dict
        :param default_base: The default base url to use for object ids.
        :type default_base: str

        :returns: The completed JSON-LD context
        :rtype: dict
        """
        obj_type = get_object_type(obj)
        context = {
            "@base": urljoin(default_base, obj["@id"]),
            "@vocab": self.get_schema_url(obj_type),
        }
        # add in defaults
        context.update(self.contexts[None])
        for t in obj["@type"]:
            if t in self.contexts:
                context.update(self.contexts[t])
        return context

    def get_experiment(self, obj_id):
        """Retrieve an ENCODE experiment as class

        This was an attempt to streamline some trying to find relationships between
        the various file objects, but is specific to ENCODE.

        :param obj_id: The @id for an experiment

        :returns: a class initialized with the JSON-LD object
        :rtype: EncodeExperiment
        """
        obj = self.get_json(obj_id)
        if obj["@type"][0] != "Experiment":
            raise ValueError("Object is not an experiment")

        return EncodeExperiment(self, obj)

    def get_json(self, obj_id, **kwargs):
        """GET an ENCODE object as JSON and return as dict

        Uses prepare_url to allow url short-cuts
        if no keyword arguments are specified it will default to adding limit=all
        Alternative keyword arguments can be passed in and will be sent to the host.

        see :py:meth:`ENCODED.get_response` for other known keyword arguments that
        can be passed to this function.

        :param obj_id: an object id url fragment, typically something like
                       "/measurement-set/accession/"
        :type obj_id: str

        :returns: the JSON retrieved from the portal
        :rtype: dict
        """
        kwargs["headers"] = self.json_headers

        response = self.get_response(obj_id, **kwargs)
        data = response.json()
        response.close()
        return data

    def get_jsonld(self, obj_id, **kwargs):
        """Get ENCODE object as JSONLD annotated with classses contexts

        see get_json for documentation about what keywords can be passed.

        :param obj_id: an object id url fragment, typically something like
                       "/measurement-set/accession/"
        :type obj_id: str

        Also accepts other arbitrary keyword objects see
        :py:meth:`ENCODED.get_response` for known arguments

        :returns: the JSON retrieved from the portal, with any additional
                  JSON-LD context values added to the object
        :rtype: dict

        """
        url = self.prepare_url(obj_id)
        json = self.get_json(obj_id, **kwargs)
        self.add_jsonld_context(json, url)
        return json

    def get_response(self, fragment, **kwargs):
        """GET an ENCODED url and return the requests request

        Uses prepare_url to allow url short-cuts
        if no keyword arguments are specified it will default to adding limit=all
        Alternative keyword arguments can be passed in and will be sent to the host.

        Known keywords are:

            limit
                (integer or 'all') how many records to return, all for all of them

            embed
                (bool) if true expands linking ids into their associated object.

            format
                text/html or application/json

        :param fragment: a url fragment
        :type fragment: str

        This function accepts an arbitrary set of other keyword arguments
        which will be passed as arguments to the database server.

        :returns: returns the Response object from requests
        :rtype: requests.Response
        """
        url = self.prepare_url(fragment)
        LOGGER.info("requesting url: {}".format(url))

        # do the request
        LOGGER.debug("username: %s, password: %s", self.username, self.password)
        arguments = {}
        if self.username and self.password:
            arguments["auth"] = self.auth

        if "stream" in kwargs:
            arguments["stream"] = kwargs["stream"]
            del kwargs["stream"]

        if "headers" in kwargs:
            arguments["headers"] = kwargs["headers"]
            del kwargs["headers"]

        response = self._session.get(
            url, params=kwargs, **arguments
        )
        if not response.status_code == requests.codes.ok:
            LOGGER.warning(
                "Error http status: {} for {}".format(response.status_code, url)
            )
            response.raise_for_status()

        return response

    def get_schema_url(self, object_type):
        """Create the ENCODED jsonschema url.

        Return the ENCODED object schema url be either
        object type name or the collection name one posts to.

        For example

        get_schema_url('experiment') and
        get_schema_url('/experiments/') both resolve to
        SERVER/profiles/experiment.json

        :param object_type: either an ENCODED object name or a collection
        :type object_type: str

        :returns: Schema URL

        """
        object_type = COLLECTION_TO_TYPE.get(object_type, object_type).lower()
        schema_name = SCHEMA_TYPE_OVERRIDES.get(object_type, object_type)

        return self.prepare_url(ENCODED_SCHEMA_ROOT + schema_name + ".json") + "#"

    def get_accession_name(self, collection):
        """Lookup common object accession name given a collection name.

        This is used by the sheet parsing code to know what column the
        objects accession is in.

        :param collection: Object at the portal are grouped by their
                           collection type
        :type collection: str

        :return: given the name of a collection object return what the primary
                 key attribute name should be. Usually it's accession, but
                 there were a few objects that used uuid.
        :type: str
        """
        collection_to_accession_name = {
            "analysis_set": "accession",
            "alignment_file": "accession",
            "/annotations/": "accession",
            "/biosamples/": "accession",
            "configuration_file": "accession",
            "curated_set": "accession",
            "document": "uuid",
            "/documents/": "uuid",
            "/experiments/": "accession",
            "/files/": "accession",
            "genome_browser_annotation_file": "accession",
            "image_file": "accession",
            "/libraries/": "accession",
            "matrix_file": "accession",
            "measurement_set": "accession",
            "model_file": "accession",
            "/mouse-donors/": "accession",
            "reference_file": "accession",
            "/replicates/": "uuid",
            "rodent_donor": "accession",
            "sequence_file": "accession",
            "signal_file": "accession",
            "tabular_file": "accession",
            "tissue": "accession",
            "multiplexed_sample": "accession",
        }

        accession_name = collection_to_accession_name.get(collection, None)
        if accession_name is None:
            raise RuntimeError(
                "Update list of collection to accession names for {}".format(collection)
            )

        return accession_name

    def _is_encoded_object(self, obj):
        """Test to see if an object is a JSON-LD object

        Some of the nested dictionaries lack the @id or @type
        information necessary to convert them.

        :param obj: An dictionary that may contain an encoded object
        :type obj: dict

        :returns: True if it is has a mapping and has an @id and @type
                  object, otherwise false.
        :rtype: bool
        """
        if not isinstance(obj, Mapping):
            return False

        if "@id" in obj and "@type" in obj:
            return True
        return False

    def patch_json(self, obj_id, changes):
        """Given a dictionary of changes push them as a HTTP patch request

        :param obj_id: the object id as a url or url fragment that
                       should be modified.
        :type obj_id: str
        :param changes: a dictionary of attributes that should be
                        updated at the portal
        :type changes: dict

        :raises: requests.HTTPError if something goes wrong with the PATCH
        """
        url = self.prepare_url(obj_id)
        LOGGER.info("PATCHing to %s", url)
        payload = json.dumps(changes)
        response = requests.patch(
            url, auth=self.auth, headers=self.json_headers, data=payload
        )
        if response.status_code != requests.codes.ok:
            LOGGER.error("Error http status: {}".format(response.status_code))
            LOGGER.error("Response: %s", response.text)
            response.raise_for_status()
        return response.json()

    def put_json(self, obj_id, new_object):
        """Create a new object at the portal using a PUT request

        :param obj_id: the object id as a url or url fragment that
                       should be modified.
        :type obj_id: str
        :param new_object: a dictionary of a metadata object ready to
                           create at the portal
        :type new_object: dict

        :returns: json object from the portal

        :raises: requests.HTTPError if something goes wrong with the PUT
        """
        url = self.prepare_url(obj_id)
        LOGGER.info("PUTing to %s", url)
        payload = json.dumps(new_object)
        response = requests.put(
            url, auth=self.auth, headers=self.json_headers, data=payload
        )
        if response.status_code != requests.codes.created:
            LOGGER.error("Error http status: {}".format(response.status_code))
            response.raise_for_status()
        return response.json()

    def post_json(self, collection_id, new_object):
        """Create a new object at the portal using a POST request

        :param obj_id: the object id as a url or url fragment that
                       should be modified.
        :type obj_id: str
        :param new_object: a dictionary of a metadata object ready to
                           create at the portal
        :type new_object: dict

        :returns: the response json object
        :rtype: dict

        :raises: requests.HTTPError if something goes wrong with the POST
        """
        url = self.prepare_url(collection_id)
        LOGGER.info("POSTing to %s", url)
        payload = json.dumps(new_object)

        response = requests.post(
            url, auth=self.auth, headers=self.json_headers, data=payload
        )
        if response.status_code != requests.codes.created:
            LOGGER.error("http status: {}".format(response.status_code))
            LOGGER.error("message: {}".format(response.content))
            response.raise_for_status()
        return response.json()

    def post_sheet(
        self, collection, sheet, dry_run=True, verbose=False, validator=None
    ):
        """Create new ENCODED objects using metadata encoded in pandas DataFrame

        The DataFrame column names need to encode the attribute names,
        and in some cases also include some additional type information.
        (see TypedColumnParser)

        :param collection: name of collection to create new objects in
        :type collection: str
        :param sheet: DataFrame with objects to create,
                      assuming the appropriate accession number is empty.
                      additional the accession number and uuid is updated if
                      the object is created.
        :type sheet: pandas.DataFrame
        :param dry_run: whether or not to skip the code to post the objects
        :type dry_run: bool
        :param verbose: print the http responses.
        :type verbose: bool

        :returns: list of created objects.
        :rtype: list[dict]

        :raises: jsonschema.ValidationError if the object doesn't validate against
                 the encoded jsonschema.
        """
        accession_name = self.get_accession_name(collection)

        to_create = self.prepare_objects_from_sheet(
            collection, sheet, validator=validator
        )

        created = []
        for i, new_object in to_create:
            if new_object:
                accession = new_object.get(accession_name)
                uuid = new_object.get("uuid")
                description = new_object.get("description")

                posted_object = self.post_object_from_row(
                    collection, i, new_object, dry_run, verbose
                )
                created.append(posted_object)

                if posted_object:
                    accession = posted_object.get("accession")
                    uuid = posted_object.get("uuid")
                    if "accession" in sheet.columns:
                        sheet.loc[i, "accession"] = accession
                    if "uuid" in sheet.columns:
                        sheet.loc[i, "uuid"] = uuid

                    description = posted_object.get("description")

                if description is None:
                    description = new_object.get("aliases")

                LOGGER.info("row {} ({}) -> {}".format((i + 2), description, accession))
                # +2 comes from python row index + 1 to convert to
                # one based indexing + 1 to account for
                # row removed by header parsing

        return created

    def prepare_objects_from_sheet(self, collection, sheet, validator=None):
        """Do a batch posting of a number of objects for a specific collection

        This is typically used to create a number of objects for a particular
        object type at the portal.

        :param collection: An endpoint that represents a collection. Typically
                           things like "sequencing_file" or "measurement_set".
        :type collection: str
        :param sheet: A spreadsheet like table with rows representing objects and
                      columns representing the various attributes in those objects.
        :type sheet: pandas.DataFrame
        :param validator: the objects created from the sheet will be passed
                          through the DCCValidator to make sure they are correct.
                          However one needs to use a shared DCCValidator object
                          to resolve aliases defined between sheets.
        :type validator: DCCValidator, optional
        """
        accession_name = self.get_accession_name(collection)
        to_create = []
        if validator is None:
            validator = DCCValidator(self)

        # delete any aliases for this sheet.
        self.remove_sheet_aliases(validator, sheet)

        for i, row in sheet.iterrows():
            new_object = {}
            for name, value in row.items():
                if pandas.notnull(value):
                    name, value = typed_column_parser(name, value)
                    if name is None:
                        continue
                    new_object[name] = value

            if len(new_object) > 0:
                try:
                    if new_object.get(accession_name) is None or (
                        not new_object.get(accession_name).lower().startswith("skip")
                    ):
                        validator.validate(new_object, collection)
                except jsonschema.ValidationError as e:
                    LOGGER.error("Validation error row %s", i)
                    raise e

                if new_object.get(accession_name) is None:
                    to_create.append((i, new_object))
                else:
                    to_create.append((i, None))
            else:
                to_create.append((i, None))

        return to_create

    def post_object_from_row(self, collection, i, new_object, dry_run=True, verbose=True):
        """Post one row of a dataframe.

        :param collection:
        :type collection: str
        :param i:
        :type i: int
        :param new_object:
        :type new_object: dict
        :param dry_run: Defaults to True, and if True it just runs the
                        validation and marks if the object would be
                        created. When false it actually posts the object
                        to the portal
        :type dry_run: bool, optional
        :param verbose: If True print the requests.Response object. This can be
                        useful for recovery if it crashes mid notebook.
        :type verbose: bool, optional

        :return: The created object from the portal
        :rtype: dict

        :raises: requests.HTTPError
        """
        accession_name = self.get_accession_name(collection)

        if not dry_run:
            response = self.post_json(collection, new_object)
            if verbose:
                print("Reponse {}".format(response))

            obj = response["@graph"][0]

            accession = obj.get(accession_name)
            if not accession:
                accession = obj.get("uuid")

            print("row {} created: {}".format(i, accession))
            return obj
        else:
            new_object[accession_name] = "would create"
            return new_object

    def prepare_url(self, request_url):
        """This attempts to provide some convienence for accessing a URL

        Given a url fragment it will default to :
        * requests over http
        * requests to self.server

        This allows fairly flexible urls. e.g:

          - prepare_url('/experiments/ENCSR000AEG')
          - prepare_url('submit.encodedcc.org/experiments/ENCSR000AEG')
          - prepare_url('http://submit.encodedcc.org/experiments/ENCSR000AEG?limit=all')

        should all return the same url

        :param request_url: some portion of a url fragment
        :type request_url: str

        :return: A complete and fully qualified URL.
        :rtype: str
        """
        # clean up potentially messy urls
        url = urlparse(request_url)._asdict()
        if not url["scheme"]:
            url["scheme"] = self.scheme
        elif url["scheme"] not in ("http", "https"):
            # this is an alias
            url["scheme"] = self.scheme
            url["path"] = request_url
        if not url["netloc"]:
            url["netloc"] = self.server
        url = urlunparse(url.values())
        return url

    def remove_sheet_aliases(self, validator, sheet):
        """Remove aliases from a validator for a particular sheet

        When repeatedly running the validator, the check to see if
        aliases are duplicated will be triggered unless we remove them
        before we re-validate the same sheet.

        :param validator: a validator object that needs to have some aliases
                          removed.
        :type validator: DCCValidator
        :param sheet: the current sheet of objects that need to have their
                      aliases removed.
        :type sheet: pandas.DataFrame

        """
        field_name = "aliases:array"
        if field_name in sheet.columns:
            for field in sheet[field_name]:
                name, aliases = typed_column_parser(field_name, field)
                for alias in aliases:
                    if alias in validator._aliases:
                        del validator._aliases[alias]

    def search_jsonld(self, **kwargs):
        """Send search request to ENCODED

        to do a general search call with

        search_jsonld(searchTerm="term")

        other known keywords allowed on search are:

        limit
             (integer or 'all') how many records to return, all for all of them
             to be kinder to the server, set an arbitrary limit of 5,000. override
             with all if all is really needed
        type
             portal object type with a leading @ sign

        """
        if len(kwargs) == 0:
            kwargs["limit"] = "5000"

        url = self.prepare_url("/search/")
        result = self.get_json(url, **kwargs)
        self.convert_search_to_jsonld(result)
        return result

    def convert_search_to_jsonld(self, result):
        """Add the context to search result

        Also remove hard to handle nested attributes
          e.g. remove object.term when we have no id
        """
        graph = result["@graph"]
        for i, obj in enumerate(graph):
            # suppress nested attributes
            graph[i] = {k: v for k, v in obj.items() if "." not in k}

        self.add_jsonld_context(result, self.prepare_url(result["@id"]))
        return result

    @property
    def user(self):
        """Return the user object that contains an access_key

        The programatic API uses an access code, which is different
        from the user ID.

        :return: The dictionary object representing the logged in user.
        :rtype: dict | None
        """
        if self._user is None:
            if self.username is None:
                return None

            access_key = self.get_json("/access_key/" + self.username)
            if access_key is None:
                return None

            self._user = self.get_json(access_key["user"])

        return self._user


class DCCValidator:
    def __init__(self, server):
        self._server = server
        self._aliases = {}
        self._http_cache = {}
        self._schemas = {}
        self._dcc_validators = {}
        self._current_request_method = "POST"

    def __getitem__(self, object_type):
        """Return customized validator for the object type

        :param object_type: One of the DCC object types like biosample or library
        :type object_type: str

        :returns: Either a cached jsonschema validator or creates a new one if needed
        """
        if object_type not in self._schemas:
            schema_url = self._server.get_schema_url(object_type)
            if not schema_url:
                raise ValueError("Unable to construct schema url")

            self._schemas[object_type] = self._server.get_json(schema_url)

        if object_type not in self._dcc_validators:
            schema = self._schemas[object_type]
            self._dcc_validators[object_type] = self.create_dcc_validator(schema)

        return self._dcc_validators[object_type]

    def get_json(self, object_id):
        """Caching get_json to speed up validation

        :param str object_id: DCC object id, can be a framgment
           a full url, or an alias
        """
        if object_id in self._aliases:
            item = self._aliases[object_id]
        elif object_id in self._http_cache:
            item = self._http_cache[object_id]
        else:
            item = self._server.get_json(object_id)
            self._http_cache[object_id] = item
        return item

    def validate(self, obj, object_type):
        """Validate an object against the ENCODED schema

        Also checks to make sure RNA-seq is spelled as RNA-seq and not any
        other capitalization

        :param obj: An object dictionary to be submitted to the portal
        :type obj: dict
        :param object_type: What portal schema should we validate this object
                            against.

        :raises: ValidationError if the object does not conform to the schema.
        """
        object_type = object_type if object_type else get_object_type(obj)

        hidden = self.strip_jsonld_attributes(obj)
        hidden = self.strip_uuid(hidden)
        self[object_type].validate(hidden)
        self.update_aliases(hidden)

        # Additional validation rules passed down from the DCC for our grant
        assay_term_name = hidden.get("assay_term_name")
        if assay_term_name is not None:
            if assay_term_name.lower() == "rna-seq":
                if assay_term_name != "RNA-seq":
                    raise jsonschema.ValidationError(
                        "Incorrect capitialization of RNA-seq"
                    )

    def strip_jsonld_attributes(self, obj):
        """Make copy of object with JSON-LD attributes

        :param obj: dictionary to POST to the DCC as an json document
        :type obj: dict

        :return: copy of object with jsonld attributes @id, and @type removed.
        :rtype: dict
        """
        hidden = obj.copy()
        if "@id" in hidden:
            del hidden["@id"]
        if "@type" in hidden:
            del hidden["@type"]
        return hidden

    def strip_uuid(self, obj):
        """Make copy of object with uuid removed

        :param obj: dictionary to POST to the DCC as an json document
        :type obj: dict

        :return: copy of object with uuid attribute removed.
        :rtype: dict
        """
        hidden = obj.copy()
        if "uuid" in hidden:
            del hidden["uuid"]
        return hidden

    def update_aliases(self, obj):
        """Save object under any listed aliases.

        This allows us to find it again for validation, even if it hasn't
        been submitted yet.

        :param obj: dictionary to POST to the DCC as an json document
        :type obj: dict

        :returns: None
        """
        # store aliases
        for a in obj.get("aliases", []):
            if a in self._aliases and self._aliases[a] != obj:
                raise DuplicateAliasError("{} was already used else where".format(a))
            self._aliases[a] = obj

    def create_dcc_validator(self, schema):
        """Return jsonschema validator customized for ENCODE DCC schemas

        The validator is also overloaded to handle validating of the
        more complex attribute types like linkTo or
        calculatedProperty.

        :param schema: dictionary containing a loaded jsonschema document
        :type schema: dict holding a jsonschema.

        :return: customized jsonschema validator for the provided schema
        :rtype: jsonschema.Draft4Validator

        """
        validator = jsonschema.validators.extend(
            jsonschema.Draft4Validator,
            validators={
                "calculatedProperty": self.calculatedPropertyValidator,
                "linkTo": self.linkToValidator,
                "linkFrom": self.linkFromValidator,
                "permission": self.permissionValidator,
                "requestMethod": self.requestMethodValidator,
            },
        )
        return validator(schema)

    def unimplementedValidator(self, *args):
        raise NotImplementedError("Unimplementated validator: %s" % (str(args)))

    def calculatedPropertyValidator(self, validator, tag, instance, schema):
        """Forbid submitting objects with a calculated property"""
        yield jsonschema.ValidationError(
            'submission of calculatedProperty "%s" is disallowed' % (tag,)
        )

    def linkFromValidator(self, validator, linkFrom, instance, schema):
        if schema["linkFrom"] == "QualityMetric.quality_metric_of":
            if "quality_metric_of" not in instance:
                yield jsonschema.ValidationError(
                    "required tag quality_metric_of missing"
                )
            object_id = instance["quality_metric_of"]
            obj = self.get_json(object_id)
            if obj is None:
                yield jsonschema.ValidationError(
                    "quality_metric_of {} was not found".format(object_id)
                )

    def linkToValidator(self, validator, linkTo, instance, schema):
        if not validator.is_type(instance, "string"):
            return

        # hack for ['Dataset']
        if isinstance(linkTo, list):
            # In release 57 one linkTo property had a single item as a list
            # the dcc is trying to get rid of it and convert it to a
            # string, but on the off chance they change and start using
            # lists again, lets flag this implementation because it'd be
            # inadequate.
            linkTo = linkTo[0]

        try:
            try:
                UUID(instance)
                object_id = instance
            except ValueError:
                # a hack to detect if we have an alias?
                if ":" in instance:
                    object_id = instance
                else:
                    collection = TYPE_TO_COLLECTION.get(linkTo)
                    object_id = urljoin(collection, instance)
            item = self.get_json(object_id)
        except HTTPError as e:
            yield jsonschema.ValidationError(
                "%s doesn't exist: %s" % (object_id, str(e))
            )

        linkEnum = schema.get("linkEnum")
        if linkEnum is not None:
            if not validator.is_type(linkEnum, "array"):
                raise Exception("Bad schema")
            if not any(enum_uuid == item["uuid"] for enum_uuid in linkEnum):
                reprs = ", ".join(repr(it) for it in linkTo)
                error = "%r is not one of %s" % (instance, reprs)
                yield jsonschema.ValidationError(error)
                return

        if schema.get("linkSubmitsFor"):
            if self._server.user is not None:
                user = self._server.user
                submits_for = [
                    self.get_json(s) for s in user.get("submits_for")
                ]
                if submits_for is not None and not any(
                    lab["uuid"] == item["uuid"] for lab in submits_for
                ):
                    name = "{} {}".format(
                        user.get("first_name", "Unknown"),
                        user.get("last_name", "User"))
                    error = "{} is not a lab in user {} can submits_for".format(
                        instance, name)
                    LOGGER.error(error)
                    return

    def permissionValidator(self, validator, permission, instance, schema):
        if not validator.is_type(permission, "string"):
            raise Exception("Bad Schema")

        admin_permissions = [
            "add_unvalidated",
            "edit_unvalidate",
            "impersonate",
            "import_items",
            "submit_for_any",
            "view_raw",
            # system permissions
            "expand",
            "index",
        ]

        if instance in admin_permissions:
            yield jsonschema.ValidationError(
                "Submitting property that requires admin permissions on %s"
                % (schema["title"])
            )

    def requestMethodValidator(self, validator, method, instance, schema):
        # see https://github.com/ENCODE-DCC/snovault src/snovault/schema_utils.py requestMethod
        # for DACCs implementation
        if isinstance(method, str):
            method = [method]

        if self._current_request_method not in method:
            yield jsonschema.ValidationError(
                "request method {} is not one of {} for {}".format(
                    self._current_request_method,
                    method,
                    schema["title"]))


class TypedColumnParser(object):
    """Functions to help parse columns annotated with types in DataFrames
    """
    @staticmethod
    def parse_sheet_array_type(value):
        """Helper function to parse :array columns in sheet"""
        if pandas.isnull(value):
            return []
        else:
            return re.split(r",\s*", value)

    @staticmethod
    def parse_sheet_integer_type(value):
        """Helper function to parse :integer columns in sheet"""
        return int(value)

    @staticmethod
    def parse_sheet_number_type(value):
        """Helper function to parse :number columns in sheet"""
        # JSON Schema seems to say number can be either
        # Integers or Floats, so lets try keeping the
        # type straight.
        try:
            value = int(value)
        except ValueError:
            value = float(value)
        return value

    @staticmethod
    def parse_sheet_boolean_type(value):
        """Helper function to parse :boolean columns in sheet"""
        return bool(value)

    @staticmethod
    def parse_sheet_timestamp_type(value):
        """Helper function to parse :date columns in sheet"""
        if isinstance(value, str):
            if "/" in value:
                parsed = datetime.strptime(value, "%m/%d/%Y")
            elif "-" in value:
                parsed = datetime.strptime(value, "%Y-%m-%d")
            else:
                raise ValueError("Unrecognized date format {}".format(value))

            LOGGER.warning("Interpreting {} as {}".format(value, parsed))
            value = parsed

        return value.strftime("%Y-%m-%d")

    @staticmethod
    def parse_sheet_string_type(value):
        """Helper function to parse :string columns in sheet (the default)"""
        return str(value)

    @staticmethod
    def parse_sheet_json_type(value):
        """Helper function to parse :json columns in sheet"""
        return json.loads(value)

    def __getitem__(self, name):
        parser = {
            "array": self.parse_sheet_array_type,
            "boolean": self.parse_sheet_boolean_type,
            "integer": self.parse_sheet_integer_type,
            "number": self.parse_sheet_number_type,
            "date": self.parse_sheet_timestamp_type,
            "json": self.parse_sheet_json_type,
            "string": self.parse_sheet_string_type,
        }.get(name)
        if parser:
            return parser
        else:
            raise RuntimeError("unrecognized column type: {}".format(name))

    def __call__(self, header, value):
        header = header.split(":")
        column_type = "string"
        if len(header) > 1:
            if header[1] == "skip":
                return None, None
            else:
                column_type = header[1]
        return header[0], self[column_type](value)


typed_column_parser = TypedColumnParser()


class Document:
    """Helper class for registering documents

    Typically used like this in a notebook.

    lysis_uuid = None
    lysis = Document(url_to_pdf, 'extraction protocol', 'Lysis Protocol')
    lysis.create_if_needed(server, lysis_uuid)

    The result of create_if_needed will include a uuid and then I copied that
    back up to replace the uuid so I could get get the same object, so after
    running this once I then have

    lysis_uuid = 'f0cc5a7f-96a5-4970-9f46-317cc8e2d6a4'

    """
    award = "U54HG006998"
    lab = "/labs/barbara-wold"

    def __init__(self, url, document_type, description, aliases=None, filename=None, server=None):
        """Construct a Document object.

        :param url: path to where the document file is located
        :type url: str or Path
        :param document_type: what is the document type from the list of known portal
                              document types
        :type document_type: str
        :param aliases: Various aliases to attach to this document
        :type aliases: list, optional
        :param filename: A filename used to override the automatic filename
                         derived from the url. Occasionally the path includes
                         values not allowed in the url fragment posted to
                         the portal.
        :type filename: str or Path
        :param server: Server object to use for querying the portal
        :type server: ENCODED
        """
        assert server is None or isinstance(server, ENCODED)

        self.url = Path(url)
        if filename is None:
            self.filename = os.path.basename(url)
        else:
            self.filename = filename
        self.document_type = document_type
        self.description = description

        self.references = []
        self.aliases = None
        if aliases:
            if isinstance(aliases, list):
                self.aliases = aliases
            else:
                raise ValueError("Aliases needs to be a list")
        self.content_type = None
        self.document = None
        self.md5sum = None
        self.urls = None
        self.uuid = None

        self.get_document(server)

    def get_document(self, server=None):
        """Get the document from either the local filesystem or the database

        It will try the local filesystem using url fragment relative
        to the current directory first, otherwise if a server is provided
        and the url starts with http or https it will request the document
        object from the server

        :param server: a configured ENCODED server object to use for API requests
        :type server: ENCODED, optional

        This will update several attributes of the current document object
        after being called.
        """
        parsed_url = urlparse(str(self.url))
        if self.url.exists():
            with open(self.url, "rb") as instream:
                if self.url.suffix in document_mime_type_default:
                    self.content_type = document_mime_type_default[self.url.suffix]
                else:
                    raise ValueError(
                        "Unrecognized filename extension {}".format(self.url.suffix))
                self.document = instream.read()
                self.md5sum = hashlib.md5(self.document)
        elif server is not None and parsed_url.scheme in ("http", "https"):
            req = server._session.get(self.url)
            if req.status_code == 200:
                self.content_type = req.headers["content-type"]
                self.document = req.content
                self.md5sum = hashlib.md5(self.document)
                self.urls = [self.url]
        else:
            raise ValueError("Unable to retrieve document {}".format(self.url))

    def create_payload(self):
        document_payload = {
            "attachment": make_attachment(
                self.url,
                mime_type=self.content_type,
                remote_filename=self.filename,
            ),
            "document_type": self.document_type,
            "description": self.description,
            "award": self.award,
            "lab": self.lab,
        }
        if self.aliases:
            document_payload["aliases"] = self.aliases
        if self.references:
            document_payload["references"] = self.references
        if self.urls:
            document_payload["urls"] = self.urls

        return document_payload

    def post(self, server, validator):
        document_payload = self.create_payload()
        validator.validate(document_payload, "document")
        return server.post_json("/documents/", document_payload)

    def save(self, filename):
        payload = self.create_payload()
        with open(filename, "w") as outstream:
            outstream.write(pformat(payload))

    def create_if_needed(self, server, uuid, validator):
        self.uuid = uuid
        if uuid is None:
            return self.post(server, validator)
        else:
            return server.get_json(uuid, embed=False)


class EncodeFileCache(Mapping):
    def __init__(self, server, dataset, json=None):
        self._server = server
        self._dataset = dataset
        self._file_cache = []

    def _update_json(self, json):
        for obj in json:
            assert obj.get("@type") == ["File", "Item"]
            assert obj.get("dataset") == self._dataset
            self._file_cache[obj["@id"]] = obj

    def __getitem__(self, key):
        if key not in self._file_cache:
            obj = self._server.get_json(key)
            self._update_json([obj])

        return self._file_cache[key]

    def __iter__(self):
        for key in self._file_cache:
            yield key

    def __len__(self):
        return len(self._file_cache)


class EncodeExperiment(Mapping):
    """Helper class for accessing ENCODED experiment objects"""

    def __init__(self, server, json=None):
        self._server = server
        self._json = json
        self._preferred_analysis = None
        self._replicates = []
        self._schema_version_check()
        self._replicate_file_map = {}
        if self._json is not None:
            self._files = EncodeFileCache(
                self._server, self._json["@id"], self._json["files"])
            self._calculate_derived_from()

    def _schema_version_check(self):
        supported_versions = {"33", "34", "35", "36", "37"}
        if self.schema_version not in supported_versions:
            LOGGER.warning(
                "New schema version {} may not be supported".format(self.schema_version)
            )

    def _calculate_derived_from(self):
        def find_source_replicate(derived_map, current_file, replicate_map):
            if current_file in replicate_map:
                return replicate_map[current_file]
            return find_source_replicate(
                derived_map, derived_map[current_file], replicate_map
            )

        file_replicate_map = {}

        # this finds replicate information attached to read files
        our_files = set((f["@id"] for f in self._json["files"]))
        for f in self._json["files"]:
            if "replicate" in f:
                file_replicate_map[f["@id"]] = f["replicate"]["@id"]

        derived_map = {}
        for file_id in self.preferred_analysis.get("files", []):
            f = self._server.get_json(file_id)
            for derived_from in f.get("derived_from", []):
                analyses = [x["@id"] for x in f.get("analyses", [])]
                if derived_from in our_files:
                    derived_map[f["@id"]] = derived_from

        for derived_file in derived_map:
            file_replicate_map[derived_file] = find_source_replicate(
                derived_map, derived_file, file_replicate_map
            )

        self._replicate_file_map = {}
        for file_id in file_replicate_map:
            self._replicate_file_map.setdefault(file_replicate_map[file_id], []).append(
                file_id
            )

        return self._replicate_file_map

    @property
    def replicates(self):
        if self._json is None:
            return

        if len(self._replicates) != len(self._json["replicates"]):
            for replicate in self._json["replicates"]:
                self._replicates.append(EncodeReplicate(replicate, self))

        return self._replicates

    @property
    def preferred_analysis(self):
        if self._preferred_analysis is None:
            default_analysis_id = self._json.get("default_analysis")
            if default_analysis_id is None:
                self._preferred_analysis = {}
            else:
                for analysis in self._json["analyses"]:
                    if analysis["@id"] == default_analysis_id:
                        self._preferred_analysis = analysis

        return self._preferred_analysis

    def __getattr__(self, key):
        if self._json is None:
            LOGGER.warn("Uninitialized Experiment object")
            return
        else:
            return self._json[key]

    def __getitem__(self, key):
        if key == "replicates":
            return self.replicates
        return self._json.get(key)

    def __iter__(self):
        for key in self._json:
            yield key

    def __len__(self):
        return len(self._json)

    def __repr__(self):
        return "<EncodeExperiment: {}>".format(self._json["@id"])


class EncodeReplicate(Mapping):
    def __init__(self, replicate, experiment):
        self._experiment = experiment
        obj_type = get_object_type(replicate)
        if obj_type != "Replicate":
            raise ValueError("Not a replicate type: {}".format(obj_type))
        self._json = replicate
        self._files = []

    @property
    def files(self):
        if self._json["@id"] not in self._experiment._replicate_file_map:
            print("id not in", self._json["@id"])
            return []

        file_ids = self._experiment._replicate_file_map[self._json["@id"]]
        known_files = [x["@id"] for x in self._files]
        for file_id in file_ids:
            if file_id not in known_files:
                f = self._experiment._server.get_json(file_id)
                LOGGER.debug("File", f["@id"], "dataset?", f["dataset"], "matches experiment", self._experiment["accession"])
                self._files.append(EncodeFile(f, self))

        return self._files

    def __getattr__(self, key):
        return self._json.get(key)

    def __getitem__(self, key):
        if key == "files":
            return self.files
        else:
            return self._json.get(key)

    def __iter__(self):
        for key in self._json:
            yield key
        # also report our extra attributes
        for key in ["files"]:
            yield key

    def __len__(self):
        return len(self._json)

    def __repr__(self):
        return "<EncodeReplicate: {}>".format(self._json["@id"])


class EncodeFile(Mapping):
    def __init__(self, json, replicate):
        self._json = json
        self._replicate = replicate
        self._response = None

    @property
    def content(self):
        if self._response is None:
            server = self._replicate._experiment._server
            url = server.prepare_url(self.href)
            self._response = requests.get(url)

        return self._requests.content

    @property
    def quality_metrics(self):
        if "quality_metrics" not in self._json:
            return []

        for metric in self._json["quality_metrics"]:
            metric_type = get_object_type(metric)
            if metric_type not in QUALITY_METRIC_PARSERS:
                raise RuntimeError("Need to implement parser for {}".format(metric_type))

            parser = QUALITY_METRIC_PARSERS[metric_type]
            value = parser(metric)
            yield value

    def __getattr__(self, key):
        return self._json.get(key)

    def __getitem__(self, key):
        if key == "quality_metrics":
            return self.quality_metrics
        else:
            return self._json.get(key)

    def __iter__(self):
        for key in self._json:
            yield key

    def __len__(self):
        return len(self._json)

    def __repr__(self):
        return "<EncodeFile: {} {}>".format(
            self._json["@id"], self._json["output_type"]
        )


def parse_pct(value):
    if value.endswith("%"):
        value = value[:-1]
    return float(value)


def parse_metric(record, column_parser):
    results = {}
    for name, value in record.items():
        if name in column_parser:
            value = column_parser[name](value)
        results[name] = value

    return results


def parse_samtools_stats(record):
    """Parse SamtoolsFlagstatsQualityMetric file"""
    columns = {
        "quality_metric_of": list,
        "duplicates": int,
        "duplicates_qc_failed": int,
        "mapped": int,
        "mapped_pct": parse_pct,
        "mapped_qc_failed": int,
        "total": int,
        "total_qc_failed": int,
    }
    return parse_metric(record, columns)


def parse_star_stats(record):
    """Pase StarQualityMetric results"""
    columns = {
        "% of chimeric reads": parse_pct,
        "% of reads mapped to multiple loci": parse_pct,
        "% of reads mapped to too many loci": parse_pct,
        "% of reads unmapped: other": parse_pct,
        "% of reads unmapped: too many mismatches": parse_pct,
        "% of reads unmapped: too short": parse_pct,
        "Average input read length": float,
        "Average mapped length": float,
        "Deletion average length": float,
        "Deletion rate per base": parse_pct,
        "Insertion average length": float,
        "Insertion rate per base": parse_pct,
        "Mapping speed, Million of reads per hour": float,
        "Mismatch rate per base, %": parse_pct,
        "Number of chimeric reads": int,
        "Number of input reads": int,
        "Number of reads mapped to multiple loci": int,
        "Number of reads mapped to too many loci": int,
        "Number of splices: AT/AC": int,
        "Number of splices: Annotated (sjdb)": int,
        "Number of splices: GC/AG": int,
        "Number of splices: GT/AG": int,
        "Number of splices: Non-canonical": int,
        "Number of splices: Total": int,
        "Uniquely mapped reads %": parse_pct,
        "Uniquely mapped reads number": int,
    }
    return parse_metric(record, columns)


def parse_gene_type_quantification(record):
    """Parse GeneTypeQuantificationQualityMetric results"""
    columns = {
        "Mt_rRNA": int,
        "antisense": int,
        "miRNA": int,
        "processed_transcript": int,
        "protein_coding": int,
        "rRNA": int,
        "ribozyme": int,
        "sRNA": int,
        "scaRNA": int,
        "sense_intronic": int,
        "sense_overlapping": int,
        "snRNA": int,
        "snoRNA": int,
        "spikein": int,
    }
    return parse_metric(record, columns)


def parse_mad_metric(record):
    """Parse MadQualityMetric results"""

    def get_href(record):
        return record["href"]

    columns = {
        "attachment": get_href,
        "SD of log ratios": float,
        "Pearson correlation": float,
        "Spearman correlation": float,
        "MAD of log ratios": float,
    }
    return parse_metric(record, columns)


def parse_genes_detected(record):
    """Parse GeneQuantificationQualityMetric results"""
    columns = {
        "quality_metric_of": list,
        "number_of_genes_detected": int,
    }
    return parse_metric(record, columns)

def parse_star_solo_quality_metric(record):
    # mode controls some attributes
    columns = {
        'barcode_rank_plot': str,
        'estimated_number_of_cells': int,
        'fraction_of_unique_reads_in_cells': float,
        'mean_UMI_per_cell': int,
        'mean_genefull_ex50pas_per_cell': int,
        'mean_reads_per_cell': int,
        'median_UMI_per_cell': int,
        'median_genefull_ex50pas_per_cell': int,
        'median_reads_per_cell': int,
        'number_of_reads': int,
        'q30_bases_in_CB_UMI': float,
        'q30_bases_in_rna_read': float,
        'reads_mapped_to_genefull_ex50pas_unique_and_multiple_gene_ex50pas': float,
        'reads_mapped_to_genefull_ex50pas_unique_genefull_ex50pas': float,
        'reads_mapped_to_genome_unique': float,
        'reads_mapped_to_genome_unique_and_multiple': float,
        'reads_with_valid_barcodes': float,
        'sequencing_saturation': float,
        'total_genefull_ex50pas_detected': int,
        'umis_in_cells': int,
        'unique_reads_in_cells_mapped_to_genefull_ex50pas': int,
    }
    return parse_metric(record, columns)


def parse_scrna_seq_counts_summary(record):
    columns = {
        'counts_violin_plot': str,
        'total_counts_vs_genes_by_count': str,
        'total_counts_vs_pct_mitochondria': str,
    }
    return parse_metric(record, columns)


QUALITY_METRIC_PARSERS = {
    "SamtoolsFlagstatsQualityMetric": parse_samtools_stats,
    "StarQualityMetric": parse_star_stats,
    "GeneTypeQuantificationQualityMetric": parse_gene_type_quantification,
    "MadQualityMetric": parse_mad_metric,
    "GeneQuantificationQualityMetric": parse_genes_detected,
    "StarSoloQualityMetric": parse_star_solo_quality_metric,
    "ScrnaSeqCountsSummaryQualityMetric": parse_scrna_seq_counts_summary,
}


def make_attachment(local_filename, mime_type=None, remote_filename=None):
    local_filename = Path(local_filename)
    if mime_type is None:
        mime_type = document_mime_type_default[local_filename.suffix]
        if mime_type is None:
            raise ValueError("Unrecognized filename extension")

    with open(local_filename, "rb") as instream:
        document = instream.read()

    if remote_filename is None:
        remote_filename = local_filename.name

    payload = {
        'download': quote_plus(str(remote_filename)),
        'type': mime_type,
        'href': 'data:{};base64,'.format(mime_type) +
                base64.b64encode(document).decode("ascii"),
        'md5sum': hashlib.md5(document).hexdigest()
    }
    return payload

def filter_aws_credentials(value, verbose=False):
    if verbose:
        return value

    value = value.copy()
    replacements = {
        "access_key": "*** AccessKey ***",
        "secret_key": "*** SecretKey ***",
        "session_token": "*** SessionToken ***",
    }
    for k in value:
        print(k)
        if k in replacements:
            value[k] = replacements[k]
            print("  ", value[k])
        elif isinstance(value[k], dict):
            value[k] = filter_aws_credentials(value[k])
        elif isinstance(value[k], list):
            newlist = []
            for x in value[k]:
                if isinstance(x, dict):
                    newlist.append(filter_aws_credentials(x, verbose=verbose))
                else:
                    newlist.append(x)
            value[k] = newlist

    print(value)
    return value


if __name__ == "__main__":
    # try it
    from rdflib import Graph

    model = Graph()
    logging.basicConfig(level=logging.DEBUG)
    encoded = ENCODED("test.encodedcc.org")
    body = encoded.get_jsonld("/experiments/ENCSR000AEC/")
    model.parse(body)
    print(model.serialize(format="turtle"))

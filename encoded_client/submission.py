"""Common submission elements
"""
import logging
import os
import re

from urllib.error import HTTPError, URLError

from rdflib import Graph, Literal, Namespace, URIRef
from rdflib.namespace import RDF

from .rdfns import (
    dafTermOntology,
    libraryOntology,
    submissionLog,
    submissionOntology,
)
from .hashfile import make_md5sum
from .fastqname import FastqName
from .daf import MetadataLookupException, ModelException, get_submission_uri

# from htsworkflow.util import opener
from xopen import xopen

LOGGER = logging.getLogger(__name__)


class Submission(object):
    def __init__(self, name, model, host):
        self.name = name
        self.model = model

        self.submissionSet = get_submission_uri(self.name)
        self.submissionSetNS = Namespace(str(self.submissionSet) + "#")
        self.libraryNS = Namespace("{0}/library/".format(host))
        self.flowcellNS = Namespace("{0}/flowcell/".format(host))

        self.__view_map = None

    def scan_submission_dirs(self, result_map):
        """Examine files in our result directory"""
        for lib_id, result_dir in result_map.items():
            LOGGER.info("Importing %s from %s" % (lib_id, result_dir))
            try:
                self.import_analysis_dir(result_dir, lib_id)
            except MetadataLookupException as e:
                LOGGER.error("Skipping %s: %s" % (lib_id, str(e)))

    def import_analysis_dir(self, analysis_dir, library_id):
        """Import a submission directories and update our model as needed"""
        # attributes = get_filename_attribute_map(paired)
        libNode = self.libraryNS[library_id + "/"]

        self._add_library_details_to_model(libNode)

        submission_files = os.listdir(analysis_dir)
        for filename in submission_files:
            pathname = os.path.abspath(os.path.join(analysis_dir, filename))
            self.construct_file_attributes(analysis_dir, libNode, pathname)

    def analysis_nodes(self, result_map):
        """Return an iterable of analysis nodes"""
        for result_dir in result_map.values():
            an_analysis = self.get_submission_node(result_dir)
            yield an_analysis

    def construct_file_attributes(self, analysis_dir, libNode, pathname):
        """Looking for the best extension
        The 'best' is the longest match

        :Args:
        filename (str): the filename whose extention we are about to examine
        """
        path, filename = os.path.split(pathname)

        LOGGER.debug("Searching for view")
        file_type = self.find_best_match(filename)
        if file_type is None:
            LOGGER.warning("Unrecognized file: {0}".format(pathname))
            return None
        if str(file_type) == str(libraryOntology["ignore"]):
            return None

        an_analysis_name = self.make_submission_name(analysis_dir)
        an_analysis = self.get_submission_node(analysis_dir)
        file_classifications = list(self.model.objects(file_type, RDF["type"]))
        if len(file_classifications) == 0:
            errmsg = "Could not find class for {0}"
            LOGGER.warning(errmsg.format(str(file_type)))
            return
        file_classification = file_classifications[0]

        self.model.add(
            (
                self.submissionSetNS[""],
                submissionOntology["has_submission"],
                an_analysis,
            )
        )
        self.model.add(
            (an_analysis, submissionOntology["name"], Literal(an_analysis_name))
        )
        self.model.add((an_analysis, RDF["type"], submissionOntology["submission"]))
        self.model.add((an_analysis, submissionOntology["library"], libNode))

        LOGGER.debug("Adding statements to {0}".format(str(an_analysis)))
        # add track specific information
        self.model.add(
            (an_analysis, dafTermOntology["paired"], Literal(self._is_paired(libNode)))
        )
        self.model.add((an_analysis, dafTermOntology["submission"], an_analysis))

        # add file specific information
        fileNode = self.make_file_node(pathname, an_analysis)
        self.add_md5s(filename, fileNode, analysis_dir)
        self.add_file_size(filename, fileNode, analysis_dir)
        self.add_read_length(filename, fileNode, analysis_dir)
        self.add_fastq_metadata(filename, fileNode)
        self.add_label(file_type, fileNode, libNode)
        self.model.add((fileNode, RDF["type"], file_type))
        self.model.add((fileNode, libraryOntology["library"], libNode))

        LOGGER.debug("Done.")

    def make_file_node(self, pathname, submissionNode):
        """Create file node and attach it to its submission."""
        # add file specific information
        path, filename = os.path.split(pathname)
        pathname = os.path.abspath(pathname)
        fileNode = URIRef("file://" + pathname)
        self.model.add((submissionNode, dafTermOntology["has_file"], fileNode))
        self.model.add((fileNode, dafTermOntology["filename"], Literal(filename)))
        self.model.add(
            (
                fileNode,
                dafTermOntology["relative_path"],
                Literal(os.path.relpath(pathname)),
            )
        )
        return fileNode

    def add_md5s(self, filename, fileNode, analysis_dir):
        LOGGER.debug("Updating file md5sum")
        submission_pathname = os.path.join(analysis_dir, filename)
        md5 = make_md5sum(submission_pathname)
        if md5 is None:
            errmsg = "Unable to produce md5sum for {0}"
            LOGGER.warning(errmsg.format(submission_pathname))
        else:
            self.model.add((fileNode, dafTermOntology["md5sum"], Literal(md5)))

    def add_file_size(self, filename, fileNode, analysis_dir):
        submission_pathname = os.path.join(analysis_dir, filename)
        file_size = os.stat(submission_pathname).st_size
        self.model.add((fileNode, dafTermOntology["file_size"], Literal(file_size)))
        LOGGER.debug("Updating file size: %d", file_size)

    def add_read_length(self, filename, fileNode, analysis_dir):
        submission_pathname = os.path.join(analysis_dir, filename)
        stream = xopen(submission_pathname, "rt")
        header = stream.readline().strip()
        sequence = stream.readline().strip()
        stream.close()
        read_length = len(sequence)
        self.model.add((fileNode, libraryOntology["read_length"], Literal(read_length)))
        LOGGER.debug("Updating read length: %d", read_length)

    def add_fastq_metadata(self, filename, fileNode):
        # How should I detect if this is actually a fastq file?
        try:
            fqname = FastqName(filename=filename)
        except ValueError:
            # currently its just ignore it if the fastq name parser fails
            return

        terms = [
            ("flowcell", libraryOntology["flowcell_id"]),
            ("lib_id", libraryOntology["library_id"]),
            ("lane", libraryOntology["lane_number"]),
            ("read", libraryOntology["read"]),
        ]
        for file_term, model_term in terms:
            value = fqname.get(file_term)
            if value is not None:
                s = (fileNode, model_term, value.toPython())
                self.model.add(s)

        if "flowcell" in fqname:
            value = self.flowcellNS[fqname["flowcell"] + "/"]
            s = (fileNode, libraryOntology["flowcell"], value)
            self.model.add(s)

    def add_label(self, file_type, file_node, lib_node):
        """Add rdfs:label to a file node"""
        from django.template import Template, loader

        # template_term = libraryOntology['label_template']
        template_term = libraryOntology["label_template"]
        label_templates = list(self.model.objects(file_type, template_term))
        if len(label_templates) > 0:
            label_template = label_templates[0]
            template = loader.get_template("submission_view_rdfs_label_metadata.sparql")
            context = {
                "library": str(lib_node.uri),
            }
            for row in self.execute_query(template, context):
                label = Template(label_template).render(row)
                s = (file_node, rdfsNS["label"], unicode(label))
                self.model.add(s)

    def _add_library_details_to_model(self, libNode):
        # attributes that can have multiple values
        set_attributes = set(
            (
                libraryOntology["has_lane"],
                libraryOntology["has_mappings"],
                dafTermOntology["has_file"],
            )
        )
        tmpmodel = Graph()
        try:
            tmpmodel.parse(source=libNode, format="rdfa")
        except (HTTPError, URLError) as e:
            LOGGER.error(e)
            return

        LOGGER.debug("Scanning %s", str(libNode))
        toadd = []
        for stmt in tmpmodel:
            s, p, o = stmt
            # always add "collections"
            if p in set_attributes:
                toadd.append(stmt)
                continue
            # don't override things we already have in the model
            targets = list(self.model.objects(s, p))
            if len(targets) == 0:
                toadd.append(stmt)

        for stmt in toadd:
            self.model.add(stmt)

        self._add_lane_details(libNode)
        self._add_flowcell_details()

    def _add_lane_details(self, libNode):
        """Import lane details"""
        query = (libNode, libraryOntology["has_lane"], None)
        lanes = []
        for lane_stmt in self.model.triples(query):
            lanes.append(lane_stmt[2])

        for lane in lanes:
            LOGGER.debug("Importing %s" % (lane,))
            self.model.parse(source=lane, format="rdfa")

    def _add_flowcell_details(self):
        template = loader.get_template("aws_flowcell.sparql")

        for r in self.execute_query(template, {}):
            flowcell = r["flowcell"]
            self.model.parse(source=flowcell, format="rdfa")

    def find_best_match(self, filename):
        """Search through potential filename matching patterns"""
        if self.__view_map is None:
            self.__view_map = self._get_filename_view_map()

        results = []
        for pattern, view in self.__view_map.items():
            if re.match(pattern, filename):
                results.append(view)

        if len(results) > 1:
            msg = "%s matched multiple views %s" % (filename, [str(x) for x in results])
            raise ModelException(msg)
        elif len(results) == 1:
            return results[0]
        else:
            return None

    def _get_filename_view_map(self):
        """Query our model for filename patterns

        return a dictionary of compiled regular expressions to view names
        """
        filename_query = (None, dafTermOntology["filename_re"], None)

        patterns = {}
        for s in self.model.triples(filename_query):
            view_name = s[0]
            literal_re = s[2].value
            LOGGER.debug("Found: %s" % (literal_re,))
            try:
                filename_re = re.compile(literal_re)
            except re.error as e:
                LOGGER.error("Unable to compile: %s" % (literal_re,))
            patterns[literal_re] = view_name
        return patterns

    def make_submission_name(self, analysis_dir):
        analysis_dir = os.path.normpath(analysis_dir)
        analysis_dir_name = os.path.split(analysis_dir)[1]
        if len(analysis_dir_name) == 0:
            raise RuntimeError(
                "Submission dir name too short: {0}".format(analysis_dir)
            )
        return analysis_dir_name

    def get_submission_node(self, analysis_dir):
        """Convert a submission directory name to a submission node"""
        submission_name = self.make_submission_name(analysis_dir)
        return self.submissionSetNS[submission_name]

    def _get_library_attribute(self, libNode, attribute):
        if not isinstance(libNode, URIRef):
            raise ValueError("libNode must be a URIRef")
        if not isinstance(attribute, URIRef):
            attribute = libraryOntology[attribute]

        targets = list(self.model.objects(libNode, attribute))
        if len(targets) > 0:
            return self._format_library_attribute(targets)
        else:
            return None

        # targets = self._search_same_as(libNode, attribute)
        # if targets is not None:
        #    return self._format_library_attribute(targets)

        # we don't know anything about this attribute
        self._add_library_details_to_model(libNode)

        targets = list(self.model.objects(libNode, attribute))
        if len(targets) > 0:
            return self._format_library_attribute(targets)

        return None

    def _format_library_attribute(self, targets):
        if len(targets) == 0:
            return None
        elif len(targets) == 1:
            return targets[0].toPython()
        elif len(targets) > 1:
            return [t.toPython() for t in targets]

    def _is_paired(self, libNode):
        """Determine if a library is paired end"""
        library_type = self._get_library_attribute(libNode, "library_type")
        if library_type is None:
            errmsg = "%s doesn't have a library type"
            raise ModelException(errmsg % (str(libNode),))

        single = [
            "CSHL (lacking last nt)",
            "Single End (non-multiplexed)",
            "Small RNA (non-multiplexed)",
        ]
        paired = [
            "Barcoded Illumina",
            "Multiplexing",
            "NEBNext Multiplexed",
            "NEBNext Small RNA",
            "Nextera",
            "Paired End (non-multiplexed)",
            "Dual Index Illumina",
        ]
        if library_type in single:
            return False
        elif library_type in paired:
            return True
        else:
            raise MetadataLookupException(
                "Unrecognized library type %s for %s" % (library_type, str(libNode))
            )

    def execute_query(self, template, context):
        """Execute the query, returning the results"""
        formatted_query = template.render(context)
        LOGGER.debug(formatted_query)
        rdfstream = self.model.query(str(formatted_query))
        results = []
        for record in rdfstream:
            d = {}
            for key, value in record.asdict().items():
                d[key] = value.toPython()
            results.append(d)
        return results


def list_submissions(model):
    """Return generator of submissions in this model."""
    query_body = """
      PREFIX subns: <http://jumpgate.caltech.edu/wiki/UcscSubmissionOntology#>

      select distinct ?submission
      where { ?submission subns:has_submission ?library_dir }
    """
    from .rdfhelp import strip_namespace

    q = (None, submissionOntology["has_submission"], None)
    submissions = set()
    for statement in model.triples(q):
        s = strip_namespace(submissionLog, statement[0])
        if s[-1] in ["#", "/", "?"]:
            s = s[:-1]
        submissions.add(s)

    return list(submissions)

"""Help collect and process results for submission
"""
from collections.abc import MutableMapping
import os
import shutil
import logging

from collections import namedtuple

LOGGER = logging.getLogger(__name__)


class ResultMap(MutableMapping):
    """Store list of results"""

    def __init__(self):
        self.results_order = []
        self.results = {}

    def __iter__(self):
        for item in self.results_order:
            yield item

    def __len__(self):
        l = len(self.results)
        assert l == len(self.results_order)
        return l

    def __setitem__(self, key, value):
        self.results_order.append(key)
        self.results[key] = value

    def __getitem__(self, key):
        return self.results[key]

    def __delitem__(self, key):
        del self.results[key]
        i = self.results_order.index(key)
        del self.results_order[i]

    def add_results_from_file(self, filename):
        pathname = os.path.abspath(filename)
        basepath, name = os.path.split(pathname)
        results = read_result_list(filename)
        for lib_id, lib_path in results:
            if not os.path.isabs(lib_path):
                lib_path = os.path.join(basepath, lib_path)
            self[lib_id] = lib_path

    def make_tree_from(self, source_path, destpath=None, link=True):
        """Create a tree using data files from source path."""
        if destpath is None:
            destpath = os.getcwd()

        LOGGER.debug("Source_path: %s", source_path)
        LOGGER.debug("Dest_path: %s", destpath)
        for lib_id in self.results_order:
            lib_path = self.results[lib_id]
            LOGGER.debug("lib_path: %s", lib_path)
            if os.path.isabs(lib_path):
                lib_path = os.path.relpath(lib_path, destpath)

            LOGGER.debug("lib_path: %s", lib_path)
            lib_destination = os.path.join(destpath, lib_path)
            if not os.path.exists(lib_destination):
                LOGGER.info("Making dir {0}".format(lib_destination))
                os.mkdir(lib_destination)

            source_rel_dir = os.path.join(source_path, lib_path)
            source_lib_dir = os.path.abspath(source_rel_dir)
            LOGGER.debug("source_lib_dir: %s", source_lib_dir)

            for filename in os.listdir(source_lib_dir):
                source_pathname = os.path.join(source_lib_dir, filename)
                target_pathname = os.path.join(lib_destination, filename)
                if not os.path.exists(source_pathname):
                    raise IOError("{0} does not exist".format(source_pathname))
                if not (
                    os.path.exists(target_pathname) or os.path.isdir(source_pathname)
                ):
                    if link:
                        os.symlink(source_pathname, target_pathname)
                    else:
                        shutil.copy(source_pathname, target_pathname)
                    LOGGER.info(
                        "LINK {0} to {1}".format(source_pathname, target_pathname)
                    )


def read_result_list(filename):
    """
    Read a file that maps library id to result directory.
    Does not support spaces in filenames.

    For example:
      10000 result/foo/bar
    """
    stream = open(filename, "r")
    results = parse_result_list(stream)
    stream.close()
    return results


def parse_result_list(stream):
    results = []
    for line in stream:
        line = line.rstrip()
        if not line.startswith("#") and len(line) > 0:
            library_id, result_dir = line.split()
            results.append((library_id, result_dir))
    return results

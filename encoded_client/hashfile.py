"""Utility to make md5sums of a file caching as a parallel file
"""
import logging
import os
from subprocess import Popen, PIPE

logger = logging.getLogger(__name__)


def make_md5sum(filename):
    """Quickly find the md5sum of a file"""
    md5_cache = os.path.join(filename + ".md5")
    if os.path.exists(md5_cache):
        logger.debug("Found md5sum in {0}".format(md5_cache))
        stream = open(md5_cache, "rt")
        lines = stream.readlines()
        md5sum = parse_md5sum_line(lines, filename)
    else:
        md5sum = make_md5sum_unix(filename, md5_cache)
    return md5sum


def make_md5sum_unix(filename, md5_cache):
    cmd = ["md5sum", filename]
    logger.debug("Running {0}".format(" ".join(cmd)))
    p = Popen(cmd, stdout=PIPE)
    stdin, stdout = p.communicate()
    retcode = p.wait()
    logger.debug("Finished {0} retcode {1}".format(" ".join(cmd), retcode))
    if retcode != 0:
        logger.error("Trouble with md5sum for {0}".format(filename))
        return None
    lines = [x.decode("utf-8") for x in stdin.splitlines()]
    md5sum = parse_md5sum_line(lines, filename)
    if md5sum is not None:
        logger.debug("Caching sum in {0}".format(md5_cache))
        stream = open(md5_cache, "wt")
        stream.write(stdin.decode("utf-8"))
        stream.close()
    return md5sum


def parse_md5sum_line(lines, filename):
    md5sum, md5sum_filename = lines[0].split()
    md5sum_filename = os.path.basename(md5sum_filename)
    filename = os.path.basename(filename)
    if md5sum_filename != filename:
        errmsg = "MD5sum and I disagre about filename. {0} != {1}"
        logger.error(errmsg.format(filename, md5sum_filename))
        return None
    return md5sum


def main(cmdline=None):
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "-v", "--verbose", action="store_true", default="False", help="verbose logging"
    )
    parser.add_argument("filenames", nargs="+", help="filenames to build md5 cache for")
    args = parser.parse_args(cmdline)

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    for filename in args.filenames:
        make_md5sum(filename)


if __name__ == "__main__":
    main()

"""Standardize reading and writing fastq submission names.
"""
import collections
import re

PAIRED_TEMPLATE = (
    "{lib_id}_{flowcell}_c{cycle}_l{lane}_r{read}.fastq{compression_extension}"
)
SINGLE_TEMPLATE = "{lib_id}_{flowcell}_c{cycle}_l{lane}.fastq{compression_extension}"

FASTQ_RE = re.compile(
    "(?P<lib_id>[^_]+)_(?P<flowcell>[^_]+)_"
    "c(?P<cycle>[\d]+)_l(?P<lane>[\d]+)(_r(?P<read>[\d]))?\.fastq(?P<compression_extension>.[\w]+)?"
)


class FastqName(collections.abc.Mapping):
    """Utility class to convert to the standardized submission fastq name."""

    def __init__(self, is_paired=None, **kwargs):
        """Create a fastq name handler.

        Takes filename or common attributes like flowcell, lib_id, lane, read, cycle
        """
        self._attributes = (
            "flowcell",
            "lib_id",
            "lane",
            "read",
            "cycle",
            "compression_extension",
        )
        self._is_paired = is_paired

        if len(kwargs) == 0:
            return
        if "filename" in kwargs:
            self._init_by_filename(**kwargs)
        else:
            self._init_by_attributes(**kwargs)

    def _init_by_attributes(self, **kwargs):
        for k in self._attributes:
            value = None
            if k in kwargs:
                value = kwargs[k]
            self[k] = value

    def _init_by_filename(self, filename):
        match = FASTQ_RE.match(filename)
        if match is None:
            raise ValueError('Is "{0}" a submission fastq?'.format(filename))

        for k in self._attributes:
            self[k] = match.group(k)

    def _get_is_paired(self):
        if self._is_paired is None:
            return getattr(self, "read", None) is not None
        else:
            return self._is_paired

    def _set_is_paired(self, value):
        self._is_paired = value

    is_paired = property(_get_is_paired, _set_is_paired)

    def _is_valid(self):
        if self.is_paired and self["read"] is None:
            return False

        for k in self.keys():
            if k == "read":
                continue
            elif k == "compression_extension":
                if self[k] not in (None, "", ".gz", ".bz2"):
                    return False
            elif getattr(self, k, None) is None:
                return False
        return True

    is_valid = property(_is_valid)

    def _get_filename(self):
        if not self.is_valid:
            raise ValueError(
                "Please set all needed variables before generating a filename"
            )

        T = PAIRED_TEMPLATE if self.is_paired else SINGLE_TEMPLATE
        attributes = {}
        for k in self:
            v = self[k]
            attributes[k] = v if v is not None else ""
        return T.format(**attributes)

    filename = property(_get_filename)

    def __iter__(self):
        return iter(self._attributes)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        if key in self._attributes:
            setattr(self, key, value)
        else:
            raise ValueError("Unrecognized key {0}".format(key))

    def __len__(self):
        return len([k for k in self if self[k] is not None])

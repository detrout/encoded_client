#!/usr/bin/env python
from __future__ import absolute_import

import os
import shutil

from unittest import TestCase, defaultTestLoader

from ..results import ResultMap
from .submission_test_common import (
    generate_sample_results_tree,
    S1_NAME,
    S2_NAME,
    S1_FILES,
    S2_FILES,
)


class TestResultMap(TestCase):
    def setUp(self):
        generate_sample_results_tree(self, "results_test")

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_dict_like(self):
        """Make sure the result map works like an ordered dictionary"""
        results = ResultMap()
        results["1000"] = "dir1000"
        results["2000"] = "dir2000"
        results["1500"] = "dir1500"

        self.assertEqual(list(results.keys()), ["1000", "2000", "1500"])
        self.assertEqual(list(results.values()), ["dir1000", "dir2000", "dir1500"])
        self.assertEqual(
            list(results.items()),
            [("1000", "dir1000"), ("2000", "dir2000"), ("1500", "dir1500")],
        )

        self.assertEqual(results["1000"], "dir1000")
        self.assertEqual(results["1500"], "dir1500")
        self.assertEqual(results["2000"], "dir2000")

        self.assertTrue(u"2000" in results)
        self.assertTrue("2000" in results)
        self.assertFalse(u"77777" in results)
        self.assertFalse("77777" in results)

    def test_make_from_absolute(self):
        """Test that make from works if ResultMap has absolute paths"""
        results = ResultMap()
        sample1_dir = os.path.join(self.resultdir, S1_NAME)
        sample2_dir = os.path.join(self.resultdir, S2_NAME)
        results["1000"] = sample1_dir
        results["2000"] = sample2_dir

        results.make_tree_from(self.sourcedir, self.resultdir)
        self.assertTrue(os.path.isdir(sample1_dir))
        self.assertTrue(os.path.isdir(sample2_dir))

        for f in S1_FILES + S2_FILES:
            self.assertTrue(os.path.islink(os.path.join(self.resultdir, f)))

    def test_make_from_filename(self):
        """Test that make from works if ResultMap has no path"""
        results = ResultMap()
        results["1000"] = S1_NAME
        results["2000"] = S2_NAME

        results.make_tree_from(self.sourcedir, self.resultdir)
        sample1_dir = os.path.join(self.resultdir, S1_NAME)
        sample2_dir = os.path.join(self.resultdir, S2_NAME)
        self.assertTrue(os.path.isdir(sample1_dir))
        self.assertTrue(os.path.isdir(sample2_dir))

        for f in S1_FILES + S2_FILES:
            self.assertTrue(os.path.islink(os.path.join(self.resultdir, f)))

    def test_make_from_shared_directory(self):
        """Split multiple datasets stored in a single directory"""
        self.skipTest("not implemented yet")
        results = ResultMap()
        results["S1"] = os.path.join(SCOMBINED_NAME, "s1*")
        results["S2"] = os.path.join(SCOMBINED_NAME, "s2*")


def suite():
    suite = defaultTestLoader.loadTestsFromTestCase(TestResultMap)
    return suite


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.DEBUG)
    from unittest import main

    main(defaultTest="suite")

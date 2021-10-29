from unittest import TestCase
import tempfile


from ..hashfile import (
    make_md5sum,
)


class TestHashfile(TestCase):
    def test_make_md5sum(self):
        with tempfile.NamedTemporaryFile(suffix="_hashfile") as tmp:
            tmp.write("hello".encode("utf-8"))
            sum1 = make_md5sum(tmp.name)

            sum2 = make_md5sum(tmp.name)

            self.assertEqual(sum1, sum2)
            self.assertEqual(sum1, "d41d8cd98f00b204e9800998ecf8427e")

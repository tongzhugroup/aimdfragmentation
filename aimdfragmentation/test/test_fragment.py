import os
import unittest
import logging

import pkg_resources

from aimdfragmentation import AIMDFragmentation


class Test_all(unittest.TestCase):
    def test_fragment(self):
        xyzfilename = 'test_fragment.xyz'
        with open(xyzfilename, 'w') as f:
            f.write(pkg_resources.resource_string(
                __name__, 'test.xyz').decode())
        af = AIMDFragmentation(cutoff=6.0, xyzfilename=xyzfilename)
        af.run()
        for filename in (af.outputenergyfile, af.outputfile):
            with open(filename) as f:
                logging.info(filename)
                print(f.read())


if __name__ == '__main__':
    unittest.main()

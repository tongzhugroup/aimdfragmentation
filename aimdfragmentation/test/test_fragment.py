import os
import unittest
import logging

import pkg_resources

from aimdfragmentation import AIMDFragmentation


class Test_all(unittest.TestCase):
    def test_fragment(self):
        folder = "testfiles"
        xyzfilename = 'test_fragment.xyz'
        if not os.path.exists(folder):
            os.makedirs(folder)
        with open(os.path.join(folder, xyzfilename), 'wb') as f:
            f.write(pkg_resources.resource_string(__name__, xyzfilename))
        af = AIMDFragmentation(cutoff=6.0, xyzfilename=xyzfilename)
        af.run()
        with open(af.outputenergyfile) as f:
            logging.info("Energy:")
            print(f.read())
        with open(af.outputenergyfile) as f:
            logging.info("Force:")
            print(f.read())


if __name__ == '__main__':
    unittest.main()

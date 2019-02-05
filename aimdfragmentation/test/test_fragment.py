import os
import logging
import tempfile

import pkg_resources

from aimdfragmentation import AIMDFragmentation


class Test_all:
    def test_fragment(self):
        folder = tempfile.mkdtemp(prefix='testfiles-', dir='.')
        logging.info(f'Folder: {folder}:')
        os.chdir(folder)
        xyzfilename = 'test.xyz'
        with open(xyzfilename, 'w') as f:
            f.write(pkg_resources.resource_string(
                __name__, xyzfilename).decode())
        af = AIMDFragmentation(cutoff=6.0, xyzfilename=xyzfilename)
        af.run()
        for filename in (af.outputenergyfile, af.outputfile):
            with open(filename) as f:
                logging.info(filename)
                print(f.read())

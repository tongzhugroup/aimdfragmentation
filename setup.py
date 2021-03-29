"""Use 'pip install .' to install."""
from os import path

from setuptools import setup, find_packages, Extension

if __name__ == '__main__':
    print(__doc__)
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'docs', 'README.md')) as f:
        long_description = f.read()

    tests_require = ['pytest-sugar', 'pytest-cov']
    setup(name='aimdfragmentation',
          description='AIMD Fragmentation Calculation',
          keywords="AIMD Fragmentation",
          url='https://github.com/njzjz/aimdfragmentation',
          author='Jinzhe Zeng',
          author_email='jzzeng@stu.ecnu.edu.cn',
          install_requires=[
              'numpy', 'gaussianrunner>=1.0.19', 'ase', 'coloredlogs',
              'openbabel-wheel',
          ],
          test_suite='aimdfragmentation.test',
          tests_require=tests_require,
          extras_require={
              "test": tests_require,
          },
          package_data={'aimdfragmentation': ['test/test.xyz']},
          packages=find_packages(),
          python_requires='~=3.6',
          use_scm_version=True,
          setup_requires=[
              'setuptools_scm',
              'pytest-runner',
              "setuptools>=18.0",
              "cython",
          ],
          long_description=long_description,
          long_description_content_type='text/markdown',
          classifiers=[
              "Natural Language :: English",
              "Operating System :: POSIX :: Linux",
              "Operating System :: Microsoft :: Windows",
              "Programming Language :: Python :: 3.6",
              "Programming Language :: Python :: 3.7",
              "Topic :: Scientific/Engineering :: Chemistry",
              "Topic :: Software Development :: Libraries :: Python Modules",
              "Topic :: Software Development :: Version Control :: Git",
          ],
          zip_safe=True,
          ext_modules=[
              Extension("aimdfragmentation.dps", sources=[
                  "aimdfragmentation/dps.pyx", "aimdfragmentation/c_stack.cpp"], language="c++"),
          ],
          )

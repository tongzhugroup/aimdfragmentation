"""Use 'pip install .' to install."""
from os import path

from setuptools import setup, find_packages

if __name__ == '__main__':
    print(__doc__)
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'docs', 'README.md')) as f:
        long_description = f.read()

    setup(name='aimdfragmentation',
          description='AIMD Fragmentation Calculation',
          keywords="AIMD Fragmentation",
          url='https://github.com/njzjz/aimdfragmentation',
          author='Jinzhe Zeng',
          author_email='jzzeng@stu.ecnu.edu.cn',
          install_requires=['numpy', 'gaussianrunner', 'ase', 'coloredlogs'],
          test_suite='aimdfragmentation.test',
          tests_require=['pytest-sugar'],
          package_data={'aimdfragmentation': ['test/test.xyz']},
          packages=find_packages(),
          python_requires='~=3.6',
          use_scm_version=True,
          setup_requires=['setuptools_scm', 'pytest-runner'],
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
          )

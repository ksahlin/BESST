from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='BESST',
    version='2.2.8',
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
      ],
    scripts=['runBESST', 'scripts/reads_to_ctg_map.py'],
    description='Scaffolder for genomic assemblies.',
    author='Kristoffer Sahlin',
    author_email='kristoffer.sahlin@scilifelab.se',
    url='https://github.com/ksahlin/BESST',
    license='GPLv3',
    long_description=long_description,
    install_requires=['pysam==0.8.4',
                      'networkx<=1.9',
                      'mathstats>=0.2.6.5',
                      'scipy==1.0.0'],
)

from setuptools import setup, find_packages
#from distutils.core import setup
setup(
    name='BESST',
    version='1.0',
    packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
    #py_modules=['BESST.py'],
    #packages=['BESST',],
    scripts = ['runBESST'],
    description='Scaffolder for genomic assemblies.',
    author='Kristoffer Sahlin',
    author_email='kristoffer.sahlin@scilifelab.se',
    url='https://github.com/ksahlin/BESST',
    license='GPLv3',
    long_description=open('README.md').read(),
    #requires=['pysam (>=0.6)','networkx (>=1.4)'],
    install_requires=[ 'pysam==0.6',
                      'networkx>=1.4'],
    #platforms=['Unix', 'Linux', 'Mac OS']
)

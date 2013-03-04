from distutils.core import setup
setup(
    name='BESST',
    version='1.0dev',
    description='Scaffolder for genomic assemblies.',
    author='Kristoffer Sahlin',
    author_email='kristoffer.sahlin@scilifelab.se',
    url='https://github.com/ksahlin/BESST',
    #py_modules=['BESST.py'],
    packages=['BESST',],
    scripts = ['runBESST'],
    license='GPLv3',
    long_description=open('README.txt').read(),
    requires=['pysam (>=0.6)','networkx (>=1.4)'],
    platforms=['Unix', 'Linux', 'Mac OS']
)

'''
Created on Aug 4, 2013

@author: ksahlin
'''

import sys


def module_exists(module_name):
    try:
        __import__(module_name)
    except ImportError:
        return False
    else:
        return True



##
# Error handling for modules 

def check_module(module):
    if not module_exists(module):
        sys.stdout.write('Uninstalled module: ' + module + '\n')
        sys.stdout.write('Please install ' + module + ' using one of the two provided commands from terminal:\n \
        pip install ' + module + '\n \
        easy_install ' + module + '\n')
        sys.exit()



def parse_check(parser, options):

    ##
    # Error handling when parsing arguments

    if not all([x == None or len(x) == len(options.bamfiles) for x in [options.stddev , options.mean , options.readlen, options.edgesupport, options.covcutoff, options.threshold, options.minsize]]):
        parser.error("Same number of arguments are required")

    if (options.mean and not options.stddev) or (not options.mean and options.stddev):
        parser.error("either both or none of -m and -s is required")
    if (options.threshold and not options.minsize) or (not options.threshold and options.minsize):
        parser.error("either both or none of -T and -k is required")
    if not options.contigfile:
        parser.error("parameter -c (a fasta contig file) is required")
    if not options.bamfiles:
        parser.error("parameter -f (BAM files) is required")

    ##
    # Error handling for files

    for file_ in options.bamfiles:
        try:
            open(file_)
        except IOError as e:
            sys.exit("couldn't find BAM file: " + file_ + " check that the path is correct and that the file exists")
        try:
            open(file_ + '.bai')
        except IOError as e:
            print "couldn't find index file: ", file_ + '.bai', " check that the path is correct and that the file exists"
            sys.exit(0)

    try:
        open(options.contigfile)
    except IOError as e:
        sys.exit("couldn't open contig file " + options.contigfile + " check that the path is correct and that the file exists")

    return()

##
# Throw warning if contig present in bam file but not seen in contig fasta file.
# This can happen if fasta file is filtered after alignment

def unknown_contig(ctg):
    sys.stdout.write('Contig: ' + ctg + ' was seen in BAM file but is not present in contig fasta file.\
    If this is expected (i.e. due to filtering after alignment), please ignore.')

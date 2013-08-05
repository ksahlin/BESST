'''
Created on Aug 4, 2013

@author: ksahlin
'''

import sys
#from argparse import ArgumentError


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



def parse_check(arg, parser):

    ##
    # Error handling for files

    for file_ in arg.bamfiles:
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
        open(arg.contigfile)
    except IOError as e:
        sys.exit("couldn't open contig file " + arg.contigfile + " check that the path is correct and that the file exists")


    ##
    # Error handling when parsing arguments

    if not all([x == None or len(x) == len(arg.bamfiles) for x in [arg.stddev , arg.mean , arg.readlen, arg.edgesupport, arg.covcutoff, arg.threshold, arg.minsize]]):
        parser.error("If any of the below following options are specified, they should have the same number of arguments as the number of BAM files.\n {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7} ".format("bamfiles", "mean", "stddev", "readlen", "edgesupport",
                          "covcutoff", "threshold", "minsize"))

    if (arg.mean and not arg.stddev) or (not arg.mean and arg.stddev):
        parser.error("either both or none of -m and -s is required")
    if (arg.threshold and not arg.minsize) or (not arg.threshold and arg.minsize):
        parser.error("either both or none of -T and -k is required")
    if not arg.contigfile:
        parser.error("parameter -c (a fasta contig file) is required")
    if not arg.bamfiles:
        parser.error("parameter -f (BAM files) is required")

    return()

##
# Throw warning if contig present in bam file but not seen in contig fasta file.
# This can happen if fasta file is filtered after alignment

def unknown_contig(ctg):
    sys.stdout.write('Contig: ' + ctg + ' was seen in BAM file but is not present in contig fasta file.\
    If this is expected (i.e. due to filtering after alignment), please ignore.')

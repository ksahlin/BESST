#!/usr/bin/env python
'''
Created on Jul 25, 2013

@author: ksahlin
'''
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from time import time

import pysam

##
# Converts a sam file to a bam file using pysam.
#
# @param sam_path Path of the .sam file.
# @param bam_path Path of the resulting .bam file.
#
def sam_to_bam(sam_path, bam_path):
    sam_file = pysam.Samfile(sam_path, "r")
    bam_file = pysam.Samfile(bam_path, "wb", template=sam_file)

    for alignment in sam_file:
        bam_file.write(alignment)

##
# Maps the given paired end reads using bwa, and writes a
# sorted .bam file in the given output file.
#
# @param pe1_path Path of the first reads.
# @param pe2_path Path of the second reads.
# @param genome_path Path to the reference genome.
# @param output_path Path of the output file without extension ".bam".
#
def bwa_sampe(pe1_path, pe2_path, genome_path, output_path):
    print 'Aligning with bwa aln/sampe'
    start = time()
    work_dir = tempfile.mkdtemp()
    genome_db = os.path.join(work_dir, "genome")
    pe1_output = os.path.join(work_dir, "pe1.sai")
    pe2_output = os.path.join(work_dir, "pe2.sai")
    bwa_output = os.path.join(work_dir, "output.sam")
    stderr_file = open(output_path+'.bwa.1','w')

    #null = open("/dev/null")
    subprocess.check_call([ "bwa", "index", "-p", genome_db, genome_path ], stderr=stderr_file)
    with open(pe1_output, "w") as pe1_file:
        subprocess.check_call([ "bwa", "aln", "-t 8", genome_db, pe1_path ], stdout=pe1_file, stderr=stderr_file)

    with open(pe2_output, "w") as pe2_file:
        subprocess.check_call([ "bwa", "aln", "-t 8", genome_db, pe2_path ], stdout=pe2_file, stderr=stderr_file)

    with open(bwa_output, "w") as bwa_file:
        subprocess.check_call([ "bwa", "sampe",
                                genome_db,
                                pe1_output, pe2_output,
                                pe1_path, pe2_path ], stdout=bwa_file , stderr=stderr_file)

    elapsed = time() - start
    print 'Time elapsed for bwa aln/sampe: ', elapsed

    sam_to_bam(bwa_output, bwa_output + ".bam")
    pysam.sort(bwa_output + ".bam", output_path)
    pysam.index(output_path + '.bam')

    shutil.rmtree(work_dir)

def bwa_mem(pe1_path, pe2_path, genome_path, threads, output_path):
    print 'Aligning with bwa mem'
    start = time()
    work_dir = tempfile.mkdtemp()
    genome_db = os.path.join(work_dir, "genome")
    pe1_output = os.path.join(work_dir, "pe1.sai")
    pe2_output = os.path.join(work_dir, "pe2.sai")
    bwa_output = os.path.join(work_dir, "output.sam")
    stderr_file = open(output_path+'.bwa.1','w')

    #null = open("/dev/null")
    subprocess.check_call([ "bwa", "index", "-p", genome_db, genome_path ], stderr=stderr_file)
    with open(bwa_output, "w") as bwa_file:
        subprocess.check_call([ "bwa", "mem", "-t", threads,
                                genome_db, pe1_path, pe2_path ],
                              stdout=bwa_file,
                              stderr=stderr_file)

    elapsed = time() - start
    print 'Time elapsed for bwa mem: ', elapsed
    sam_to_bam(bwa_output, bwa_output + ".bam")
    pysam.sort(bwa_output + ".bam", output_path)
    pysam.index(output_path + '.bam')
   
    shutil.rmtree(work_dir)

def map_single_reads(pe_path, genome_path, output_path):
    work_dir = tempfile.mkdtemp()
    genome_db = os.path.join(work_dir, "genome")
    pe_output = os.path.join(work_dir, "pe.sai")
    bwa_output = os.path.join(work_dir, "output.sam")

    null = open("/dev/null")
    subprocess.check_call([ "bwa", "index", "-p", genome_db, genome_path ], stderr=null)
    with open(pe_output, "w") as pe_file:
        subprocess.check_call([ "bwa", "aln", genome_db, pe_path ], stdout=pe_file, stderr=null)

    with open(bwa_output, "w") as bwa_file:
        subprocess.check_call([ "bwa", "samse",
                                "-r", "@RG\tID:ILLUMINA\tSM:48_2\tPL:ILLUMINA\tLB:LIB1",
                                genome_db,
                                pe_output,
                                pe_path ], stdout=bwa_file, stderr=null)

    shutil.move(bwa_output, output_path)
    shutil.rmtree(work_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Maps the given reads with bwa.")
    parser.add_argument('pe1_path', type=str, help='Path to the first reads in a read pair (if paired reads). Just the reads if single reads')
    parser.add_argument('pe2_path', type=str, nargs='?', default=False, help='Path to the second pairs. Leave unspecified if single reads.')
    parser.add_argument('genome_path', type=str, help='Path to the reference genome/contigs.')
    parser.add_argument('output_path', type=str, help='Output path of resulting .bam and .bai file.')
    parser.add_argument('--threads', type=str, default='8', required=False, help='Number of threads for bwa mem.')
    parser.add_argument('--nomem', action="store_true", required=False,
                        help='bwa mem default, If flag specified thescript uses old bwa algorithm with "aln" and "sampe". ')

    args = parser.parse_args()


    if args.pe2_path and args.nomem:
        bwa_sampe(args.pe1_path, args.pe2_path, args.genome_path, args.output_path)
    elif args.pe2_path:
        bwa_mem(args.pe1_path, args.pe2_path, args.genome_path, args.threads, args.output_path)
    else:
        map_single_reads(args.pe1_path, args.genome_path, args.output_path)

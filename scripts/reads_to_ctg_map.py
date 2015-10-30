#!/usr/bin/env python
'''
Created on Jul 25, 2013
Modified on Oct 30, 2015

@author: ksahlin
@author: sidorov-si
'''
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime 
from sys import stdout
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
def bwa_sampe(pe1_path, pe2_path, genome_path, output_path, tmp_dir, bwa_path, clear_tmp):
    print
    print 'Aligning with bwa aln/sampe.'
    stdout.flush()
    start = datetime.now()
    work_dir = tempfile.mkdtemp() if tmp_dir == None else tmp_dir
    genome_db = os.path.join(work_dir, "genome")
    pe1_output = os.path.join(work_dir, "pe1.sai")
    pe2_output = os.path.join(work_dir, "pe2.sai")
    bwa_output = os.path.join(work_dir, "output.sam")
    stderr_file = open(output_path + '.bwa.1', 'w')
    print 'Temp directory:', work_dir
    print 'Output path:   ', output_path
    print 'Stderr file:   ', output_path + '.bwa.1' 
    stdout.flush()

    print 'Make bwa index...',
    stdout.flush()
    subprocess.check_call([ bwa_path, "index", "-p", genome_db, genome_path ], stderr=stderr_file)
    print 'Done.'
    stdout.flush()
    with open(pe1_output, "w") as pe1_file:
        print 'Align forward reads with bwa aln...',
        stdout.flush()
        subprocess.check_call([ bwa_path, "aln", "-t 8", genome_db, pe1_path ], 
                              stdout=pe1_file, stderr=stderr_file)
        print 'Done.'
        stdout.flush()

    with open(pe2_output, "w") as pe2_file:
        print 'Align reverse reads with bwa aln...',
        stdout.flush()
        subprocess.check_call([ bwa_path, "aln", "-t 8", genome_db, pe2_path ], 
                              stdout=pe2_file, stderr=stderr_file)
        print 'Done.'
        stdout.flush()

    with open(bwa_output, "w") as bwa_file:
        print 'Start bwa sampe...',
        stdout.flush()
        subprocess.check_call([ bwa_path, "sampe",
                                genome_db,
                                pe1_output, pe2_output,
                                pe1_path, pe2_path ], stdout=bwa_file, stderr=stderr_file)
        print 'Done.'
        stdout.flush()

    elapsed = datetime.now() - start
    print 'Time elapsed for bwa index and aln/sampe:', elapsed

    print
    print 'Convert SAM to BAM...',
    stdout.flush()
    start = datetime.now()
    sam_to_bam(bwa_output, bwa_output + ".bam")
    if clear_tmp:
        os.remove(bwa_output)
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for SAM to BAM conversion:', elapsed

    print
    print 'Sort BAM...',
    stdout.flush()
    start = datetime.now()
    pysam.sort(bwa_output + ".bam", output_path)
    if clear_tmp:
        os.remove(bwa_output + ".bam")
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for BAM sorting:', elapsed

    print
    print 'Index BAM...',
    stdout.flush()
    start = datetime.now()
    pysam.index(output_path + '.bam')
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for BAM indexing:', elapsed

    print
    print 'Remove temp files...',
    stdout.flush()
    start = datetime.now()
    shutil.rmtree(work_dir)
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for temp files removing:', elapsed
    stdout.flush()


def bwa_mem(pe1_path, pe2_path, genome_path, threads, output_path, tmp_dir, bwa_path, clear_tmp):
    print
    print 'Aligning with bwa mem.'
    stdout.flush()
    start = datetime.now()
    work_dir = tempfile.mkdtemp() if tmp_dir == None else tmp_dir
    genome_db = os.path.join(work_dir, "genome")
    pe1_output = os.path.join(work_dir, "pe1.sai")
    pe2_output = os.path.join(work_dir, "pe2.sai")
    bwa_output = os.path.join(work_dir, "output.sam")
    stderr_file = open(output_path + '.bwa.1', 'w')
    print 'Temp directory:', work_dir
    print 'Output path:   ', output_path
    print 'Stderr file:   ', output_path + '.bwa.1'
    stdout.flush()

    print 'Make bwa index...',
    stdout.flush()
    subprocess.check_call([ bwa_path, "index", "-p", genome_db, genome_path ], stderr=stderr_file)
    print 'Done.'
    stdout.flush()
    with open(bwa_output, "w") as bwa_file:
        print 'Align with bwa mem...',
        stdout.flush()
        subprocess.check_call([ bwa_path, "mem", "-t", threads,
                                genome_db, pe1_path, pe2_path ],
                                stdout=bwa_file,
                                stderr=stderr_file)
        print 'Done.'
        stdout.flush()

    elapsed = datetime.now() - start
    print 'Time elapsed for bwa index and mem: ', elapsed

    print
    print 'Convert SAM to BAM...',
    stdout.flush()
    start = datetime.now()
    sam_to_bam(bwa_output, bwa_output + ".bam")
    if clear_tmp:
        os.remove(bwa_output)
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for SAM to BAM conversion:', elapsed

    print
    print 'Sort BAM...',
    stdout.flush()
    start = datetime.now()
    pysam.sort(bwa_output + ".bam", output_path)
    if clear_tmp:
        os.remove(bwa_output + ".bam")
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for BAM sorting:', elapsed

    print
    print 'Index BAM...',
    stdout.flush()
    start = datetime.now()
    pysam.index(output_path + '.bam')
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for BAM indexing:', elapsed
   
    print
    print 'Remove temp files...',
    stdout.flush()
    start = datetime.now()
    shutil.rmtree(work_dir)
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for temp files removing:', elapsed
    stdout.flush()


def map_single_reads(pe_path, genome_path, output_path, bwa_path):
    print
    print 'Aligning with bwa aln/samse.'
    stdout.flush()
    start = datetime.now()
    work_dir = tempfile.mkdtemp()
    genome_db = os.path.join(work_dir, "genome")
    pe_output = os.path.join(work_dir, "pe.sai")
    bwa_output = os.path.join(work_dir, "output.sam")
    print 'Temp directory:', work_dir
    print 'Output path:   ', output_path
    print 'Stderr file:    /dev/null'
    stdout.flush()

    null = open("/dev/null")
    print 'Make bwa index...',
    stdout.flush()
    subprocess.check_call([ bwa_path, "index", "-p", genome_db, genome_path ], stderr=null)
    print 'Done.'
    stdout.flush()
    with open(pe_output, "w") as pe_file:
        print 'Align with bwa aln...',
        stdout.flush()
        subprocess.check_call([ bwa_path, "aln", genome_db, pe_path ], stdout=pe_file, stderr=null)
        print 'Done.'
        stdout.flush()

    with open(bwa_output, "w") as bwa_file:
        print 'Start bwa samse...',
        stdout.flush()
        subprocess.check_call([ bwa_path, "samse",
                                "-r", "@RG\tID:ILLUMINA\tSM:48_2\tPL:ILLUMINA\tLB:LIB1",
                                genome_db,
                                pe_output,
                                pe_path ], stdout=bwa_file, stderr=null)
        print 'Done.'
        stdout.flush()

    elapsed = datetime.now() - start
    print 'Time elapsed for bwa index and aln/samse:', elapsed

    print
    print 'Copy the result to the output directory and remove temp files...',
    stdout.flush()
    start = datetime.now()
    shutil.move(bwa_output, output_path)
    shutil.rmtree(work_dir)
    elapsed = datetime.now() - start
    print 'Done.'
    print 'Time elapsed for copying result to the output directory and removing temp files:', elapsed
    stdout.flush()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Maps the given reads with bwa.")
    parser.add_argument('pe1_path', type=str, help='Path to the first reads in a read pair (if paired reads). Just the reads if single reads')
    parser.add_argument('pe2_path', type=str, nargs='?', default=False, help='Path to the second pairs. Leave unspecified if single reads.')
    parser.add_argument('genome_path', type=str, help='Path to the reference genome/contigs.')
    parser.add_argument('output_path', type=str, help='Output path of resulting .bam and .bai file.')
    parser.add_argument('--tmp_path', type=str, required=False, help='Output path of temporary files.')
    parser.add_argument('--bwa_path', type=str, default='bwa', required=False, help='Path to bwa binary with bwa binary name at the end.')
    parser.add_argument('--threads', type=str, default='8', required=False, help='Number of threads for bwa mem.')
    parser.add_argument('--clear', action="store_true", required=False, 
                        help='Remove SAM file when BAM is already created, and remove BAM file when sorted BAM is already created.')
    parser.add_argument('--nomem', action="store_true", required=False,
                        help='bwa mem default, If flag specified the script uses old bwa algorithm with "aln" and "sampe". ')

    args = parser.parse_args()

    tmp_path = args.tmp_path
    if tmp_path != None and not os.path.exists(tmp_path):
        os.makedirs(tmp_path)

    output_path = args.output_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    print
    print 'pe1_path:', args.pe1_path
    print 'pe2_path:', args.pe2_path
    print 'genome_path:', args.genome_path
    print 'output_path:', output_path
    if tmp_path != None:  
        print "tmp_path:", tmp_path
    print 'bwa path:', args.bwa_path
    print 'number of threads:', args.threads
    print 'Remove temp SAM and BAM files:',
    if args.clear:
        print 'Yes'
    else: 
        print 'No'
    print 'Use bwa aln and sampe instead of bwa mem:',
    if args.nomem:
        print 'Yes'
    else: 
        print 'No'
    stdout.flush()

    print 
    print 'Start processing.'
    stdout.flush()
    if args.pe2_path and args.nomem:
        bwa_sampe(args.pe1_path, args.pe2_path, args.genome_path, output_path, tmp_path, 
                  args.bwa_path, args.clear)
    elif args.pe2_path:
        bwa_mem(args.pe1_path, args.pe2_path, args.genome_path, args.threads, output_path, 
                tmp_path, args.bwa_path, args.clear)
    else:
        map_single_reads(args.pe1_path, args.genome_path, args.output_path, args.bwa_path)
    print
    print 'Processing is finished.'
    stdout.flush()



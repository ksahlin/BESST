'''
Created on Mar 4, 2013

@author: ksahlin
'''
'''
    Created on Sep 29, 2011
    
    @author: ksahlin
    
    This file is part of BESST.
    
    BESST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    BESST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with BESST.  If not, see <http://www.gnu.org/licenses/>.
    '''


import os
from time import time

import Parameter
import CreateGraph as CG
import MakeScaffolds as MS
import GenerateOutput as GO

def ReadInContigseqs(contigfile):
    cont_dict = {}
    k = 0
    temp = ''
    accession = ''
    for line in contigfile:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            cont_dict[accession] = ''
            k += 1
        elif line[0] == '>':
            cont_dict[accession] = temp
            temp = ''
            accession = line[1:].strip().split()[0]
        else:
            temp += line.strip()
    cont_dict[accession] = temp
    return(cont_dict)

def main(options):
    Parameter.MemoryUsage()
    tot_start = time()
    Contigs = {} # contig dict that stores contig objects 
    Scaffolds = {}     # scaffold dict with contig objects for easy fetching of all contigs in a scaffold
    small_contigs = {}
    small_scaffolds = {}

    param = Parameter.parameter() # object containing all parameters (user specified, defaulted and comuted along the way.)
    param.scaffold_indexer = 1 # global indicator for scaffolds, used to index scaffolds when they are created
    param.rel_weight = options.relweight
    param.multiprocess = options.multiprocess
    param.development = options.development
    options.output = options.output + '/BESST_output'
    os.popen('mkdir -p ' + options.output)
    param.information_file = open(options.output + '/Statistics.txt', 'w')
    Information = param.information_file
    open(options.output + '/haplotypes.fa', 'w')
    #Read in the sequences of the contigs in memory
    C_dict = ReadInContigseqs(open(options.contigfile, 'r'))

    #iterate over libraries
    param.first_lib = True
    for i in xrange(len(options.bamfiles)):
        start = time()
        param.bamfile = options.bamfiles[i]
        param.mean_ins_size = options.mean[i] if options.mean != None else None
        param.ins_size_threshold = options.threshold[i] if options.threshold != None else None
        param.edgesupport = options.edgesupport[i] if options.edgesupport != None else 5
        param.read_len = options.readlen[i] if options.readlen != None else None
        param.output_directory = options.output
        param.std_dev_ins_size = options.stddev[i] if options.stddev != None else None
        param.contig_threshold = options.minsize[i] if options.minsize != None else None
        param.cov_cutoff = options.covcutoff[i] if options.covcutoff != None else None
        param.hapl_ratio = options.haplratio
        param.hapl_threshold = options.haplthreshold
        param.detect_haplotype = options.haplotype
        param.detect_duplicate = options.duplicate
        param.fosmidpool = options.fosmidpool
        param.extend_paths = options.extendpaths
        print >> Information, '\nPASS ' + str(i + 1) + '\n\n'
        print 'Creating contig graph with library: ', param.bamfile
        (G, G_prime, param) = CG.PE(Contigs, Scaffolds, Information, param.output_directory, C_dict, param, small_contigs, small_scaffolds)      #Create graph, single out too short contigs/scaffolds and store them in F
        param.first_lib = False   #not the first lib any more
        if G == None:
            print '0 contigs/super-contigs passed the length criteria of this step. Exiting and printing results.. '
            break
        elapsed = time() - start
        print >> Information, 'Time elapsed for creating graph, iteration ' + str(i) + ': ' + str(elapsed) + '\n'
        start = time()

        Parameter.MemoryUsage()
        print 'Constructed contig graph. Start BESST algorithm for creating scaffolds. ',
        MS.Algorithm(G, G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, Information, param)   # Make scaffolds, store the complex areas (consisting of contig/scaffold) in F, store the created scaffolds in Scaffolds, update Contigs
        elapsed = time() - start

        Parameter.MemoryUsage()

        print >> Information, 'Time elapsed for making scaffolds, iteration ' + str(i) + ': ' + str(elapsed) + '\n'
        print 'Writing out scaffolding results for step', i + 1, ' ...'

        F = [] #list of (ordered) lists of tuples containing (contig_name, direction, position, length, sequence). The tuple is a contig within a scaffold and the list of tuples is the scaffold.
        for scaffold_ in small_scaffolds:
            S_obj = small_scaffolds[scaffold_]
            list_of_contigs = S_obj.contigs   #list of contig objects contained in scaffold object
            F = GO.WriteToF(F, small_contigs, list_of_contigs)
        for scaffold_ in Scaffolds.keys(): #iterate over keys in hash, so that we can remove keys while iterating over it
            ###  Go to function and print to F
            ### Remove Scaf_obj from Scaffolds and Contig_obj from contigs
            S_obj = Scaffolds[scaffold_]
            list_of_contigs = S_obj.contigs   # List of contig objects contained in scaffold object
            F = GO.WriteToF(F, Contigs, list_of_contigs)


        GO.PrintOutput(F, Information, param.output_directory, param, i + 1)
        Parameter.MemoryUsage()

    ### Calculate stats for last scaffolding step    
    scaf_lengths = [Scaffolds[scaffold_].s_length for scaffold_ in Scaffolds.keys()]
    sorted_lengths = sorted(scaf_lengths, reverse=True)
    scaf_lengths_small = [small_scaffolds[scaffold_].s_length for scaffold_ in small_scaffolds.keys()]
    sorted_lengths_small = sorted(scaf_lengths_small, reverse=True)
    NG50, LG50 = CG.CalculateStats(sorted_lengths, sorted_lengths_small, param, Information)
    param.current_LG50 = LG50
    param.current_NG50 = NG50

    elapsed = time() - tot_start
    print >> Information, 'Total time for scaffolding: ' + str(elapsed) + '\n'
    print 'Finished\n\n '
    Parameter.MemoryUsage()

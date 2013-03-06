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


import sys
from collections import defaultdict
from math import pi
from time import time

import networkx as nx
import pysam

import Contig
import Scaffold
from Parameter import counters
from Norm import normcdf, normpdf, erf
import GenerateOutput as GO
import GapCalculator as GC

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    pass

def PE(Contigs, Scaffolds, Information, output_dest, C_dict, param, small_contigs, small_scaffolds):
    G = nx.Graph()
    G_prime = nx.Graph()  # If we want to do path extension with small contigs
    print >> Information, 'Parsing BAM file...'

    with pysam.Samfile(param.bamfile, 'rb') as bam_file:

        #### Get parameters -r, -m, -s, -T, -t for library ####

        print >> Information, 'Computing parameters not set by user...'
        GetParams(bam_file, param, Scaffolds, small_scaffolds, Contigs, small_contigs, Information)


        ##### Initialize contig and scaffold objects ######

        if param.first_lib:
            start_init = time()
            InitializeObjects(bam_file, Contigs, Scaffolds, param, Information, G_prime, small_contigs, small_scaffolds, C_dict)
            print >> Information, 'Time initializing BESST objects: ', time() - start_init

        else:
            #Clean contig_library/scaffold_library
            start_clean = time()
            CleanObjects(Contigs, Scaffolds, param, Information, small_contigs, small_scaffolds)
            print >> Information, 'Time cleaning BESST objects for next library: ', time() - start_clean

        print >> Information, 'Nr of contigs/scaffolds included in scaffolding: ' + str(len(Scaffolds)) #,Scaffolds.keys()
        if len(Scaffolds) == 0:
            return(G, G_prime, param)

        ### initialize graph objects two nodes per contig "left" and "right" node. ###    
        tot_start = time()
        if param.extend_paths:
            InitializeGraph(Scaffolds, G, Information)

            #small contig graph contains all scaffolds
            InitializeGraph(small_scaffolds, G_prime, Information)
            InitializeGraph(Scaffolds, G_prime, Information)
        else:
            InitializeGraph(Scaffolds, G, Information)
        print >> Information, 'Total time elapsed for initializing Graph: ', time() - tot_start

        #for coverage computation
        cont_aligned_len = {}
        for contig in Contigs:
            cont_aligned_len[contig] = [0, Contigs[contig].length]

        #if param.extend_paths:
        for contig in small_contigs:
            cont_aligned_len[contig] = [0, small_contigs[contig].length]

        #initialize counters for library
        counter = counters(0, 0, 0, 0, -1, -1, 0)
        fishy_edges = defaultdict(int)
        ctr = 0
        # Create the link edges in the graph by fetching info from bam file  
        print >> Information, 'Reading bam file and creating scaffold graph...'
        staart = time()
        for alignedread in bam_file:
            try: #check that read is aligned OBS: not with is_unmapped since this flag is fishy for e.g. BWA
                contig1 = bam_file.getrname(alignedread.rname)
                contig2 = bam_file.getrname(alignedread.mrnm)
            except ValueError:
                continue
            #TODO:Repeats (and haplotypes) may have been singled out, we need this statement (or a smarter version of it)
            if (contig1 in Contigs or contig1 in small_contigs) and (contig2 in Contigs or contig2 in small_contigs):
                pass
            else:
                continue

            #TODO: this if-statement is an ad hoc implementation to deal with BWA's buggy SAM-flag reporting
            #if BWA fixes this -> remove this statement. If the links in fishy edges is equal to or ore than
            #the links in the graph G or G'. The edge will be removed.
            if alignedread.is_unmapped and alignedread.is_read1: # and contig1 != contig2: 
                #Some BWA error in mappings can still slip through, these edges are caracterized by very few links                 
                try:
                    cont_obj1 = small_contigs[contig1]
                    scaf_obj1 = small_scaffolds[cont_obj1.scaffold]
                except KeyError:
                    cont_obj1 = Contigs[contig1]
                    scaf_obj1 = Scaffolds[cont_obj1.scaffold]
                try:
                    cont_obj2 = small_contigs[contig2]
                    scaf_obj2 = small_scaffolds[cont_obj2.scaffold]
                except KeyError:
                    cont_obj2 = Contigs[contig2]
                    scaf_obj2 = Scaffolds[cont_obj2.scaffold]

                if scaf_obj2.name != scaf_obj1.name:
                    (side1, side2) = CheckDir(cont_obj1, cont_obj2, alignedread)
                    #get scaffold name for contig
                    s1 = Contigs[contig1].scaffold if contig1 in Contigs else small_contigs[contig1].scaffold
                    s2 = Contigs[contig2].scaffold if contig2 in Contigs else small_contigs[contig2].scaffold
                    fishy_edges[((s1, side1), (s2, side2))] += 1
                    fishy_edges[((s2, side2), (s1, side1))] += 1
                    ctr += 1

            ## add to coverage computation if contig is still in the list of considered contigs

            cont_aligned_len[contig1][0] += alignedread.rlen

            if contig1 != contig2 and alignedread.mapq == 0:
                counter.non_unique += 1
            if contig1 != contig2 and alignedread.is_read2 and not alignedread.is_unmapped and alignedread.mapq > 10:
                #check how many non unique reads out of the useful ones (mapping to two different contigs)
                if contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold: # and alignedread.tags[0][1] == 'U':
                    cont_obj1 = Contigs[contig1]
                    cont_obj2 = Contigs[contig2]
                    scaf_obj1 = Scaffolds[cont_obj1.scaffold]
                    scaf_obj2 = Scaffolds[cont_obj2.scaffold]
                    is_dupl = CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G, param, alignedread, counter, contig1, contig2)
                    if param.extend_paths and not is_dupl:
                        counter.prev_obs1 = -1
                        counter.prev_obs2 = -1
                        CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2)
                elif param.extend_paths:
                    try:
                        cont_obj1 = small_contigs[contig1]
                        scaf_obj1 = small_scaffolds[cont_obj1.scaffold]
                        value1 = 1
                    except KeyError:
                        cont_obj1 = Contigs[contig1]
                        scaf_obj1 = Scaffolds[cont_obj1.scaffold]
                        value1 = 0
                    try:
                        cont_obj2 = small_contigs[contig2]
                        scaf_obj2 = small_scaffolds[cont_obj2.scaffold]
                        value2 = 1
                    except KeyError:
                        cont_obj2 = Contigs[contig2]
                        scaf_obj2 = Scaffolds[cont_obj2.scaffold]
                        value2 = 0
                    if value1 and value2 and scaf_obj1 != scaf_obj2:
                        CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2)
                    elif value1 and not value2:
                        CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2)
                    elif not value1 and value2:
                        CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2)


                elif contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold:
########################Use to validate scaffold in previous step here ############
                    pass
        print >> Information, 'ELAPSED reading file:', time() - staart
        print >> Information, 'NR OF FISHY READ LINKS: ', ctr #,len(fishy_edges)/2.0, fishy_edges

        print >> Information, 'USEFUL READS (reads mapping to different contigs): ', counter.count
        print >> Information, 'Non unique portion out of "USEFUL READS"  (filtered out from scaffolding): ', counter.non_unique
        #print 'Non unique used for scaf: ', non_unique_for_scaf
        print >> Information, 'Reads with too large insert size from "USEFUL READS" (filtered out): ', counter.reads_with_too_long_insert
        if param.detect_duplicate:
            print >> Information, 'Number of duplicated reads indicated and removed: ', counter.nr_of_duplicates



    ##### Calc coverage for all contigs with current lib here #####
        sum_x = 0
        sum_x_sq = 0
        n = 0
        for contig in cont_aligned_len:
            cont_coverage = cont_aligned_len[contig][0] / float(cont_aligned_len[contig][1])
            try:
                Contigs[contig].coverage = cont_coverage
            except KeyError:
                small_contigs[contig].coverage = cont_coverage
#            else:
#                try:
#                    Contigs[contig].coverage=cont_coverage
#                except KeyError:
#                    pass
            sum_x += cont_coverage
            sum_x_sq += cont_coverage ** 2
            n += 1

        mean_cov, std_dev_cov = CalculateMeanCoverage(Contigs, param.first_lib, output_dest, param.bamfile, Information)
        param.mean_coverage = mean_cov
        param.std_dev_coverage = std_dev_cov
        if param.first_lib:
            Contigs, Scaffolds, G = RepeatDetector(Contigs, Scaffolds, G, param, G_prime, small_contigs, small_scaffolds, Information)


        ### Remove edges created by false reporting of BWA ###
        RemoveBugEdges(G, G_prime, fishy_edges, param, Information)


        #### temp check mean and std_dev and sign value of all links of all edges ###
        GiveScoreOnEdges(G, Scaffolds, small_scaffolds, Contigs, param, Information)

        GiveScoreOnEdges(G_prime, Scaffolds, small_scaffolds, Contigs, param, Information)


        #Remove all edges with link support less than 3 to be able to compute statistics: 
        cntr_sp = 0
        for edge in G_prime.edges():
            if G_prime[edge[0]][edge[1]]['nr_links'] != None:
                if G_prime[edge[0]][edge[1]]['nr_links'] < 3:
                    G_prime.remove_edge(edge[0], edge[1])
                    cntr_sp += 1
        print >> Information, 'Number of fishy edges in G_prime', cntr_sp



    return(G, G_prime, param)


def GiveScoreOnEdges(G, Scaffolds, small_scaffolds, Contigs, param, Information):

    cnt_sign = 0

    for edge in G.edges():
        mean_ = 0
        std_dev = 0
        if G[edge[0]][edge[1]]['nr_links'] != None:
            n = G[edge[0]][edge[1]]['nr_links']
            obs_squ = G[edge[0]][edge[1]]['obs_sq']
            mean_ = G[edge[0]][edge[1]]['obs'] / float(n)
            data_observation = (n * param.mean_ins_size - G[edge[0]][edge[1]]['obs']) / float(n)
            try:
                len1 = Scaffolds[ edge[0][0] ].s_length
            except KeyError:
                len1 = small_scaffolds[ edge[0][0] ].s_length
            try:
                len2 = Scaffolds[ edge[1][0] ].s_length
            except KeyError:
                len2 = small_scaffolds[ edge[1][0] ].s_length
            if 2 * param.std_dev_ins_size < len1 and 2 * param.std_dev_ins_size < len2:
                gap = GC.GapEstimator(param.mean_ins_size, param.std_dev_ins_size, param.read_len, data_observation, len1, len2)
            else:
                gap = data_observation

            G[edge[0]][edge[1]]['gap'] = int(gap)
            if -gap > len1 or -gap > len2:
                G[edge[0]][edge[1]]['score'] = 0
                continue
            if 2 * param.std_dev_ins_size < len1 and 2 * param.std_dev_ins_size < len2:
                std_dev_d_eq_0 = CalcTheoreticalStdDev(param.mean_ins_size, param.std_dev_ins_size, param.read_len, len1, len2, gap)
            else:
                std_dev_d_eq_0 = 2 ** 32
            try:
                std_dev = ((obs_squ - n * mean_ ** 2) / (n - 1)) ** 0.5
                chi_sq = (n - 1) * (std_dev ** 2 / std_dev_d_eq_0 ** 2)
            except ZeroDivisionError:
                std_dev = 2 ** 32
                chi_sq = 0


            k = MaxObsDistr(n)
            if 2 * param.read_len < len1 and 2 * param.read_len < len2:
                span_max1 = min(param.mean_ins_size + k * param.std_dev_ins_size - 2 * param.read_len, len1 - param.read_len + max(0, gap))
                span_max2 = min(param.mean_ins_size + k * param.std_dev_ins_size - 2 * param.read_len, len2 - param.read_len + max(0, gap))
                try:
                    span_obs1 = Scaffolds[ edge[0][0] ].upper_right_nbrs_obs[edge[1]] - Scaffolds[ edge[0][0] ].lower_right_nbrs_obs[edge[1]] if edge[0][1] == 'R' else Scaffolds[ edge[0][0] ].upper_left_nbrs_obs[edge[1]] - Scaffolds[ edge[0][0] ].lower_left_nbrs_obs[edge[1]]
                except KeyError:
                    span_obs1 = small_scaffolds[ edge[0][0] ].upper_right_nbrs_obs[edge[1]] - small_scaffolds[ edge[0][0] ].lower_right_nbrs_obs[edge[1]] if edge[0][1] == 'R' else small_scaffolds[ edge[0][0] ].upper_left_nbrs_obs[edge[1]] - small_scaffolds[ edge[0][0] ].lower_left_nbrs_obs[edge[1]]
                try:
                    span_obs2 = Scaffolds[ edge[1][0] ].upper_right_nbrs_obs[edge[0]] - Scaffolds[ edge[1][0] ].lower_right_nbrs_obs[edge[0]] if edge[1][1] == 'R' else Scaffolds[ edge[1][0] ].upper_left_nbrs_obs[edge[0]] - Scaffolds[ edge[1][0] ].lower_left_nbrs_obs[edge[0]]
                except KeyError:
                    span_obs2 = small_scaffolds[ edge[1][0] ].upper_right_nbrs_obs[edge[0]] - small_scaffolds[ edge[1][0] ].lower_right_nbrs_obs[edge[0]] if edge[1][1] == 'R' else small_scaffolds[ edge[1][0] ].upper_left_nbrs_obs[edge[0]] - small_scaffolds[ edge[1][0] ].lower_left_nbrs_obs[edge[0]]


                span_score1 = min((max(0, gap) + 2 * param.read_len + span_obs1) / float(span_max1) , float(span_max1) / (max(0, gap) + 2 * param.read_len + span_obs1)) if span_obs1 > 0 else 0
                span_score2 = min((max(0, gap) + 2 * param.read_len + span_obs2) / float(span_max2) , float(span_max2) / (max(0, gap) + 2 * param.read_len + span_obs2)) if span_obs2 > 0 else 0
                span_score = min(span_score1, span_score2)
            else:
                span_score = 0


            try:
                std_dev_score = min(std_dev / std_dev_d_eq_0, std_dev_d_eq_0 / std_dev) #+ span_score #+ min(n/E_links, E_links/float(n))
            except ZeroDivisionError:
                std_dev_score = 0
                sys.stderr.write(str(std_dev) + ' ' + str(std_dev_d_eq_0) + ' ' + str(span_score) + '\n')

            G[edge[0]][edge[1]]['score'] = std_dev_score + span_score if std_dev_score > 0.5 and span_score > 0.5 else 0

    for edge in G.edges():
        if G[edge[0]][edge[1]]['nr_links'] != None:
            try:
                G[edge[0]][edge[1]]['score']
            except KeyError:
                sys.stderr.write(str(G[edge[0]][edge[1]]) + ' ' + str(Scaffolds[edge[0][0]].s_length) + ' ' + str(Scaffolds[edge[1][0]].s_length))
    print >> Information, 'Number of significantly spurious edges:', cnt_sign

    return()




def CalcTheoreticalStdDev(mean, stdDev, readLen, c1Len, c2Len, d):
    def E_O(d, mean, stdDev, c1Len, c2Len, readLen):

        def CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen):

            term1 = -0.5 * (erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) + erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))))
            term2 = +0.5 * (erf((d + c_max + readLen - mean) / (2 ** 0.5 * float(stdDev))) + erf((d + c_min + readLen - mean) / (2 ** 0.5 * float(stdDev))))
            g_prime_d = term1 + term2
            return g_prime_d

        def CalcGd(d, c1Len, c2Len, readLen):
            #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
            term1 = (c_min - readLen + 1) / 2.0 * (erf((c_max + d + readLen - mean) / ((2 ** 0.5) * stdDev)) - erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

            term2 = (c_min + c_max + d - mean + 1) / 2.0 * (erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) - erf((c_max + readLen + d - mean) / (2 ** 0.5 * float(stdDev))))

            term3 = (d + 2 * readLen - mean - 1) / 2.0 * (erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))) - erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

            term4 = stdDev / ((2 * pi) ** 0.5) * (2.718 ** (-((c_min + c_max + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))) + 2.718 ** (-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2))))

            term5 = -stdDev / ((2 * pi) ** 0.5) * (2.718 ** (-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))) + 2.718 ** (-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))))

            g_d = term1 + term2 + term3 + term4 + term5
            return g_d

        c_min = min(c1Len, c2Len)
        c_max = max(c1Len, c2Len)
        g_prime_d = CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen)
        g_d = CalcGd(d, c1Len, c2Len, readLen)
        a = stdDev ** 2 * g_prime_d + mean * g_d
        E_o = a / g_d - d
        return E_o

    def E_O_square(d, mean, stdDev, c1Len, c2Len, readLen):

        def CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen):

            term1 = -0.5 * (erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) + erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))))
            term2 = +0.5 * (erf((d + c_max + readLen - mean) / (2 ** 0.5 * float(stdDev))) + erf((d + c_min + readLen - mean) / (2 ** 0.5 * float(stdDev))))
            g_prime_d = term1 + term2
            return g_prime_d

        def CalcGd(d, c1Len, c2Len, readLen):
            #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
            term1 = (c_min - readLen + 1) / 2.0 * (erf((c_max + d + readLen - mean) / ((2 ** 0.5) * stdDev)) - erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

            term2 = (c_min + c_max + d - mean + 1) / 2.0 * (erf((c_min + c_max + d + 1 - mean) / (2 ** 0.5 * float(stdDev))) - erf((c_max + readLen + d - mean) / (2 ** 0.5 * float(stdDev))))

            term3 = (d + 2 * readLen - mean - 1) / 2.0 * (erf((d + 2 * readLen - 1 - mean) / (2 ** 0.5 * float(stdDev))) - erf((c_min + d + readLen - mean) / ((2 ** 0.5) * stdDev)))

            term4 = stdDev / ((2 * pi) ** 0.5) * (2.718 ** (-((c_min + c_max + d + 1 - mean) ** 2) / (float(2 * stdDev ** 2))) + 2.718 ** (-((d + 2 * readLen - 1 - mean) ** 2) / (float(2 * stdDev ** 2))))

            term5 = -stdDev / ((2 * pi) ** 0.5) * (2.718 ** (-((c_max + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))) + 2.718 ** (-((c_min + readLen + d - mean) ** 2) / (float(2 * stdDev ** 2))))
            g_d = term1 + term2 + term3 + term4 + term5
            return g_d

        def CalcB(d, c_min, c_max, c1Len, c2Len, readLen, g_d, g_prime_d):
            c1 = 0
            c2 = normpdf(d + 2 * readLen - 1, mean, stdDev) - normpdf(c_min + d + readLen , mean, stdDev) - normpdf(c_max + d + readLen , mean, stdDev) + normpdf(c_min + c_max + d + 1 , mean, stdDev)

            b = stdDev ** 4 * (c1 + c2) + mean ** 2 * g_d + stdDev ** 2 * g_d + 2 * stdDev ** 2 * mean * g_prime_d
            return b, c1, c2

        c_min = min(c1Len, c2Len)
        c_max = max(c1Len, c2Len)
        g_prime_d = CalcG_prime_d(d, c_min, c_max, c1Len, c2Len, readLen)
        g_d = CalcGd(d, c1Len, c2Len, readLen)
        a = stdDev ** 2 * g_prime_d + mean * g_d
        b, c1, c2 = CalcB(d, c_min, c_max, c1Len, c2Len, readLen, g_d, g_prime_d)
        E_o_square = b / g_d - 2 * d * (stdDev ** 2 * (g_prime_d / g_d) + mean) + d ** 2
        return E_o_square



    e_o = E_O(d, mean, stdDev, c1Len, c2Len, readLen)
    e_o_square = E_O_square(d, mean, stdDev, c1Len, c2Len, readLen)
    if e_o_square - e_o ** 2 < 0:
        sys.stderr.write('Std_dev(O|d =' + str(d) + ')=' + str(e_o_square) + '  ' + str(e_o ** 2) + str(c1Len) + ' ' + str(c2Len))
        std_dev = 0
    else:
        std_dev = (e_o_square - e_o ** 2) ** 0.5
    return(std_dev)

def CheckDir(cont_obj1, cont_obj2, alignedread):
    (read_dir, mate_dir) = (not alignedread.is_reverse, not alignedread.mate_is_reverse)
    cont_dir1 = cont_obj1.direction  #if pos : L if neg: R
    #position2 cont2/scaf2                        
    cont_dir2 = cont_obj2.direction
    (obs1, obs2, scaf_side1, scaf_side2) = PosDirCalculatorPE(cont_dir1, read_dir, 0, 0, 0, 0, cont_dir2, mate_dir, 0, 0, 0, 0, 0)
    return(scaf_side1, scaf_side2)

def RemoveBugEdges(G, G_prime, fishy_edges, param, Information):
    edges_removed = 0
    for edge_tuple, nr_links in fishy_edges.items():
        if param.extend_paths:
            if edge_tuple[1] in G_prime and edge_tuple[0] in G_prime[edge_tuple[1]]:
                if nr_links >= G_prime[edge_tuple[0]][edge_tuple[1]]['nr_links']:
                    G_prime.remove_edge(edge_tuple[0], edge_tuple[1])
                    edges_removed += 1
            if edge_tuple[1] in G and edge_tuple[0] in G[edge_tuple[1]]:
                if nr_links >= G[edge_tuple[0]][edge_tuple[1]]['nr_links']:
                    G.remove_edge(edge_tuple[0], edge_tuple[1])
        else:
            if edge_tuple[1] in G and edge_tuple[0] in G[edge_tuple[1]]:
                if nr_links >= G[edge_tuple[0]][edge_tuple[1]]['nr_links']:
                    G.remove_edge(edge_tuple[0], edge_tuple[1])
                    edges_removed += 1
    print >> Information, 'Number of BWA buggy edges removed: ', edges_removed
    return()

def InitializeGraph(dict_with_scaffolds, graph, Information):
    cnt = 1
    start1 = time()
    for scaffold_ in dict_with_scaffolds:
        graph.add_edge((scaffold_, 'L'), (scaffold_, 'R'), nr_links=None)    #this is a scaffold object but can be both a single contig or a scaffold.
        graph.node[(scaffold_, 'L')]['length'] = dict_with_scaffolds[scaffold_].s_length
        graph.node[(scaffold_, 'R')]['length'] = dict_with_scaffolds[scaffold_].s_length
        #dict_with_scaffolds[ scaffold_ ].scaffold_left_nbrs = {}
        #dict_with_scaffolds[ scaffold_ ].scaffold_right_nbrs  = {}  
        if cnt % 100000 == 0 and cnt > 0:
            elapsed = time() - start1
            print >> Information, 'Total nr of keys added: ', cnt, 'Time for adding last 100 000 keys: ', elapsed
            start1 = time()
        cnt += 1
    return()

def constant_large():
    return 2 ** 32
def constant_small():
    return -1

def InitializeObjects(bam_file, Contigs, Scaffolds, param, Information, G_prime, small_contigs, small_scaffolds, C_dict):
    singeled_out = 0
    contig_threshold = param.contig_threshold
    cont_lengths = bam_file.lengths
    cont_lengths = [int(nr) for nr in cont_lengths]  #convert long to int object
    cont_names = bam_file.references

    #Calculate NG50 and LG 50
    param.tot_assembly_length = sum(cont_lengths)
    sorted_lengths = sorted(cont_lengths, reverse=True)
    NG50, LG50 = CalculateStats(sorted_lengths, [], param, Information)
    param.current_LG50 = LG50
    param.current_NG50 = NG50
    #extend_paths = param.extend_paths
    counter = 0
    start = time()
    for i in range(0, len(cont_names)):
        counter += 1
        if counter % 100000 == 0:
            print >> Information, 'Time adding 100k keys', time() - start
            start = time()
        if cont_lengths[i] >= contig_threshold:
            C = Contig.contig(cont_names[i])   # Create object contig
            C.length = cont_lengths[i]
            C.sequence = C_dict[cont_names[i]]
            del C_dict[cont_names[i]]
            scaf_length = C.length        # Initially, scaffold consists of only this contig    
            C.direction = True              # always in same direction first, False=reverse
            C.position = 0                  #position always 0
            #C.links = {}
            Contigs[C.name] = C              # Create a dict with name as key and the object container as value
            S = Scaffold.scaffold(param.scaffold_indexer, [C], scaf_length, defaultdict(constant_large), defaultdict(constant_large), defaultdict(constant_small), defaultdict(constant_small))  # Create object scaffold
            Scaffolds[S.name] = S
            C.scaffold = S.name
            param.scaffold_indexer += 1
        else:
            if cont_lengths[i] > 0: #In case of contigs with size 0 (due to some error in fasta file)
                C = Contig.contig(cont_names[i])   # Create object contig
                C.length = cont_lengths[i]
                C.sequence = C_dict[cont_names[i]]
                del C_dict[cont_names[i]]
                scaf_length = C.length        # Initially, scaffold consists of only this contig    
                C.direction = True              # always in same direction first, False=reverse
                C.position = 0                  #position always 0
                small_contigs[C.name] = C              # Create a dict with name as key and the object container as value
                S = Scaffold.scaffold(param.scaffold_indexer, [C], scaf_length, defaultdict(constant_large), defaultdict(constant_large), defaultdict(constant_small), defaultdict(constant_small))  # Create object scaffold
                small_scaffolds[S.name] = S
                C.scaffold = S.name
                param.scaffold_indexer += 1
                singeled_out += 1
    del C_dict


    print >> Information, 'Nr of contigs that was singeled out due to length constraints ' + str(singeled_out)
    return()

def CleanObjects(Contigs, Scaffolds, param, Information, small_contigs, small_scaffolds):
    singeled_out = 0
    scaf_lengths = [Scaffolds[scaffold_].s_length for scaffold_ in Scaffolds.keys()]
    sorted_lengths = sorted(scaf_lengths, reverse=True)
    scaf_lengths_small = [small_scaffolds[scaffold_].s_length for scaffold_ in small_scaffolds.keys()]
    sorted_lengths_small = sorted(scaf_lengths_small, reverse=True)
    NG50, LG50 = CalculateStats(sorted_lengths, sorted_lengths_small, param, Information)
    param.current_LG50 = LG50
    param.current_NG50 = NG50
    for scaffold_ in Scaffolds.keys(): #iterate over keys in hash, so that we can remove keys while iterating over it
        if Scaffolds[scaffold_].s_length < param.contig_threshold:
            ###  Switch from Scaffolds to small_scaffolds (they can still be used in the path extension)
            ### Remove Scaf_obj from Scaffolds and Contig_obj from contigs
            S_obj = Scaffolds[scaffold_]
            list_of_contigs = S_obj.contigs   #list of contig objects contained in scaffold object
            GO.ChangeToSmallContigs(Contigs, list_of_contigs, small_contigs)
            scaf_obj = Scaffolds[scaffold_]
            small_scaffolds[scaffold_] = scaf_obj
            del Scaffolds[scaffold_]
            singeled_out += 1

    print >> Information, 'Nr of contigs/scaffolds that was singeled out due to length constraints ' + str(singeled_out)
    return()

def CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G, param, alignedread, counter, contig1, contig2):
    if alignedread.mapq == 0:
        counter.non_unique_for_scaf += 1
    counter.count += 1
    (read_dir, mate_dir) = (not alignedread.is_reverse, not alignedread.mate_is_reverse)
    #Calculate actual position on scaffold here
    #position1 cont/scaf1
    cont_dir1 = cont_obj1.direction  #if pos : L if neg: R
    cont1_pos = cont_obj1.position
    readpos = alignedread.pos
    cont1_len = cont_obj1.length
    s1len = scaf_obj1.s_length
    #position2 cont2/scaf2                        
    cont_dir2 = cont_obj2.direction
    cont2_pos = cont_obj2.position
    matepos = alignedread.mpos
    cont2_len = cont_obj2.length
    s2len = scaf_obj2.s_length
    (obs1, obs2, scaf_side1, scaf_side2) = PosDirCalculatorPE(cont_dir1, read_dir, cont1_pos, readpos, s1len, cont1_len, cont_dir2, mate_dir, cont2_pos, matepos, s2len, cont2_len, param.read_len)
    if obs1 == counter.prev_obs1 and obs2 == counter.prev_obs2:
        counter.nr_of_duplicates += 1
        if param.detect_duplicate:
            return(1)

    if obs1 + obs2 < param.mean_ins_size + 6 * param.std_dev_ins_size and obs1 > 25 and obs2 > 25:# 2**32: #param.ins_size_threshold: 
        if scaf_side1 == 'R':
            scaf_obj1.lower_right_nbrs_obs[(scaf_obj2.name, scaf_side2)] = obs1 if obs1 < scaf_obj1.lower_right_nbrs_obs[(scaf_obj2.name, scaf_side2)] and scaf_obj1.lower_right_nbrs_obs[(scaf_obj2.name, scaf_side2)] > 0 else scaf_obj1.lower_right_nbrs_obs[(scaf_obj2.name, scaf_side2)]
            scaf_obj1.upper_right_nbrs_obs[(scaf_obj2.name, scaf_side2)] = obs1 if obs1 > scaf_obj1.upper_right_nbrs_obs[(scaf_obj2.name, scaf_side2)] else scaf_obj1.upper_right_nbrs_obs[(scaf_obj2.name, scaf_side2)]
        if scaf_side1 == 'L':
            scaf_obj1.lower_left_nbrs_obs[(scaf_obj2.name, scaf_side2)] = obs1 if obs1 < scaf_obj1.lower_left_nbrs_obs[(scaf_obj2.name, scaf_side2)] and scaf_obj1.lower_left_nbrs_obs[(scaf_obj2.name, scaf_side2)] > 0 else scaf_obj1.lower_left_nbrs_obs[(scaf_obj2.name, scaf_side2)]
            scaf_obj1.upper_left_nbrs_obs[(scaf_obj2.name, scaf_side2)] = obs1 if obs1 > scaf_obj1.upper_left_nbrs_obs[(scaf_obj2.name, scaf_side2)] else scaf_obj1.upper_left_nbrs_obs[(scaf_obj2.name, scaf_side2)]
        if scaf_side2 == 'R':
            scaf_obj2.lower_right_nbrs_obs[(scaf_obj1.name, scaf_side1)] = obs2 if obs2 < scaf_obj2.lower_right_nbrs_obs[(scaf_obj1.name, scaf_side1)] and scaf_obj2.lower_right_nbrs_obs[(scaf_obj1.name, scaf_side1)] > 0 else scaf_obj2.lower_right_nbrs_obs[(scaf_obj1.name, scaf_side1)]
            scaf_obj2.upper_right_nbrs_obs[(scaf_obj1.name, scaf_side1)] = obs2 if obs2 > scaf_obj2.upper_right_nbrs_obs[(scaf_obj1.name, scaf_side1)] else scaf_obj2.upper_right_nbrs_obs[(scaf_obj1.name, scaf_side1)]
        if scaf_side2 == 'L':
            scaf_obj2.lower_left_nbrs_obs[(scaf_obj1.name, scaf_side1)] = obs2 if obs2 < scaf_obj2.lower_left_nbrs_obs[(scaf_obj1.name, scaf_side1)] and scaf_obj2.lower_left_nbrs_obs[(scaf_obj1.name, scaf_side1)] > 0 else scaf_obj2.lower_left_nbrs_obs[(scaf_obj1.name, scaf_side1)]
            scaf_obj2.upper_left_nbrs_obs[(scaf_obj1.name, scaf_side1)] = obs2 if obs2 > scaf_obj2.upper_left_nbrs_obs[(scaf_obj1.name, scaf_side1)] else scaf_obj2.upper_left_nbrs_obs[(scaf_obj1.name, scaf_side1)]
        if (scaf_obj2.name, scaf_side2) not in G[(scaf_obj1.name, scaf_side1)]:
            G.add_edge((scaf_obj2.name, scaf_side2), (scaf_obj1.name, scaf_side1), nr_links=1, obs=obs1 + obs2)
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['obs_sq'] = (obs1 + obs2) ** 2
        else:
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['nr_links'] += 1
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['obs'] += obs1 + obs2
            #G.edge[(scaf_obj1.name,scaf_side1)][(scaf_obj2.name,scaf_side2)]['obs'] += obs1+obs2
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['obs_sq'] += (obs1 + obs2) ** 2
    else:
        counter.reads_with_too_long_insert += 1

        ## add to haplotype graph here!!

    counter.prev_obs1 = obs1
    counter.prev_obs2 = obs2
    return(0)


def CalculateStats(sorted_contig_lengths, sorted_contig_lengths_small, param, Information):
    cur_length = 0
    nr_conts = 0
    LG50 = 0
    NG50 = 0
    for contig_length in sorted_contig_lengths:
        cur_length += contig_length
        nr_conts += 1
        if cur_length >= param.tot_assembly_length / 2.0:
            LG50 = contig_length
            NG50 = nr_conts
            break
    if LG50 == 0:
        for contig_length in sorted_contig_lengths_small:
            cur_length += contig_length
            nr_conts += 1
            if cur_length >= param.tot_assembly_length / 2.0:
                LG50 = contig_length
                NG50 = nr_conts
                break
    print >> Information, 'LG50: ', LG50, 'NG50: ', NG50, 'Initial contig assembly length: ', param.tot_assembly_length
    return(NG50, LG50)

def CalculateMeanCoverage(Contigs, first_lib, output_dest, bamfile, Information):
    # tuples like (cont lenght, contig name)
    list_of_cont_tuples = [(Contigs[contig].length, contig) for contig in Contigs]
    #sorted as longest first
    list_of_cont_tuples = sorted(list_of_cont_tuples, key=lambda tuple: tuple[0], reverse=True)
    #coverages of longest contigs
    longest_contigs = list_of_cont_tuples[:1000]
    cov_of_longest_contigs = [Contigs[contig[1]].coverage for contig in longest_contigs]
    #Calculate mean coverage from the 1000 longest contigs
    n = float(len(cov_of_longest_contigs))
    mean_cov = sum(cov_of_longest_contigs) / n
    std_dev = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_cov + mean_cov ** 2), cov_of_longest_contigs))) / (n - 1)) ** 0.5
    extreme_obs_occur = True
    print >> Information, 'Mean coverage before filtering out extreme observations = ', mean_cov
    print >> Information, 'Std dev of coverage before filtering out extreme observations= ', std_dev

    ## SMOOTH OUT THE MEAN HERE by removing extreme observations ## 
    while extreme_obs_occur:
        extreme_obs_occur, filtered_list = RemoveOutliers(mean_cov, std_dev, cov_of_longest_contigs)
        n = float(len(filtered_list))
        try:
            mean_cov = sum(filtered_list) / n
        except ZeroDivisionError:
            break
        std_dev = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_cov + mean_cov ** 2), filtered_list))) / (n - 1)) ** 0.5
        cov_of_longest_contigs = filtered_list

    print >> Information, 'Mean coverage after filtering = ', mean_cov
    print >> Information, 'Std coverage after filtering = ', std_dev
    print >> Information, 'Length of longest contig in calc of coverage: ', longest_contigs[0][0]
    print >> Information, 'Length of shortest contig in calc of coverage: ', longest_contigs[-1][0]

    try:
        import matplotlib
        plt.hist(cov_of_longest_contigs, bins=50)
        library = bamfile.split('/')[-1]
        plt.savefig(output_dest + "/BESST_cov_1000_longest_cont" + library + ".png")
        plt.clf()
    except ImportError:
        pass
    return(mean_cov, std_dev)

def RemoveOutliers(mean_cov, std_dev, cov_list):
    k = MaxObsDistr(len(cov_list))
    filtered_list = list(filter((lambda x : (x < mean_cov + k * std_dev and x < 2 * mean_cov)), cov_list))
    if len(cov_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)

def RepeatDetector(Contigs, Scaffolds, G, param, G_prime, small_contigs, small_scaffolds, Information):
    output_dest = param.output_directory
    mean_cov = param.mean_coverage
    std_dev = param.std_dev_coverage
    cov_cutoff = param.cov_cutoff
    Repeats = []
    count_repeats = 0
    count_hapl = 0
    nr_of_contigs = len(Contigs)
    k = MaxObsDistr(nr_of_contigs)
    repeat_thresh = max(mean_cov + k * std_dev, 2 * mean_cov - std_dev)
    if cov_cutoff:
        repeat_thresh = cov_cutoff
    print >> Information, 'Detecting repeats..'
    for contig in Contigs:
        if Contigs[contig].coverage > repeat_thresh:
            count_repeats += 1
            cont_obj_ref = Contigs[contig]
            Repeats.append(cont_obj_ref)
            scaf_ = Contigs[contig].scaffold
            del Scaffolds[scaf_]
            G.remove_nodes_from([(scaf_, 'L'), (scaf_, 'R')])
            if param.extend_paths:
                G_prime.remove_nodes_from([(scaf_, 'L'), (scaf_, 'R')])
        #TEST DETECTOR
        if param.detect_haplotype and Contigs[contig].coverage < mean_cov / 2.0 + param.hapl_threshold * std_dev: # < mean_cov - 2.5*std_dev:
            count_hapl += 1
            #have indicator of possible haplotype
            Contigs[contig].is_haplotype = True


#    if param.extend_paths:
    for contig in small_contigs:
        if small_contigs[contig].coverage > repeat_thresh:
            count_repeats += 1
            cont_obj_ref = small_contigs[contig]
            Repeats.append(cont_obj_ref)
            scaf_ = small_contigs[contig].scaffold
            del small_scaffolds[scaf_]
            G_prime.remove_nodes_from([(scaf_, 'L'), (scaf_, 'R')])
        #TEST DETECTOR
        if param.detect_haplotype and small_contigs[contig].coverage < mean_cov / 2.0 + param.hapl_threshold * std_dev: # < mean_cov - 2.5*std_dev:
            count_hapl += 1
            #have indicator of possible haplotype
            small_contigs[contig].is_haplotype = True


    GO.PrintOutRepeats(Repeats, Contigs, output_dest, small_contigs)
    print >> Information, 'Removed a total of: ', count_repeats, ' repeats.'
    if param.detect_haplotype:
        print >> Information, 'Marked a total of: ', count_hapl, ' potential haplotypes.'
    return(Contigs, Scaffolds, G)

def MaxObsDistr(nr_of_obs):
    # Here, we choose the quantile that separates "normal" contigs from repeats. We want to find 
    # the k s.t. the probability of marking one or more normal contigs (from the hole set of contigs) 
    # as repeats (false positives) is less than p=0.05. This is equivalent to: Choosing p in Bin(n,p)
    # (where n = nr of contigs) s.t. P(k=0)>= 0.95 (no successes if a success is a false positive).
    # We get P(k=0)= choose(n,k)*p**n => p = 1 - (0.95*n)**n. With this value of p, if X~N(mean,sigma),
    # we want to derive k from P_x(x < mean + k*sigma) = 1-p. This gives the k that is returned from this function. 
    #from scipy.stats import norm
    import math
    def rational_approximation(t):

        # Abramowitz and Stegun formula 26.2.23.
        # The absolute value of the error should be less than 4.5 e-4.
        c = [2.515517, 0.802853, 0.010328]
        d = [1.432788, 0.189269, 0.001308]
        numerator = (c[2] * t + c[1]) * t + c[0]
        denominator = ((d[2] * t + d[1]) * t + d[0]) * t + 1.0
        return t - numerator / denominator

    def normal_CDF_inverse(p):

        assert p > 0.0 and p < 1

        # See article above for explanation of this section.
        if p < 0.5:
            # F^-1(p) = - G^-1(p)
            return -rational_approximation(math.sqrt(-2.0 * math.log(p)))
        else:
            # F^-1(p) = G^-1(1-p)
            return rational_approximation(math.sqrt(-2.0 * math.log(1.0 - p)))

    p = 1 - (0.95) ** (1 / float(nr_of_obs))
    k = normal_CDF_inverse(1 - p)

    return(k)


def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list))
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)

def GetParams(bam_file, param, Scaffolds, small_scaffolds, Contigs, small_contigs, Information):
    import sys
    informative_pair = set([147, 163]) #161,145,129,177,
    cont_names = bam_file.references
    cont_lengths = bam_file.lengths
    #cont_lengths=[int(nr) for nr in cont_lengths]  #convert long to int object
    cont_lengths_list = list(cont_lengths)
    indexes = [i for i in range(0, len(cont_lengths_list))]
    from heapq import nlargest
    largest_contigs_indexes = nlargest(1000, indexes, key=lambda i: cont_lengths_list[i]) #get indexes of the 1000 longest contigs

    if not param.read_len: # user has not specified read len  
        #get read length
        try:
            iter = bam_file.fetch(cont_names[largest_contigs_indexes[0]])
        except ValueError:
            sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
            sys.exit(0)
        nr_reads = 0
        tot_read_len = 0
        for read in iter:
            if read.rlen != 0:
                tot_read_len += read.rlen
                nr_reads += 1
            else:
                tot_read_len += read.alen
                nr_reads += 1
        param.read_len = tot_read_len / float(nr_reads)


    if param.mean_ins_size and param.std_dev_ins_size and not param.ins_size_threshold: # user has specified mean and std dev but no thresholds
        param.ins_size_threshold = param.mean_ins_size + 4 * param.std_dev_ins_size
        if param.extend_paths:
            param.contig_threshold = param.ins_size_threshold
        else:
            param.contig_threshold = param.mean_ins_size + (param.std_dev_ins_size / float(param.mean_ins_size)) * param.std_dev_ins_size
        print >> Information, '-T', param.ins_size_threshold, '-t', param.contig_threshold

    if not param.mean_ins_size: # user has not specified mean and std dev. (and no thresholds)
        counter = 1
        ins_size_reads = []
        for index in largest_contigs_indexes:
            try:
                iter = bam_file.fetch(cont_names[index])
            except ValueError:
                sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
                sys.exit(0)
            for read in iter:
                if read.flag in informative_pair and read.rname == read.mrnm  and not read.mate_is_unmapped and not read.is_unmapped and read.mapq > 10: ##read.tid == read.rnext and not read.mate_is_unmapped and not read.is_unmapped: #
                    ins_size_reads.append(abs(read.tlen))
                    counter += 1
                if counter > 1000000:
                    break
            if counter > 1000000:
                break

        #get mean and std dev here. 
        #Assure that there were enough reads  for computation of mean and variance
        if len(ins_size_reads) <= 1000:
            sys.stderr.write('To few valid read alignments exists to compute mean and variance of library (need at least 1000 observations). Got only ' + str(len(ins_size_reads)) + ' valid alignments. Please specify -m and -s to the program. \nPrinting out scaffolds produced in earlier steps...')
            sys.stderr.write('\nterminating...\n')
            sys.exit(0)

        ## SMOOTH OUT THE MEAN HERE by removing extreme observations## 
        n = float(len(ins_size_reads))
        mean_isize = sum(ins_size_reads) / n
        std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5
        print >> Information, 'Mean before filtering :', mean_isize
        print >> Information, 'Std_est  before filtering: ', std_dev_isize
        extreme_obs_occur = True
        while extreme_obs_occur:
            extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, ins_size_reads)
            n = float(len(filtered_list))
            mean_isize = sum(filtered_list) / n
            std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n - 1)) ** 0.5
            ins_size_reads = filtered_list

        n = float(len(ins_size_reads))
        mean_isize = sum(ins_size_reads) / n
        std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5

        print >> Information, 'Mean converged:', mean_isize
        print >> Information, 'Std_est converged: ', std_dev_isize

        param.mean_ins_size = mean_isize
        param.std_dev_ins_size = std_dev_isize
    if not param.ins_size_threshold:
        param.ins_size_threshold = param.mean_ins_size + 4 * param.std_dev_ins_size
        if param.extend_paths:
            param.contig_threshold = param.ins_size_threshold
        else:
            param.contig_threshold = param.mean_ins_size + (param.std_dev_ins_size / float(param.mean_ins_size)) * param.std_dev_ins_size

    print >> Information, ''
    print >> Information, 'Mean of library set to:', param.mean_ins_size
    print >> Information, 'Standard deviation of library set to: ', param.std_dev_ins_size
    print >> Information, '-T (library insert size threshold) set to: ', param.ins_size_threshold
    print >> Information, '-k set to (Scaffolding with contigs larger than): ', param.contig_threshold
    print >> Information, 'Number of links required to create an edge: ', param.edgesupport
    print >> Information, 'Read length set to: ', param.read_len
    print >> Information, 'Relative weight of dominating link set to (default=3): ', param.rel_weight
    print >> Information, ''
    return()

def PosDirCalculatorPE(cont_dir1, read_dir, cont1pos, readpos, s1len, cont1_len, cont_dir2, mate_dir, cont2pos, matepos, s2len, cont2_len, read_len):
    if cont_dir1 and read_dir:
        obs1 = s1len - cont1pos - readpos
        read_side1 = 'R'
    if cont_dir2 and mate_dir:
        obs2 = s2len - cont2pos - matepos
        read_side2 = 'R'
    if (not cont_dir1) and read_dir:
        obs1 = cont1pos + (cont1_len - readpos)
        read_side1 = 'L'
    if (not cont_dir2) and mate_dir:
        obs2 = cont2pos + (cont2_len - matepos)
        read_side2 = 'L'
    if cont_dir1 and not read_dir:
        obs1 = cont1pos + readpos + read_len
        read_side1 = 'L'
    if cont_dir2 and not mate_dir:
        obs2 = cont2pos + matepos + read_len
        read_side2 = 'L'
    if not cont_dir1 and not read_dir:
        obs1 = s1len - cont1pos - (cont1_len - readpos - read_len)
        read_side1 = 'R'
    if not cont_dir2 and not mate_dir:
        obs2 = s2len - cont2pos - (cont2_len - matepos - read_len)
        read_side2 = 'R'

    if read_side1 == 'L':
        scaf_side1 = 'L'
    if read_side2 == 'L':
        scaf_side2 = 'L'
    if read_side1 == 'R':
        scaf_side1 = 'R'
    if read_side2 == 'R':
        scaf_side2 = 'R'
    return(int(obs1), int(obs2), scaf_side1, scaf_side2)








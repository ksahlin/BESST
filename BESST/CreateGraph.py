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
import os
from collections import defaultdict, Counter
from math import pi, sqrt
from time import time

from scipy.stats import ks_2samp
import networkx as nx

import Contig
import Scaffold
from Parameter import counters
from mathstats.normaldist import normal
import GenerateOutput as GO
from mathstats.normaldist.truncatedskewed import param_est
import mathstats.log_normal_param_est as lnpe
import errorhandle
import plots
import e_nr_links




def PE(Contigs, Scaffolds, Information, C_dict, param, small_contigs, small_scaffolds, bam_file):
    G = nx.Graph()
    G_prime = nx.Graph()  # If we want to do path extension with small contigs
    print >> Information, 'Parsing BAM file...'

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

    if len(Scaffolds) == 0:
        if not os.path.isfile(param.output_directory + '/repeats.fa'):
            repeat_file = open(param.output_directory + '/repeats.fa', 'w')
        return(G, G_prime)

    ### initialize graph objects two nodes per contig "left" and "right" node. ###    
    tot_start = time()

    if param.no_score:
        # Only do path search
        #small contig graph contains all scaffolds
        InitializeGraph(small_scaffolds, G_prime, Information)
        InitializeGraph(Scaffolds, G_prime, Information)        
    elif param.extend_paths:
        # do scoring and extend paths
        InitializeGraph(Scaffolds, G, Information)
        #small contig graph contains all scaffolds
        InitializeGraph(small_scaffolds, G_prime, Information)
        InitializeGraph(Scaffolds, G_prime, Information)
    else:
        # only do scoring
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

    if param.development:
        import guppy
        h = guppy.hpy()
        before_ctg_graph = h.heap()
        print 'Just before creating contig graph:\n{0}\n'.format(before_ctg_graph)

    for alignedread in bam_file:
        # try: #check that read is aligned OBS: not with is_unmapped since this flag is fishy for e.g. BWA
        #     contig1 = bam_file.getrname(alignedread.rname)
        #     contig2 = bam_file.getrname(alignedread.mrnm)
        # except ValueError:
        #     continue
        
        try:
            contig1 = param.contig_index[alignedread.rname]
            contig2 = param.contig_index[alignedread.mrnm]
            #assert contig1 == contig1_old
            #assert contig2 == contig2_old
        except KeyError:
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
                (side1, side2) = CheckDir(cont_obj1, cont_obj2, alignedread,param)
                #get scaffold name for contig
                s1 = Contigs[contig1].scaffold if contig1 in Contigs else small_contigs[contig1].scaffold
                s2 = Contigs[contig2].scaffold if contig2 in Contigs else small_contigs[contig2].scaffold
                fishy_edges[((s1, side1), (s2, side2))] += 1
                fishy_edges[((s2, side2), (s1, side1))] += 1
                ctr += 1

        ## add to coverage computation if contig is still in the list of considered contigs
        #print contig1, contig2, alignedread.is_read2
        cont_aligned_len[contig1][0] += alignedread.qlen
        if contig1 != contig2 and alignedread.mapq == 0:
            counter.non_unique += 1  # check how many non unique reads out of the useful ones (mapping to two different contigs)

        if contig1 != contig2 and alignedread.is_read2 and not alignedread.is_unmapped and alignedread.mapq >= param.min_mapq:
            if contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold:
                cont_obj1 = Contigs[contig1]
                cont_obj2 = Contigs[contig2]
                scaf_obj1 = Scaffolds[cont_obj1.scaffold]
                scaf_obj2 = Scaffolds[cont_obj2.scaffold]
                if not param.no_score:
                    is_dupl = CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G, param, alignedread, counter, contig1, contig2)
                else:
                    is_dupl = CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2, save_obs=False)

                if param.extend_paths and not is_dupl and not param.no_score:
                    counter.prev_obs1 = -1
                    counter.prev_obs2 = -1
                    CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2, save_obs=False)
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
                    CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2, save_obs=False)
                elif value1 and not value2:
                    CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2,save_obs=False)
                elif not value1 and value2:
                    CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G_prime, param, alignedread, counter, contig1, contig2,save_obs=False)


            elif contig1 in Contigs and contig2 in Contigs and Contigs[contig2].scaffold != Contigs[contig1].scaffold:
########################Use to validate scaffold in previous step here ############
                pass
    print >> Information, 'ELAPSED reading file:', time() - staart
    print >> Information, 'NR OF FISHY READ LINKS: ', ctr

    print >> Information, 'Number of USEFUL READS (reads mapping to different contigs uniquly): ', counter.count
    print >> Information, 'Number of non unique reads (at least one read non-unique in read pair) that maps to different contigs (filtered out from scaffolding): ', counter.non_unique
    print >> Information, 'Reads with too large insert size from "USEFUL READS" (filtered out): ', counter.reads_with_too_long_insert
    print >> Information, 'Initial number of edges in G (the graph with large contigs): ', len(G.edges())
    print >> Information, 'Initial number of edges in G_prime (the full graph of all contigs before removal of repats): ', len(G_prime.edges())
    
    if param.detect_duplicate:
        print >> Information, 'Number of duplicated reads indicated and removed: ', counter.nr_of_duplicates

    if param.development:
        h = guppy.hpy()
        after_ctg_graph = h.heap()
        print 'Just after creating contig graph:\n{0}\n'.format(after_ctg_graph)


##### Calc coverage for all contigs with current lib here #####
    sum_x = 0
    sum_x_sq = 0
    n = 0
    cov = []
    leng = []
    for contig in cont_aligned_len:
        cont_coverage = cont_aligned_len[contig][0] / float(cont_aligned_len[contig][1])
        try:
            Contigs[contig].coverage = cont_coverage
            cov.append(cont_coverage)
            leng.append(Contigs[contig].length)
        except KeyError:
            small_contigs[contig].coverage = cont_coverage
            cov.append(cont_coverage)
            leng.append(small_contigs[contig].length)


        sum_x += cont_coverage
        sum_x_sq += cont_coverage ** 2
        n += 1

    del cont_aligned_len



    mean_cov, std_dev_cov = CalculateMeanCoverage(Contigs, Information, param)
    param.mean_coverage = mean_cov
    param.std_dev_coverage = std_dev_cov
    if param.first_lib:
        Contigs, Scaffolds, G = RepeatDetector(Contigs, Scaffolds, G, param, G_prime, small_contigs, small_scaffolds, Information)
        if param.lower_cov_cutoff:
            filter_low_coverage_contigs(Contigs, Scaffolds, G, param, G_prime, small_contigs, small_scaffolds, Information)

    print >> Information, 'Number of edges in G (after repeat removal): ', len(G.edges())
    print >> Information, 'Number of edges in G_prime (after repeat removal): ', len(G_prime.edges())

    ### Remove edges created by false reporting of BWA ###
    RemoveBugEdges(G, G_prime, fishy_edges, param, Information)

    print >> Information, 'Number of edges in G (after filtering for buggy flag stats reporting): ', len(G.edges())
    print >> Information, 'Number of edges in G_prime  (after filtering for buggy flag stats reporting): ', len(G_prime.edges())

    if param.development:
        h = guppy.hpy()
        after_coverage = h.heap()
        print 'Just after calc coverage:\n{0}\n'.format(after_coverage)


    #Remove all edges with link support less than inferred parameter of param.edgesupport (unless specified by user) to be able to compute statistics.
    infer_spurious_link_count_threshold(G_prime, param)
    if not param.edgesupport:
        # default value of 3
        param.edgesupport = 3 
        print >> Information, 'Letting -e be {0} for this library.'.format(param.edgesupport)
    else:
        print >> Information, 'User has set -e to be {0} for this library.'.format(param.edgesupport)

    counter_low_support = 0
    for edge in G.edges():
        if G[edge[0]][edge[1]]['nr_links'] != None:
            if G[edge[0]][edge[1]]['nr_links'] < 3:
                G.remove_edge(edge[0], edge[1])
                counter_low_support += 1

    # counter_low_support = remove_edges_below_threshold(G, param)

    print >> Information, 'Removed {0} edges from graph G of border contigs.'.format(counter_low_support)
    remove_edges_below_threshold(G_prime, param)


    ## Score edges in graph
    plot = 'G'
    if not param.no_score:
        GiveScoreOnEdges(G, Scaffolds, small_scaffolds, Contigs, param, Information, plot)
    #plot = 'G_prime'
    #GiveScoreOnEdges(G_prime, Scaffolds, small_scaffolds, Contigs, param, Information, plot)

    print >> Information, 'Number of edges in G_prime  (after removing edges under -e threshold (if not specified, default is -e 3): ', len(G_prime.edges())

    if param.development:
        h = guppy.hpy()
        after_scoring = h.heap()
        print 'After scoring:\n{0}\n'.format(after_scoring)

    print >> Information, "\n -------------------------------------------------------------\n"
    print >> Information, 'Nr of contigs/scaffolds included in this pass: ' + str(len(Scaffolds) + len(small_scaffolds))
    print >> Information, 'Out of which {0} acts as border contigs.'.format(len(Scaffolds))
    return(G, G_prime)

def infer_spurious_link_count_threshold(G_prime, param):
    # Finding a base value estimation  for the cutoff be expected number of links that spans the mean (from a normal distribution)
    # this is probably an underestimation for a positive skev library and definitely
    # an underestimation for a negative skew library, however it's better with an undersetimation rather
    # than removing to many edges..

    nr_nodes = nx.number_of_nodes(G_prime)/2
    contamination_ratio = param.contamination_ratio if param.contamination_ratio else 0
    cov = param.mean_coverage *(1 - contamination_ratio )
    link_params = e_nr_links.Param(param.mean_ins_size, param.std_dev_ins_size, cov, param.read_len, 0)
    gap = param.mean_ins_size +  param.std_dev_ins_size - 2* param.read_len
    # print cov, gap, param.mean_ins_size, param.std_dev_ins_size, cov, param.read_len, 0
    expected_links_over_mean_plus_stddev = e_nr_links.ExpectedLinks(100000, 100000, gap, link_params)
    link_stats = []
    for edge in G_prime.edges():
        if G_prime[edge[0]][edge[1]]['nr_links'] != None:
            link_count = G_prime[edge[0]][edge[1]]['nr_links']
            link_stats.append(link_count)
    link_counter = Counter(link_stats)
    total_included_edges = 0
    for link_number in sorted(link_counter.keys(), reverse=True):
        # print >> param.information_file, 'Number of links: {0}\tFound in {1} edges.'.format(link_number, link_counter[link_number])
        total_included_edges += link_counter[link_number]
        print >> param.information_file, 'Nodes: {0}.\t Total edges with over {1} links:{2}. \tAverage density: {3}'.format(nr_nodes, link_number, total_included_edges, total_included_edges/float(nr_nodes))

    # print "expected_links_over_mean_plus_stddev:", expected_links_over_mean_plus_stddev
    if expected_links_over_mean_plus_stddev < 3:
        param.expected_links_over_mean_plus_stddev = 3
    else:
        param.expected_links_over_mean_plus_stddev = int(expected_links_over_mean_plus_stddev)
    print >> param.information_file, 'Letting filtering threshold in high complexity regions be {0} for this library.'.format(param.expected_links_over_mean_plus_stddev)

def remove_edges_below_threshold(graph, param):
    #Remove all edges with link support less than param.edgesupport in high complexity areas

    print >> param.information_file, 'Remove edges in high complexity areas.'
    counter = 0
    nr_nbrs_counter = []
    # print "NR EDGES", len(graph.edges())
    for node1, node2 in sorted(graph.edges(), key=lambda x: graph[x[0]][x[1]]["nr_links"]):
        if  graph[node1][node2]["nr_links"] == None:
            continue
        elif graph[node1][node2]["nr_links"] < param.expected_links_over_mean_plus_stddev:
            # print graph[node1][node2]["nr_links"]
            nbrs1 = graph.neighbors(node1)
            nbrs2 = graph.neighbors(node2)
            if len(nbrs1) > 4 and len(nbrs2) > 4:
                graph.remove_edge(node1, node2)
                # nr_nbrs_counter.append(len(nbrs))
                counter += 1
        else:
            break




    # for node in graph.nodes():
    #     nbrs = graph.neighbors(node)
    #     if len(nbrs) > 3:
    #         nr_nbrs_counter.append(len(nbrs))
    #         for nbr in nbrs:
    #             if node[0] != nbr[0] and graph[node][nbr]['nr_links'] < param.edgesupport:
    #                 graph.remove_edge(node, nbr)
    #                 counter += 1

    # density_nbr_counter = Counter(nr_nbrs_counter)
    print >> param.information_file, 'Removed total of {0} edges in high density areas.'.format(counter)
    # print >> param.information_file, "DENSITY PROFILE:"
    # for nbr_number in sorted(density_nbr_counter.keys(), reverse=True):
    #     print >> param.information_file, 'Number of nbrs: {0}\t found in {1} edges.'.format(nbr_number, density_nbr_counter[nbr_number])


    # remove edges below 3 links (independent of region complexity) to be able to compute statistics
    counter_low_support = 0
    for edge in graph.edges():
        if graph[edge[0]][edge[1]]['nr_links'] != None:
            if graph[edge[0]][edge[1]]['nr_links'] < param.edgesupport:
                graph.remove_edge(edge[0], edge[1])
                counter_low_support += 1
    print >> param.information_file, 'Removed an additional of {0} edges with low support from full graph G_prime of all contigs.'.format(counter_low_support)


def filter_low_coverage_contigs(Contigs, Scaffolds, G, param, G_prime, small_contigs, small_scaffolds, Information):
    low_coverage_contigs = []
    count_low_coverage_contigs = 0
    print >> Information, 'Removing low coverage contigs if -z_min specified..'
    for contig in Contigs:
        if Contigs[contig].coverage < param.lower_cov_cutoff:
            count_low_coverage_contigs += 1
            cont_obj_ref = Contigs[contig]
            low_coverage_contigs.append(cont_obj_ref)
            scaf_ = Contigs[contig].scaffold
            del Scaffolds[scaf_]
            G.remove_nodes_from([(scaf_, 'L'), (scaf_, 'R')])
            if param.extend_paths:
                G_prime.remove_nodes_from([(scaf_, 'L'), (scaf_, 'R')])

    for contig in small_contigs:
        if small_contigs[contig].coverage <  param.lower_cov_cutoff:
            count_low_coverage_contigs += 1
            cont_obj_ref = small_contigs[contig]
            low_coverage_contigs.append(cont_obj_ref)
            scaf_ = small_contigs[contig].scaffold
            del small_scaffolds[scaf_]
            G_prime.remove_nodes_from([(scaf_, 'L'), (scaf_, 'R')])


    GO.PrintOut_low_cowerage_contigs(low_coverage_contigs, Contigs, param.output_directory, small_contigs)
    print >> Information, 'Removed a total of: ', count_low_coverage_contigs, ' low coverage contigs. With coverage lower than ', param.lower_cov_cutoff




def get_conditional_stddevs(steps, empirical_isize_distr, max_isize):
    # create list of contiditional stddevs
    expected_stddev = []

    previous_gap = 0
    for gap in steps:
        # calculate conditional distribution given a gap
        conditional_density = [0]*(max_isize+1)

        for x in empirical_isize_distr.keys():
            w_x = max(0, x-gap+1)
            f_x = empirical_isize_distr[x]
            if w_x > 0:
                conditional_density[x] = f_x * w_x

        # get conditional stddev from distribution
        tot_density = float(sum(conditional_density))
        conditional_mu = sum(map(lambda (i, f_x): i*f_x, enumerate(conditional_density)))/tot_density
        conditional_sigma = sqrt(sum(map(lambda (i, f_x): (i-conditional_mu)**2 * f_x, enumerate(conditional_density))) / tot_density)
        

        # add to our conditional stddev vector
        if gap == 0:
            expected_stddev.append([conditional_sigma])
        else:
            expected_stddev.append([conditional_sigma]*(gap - previous_gap))
        previous_gap = gap
        # print gap, conditional_mu, conditional_sigma
    # print expected_stddev
    conditional_stddevs = [item for sublist in expected_stddev for item in sublist] # flatten list

    return conditional_stddevs



def GiveScoreOnEdges(G, Scaffolds, small_scaffolds, Contigs, param, Information, plot):

    span_score_obs = []
    std_dev_score_obs = []
    gap_obs = []
    nr_link_obs = []
    cnt_sign = 0

    if param.print_scores:
        score_file = open(os.path.join(param.output_directory,"score_file_pass_{0}.tsv".format(param.pass_number)), "w")
        print >>score_file, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format("scf1/ctg1", "o1", "scf2/ctg2", "o2", "gap", "link_variation_score", "link_dispersity_score", "number_of_links")

    if param.lognormal:
        emp_distr = param.empirical_distribution
        #print emp_distr
        isize_range = sorted(emp_distr.keys())
        max_isize = isize_range[-1]
        steps = range(0, int(max_isize*0.8), max_isize/50)  # isize_range[:int(0.8*len(isize_range)):len(isize_range)/50] # not the upper most values
        #print steps
        conditional_stddevs = get_conditional_stddevs(steps, emp_distr, max_isize)
        # print conditional_stddevs


    for edge in G.edges():
        mean_ = 0
        std_dev = 0
        if G[edge[0]][edge[1]]['nr_links'] != None:

            n = G[edge[0]][edge[1]]['nr_links']
            obs_squ = G[edge[0]][edge[1]]['obs_sq']
            mean_ = G[edge[0]][edge[1]]['obs'] / float(n)
            samples = G[edge[0]][edge[1]]['observations']
            #print samples
            # mu = sum(samples)/len(samples)
            # sigma = (sum(list(map((lambda x: x ** 2 - 2 * x * mu + mu ** 2), samples))) / (n - 1)) ** 0.5

            data_observation = (n * param.mean_ins_size - G[edge[0]][edge[1]]['obs']) / float(n)
            try:
                len1 = Scaffolds[ edge[0][0] ].s_length
            except KeyError:
                len1 = small_scaffolds[ edge[0][0] ].s_length
            try:
                len2 = Scaffolds[ edge[1][0] ].s_length
            except KeyError:
                len2 = small_scaffolds[ edge[1][0] ].s_length


            # Distribution has skew
            if param.lognormal:
                if 2 * param.std_dev_ins_size < len1 and 2 * param.std_dev_ins_size < len2:
                    samples = G[edge[0]][edge[1]]['observations']
                    gap = lnpe.GapEstimator(param.lognormal_mean, param.lognormal_sigma, param.read_len, samples, len1, c2_len=len2)
                    # print 'GAP', gap
                else:
                    gap = data_observation

            # Normal-ish distribution
            else:

                if 2 * param.std_dev_ins_size < len1 and 2 * param.std_dev_ins_size < len2:
                    gap = param_est.GapEstimator(param.mean_ins_size, param.std_dev_ins_size, param.read_len, mean_, len1, len2)
                else:
                    gap = data_observation

            G[edge[0]][edge[1]]['gap'] = int(gap)
            if -gap > len1 or -gap > len2:
                G[edge[0]][edge[1]]['score'] = 0
                continue

            #std_dev_d_eq_0 = param_est.tr_sk_std_dev(param.mean_ins_size, param.std_dev_ins_size, param.read_len, len1, len2, gap)

            if 2 * param.std_dev_ins_size < len1 and 2 * param.std_dev_ins_size < len2:
                if param.lognormal:
                    if gap > 0:
                        std_dev_d_eq_0 = conditional_stddevs[int(gap)]
                    else:
                        std_dev_d_eq_0 = conditional_stddevs[0]
                else:
                    std_dev_d_eq_0 = param_est.tr_sk_std_dev(param.mean_ins_size, param.std_dev_ins_size, param.read_len, len1, len2, gap)

            else:
                std_dev_d_eq_0 = 2 ** 32

            try:
                std_dev = ((obs_squ - n * mean_ ** 2) / (n - 1)) ** 0.5
                #chi_sq = (n - 1) * (std_dev ** 2 / std_dev_d_eq_0 ** 2)
            except ZeroDivisionError:
                std_dev = 2 ** 32
                #chi_sq = 0

            # print std_dev_d_eq_0, std_dev #, 'New sigma:', sigma 
            try:
                l1 = G[edge[0]][edge[1]][Scaffolds[edge[0][0]].name]
                del G[edge[0]][edge[1]][Scaffolds[edge[0][0]].name]
            except KeyError:
                l1 = G[edge[0]][edge[1]][small_scaffolds[edge[0][0]].name]
                del G[edge[0]][edge[1]][small_scaffolds[edge[0][0]].name]
            try:
                l2 = G[edge[0]][edge[1]][Scaffolds[edge[1][0]].name]
                del G[edge[0]][edge[1]][Scaffolds[edge[1][0]].name]
            except KeyError:
                l2 = G[edge[0]][edge[1]][small_scaffolds[edge[1][0]].name]
                del G[edge[0]][edge[1]][small_scaffolds[edge[1][0]].name]


            l1.sort()
            n_obs = len(l1)
            l1_mean = sum(l1) / float(n_obs)
            #l1_median = l1[len(l1) / 2]
            l1 = map(lambda x: x - l1_mean, l1)
            #l1 = map(lambda x: x - l1_median, l1)
            max_obs2 = max(l2)
            l2.sort(reverse=True)
            l2 = map(lambda x: abs(x - max_obs2), l2)
            l2_mean = sum(l2) / float(n_obs)
            #l2_median = l2[len(l2) / 2]
            l2 = map(lambda x: x - l2_mean, l2)
            #l2 = map(lambda x: x - l2_median, l2)
            KS_statistic, p_value = ks_2samp(l1, l2)


            #M_W_statistic, p_val = mannwhitneyu(l1, l2)

            #diff = map(lambda x: abs(abs(x[1]) - abs(x[0])), zip(l1, l2))
            #sc = sum(diff) / len(diff)

            if len(l1) < 3:
                span_score = 0
            else:
                span_score = 1 - KS_statistic

            try:
                std_dev_score = min(std_dev / std_dev_d_eq_0, std_dev_d_eq_0 / std_dev) #+ span_score #+ min(n/E_links, E_links/float(n))
            except ZeroDivisionError:
                std_dev_score = 0
                sys.stderr.write(str(std_dev) + ' ' + str(std_dev_d_eq_0) + ' ' + str(span_score) + '\n')

            G[edge[0]][edge[1]]['score'] = std_dev_score + span_score if std_dev_score > 0.5 and span_score > 0.5 else 0

            if param.plots:
                span_score_obs.append(span_score)
                std_dev_score_obs.append(std_dev_score)
                gap_obs.append(gap)
                nr_link_obs.append(n_obs)

            if param.print_scores:
                if edge[0][1] == 'R':
                    contig_objects = Scaffolds[edge[0][0]].contigs
                    contig_names = map(lambda x: x.name, contig_objects)
                    scf1 = ";".join(contig_names)
                    contig_directions_bool = map(lambda x: x.direction, contig_objects)
                    contig_directions1 = ";".join(map(lambda x: '+' if True else '-', contig_directions_bool))

                else:
                    contig_objects = Scaffolds[edge[0][0]].contigs[::-1]
                    contig_names = map(lambda x: x.name, contig_objects)
                    scf1 =";".join(contig_names)
                    contig_directions_bool = map(lambda x: x.direction, contig_objects)
                    contig_directions1 = ";".join(map(lambda x: '+' if False else '-', contig_directions_bool))

                if edge[1][1] == 'L':
                    contig_objects = Scaffolds[edge[1][0]].contigs
                    contig_names = map(lambda x: x.name, contig_objects)
                    scf2 =";".join(contig_names)
                    contig_directions_bool = map(lambda x: x.direction, contig_objects)
                    contig_directions2 = ";".join(map(lambda x: '+' if True else '-', contig_directions_bool))

                else:
                    contig_objects = Scaffolds[edge[1][0]].contigs[::-1]
                    contig_names = map(lambda x: x.name, contig_objects)
                    scf2 =";".join(contig_names)
                    contig_directions_bool = map(lambda x: x.direction, contig_objects)
                    contig_directions2 = ";".join(map(lambda x: '+' if False else '-', contig_directions_bool))

                print >>score_file, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(scf1, contig_directions1, scf2, contig_directions2, gap, std_dev_score, span_score, n_obs)

    if param.print_scores:
        score_file.close()

    if param.plots:
        plots.histogram(span_score_obs, param, bins=20, x_label='score', y_label='frequency', title='Dispersity_score_distribuion' + plot + '.' + param.bamfile.split('/')[-1])
        plots.histogram(std_dev_score_obs, param, bins=20, x_label='score', y_label='frequency', title='Standard_deviation_score_distribuion' + plot + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(std_dev_score_obs, span_score_obs, param, x_label='std_dev_score_obs', y_label='span_score_obs', title='Score_correlation' + plot + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(std_dev_score_obs, gap_obs, param, x_label='std_dev_score_obs', y_label='estimated gap size', title='Gap_to_sigma' + plot + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(span_score_obs, gap_obs, param, x_label='span_score_obs', y_label='estimated gap size', title='Gap_to_span' + plot + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(span_score_obs, nr_link_obs, param, x_label='span_score_obs', y_label='Number links', title='Obs_to_span' + plot + '.' + param.bamfile.split('/')[-1])

    for edge in G.edges():
        if G[edge[0]][edge[1]]['nr_links'] != None:
            try:
                G[edge[0]][edge[1]]['score']
            except KeyError:
                sys.stderr.write(str(G[edge[0]][edge[1]]) + ' ' + str(Scaffolds[edge[0][0]].s_length) + ' ' + str(Scaffolds[edge[1][0]].s_length))
    print >> Information, 'Number of significantly spurious edges:', cnt_sign

    return()





def CheckDir(cont_obj1, cont_obj2, alignedread, param):
    (read_dir, mate_dir) = (not alignedread.is_reverse, not alignedread.mate_is_reverse)
    cont_dir1 = cont_obj1.direction  #if pos : L if neg: R
    #position2 cont2/scaf2                        
    cont_dir2 = cont_obj2.direction
    if param.orientation == 'fr':
        (obs1, obs2, scaf_side1, scaf_side2) = PosDirCalculatorPE(cont_dir1, read_dir, 0, 0, 0, 0, cont_dir2, mate_dir, 0, 0, 0, 0, 0)
    else:
        (obs1, obs2, scaf_side1, scaf_side2) = PosDirCalculatorMP(cont_dir1, read_dir, 0, 0, 0, 0, cont_dir2, mate_dir, 0, 0, 0, 0, 0)

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
    del fishy_edges
    return()

def InitializeGraph(dict_with_scaffolds, graph, Information):
    cnt = 1
    start1 = time()
    for scaffold_ in dict_with_scaffolds:
        graph.add_edge((scaffold_, 'L'), (scaffold_, 'R'), nr_links=None)    #this is a scaffold object but can be both a single contig or a scaffold.
        graph.node[(scaffold_, 'L')]['length'] = dict_with_scaffolds[scaffold_].s_length
        graph.node[(scaffold_, 'R')]['length'] = dict_with_scaffolds[scaffold_].s_length
        if cnt % 100000 == 0 and cnt > 0:
            elapsed = time() - start1
            print >> Information, 'Total nr of keys added: ', cnt, 'Time for adding last 100 000 keys: ', elapsed
            start1 = time()
        cnt += 1
    return()

# def constant_large():
#     return 2 ** 32
# def constant_small():
#     return -1

def InitializeObjects(bam_file, Contigs, Scaffolds, param, Information, G_prime, small_contigs, small_scaffolds, C_dict):
    singeled_out = 0
    contig_threshold = param.contig_threshold
    cont_lengths = bam_file.lengths
    cont_lengths = [int(nr) for nr in cont_lengths]  #convert long to int object
    cont_names = bam_file.references

    #Calculate NG50 and LG 50
    param.tot_assembly_length = sum(cont_lengths)
    sorted_lengths = sorted(cont_lengths, reverse=True)
    N50, L50 = CalculateStats(sorted_lengths, [], param, Information)
    param.current_L50 = L50
    param.current_N50 = N50
    #extend_paths = param.extend_paths
    counter = 0
    start = time()
    for i in range(0, len(cont_names)):
        counter += 1
        if counter % 100000 == 0:
            print >> Information, 'Time adding 100k keys', time() - start
            start = time()
        if cont_names[i] not in  C_dict:
            #errorhandle.unknown_contig(cont_names[i])
            continue

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
            S = Scaffold.scaffold(param.scaffold_indexer, [C], scaf_length)  # Create object scaffold
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
                S = Scaffold.scaffold(param.scaffold_indexer, [C], scaf_length)  # Create object scaffold
                small_scaffolds[S.name] = S
                C.scaffold = S.name
                param.scaffold_indexer += 1
                singeled_out += 1
    del C_dict

    return()

def CleanObjects(Contigs, Scaffolds, param, Information, small_contigs, small_scaffolds):
    singeled_out = 0
    scaf_lengths = [Scaffolds[scaffold_].s_length for scaffold_ in Scaffolds.keys()]
    sorted_lengths = sorted(scaf_lengths, reverse=True)
    scaf_lengths_small = [small_scaffolds[scaffold_].s_length for scaffold_ in small_scaffolds.keys()]
    sorted_lengths_small = sorted(scaf_lengths_small, reverse=True)
    N50, L50 = CalculateStats(sorted_lengths, sorted_lengths_small, param, Information)
    param.current_L50 = L50
    param.current_N50 = N50
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

def CreateEdge(cont_obj1, cont_obj2, scaf_obj1, scaf_obj2, G, param, alignedread, counter, contig1, contig2, save_obs=True):
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
    if param.orientation == 'fr':
        (obs1, obs2, scaf_side1, scaf_side2) = PosDirCalculatorPE(cont_dir1, read_dir, cont1_pos, readpos, s1len, cont1_len, cont_dir2, mate_dir, cont2_pos, matepos, s2len, cont2_len, param.read_len)
    else:
        (obs1, obs2, scaf_side1, scaf_side2) = PosDirCalculatorMP(cont_dir1, read_dir, cont1_pos, readpos, s1len, cont1_len, cont_dir2, mate_dir, cont2_pos, matepos, s2len, cont2_len, param.read_len)

    if obs1 == counter.prev_obs1 and obs2 == counter.prev_obs2:
        counter.nr_of_duplicates += 1
        if param.detect_duplicate:
            return(1)

    if obs1 + obs2 < param.ins_size_threshold and obs1 > 25 and obs2 > 25:
        if (scaf_obj2.name, scaf_side2) not in G[(scaf_obj1.name, scaf_side1)]:
            G.add_edge((scaf_obj2.name, scaf_side2), (scaf_obj1.name, scaf_side1), nr_links=1, obs=obs1 + obs2)
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['obs_sq'] = (obs1 + obs2) ** 2
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['observations'] = [obs1+obs2]

            if save_obs:
                G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)][scaf_obj1.name] = [obs1]
                G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)][scaf_obj2.name] = [obs2]
        elif save_obs:
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['nr_links'] += 1
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['obs'] += obs1 + obs2
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)][scaf_obj1.name].append(obs1)
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)][scaf_obj2.name].append(obs2)
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['obs_sq'] += (obs1 + obs2) ** 2
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['observations'].append(obs1+obs2)

        else:
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['nr_links'] += 1
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['obs'] += obs1 + obs2
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['obs_sq'] += (obs1 + obs2) ** 2
            G.edge[(scaf_obj1.name, scaf_side1)][(scaf_obj2.name, scaf_side2)]['observations'].append(obs1+obs2)


    else:
        counter.reads_with_too_long_insert += 1

        ## add to haplotype graph here!!

    counter.prev_obs1 = obs1
    counter.prev_obs2 = obs2
    return(0)


def CalculateStats(sorted_contig_lengths, sorted_contig_lengths_small, param, Information):
    cur_length = 0
    nr_conts = 0
    L50 = 0
    N50 = 0
    for contig_length in sorted_contig_lengths:
        cur_length += contig_length
        nr_conts += 1
        if cur_length >= param.tot_assembly_length / 2.0:
            N50 = contig_length
            L50 = nr_conts
            break
    if N50 == 0:
        for contig_length in sorted_contig_lengths_small:
            cur_length += contig_length
            nr_conts += 1
            if cur_length >= param.tot_assembly_length / 2.0:
                N50 = contig_length
                L50 = nr_conts
                break
    print >> Information, 'L50: ', L50, 'N50: ', N50, 'Initial contig assembly length: ', param.tot_assembly_length
    return(N50, L50)

def CalculateMeanCoverage(Contigs, Information, param):
    # tuples like (cont lenght, contig name)
    list_of_cont_tuples = [(Contigs[contig].length, contig) for contig in Contigs]
    #sorted as longest first
    list_of_cont_tuples = sorted(list_of_cont_tuples, key=lambda tuple: tuple[0], reverse=True)
    #coverages of longest contigs
    longest_contigs = list_of_cont_tuples[:50000]
    cov_of_longest_contigs = [Contigs[contig[1]].coverage for contig in longest_contigs]
    #Calculate mean coverage from the 1000 longest contigs
    n = float(len(cov_of_longest_contigs))
    mean_cov = sum(cov_of_longest_contigs) / n
    # If there is only one contig above the size threshold, n can be 1
    if n==1:
        n+=1

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


    if param.plots:
        plots.histogram(cov_of_longest_contigs, param, bins=50, x_label='coverage' , y_label='frequency' , title='BESST_cov_1000_longest_cont' + param.bamfile.split('/')[-1])

    return(mean_cov, std_dev)

def RemoveOutliers(mean_cov, std_dev, cov_list):
    k = normal.MaxObsDistr(len(cov_list), 0.95)
    filtered_list = list(filter((lambda x : (x < mean_cov + k * std_dev and x < 2 * mean_cov)), cov_list))
    if len(cov_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)

def RepeatDetector(Contigs, Scaffolds, G, param, G_prime, small_contigs, small_scaffolds, Information):
    mean_cov = param.mean_coverage
    std_dev = param.std_dev_coverage
    Repeats = []
    count_repeats = 0
    count_hapl = 0
    nr_of_contigs = len(Contigs)
    k = normal.MaxObsDistr(nr_of_contigs, 0.95)
    if param.cov_cutoff:
        repeat_thresh = param.cov_cutoff
    else:
        repeat_thresh = max(mean_cov + k * std_dev, 2 * mean_cov - 3*std_dev)
    print >> Information, 'Detecting repeats..'
    plot_cov_list = []
    for contig in Contigs:
        plot_cov_list.append(Contigs[contig].coverage)
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
        plot_cov_list.append(small_contigs[contig].coverage)
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


    if param.plots:
        plotting_ = filter(lambda x: x < 1000, plot_cov_list)
        plots.histogram(plotting_, param, bins=100, x_label='coverage' , y_label='frequency' , title='contigs_coverage' + param.bamfile.split('/')[-1])

    GO.repeat_contigs_logger(Repeats, Contigs, param.output_directory, small_contigs, param)
    GO.PrintOutRepeats(Repeats, Contigs, param.output_directory, small_contigs)
    print >> Information, 'Removed a total of: ', count_repeats, ' repeats. With coverage larger than ', repeat_thresh
    #sys.exit()
    if param.detect_haplotype:
        print >> Information, 'Marked a total of: ', count_hapl, ' potential haplotypes.'
    return(Contigs, Scaffolds, G)





def PosDirCalculatorPE(cont_dir1, read_dir, cont1pos, readpos, s1len, cont1_len, cont_dir2, mate_dir, cont2pos, matepos, s2len, cont2_len, read_len):
    if cont_dir1 and read_dir:
        obs1 = s1len - cont1pos - readpos
        scaf_side1 = 'R'
    if cont_dir2 and mate_dir:
        obs2 = s2len - cont2pos - matepos
        scaf_side2 = 'R'
    if (not cont_dir1) and read_dir:
        obs1 = cont1pos + (cont1_len - readpos)
        scaf_side1 = 'L'
    if (not cont_dir2) and mate_dir:
        obs2 = cont2pos + (cont2_len - matepos)
        scaf_side2 = 'L'
    if cont_dir1 and not read_dir:
        obs1 = cont1pos + readpos + read_len
        scaf_side1 = 'L'
    if cont_dir2 and not mate_dir:
        obs2 = cont2pos + matepos + read_len
        scaf_side2 = 'L'
    if not cont_dir1 and not read_dir:
        obs1 = s1len - cont1pos - (cont1_len - readpos - read_len)
        scaf_side1 = 'R'
    if not cont_dir2 and not mate_dir:
        obs2 = s2len - cont2pos - (cont2_len - matepos - read_len)
        scaf_side2 = 'R'
    return(int(obs1), int(obs2), scaf_side1, scaf_side2)

def PosDirCalculatorMP(cont_dir1, read_dir, cont1pos, readpos, s1len, cont1_len, cont_dir2, mate_dir, cont2pos, matepos, s2len, cont2_len, read_len):
    if cont_dir1 and not read_dir:
        obs1 = s1len - cont1pos - readpos
        scaf_side1 = 'R'
    if cont_dir2 and not mate_dir:
        obs2 = s2len - cont2pos - matepos
        scaf_side2 = 'R'
    if (not cont_dir1) and not read_dir:
        obs1 = cont1pos + (cont1_len - readpos)
        scaf_side1 = 'L'
    if (not cont_dir2) and not mate_dir:
        obs2 = cont2pos + (cont2_len - matepos)
        scaf_side2 = 'L'
    if cont_dir1 and read_dir:
        obs1 = cont1pos + readpos + read_len
        scaf_side1 = 'L'
    if cont_dir2 and mate_dir:
        obs2 = cont2pos + matepos + read_len
        scaf_side2 = 'L'
    if not cont_dir1 and read_dir:
        obs1 = s1len - cont1pos - (cont1_len - readpos - read_len)
        scaf_side1 = 'R'
    if not cont_dir2 and mate_dir:
        obs2 = s2len - cont2pos - (cont2_len - matepos - read_len)
        scaf_side2 = 'R'
    return(int(obs1), int(obs2), scaf_side1, scaf_side2)



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
import time
import copy

import networkx as nx
from networkx import algorithms
import multiprocessing

import Contig, Scaffold, Parameter
import GenerateOutput as GO
from mathstats.normaldist.truncatedskewed import param_est as GC
import mathstats.log_normal_param_est as lnpe
from mathstats.normaldist import normal
import ExtendLargeScaffolds as ELS
import haplotypes as HR
import plots
import order_contigs


# def constant_large():
#     return 2 ** 32

# def constant_small():
#     return -1


def Algorithm(G, G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, Information, param):
    #search for linear streches in the graph, remove all cliques >2 and all contigs having more than three neighbors
    nr_edges = 0
    for edge in G.edges_iter():
        if G[edge[0]][edge[1]]['nr_links']:
            nr_edges += 1

    print >> Information, str(nr_edges) + ' link edges created.'
    print >> Information, 'Perform inference on scaffold graph...'

    #VizualizeGraph(G,param,Information)

    if param.detect_haplotype:
        HR.HaplotypicRegions(G, G_prime, Contigs, Scaffolds, param, Information)
    #save graph in dot format to file here

    ##If sigma specified. Pre calculate a look up table for every possible gap estimate in the common case
    ##where we have two long contigs that are linked. We do this one time per library to save computation time.
    if param.std_dev_ins_size:
        dValuesTable = GC.PreCalcMLvaluesOfdLongContigs(param.mean_ins_size, param.std_dev_ins_size, param.read_len)
    else :
        dValuesTable = None
    already_visited = set()
    if param.extend_paths:
        for node in G:
            already_visited.add(node)


    ##### Here is the scaffolding algorithm #######

    G = RemoveIsolatedContigs(G, Information)     #step1
    plot = 'G'

    if not param.no_score:
        RemoveAmbiguousRegionsUsingScore(G, G_prime, Information, param, plot) #step2
        G = RemoveIsolatedContigs(G, Information) #there are probably new isolated nodes created from step 2
        G, Contigs, Scaffolds = RemoveLoops(G, G_prime, Scaffolds, Contigs, Information, param)    #step4    
        #The contigs that made it to proper scaffolds
        (Contigs, Scaffolds, param) = NewContigsScaffolds(G, G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, Information, dValuesTable, param, already_visited)   #step5
        ##Here PathExtension algorithm between created scaffolds is called if PRO is activated

    if param.extend_paths:
        print >> Information, '\n\n\n Searching for paths BETWEEN scaffolds\n\n\n'
        PROBetweenScaf(G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, param, dValuesTable, Information)
        print >> Information, 'Nr of contigs left: ', len(G_prime.nodes()) / 2.0, 'Nr of linking edges left:', len(G_prime.edges()) - len(G_prime.nodes()) / 2.0
        print >> Information, 'Number of gaps estimated by GapEst-LP module order_contigs in this step is: {0}'.format(param.path_gaps_estimated)

    if param.plots:
        if len(param.gap_estimations) > 0:
            plots.histogram(param.gap_estimations, param, bins=50, x_label='gap size', y_label='Frequency', title='GapPredictions', nr_obs = len(param.gap_estimations))

 ## TODO: Look if we can use the new pathgaps theory to include these isolated 
 ## regions of small contigs in the scaffolding

#
#        for node in G_prime.nodes():
#            nbrs = G_prime.neighbors(node)
#            for nbr in nbrs:
#                if G_prime[node][nbr]['nr_links'] and 'score' not in G_prime[node][nbr]:
#                    G_prime.remove_edge(node, nbr)
#
#        plot = 'G_prime'
#        #print 'GOING IN:', len(G_prime.edges()), len(G_prime.nodes())
#        RemoveAmbiguousRegionsUsingScore(G_prime, G_prime, Information, param, plot)
#        G_prime, Contigs, Scaffolds = RemoveLoops(G_prime, G_prime, Scaffolds, Contigs, Information, param)
#        for node in G_prime:
#            if node[0] not in Scaffolds:
#                scaf_obj = small_scaffolds[node[0]]
#                Scaffolds[node[0]] = scaf_obj
#                cont_objects = scaf_obj.contigs
#                for obj_ in cont_objects:
#                    ctg_name = obj_.name
#                    Contigs[ctg_name] = obj_
#                    del small_contigs[ctg_name]
#                del small_scaffolds[node[0]]
#        (Contigs, Scaffolds, param) = NewContigsScaffolds(G_prime, G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, Information, dValuesTable, param, already_visited)


    ####### End of algorithm #####################


    return()



def RemoveIsolatedContigs(G, Information):
    print >> Information, 'Remove isolated nodes.'
    counter = 0
    for node in G.nodes():
        if node in G:
            nbr = G.neighbors(node)[0]
            if len(G.neighbors(node)) == 1 and len(G.neighbors(nbr)) == 1:
                counter += 1
                G.remove_nodes_from([node, nbr])
    print >> Information, str(counter) + ' isolated contigs removed from graph.'
    return(G)


def partition(pred, iterable):
    trues = []
    falses = []
    for item in iterable:
        if pred(item):
            trues.append(item)
        else:
            falses.append(item)
    return trues, falses

def remove_edges(G, G_prime, Information, param, node, score_chosen_obs, non_zero_score_removed_obs, edge_score_to_zero, record_decision):
    nbrs = G.neighbors(node)


    #Remove ambiguous edges
    filter_nbrs = filter(lambda nbr: 0 < G[node][nbr]['nr_links'] , nbrs) # remove other side of contig (has no "links")
    score_list_temp = sorted(map(lambda nbr: (G[node][nbr]['score'], nbr), filter_nbrs)) # sort list of scores for neighbors
    non_zero_score_edges, zero_score_edges = partition(lambda nbr: 0 < G[node][nbr[1]]['score'], score_list_temp) # Split list into 0-scoring and non-0-scoring edges


    #remove all zero score edges
    remove_zero_score_edges = map(lambda item: (node, item[1]) , zero_score_edges)
    G.remove_edges_from(remove_zero_score_edges)
    if param.plots:
        for item in zero_score_edges:
            edge_score_to_zero.append(item[0])

    if param.extend_paths:
        try: #we might have been removed this edge from G_prime when we did individual filtering of G_prime in CreateGraph module
            G_prime.remove_edges_from(remove_zero_score_edges)
        except nx.exception.NetworkXError:
            pass

    # Remove lower non zero score edges
    if len(non_zero_score_edges) > 1:
        if non_zero_score_edges[-2][0] / non_zero_score_edges[-1][0] > 0.8: # lowered threshold due to the unceartinaty in scores
            print >> Information, 'SCORES AMBVIVALENT', non_zero_score_edges[-1][0], non_zero_score_edges[-2][0]
            remove_non_zero_edges = map(lambda item: (node, item[1]) , non_zero_score_edges) # Edges that does not have a score of 0 but are not the highest scoring one in an ambigous region
        else:
            remove_non_zero_edges = map(lambda item: (node, item[1]) , non_zero_score_edges[:-1]) # Edges that does not have a score of 0 but are not the highest scoring one in an ambigous region

        G.remove_edges_from(remove_non_zero_edges) # remove the lower scoring edges
        if param.extend_paths:
            try: #we might have been removed this edge from G_prime when we did individual filtering of G_prime in CreateGraph module
                G_prime.remove_edges_from(remove_non_zero_edges)
            except nx.exception.NetworkXError:
                pass

    if param.plots:
        if len(non_zero_score_edges) > 0:
            if record_decision:
                score_chosen_obs.append(non_zero_score_edges[-1][0])
            for item in non_zero_score_edges[:-1]:
                non_zero_score_removed_obs.append(item[0])


    return()

def RemoveAmbiguousRegionsUsingScore(G, G_prime, Information, param, plot):
    score_chosen_obs = []
    non_zero_score_removed_obs = []
    edge_score_to_zero = []

    nr_edges_before = len(G.edges())
    print >> Information, 'Remove edges from node if more than two edges'
    #counter1 = 0

    link_edges = filter(lambda x: 'score' in x[2], G.edges(data=True))
    edge_scores_sorted = sorted(link_edges, key=lambda x: x[2]['score'], reverse=True)

    for edge in edge_scores_sorted:
        remove_edges(G, G_prime, Information, param, edge[0], score_chosen_obs, non_zero_score_removed_obs, edge_score_to_zero, True)
        remove_edges(G, G_prime, Information, param, edge[1], score_chosen_obs, non_zero_score_removed_obs, edge_score_to_zero, False)

    nr_edges_after = len(G.edges())
    #print >> Information, str(counter1) + ' ambiguous regions in graph ( a contig with more than 2 neighbors).'
    print >> Information, ' Number of edges in G before:', nr_edges_before
    print >> Information, ' Number of edges in G after:', nr_edges_after
    try:
        print >> Information, ' %-age removed edges:', 100 * (1 - (nr_edges_after / float(nr_edges_before)))
    except ZeroDivisionError:
        pass

    if param.plots:
        print len(score_chosen_obs)
        list_of_datasets = [edge_score_to_zero, score_chosen_obs , non_zero_score_removed_obs]
        plots.multiple_histogram(list_of_datasets, param, 'score', 'frequency', title='Besst_decision_scores' + plot + '.' + param.bamfile.split('/')[-1])

    return()






def RemoveLoops(G, G_prime, Scaffolds, Contigs, Information, param):
#### After the proceure above, we hope that the graph is almost perfectly linear but we can still be encountering cycles (because of repeats or haplotypic contigs that has slipped through our conditions). Thus we finally search for loops
    print >> Information, 'Contigs/scaffolds left:', len(G.nodes()) / 2
    print >> Information, 'Remove remaining cycles...'
    #graphs = nx.connected_component_subgraphs(G)
    #print 'Nr connected components',len(graphs)
    counter = 0
    #for graph in graphs:
    list_of_cycles = algorithms.cycles.cycle_basis(G)
    for cycle in list_of_cycles:
        print >> Information, 'A cycle in the scaffold graph: ' + str(cycle) + '\n'
        print >> Information, 'A cycle in the scaffold graph: ' + str(cycle)
        counter += 1
        for node in cycle:
            if node in G:
                #we split up the whole cycle into separate contigs and send them to F
                scaffold_ = node[0]
                G.remove_nodes_from([(scaffold_, 'L'), (scaffold_, 'R')])
                if param.extend_paths:
                    G_prime.remove_nodes_from([(scaffold_, 'L'), (scaffold_, 'R')])
#                    S_obj=Scaffolds[scaffold_]
#                    list_of_contigs=S_obj.contigs   #list of contig objects contained in scaffold object
#                    Contigs, F = GO.WriteToF(F,Contigs,list_of_contigs)
#                    del Scaffolds[scaffold_]
    print >> Information, str(counter) + ' cycles removed from graph.'
    return(G, Contigs, Scaffolds)


def NewContigsScaffolds(G, G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, Information, dValuesTable, param, already_visited):
### Remaining scaffolds are true sensible scaffolds, we must now update both the library of scaffold objects and the library of contig objects
    new_scaffolds_ = [ G.subgraph(c) for c in nx.connected_components(G)]
    #print 'INITIAL SCAFFOLD INDEX:',param.scaffold_indexer
    print >> Information, 'Nr of new scaffolds created in this step: ' + str(len(new_scaffolds_))
    for new_scaffold_ in new_scaffolds_:
        param.scaffold_indexer += 1
        #scaf_size=len(new_scaffold_)
        scaffold_length = 0
        contig_list = []

        ##### Here PathExtension algorithm is called if PRO is activated #####
        if param.extend_paths:
            #print 'SCAFFOLD INDEX:',param.scaffold_indexer
            PROWithinScaf(G, G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, param, new_scaffold_, dValuesTable, already_visited)

        for node in new_scaffold_:
            if len(G.neighbors(node)) == 1:
                start = node
                break
        for node in new_scaffold_:
            if len(G.neighbors(node)) == 1 and node != start:
                end = node
        #Create info to new scaffold object such as total length and the contig objects included

        prev_node = ('', '')
        pos = 0
        values = UpdateInfo(G, Contigs, small_contigs, Scaffolds, small_scaffolds, start, prev_node, pos, contig_list, scaffold_length, dValuesTable, param)
        #print values
        while len(values) !=2:
            values = UpdateInfo(G,*values)
            #print len(values)
        #(G, contig_list, scaffold_length) = UpdateInfo(G, Contigs, small_contigs, Scaffolds, small_scaffolds, start, prev_node, pos, contig_list, scaffold_length, dValuesTable, param)
        contig_list, scaffold_length = values[0],values[1]
        S = Scaffold.scaffold(param.scaffold_indexer, contig_list, scaffold_length)  #Create the new scaffold object 

        Scaffolds[S.name] = S        #include in scaffold library

        if param.extend_paths:
            # Find the ends of the old subgraph new_scaffold_. We want them to be able to relabel these end nodes as the new sides on the new scaffold object created
            #only these ends are allowed to have links because they are of size mean+ 4*sigma so nothing is supposed to span over.

            #add the new scaffold object to G_prime

            G_prime.add_node((S.name, 'L'))  #start node
            G_prime.add_node((S.name, 'R'))  # end node
            G_prime.add_edge((S.name, 'L'), (S.name, 'R'), nr_links=None)
            try:
                for nbr in G_prime.neighbors(start):
                    nr_links_ = G_prime[start][nbr]['nr_links']
                    if nr_links_:
                        obs_ = G_prime[start][nbr]['obs']
                        observations_ = G_prime[start][nbr]['observations']
                        # G_prime.add_edge((S.name, 'L'), nbr, nr_links=nr_links_, obs=obs_)
                        obs_sq_ = G_prime[start][nbr]['obs_sq']
                        G_prime.add_edge((S.name, 'L'), nbr, nr_links=nr_links_, obs=obs_, obs_sq=obs_sq_, observations=observations_)

                for nbr in G_prime.neighbors(end):
                    nr_links_ = G_prime[end][nbr]['nr_links']
                    if nr_links_:
                        obs_ = G_prime[end][nbr]['obs']
                        observations_ = G_prime[end][nbr]['observations']
                        # G_prime.add_edge((S.name, 'R'), nbr, nr_links=nr_links_, obs=obs_)
                        obs_sq_ = G_prime[end][nbr]['obs_sq']
                        G_prime.add_edge((S.name, 'R'), nbr, nr_links=nr_links_, obs=obs_, obs_sq=obs_sq_, observations=observations_)

                #remove the old scaffold objects from G_prime
                G_prime.remove_nodes_from(new_scaffold_)
            except nx.exception.NetworkXError:
                pass

    return(Contigs, Scaffolds, param)


def UpdateInfo(G, Contigs, small_contigs, Scaffolds, small_scaffolds, node, prev_node, pos, contig_list, scaffold_length, dValuesTable, param):
    scaf = node[0]
    side = node[1]
    prev_scaf = prev_node[0]
    if len(G.neighbors((scaf, side))) == 0:  #reached end of scaffol
        #find the contig with the largest position
        object_with_largest_pos_in_scaffold = max(contig_list, key=lambda object: object.position + object.length)
        scaffold_length = object_with_largest_pos_in_scaffold.position + object_with_largest_pos_in_scaffold.length
#        try:
        del Scaffolds[scaf] #finally, delete the old scaffold object
        G.remove_node((scaf, side))
#        except KeyError:
#            del small_scaffolds[scaf] #finally, delete the old scaffold object
        return(contig_list, scaffold_length)
    else:
        nbr_node = G.neighbors((scaf, side))
        nbr_scaf = nbr_node[0][0]
        nbr_side = nbr_node[0][1]
        if scaf != prev_scaf:
            if side == 'L':    #Contig/scaffold still has same orientation as in previous iteration, just update position in scaffold                                           
                #want to assign nr of links to contig object, note that in case of a "multiple contigs"-scaffold object, only the outermost contig holds the information of the total nr of links between the two scaffold objects
                #try:
                contig_objects = Scaffolds[scaf].contigs #list of contig objects
                #except KeyError:
                #    contig_objects=small_scaffolds[scaf].contigs #list of contig objects
                #Update just update position in scaffold 
                for contig in contig_objects:
                    contig.scaffold = param.scaffold_indexer
                    contig.position += pos
                    #direction unchanged
                    contig_list.append(contig)
                G.remove_node((scaf, side))
                prev_node = node
                node = (nbr_scaf, nbr_side)
                #try:                 
                pos += Scaffolds[scaf].s_length  #update position before sending it to next scaffold
                #except KeyError:
                #    pos+=small_scaffolds[scaf].s_length  #update position before sending it to next scaffold

                return Contigs, small_contigs, Scaffolds, small_scaffolds, node, prev_node, pos, contig_list, scaffold_length, dValuesTable, param

                #G, contig_list, scaffold_length = UpdateInfo(G, Contigs, small_contigs, Scaffolds, small_scaffolds, node, prev_node, pos, contig_list, scaffold_length, dValuesTable, param)

            else:  #Contig/scaffold need to change orientation as well as modify orientation in this case
                #try:
                contig_objects = Scaffolds[scaf].contigs #list of contig objects
                #except KeyError:
                #    contig_objects=small_scaffolds[scaf].contigs #list of contig objects

                for contig in contig_objects:
                    contig.scaffold = param.scaffold_indexer
                    #try:
                    curr_scaf_length = Scaffolds[scaf].s_length
                    #except KeyError:
                    #    curr_scaf_length=small_scaffolds[scaf].s_length

                    curr_pos_within_scaf = contig.position
                    contig_length = contig.length
                    contig.position = pos + (curr_scaf_length - curr_pos_within_scaf) - contig_length #updates the position within scaf
                    contig.direction = bool(True -contig.direction) #changes the direction
                    contig_list.append(contig)

                G.remove_node((scaf, side))
                prev_node = node
                node = (nbr_scaf, nbr_side)
                #try:                    
                pos += Scaffolds[scaf].s_length  #update position before sending it to next scaffold
                #except KeyError:
                #pos+=small_scaffolds[scaf].s_length  #update position before sending it to next scaffold

                return Contigs, small_contigs, Scaffolds, small_scaffolds, node, prev_node, pos, contig_list, scaffold_length, dValuesTable, param

#                G, contig_list, scaffold_length = UpdateInfo(G, Contigs, small_contigs, Scaffolds, small_scaffolds, node, prev_node, pos, contig_list, scaffold_length, dValuesTable, param)
        else:
            if  'avg_gap' not in  G[(scaf, side)][(nbr_scaf, nbr_side)]:
                #calculate gap to next scaffold
                c1_len = Scaffolds[scaf].s_length
                c2_len = Scaffolds[nbr_scaf].s_length
                if param.lognormal:
                    # print G[(scaf, side)][(nbr_scaf, nbr_side)]
                    # print G[(scaf, side)]
                    samples = G[(scaf, side)][(nbr_scaf, nbr_side)]['observations']
                    avg_gap = lnpe.GapEstimator(param.lognormal_mean, param.lognormal_sigma, param.read_len, samples, c1_len, c2_len=c2_len)
                else:

                    sum_obs = G[(scaf, side)][(nbr_scaf, nbr_side)]['obs']
                    nr_links = G[(scaf, side)][(nbr_scaf, nbr_side)]['nr_links']
                    data_observation = (nr_links * param.mean_ins_size - sum_obs) / float(nr_links)
                    mean_obs = sum_obs / float(nr_links)
                    #try:
                    c1_len = Scaffolds[scaf].s_length
                    #except KeyError:
                    #    c1_len=small_scaffolds[scaf].s_length
                    #try:
                    c2_len = Scaffolds[nbr_scaf].s_length
                    #except KeyError:
                    #   c2_len=small_scaffolds[nbr_scaf].s_length
                    #do fancy gap estimation by the bias estimator formula
                    if param.std_dev_ins_size and nr_links >= 5:
                        #pre calculated value in lookup table 
                        if c1_len > param.mean_ins_size + 4 * param.std_dev_ins_size and c2_len > param.mean_ins_size + 4 * param.std_dev_ins_size:
                            #(heuristic scale down of table to gaps of at most 2 stddevs away from mean)
                            try:
                                avg_gap = dValuesTable[int(round(data_observation, 0))]
                            except KeyError:
                                avg_gap = GC.GapEstimator(param.mean_ins_size, param.std_dev_ins_size, param.read_len, mean_obs, c1_len, c2_len)
                                #print 'Gap estimate was outside the boundary of the precalculated table, obs were: ', data_observation, 'binary search gave: ', avg_gap
                        #Do binary search for ML estimate of gap if contigs is larger than 3 times the std_dev
                        elif c1_len > param.std_dev_ins_size + param.read_len and c2_len > param.std_dev_ins_size + param.read_len:
                            avg_gap = GC.GapEstimator(param.mean_ins_size, param.std_dev_ins_size, param.read_len, mean_obs, c1_len, c2_len)
                        else:
                            #print 'now', 2 * param.std_dev_ins_size + param.read_len
                            avg_gap = int(data_observation)
                            param.gap_estimations.append( avg_gap )
                        #print 'Gapest if used:' + str(int(avg_gap)), 'Naive: ' + str(int(data_observation)), c1_len, c2_len, Scaffolds[scaf].contigs[0].name, Scaffolds[nbr_scaf].contigs[0].name
                        #See if the two contigs are in fact negatively overlapped in the delta file, , then abyss produses
                        #contigs contained in other contigs
                    #do naive gap estimation
                    else:
                        avg_gap = int(data_observation)
                        param.gap_estimations.append( avg_gap )
            else:
                avg_gap = G[(scaf, side)][(nbr_scaf, nbr_side)]['avg_gap']
                param.gap_estimations.append( avg_gap )

            if avg_gap <= 1:
                avg_gap = 1

            pos += int(avg_gap)
            G.remove_node((scaf, side))
            prev_node = node
            node = (nbr_scaf, nbr_side)
            #try:
            del Scaffolds[scaf] #finally, delete the old scaffold object
            #except KeyError:
            #    del small_scaffolds[scaf] #finally, delete the old scaffold object
            return Contigs, small_contigs, Scaffolds, small_scaffolds, node, prev_node, pos, contig_list, scaffold_length, dValuesTable, param
            #G, contig_list, scaffold_length = UpdateInfo(G, Contigs, small_contigs, Scaffolds, small_scaffolds, node, prev_node, pos, contig_list, scaffold_length, dValuesTable, param)
    return( contig_list, scaffold_length)

def get_total_length(Scaffolds, small_scaffolds, path):
    total_within_path_length = 0
    for i,node in enumerate(path[1:-1]):
        if i % 2 ==0:
            try:
                total_within_path_length += small_scaffolds[node[0]].s_length
            except KeyError:
                print 'Contig/Scaffold {0} is not in small_scaffolds'.format(node[0])
                return -1
    return total_within_path_length





def PROWithinScaf(G, G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, param, new_scaffold_, dValuesTable, already_visited):
    #loc_count = 0
    for edge in new_scaffold_.edges_iter():
        nr_links_ = G[edge[0]][edge[1]]['nr_links']
        if nr_links_:
            start = edge[0]
            end = edge[1]
            c1_len = Scaffolds[start[0]].s_length
            c2_len = Scaffolds[end[0]].s_length
            sum_obs = G[start][end]['obs']
            data_observation = (nr_links_ * param.mean_ins_size - sum_obs) / float(nr_links_)
            #avg_gap = GC.GapEstimator(param.mean_ins_size, param.std_dev_ins_size, param.read_len, data_observation, c1_len, c2_len)

            #high_score_path, bad_links, score, path_len = ELS.WithinScaffolds(G, G_prime, start, end, already_visited, param.ins_size_threshold, param)
            all_paths_sorted_wrt_score = ELS.WithinScaffolds(G, G_prime, start, end, already_visited, param.ins_size_threshold, param)

            high_score_path = False
            for path_info in reversed(all_paths_sorted_wrt_score):
                path = path_info[2]
                length_of_path = get_total_length(Scaffolds, small_scaffolds, path)
                if 0 < length_of_path < param.ins_size_threshold:
                    high_score_path, score = path, path_info[0]
                    break
            if not all_paths_sorted_wrt_score:
                continue

            # if all_paths_sorted_wrt_score[-1][1] != 0:
            #     sublist = all_paths_sorted_wrt_score[-1]
            #     bad_links = sublist[1]
            #     score = sublist[0]
            #     path_len = sublist[3]
            #     print 'Path: path length: {0}, nr bad links: {1}, score: {2} '.format((path_len - 2) / 2.0, bad_links, score)

            # # for sublist in reversed(all_paths_sorted_wrt_score):
            # #     path = sublist[2]
            # #     bad_links = sublist[1]
            # #     score = sublist[0]
            # #     path_len = sublist[3]
            # #     print 'Path: path length: {0}, nr bad links: {1}, score: {2} '.format((path_len - 2) / 2.0, bad_links, score)

            #high_score_path, score = all_paths_sorted_wrt_score[-1][2], all_paths_sorted_wrt_score[-1][0]
            #print 'Highest scoring path:{0}'.format(high_score_path)
            #print 'Score: {0}'.format(score)

            if high_score_path and score >= param.score_cutoff:

                ##################### v1.0.4.5 
                ## modified improved path gap estimation here!!

                high_score_path_copy = copy.deepcopy(high_score_path)
                G_, path = estimate_path_gaps(Contigs, small_contigs, high_score_path_copy,Scaffolds,small_scaffolds, G_prime,param)
                del path[0]
                del path[-1]
                G.remove_edge(start, end)
                G.add_edges_from(G_.edges(data=True))
                G_prime.remove_nodes_from(path[1:-1])

                # also remove the edges from inner ends of the large scaffolds so they can't be
                # involved in creating paths in "between scaffolds"

                G_prime.remove_node(path[0])
                G_prime.remove_node(path[-1])
                if path[0][1] == 'L':
                    G_prime.add_edge(path[0],(path[0][0],'R'),nr_links=None)
                else:
                    G_prime.add_edge(path[0],(path[0][0],'L'),nr_links=None)
                if path[-1][1] == 'L':
                    G_prime.add_edge(path[-1],(path[-1][0],'R'),nr_links=None)
                else:
                    G_prime.add_edge(path[-1],(path[-1][0],'L'),nr_links=None)

                # move all contig and scaffold objects from "small" structure to large structure to fit with UpdateInfo structure
                small_scafs = map(lambda i: path[i], filter(lambda i: i % 2 == 1, range(len(path) - 1)))
                for item in small_scafs:
                    scaf_obj = small_scaffolds[item[0]]
                    Scaffolds[item[0]] = scaf_obj
                    cont_objects = scaf_obj.contigs
                    for obj_ in cont_objects:
                        ctg_name = obj_.name
                        Contigs[ctg_name] = obj_
                        del small_contigs[ctg_name]
                    del small_scaffolds[item[0]]
                #sys.exit()
                ######################################



                # #loc_count += 1
                # #remove edge in G to fill in the small scaffolds
                # G.remove_edge(start, end)
                # #add small scaffolds to G
                # for i in range(0, len(high_score_path) - 1):
                #     nr_lin = G_prime[high_score_path[i]][high_score_path[i + 1]]['nr_links']
                #     try:
                #         total_dist = G_prime[high_score_path[i]][high_score_path[i + 1]]['obs']
                #         G.add_edge(high_score_path[i], high_score_path[i + 1], nr_links=nr_lin, obs=total_dist)
                #     except KeyError:
                #         G.add_edge(high_score_path[i], high_score_path[i + 1], nr_links=nr_lin)
                # #remove the small contigs from G_prime
                # G_prime.remove_nodes_from(high_score_path[1:-1])
                # # move all contig and scaffold objects from "small" structure to large structure to fit with UpdateInfo structure
                # small_scafs = map(lambda i: high_score_path[i], filter(lambda i: i % 2 == 1, range(len(high_score_path) - 1)))
                # for item in small_scafs:
                #     scaf_obj = small_scaffolds[item[0]]
                #     Scaffolds[item[0]] = scaf_obj
                #     cont_objects = scaf_obj.contigs
                #     for obj_ in cont_objects:
                #         ctg_name = obj_.name
                #         Contigs[ctg_name] = obj_
                #         del small_contigs[ctg_name]
                #     del small_scaffolds[item[0]]

    ####################################################################
    return()


def other_end(node):
    if node[1] == 'L':
        return (node[0],'R')
    elif node[1] == 'R':
        return (node[0],'L')
    else:
        print 'Node is not properly declared (has not right or left end)'
        sys.exit()



def permute_path(path, ctg_to_move, contig_after):
    pos_ctg_to_move_left = path.index((ctg_to_move,'L'))
    pos_ctg_to_move_right = path.index((ctg_to_move,'R'))

    path = filter(lambda x: x[0] != ctg_to_move, path)
    pos_of_ctg_after = map(lambda x: x[0],path).index(contig_after)
    if pos_ctg_to_move_right < pos_ctg_to_move_left:
        path.insert(pos_of_ctg_after, (ctg_to_move,'L') )
        path.insert(pos_of_ctg_after, (ctg_to_move,'R') )
    else:
        path.insert(pos_of_ctg_after, (ctg_to_move,'R') )
        path.insert(pos_of_ctg_after, (ctg_to_move,'L') )

    return path

# def permutation_ok(path,):

#         self.mean = param.mean_ins_size
#         self.stddev = param.std_dev_ins_size
#         self.read_len = param.read_len
#         self.contamination_ratio = param.contamination_ratio
#         self.contamination_mean = param.contamination_mean
#         self.contamination_stddev = param.contamination_stddev
#         self.ctgs = []
#         for i,length in enumerate(ctg_lengths):
#             self.ctgs.append(Contig(i, length))
#         self.ctgs = tuple(self.ctgs)
#         self.gaps = [0]*(len(ctg_lengths)-1) # n contigs has n-1 gaps between them, start with gap size 0  
        
#         # get positions for when all gaps are 0
#         self.update_positions()

#         # let us cheat here! Instead of calculating likeliooods of thousands of
#         # onservations we calculate the ikelihood for them average (mean) of the
#         # observations and weight it with the number of observations
#         self.mp_links = 0.0
#         self.pe_links = 0.0
#         obs_dict = {}

#         for c1,c2,is_PE_link in observations:
#             #nr_obs = len(observations[(c1,c2,is_PE_link)])
#             mean_obs, nr_obs, stddev_obs = observations[(c1,c2,is_PE_link)]
#             if is_PE_link:
#                 mean_PE_obs = self.ctgs[c1].length + self.ctgs[c2].length - observations[(c1,c2,is_PE_link)][0] + 2*param.read_len 
#                 #PE_obs = map(lambda x: self.ctgs[c1].length + self.ctgs[c2].length - x + 2*param.read_len ,observations[(c1,c2,is_PE_link)])
#                 #mean_obs = sum( PE_obs)/nr_obs
#                 obs_dict[(c1, c2, is_PE_link)] = (mean_PE_obs, nr_obs, stddev_obs)
#                 self.pe_links += nr_obs
#                 if mean_PE_obs > self.contamination_mean + 6 * self.contamination_stddev and not initial_path:
#                     self.observations = None
#                     return None
#             else:
#                 #mean_obs = sum(observations[(c1,c2,is_PE_link)])/nr_obs
#                 obs_dict[(c1, c2, is_PE_link)] = (mean_obs, nr_obs, stddev_obs)
#                 self.mp_links += nr_obs
            



def calculate_path_LP(current_path, Scaffolds, small_scaffolds, observations, param, moved_forward=None, shifted=None):
    contigs_to_indexes = {}
    indexes_to_contigs = {}
    index = 0
    #print 'CURNNNT PATH_',current_path
    for ctg in current_path:
        if ctg[0] in contigs_to_indexes:
            continue
        else:
            contigs_to_indexes[ctg[0]] = index
            indexes_to_contigs[index] =  ctg[0]

            index += 1

    contig_lengths = []
    for ctg in contigs_to_indexes:
        try:
            contig_lengths.append((contigs_to_indexes[ctg], Scaffolds[ctg].s_length))
        except KeyError:
            contig_lengths.append((contigs_to_indexes[ctg], small_scaffolds[ctg].s_length))

    #contig_lengths2 = map(lambda x: (contigs_to_indexes[x], Scaffolds[x].s_length), contigs_to_indexes) 

    #print contig_lengths
    ctg_lengths_index_ordered = map(lambda x: x[1], sorted(contig_lengths, key=lambda x: x[0]))
    #print ctg_lengths_index_ordered
    tot_links = 0
    mp_links = 0
    index_observations = {}  
    permuted_contigs = set([moved_forward, shifted])
   
    for c1,c2 in observations:  
        tot_links += observations[(c1,c2)][1]
        if current_path.index(c1) < current_path.index(c2) and current_path.index(c1) % 2 == 1 and current_path.index(c2) % 2 == 0 and param.contamination_ratio:
            PE = 1
            #print 'PE link!!',c1,c2
        elif current_path.index(c1) > current_path.index(c2) and current_path.index(c1) % 2 == 0 and current_path.index(c2) % 2 == 1 and param.contamination_ratio:
            PE = 1
            #print 'PE link!!',c1,c2
        elif current_path.index(c1) < current_path.index(c2) and current_path.index(c1) % 2 == 0 and current_path.index(c2) % 2 == 1:
            PE = 0
            mp_links += observations[(c1,c2)][1]
            #print 'MP link!!',c1,c2
        elif current_path.index(c1) > current_path.index(c2) and current_path.index(c1) % 2 == 1 and current_path.index(c2) % 2 == 0:
            mp_links += observations[(c1,c2)][1]
            PE = 0
            #print 'MP link!!',c1,c2           
        else:
            PE = -1
            #print 'Spurious link!!',c1,c2

        if PE == -1:
            # we dont want to have a spurious link to contiribute to the LP objective value
            continue
        else:
            i1, i2 = min(contigs_to_indexes[c1[0]], contigs_to_indexes[c2[0]]), max(contigs_to_indexes[c1[0]], contigs_to_indexes[c2[0]])
            index_observations[(i1,i2,PE)] = observations[(c1,c2)]

        if c1[0] in permuted_contigs and c2[0] in permuted_contigs:
            index_moved_forward = i1
            index_shifted = i2

    if moved_forward and shifted:
        is_PE_link = 1
        mean_obs, nr_obs, stddev_obs, list_of_obs = index_observations[(index_moved_forward, index_shifted, is_PE_link)]
        mean_PE_obs = ctg_lengths_index_ordered[index_moved_forward] + ctg_lengths_index_ordered[index_shifted] - mean_obs + 2*param.read_len
        if mean_PE_obs > param.contamination_mean + 6 * param.contamination_stddev:
            #print 'permuted link invalid. Continuing..'
            if all(length in ctg_lengths_index_ordered for length in [670, 2093]) or all(length in ctg_lengths_index_ordered for length in [900, 3810]) or all(length in ctg_lengths_index_ordered for length in [2528, 591]) or all(length in ctg_lengths_index_ordered for length in [734, 257, 1548]):
                print >> param.information_file, 'permuted link invalid. Continuing..', 'moved forward:', index_moved_forward, 'shifted:', index_shifted
            return None, None, None, None


        #print i1,i2 #'OBSLIST_', observations[(c1,c2)]
    #observations =  map(lambda x: (contigs_to_indexes[x[0][0]], contigs_to_indexes[x[1][0]]) = observations[x], observations)   

    #print 'MP LINK RATIO:{0}, tot_links:{1}'.format(mp_links/float(tot_links), tot_links)
    #print index_observations
    #print ctg_lengths_index_ordered
    #print index_observations

    ## 2 Get optimal LP solution for given path order
    result_path = order_contigs.main(ctg_lengths_index_ordered, index_observations, param)
    if not result_path:
        print 'We should never end up here!'
        return None,None,None,None

    return result_path, contigs_to_indexes, indexes_to_contigs, index_observations


def path_permutations_iterator(path, i, original_path):
    """
        Iterator function to return all "valid" permutations of contigs given our
        assumptions of links. This is a very limited set of permutations compared to the m!
        if m contigs
        blowup of a general permutation set. For small instances (m < 5), this will return 2^m
        permutations. If m is larger, we return 2^5 + 2^5 + ... +2^(rest < 5) permutations where the sum
        contians m/5 +1 terms. This is a huristic but well motivated by the fact that far away contigs in initial 
        solution is unlikely to be permuted so that they switch order with each other. (This depends on the isize of course)

        Another heuristic could be to instead of chosing 5, we scale down 2^m for large m by
        giving back all permutations of sets of contigs that have sum of lenghts larger than, say, mu+3sigma. This
        is another way of separating and scaling down the solution space by using the fact that far away contigs 
        in initial solution is unlikely to be permuted so that they switch order with each other. 

        In original_path : switch is always performed on contig i with contig i-1. Since each
        contig except start and en here has 2 nodes in the path object, we will always permute
        items [i,i+1] with items [i-2,i-1] in original path, where i =3, 5, 7,..., len(path) - 3.
        Object on index len(path)-1 in the list (last object) is the right border contig.  
    """
        #Path needs to be a deepcopy of original_path !!!!
    if 3<= i <= len(original_path)-3:
        # no permutation
        #yield path
        #print '1lol',i,path
        for p in path_permutations_iterator(copy.deepcopy(path), i+2, original_path):
            yield p #path_permutations_iterator(path, i+2, original_path)

        ### permuting ###
        #contigs to permute
        c1_x, c1_y, c2_x, c2_y = original_path[i-2],original_path[i-1],original_path[i],original_path[i+1]
        # current indexes of them in path
        i_c1_x, i_c1_y, i_c2_x, i_c2_y = path.index(c1_x), path.index(c1_y), path.index(c2_x), path.index(c2_y)

        # permutation is that i_c2_x --> i_c1_x , i_c2_y --> i_c1_y and all contigs after that
        # is shited 2 indicies up  
        path.insert(i_c1_x, path.pop(i_c2_x))
        path.insert(i_c1_x+1, path.pop(i_c2_y))

        # path.insert(pos_of_ctg_after, (ctg_to_move,'L') )
        # path.insert(pos_of_ctg_after, (ctg_to_move,'R') )
        # path.insert(pos_of_ctg_after, (ctg_to_move,'R') )
        # path.insert(pos_of_ctg_after, (ctg_to_move,'L') )


        # path[i_c1_x],path[i_c1_y],path[i_c2_x],path[i_c2_y] = path[i_c2_x],path[i_c2_y],path[i_c1_x],path[i_c1_y]
        #yield path
        #print '2lol',i,path
        for p in path_permutations_iterator(copy.deepcopy(path), i+2, original_path):
            yield p #path_permutations_iterator(path, i+2, original_path)
    else:
        yield path



def path_permutations_with_overlap_constraints(original_path, constraint_dict):
    """
        Iterator function to return all "valid" permutations of contigs given overlap constraints.
        This is a very limited set of permutations compared to the m!
        if m contigs.
    """
    # # Adjacancies of valid permutations according to our model
    # path = copy.deepcopy(original_path)
    # valid_permutations = set()
    # for i in range(3, len(path) - 1, 2):
    #     # switch positions of two contigs
    #     ctg_end = original_path[i]
    #     contig_end = original_path[i-2]    
    # # 

    # order original path with the constraint pairs
    path = copy.deepcopy(original_path)
    print "before constraints", path

    long_end = path[0]
    long_beginning = path[-1]

    for beginning in constraint_dict:
        c2_x = beginning
        c1_y = constraint_dict[beginning][0]

        if c2_x == long_beginning and c1_y == long_end:
            continue

        #print c1_y,c2_x 
        i_c2_x = path.index(c2_x)

        i_c1_y =  path.index(c1_y)

        #print 'HHEHER', i_c1_y, i_c2_x


        if c2_x == long_beginning:
            # print '1', path
            i_c1_x = i_c1_y - 1
            c1_x = path[i_c1_x]
            path.insert(-1, c1_x)
            path.insert(-1, c1_y)
            path.pop(i_c1_y)
            path.pop(i_c1_x)
            # print '5', path
            continue



        if c1_y == long_end:
            # print '1', path
            i_c2_y = i_c2_x+1
            c2_y = path[i_c2_y]
            path.insert(1, c2_y)
            path.insert(1, c2_x)
            path.pop(i_c2_y+2)
            path.pop(i_c2_x+2)
            # print '5', path
            continue

        i_c1_x = i_c1_y - 1
        c1_x = path[i_c1_x]
        i_c2_y = i_c2_x+1
        c2_y = path[i_c2_y]

        if i_c1_y < i_c2_x:
            # print '1', path
            path.insert(i_c2_x, c1_y)
            # print '2', path
            path.pop(i_c1_y)
            # print '3', path
            path.insert(i_c2_x-1, c1_x)
            # print '4', path
            path.pop(i_c1_x)
            # print '5', path

        elif i_c2_x < i_c1_y:
            # print '1lolz', path
            path.insert(i_c1_y, path.pop(i_c2_x))
            path.insert(i_c1_y, path.pop(i_c2_x))
            # print '5', path

    print "after constraints", path

    return path




rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

def rev_comp(string, rev_nuc):
    #rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def check_kmer_overlap(end1,end2):
    i = len(end1)
    while i > 0:
        if end1[-i:] == end2[:i]:
            return i
        i -= 1
    return i

def kmer_overlaps(path, Scaffolds, small_scaffolds, Contigs, small_contigs):

    # get sequences of contigs
    seq_dict = {} # "ctg_sequence" : "scaf_index"
    label_dict_beginnings = {}
    label_dict_ends = {}
    # print path
    for label in path:
        # get each end of 200bp

        # scaffold is forward oriented and beginning
        if (label[1] == 'L' and path.index(label) % 2 == 1):
            #print Scaffolds[label[0]].contigs, Scaffolds[label[0]].contigs[0].direction
            ctg = Scaffolds[label[0]].contigs[0]
            #print ctg.sequence
            if ctg.direction:
                ctg_seq = ctg.sequence[:200]
            else:
                ctg_seq = rev_comp(ctg.sequence[-200:], rev_nuc)
            label_dict_beginnings[label] = ctg_seq

        # scaffold is forward oriented and end
        elif (label[1] == 'R' and path.index(label) % 2 == 0):
            ctg = Scaffolds[label[0]].contigs[-1]
            #print ctg.sequence
            if ctg.direction:
                ctg_seq = ctg.sequence[-200:]
            else:
                ctg_seq = rev_comp(ctg.sequence[:200], rev_nuc)
            label_dict_ends[label] = ctg_seq

        # scaffold is reverse oriented and beginning
        elif (label[1] == 'R' and path.index(label) % 2 == 1):
            ctg = Scaffolds[label[0]].contigs[-1]
            if ctg.direction:
                ctg_seq = rev_comp(ctg.sequence[-200:], rev_nuc)
            else:
                ctg_seq = ctg.sequence[:200]
            label_dict_beginnings[label] = ctg_seq

        # scaffold is reverse oriented and end
        elif (label[1] == 'L' and path.index(label) % 2 == 0):
            ctg = Scaffolds[label[0]].contigs[0]
            if ctg.direction:
                # print "LLLOOOOL"
                ctg_seq = rev_comp(ctg.sequence[-200:], rev_nuc) #[::-1]
            else:
                ctg_seq = ctg.sequence[-200:]
            label_dict_ends[label] = ctg_seq

        else:
            print 'BUG!'

        #label_dict[label] = ctg_seq

    already_added = set()
    constraint_pairs = {}
    constraint_dict = {}

    for label1 in label_dict_ends:
        for label2 in label_dict_beginnings:
            if label1[0] != label2[0]:
                overlap = check_kmer_overlap(label_dict_ends[label1], label_dict_beginnings[label2])
                if overlap > 20:
                    print overlap
                    # print overlap, label1, label2
                    # print label_dict[label1], label_dict[label2]
                    # print
                    if label1 not in already_added and label2 not in already_added:
                        #print 'YEAH', already_added
                        constraint_pairs[label1] = (label2, overlap)
                        constraint_dict[label2] = (label1, overlap)
                        already_added.add(label1)
                        already_added.add(label2)

                    elif label1 in already_added and label2 not in already_added:
                        old_label2 = constraint_pairs[label1][0]
                        if overlap > constraint_pairs[label1][1]:
                            del constraint_pairs[label1]
                            del constraint_dict[old_label2]
                            constraint_pairs[label1] = (label2, overlap)
                            constraint_dict[label2] = (label1, overlap)
                            already_added.add(label1)
                            already_added.add(label2)
                            already_added.remove(old_label2)

                    elif label1 not in already_added and label2 in already_added:
                        old_label1 = constraint_dict[label2][0]
                        if overlap > constraint_dict[label2][1]:
                            del constraint_dict[label2]
                            del constraint_pairs[old_label1]
                            constraint_dict[label2] = (label1, overlap)
                            constraint_pairs[label1] = (label2, overlap)
                            already_added.add(label1)
                            already_added.add(label2)
                            already_added.remove(old_label1)

                    elif label1 in already_added and label2 in already_added:
                        old_label1 = constraint_dict[label2][0]
                        old_label2 = constraint_pairs[label1][0]
                        if overlap > constraint_dict[label2][1] and overlap > constraint_pairs[label1][1]:
                            del constraint_dict[label2]
                            del constraint_pairs[old_label1]
                            del constraint_pairs[label1]
                            del constraint_dict[old_label2]
                            constraint_dict[label2] = (label1, overlap)
                            constraint_pairs[label1] = (label2, overlap)
                            already_added.add(label1)
                            already_added.add(label2)
                            already_added.remove(old_label1)
                            already_added.remove(old_label2)

                    # else:
                    #     constraint_pairs[label1] = (label2, overlap)


    #print already_added

            #''.join(map(lambda ctg: ctg.sequence, Scaffolds[label[0]].contigs[0]))

    # find longest overlaps

    # generate and return constraints on permutations

    return constraint_dict


def permute_path_segment(path, ctg_to_move, ctg_after, valid_connections):
    #print 'path before permute', path
    segment_upper_index = max(path.index((ctg_to_move, 'L')), path.index((ctg_to_move, 'R')))
    segment_min_index = max(path.index((ctg_after, 'L')), path.index((ctg_after, 'R'))) + 1
    segment_lower_index = segment_min_index
    for index in range(segment_upper_index, segment_min_index, -1):
        if path[index][0] == path[index-1][0] or (path[index][0], path[index-1][0]) in valid_connections or (path[index-1][0], path[index][0]) in valid_connections:
            continue
        else:
            segment_lower_index = index

    segment_to_move = path[segment_lower_index: segment_upper_index + 1]
    #print "segment to move", segment_to_move
    path[segment_lower_index: segment_upper_index + 1] = []
    permuted_path = path[0:segment_min_index-2] + segment_to_move + path[segment_min_index-2:]
    #print "path after permute", permuted_path
    return permuted_path


def calculate_min_objective(path, Scaffolds, small_scaffolds, observations, param, valid_permutations, valid_connections):
    current_min_path = copy.deepcopy(path)
    current_min_path_instance, current_min_contigs_to_indexes, current_min_indexes_to_contigs, current_min_index_observations = calculate_path_LP(current_min_path, Scaffolds, small_scaffolds, observations, param)
    ctgs_perm_chosen = None
    #print 'WORK IS DONE'

    for edge_end1, edge_end2 in valid_permutations:
        # switch positions of two contigs
        #print "To permute:", edge_end1, edge_end2
        current_path = copy.deepcopy(path)
        ctg_to_move = edge_end1[0]
        ctg_after = edge_end2[0]
        current_path = permute_path_segment(current_path, ctg_to_move, ctg_after, valid_connections)

        #print 'Current current_path:',current_path
    ## 1 Get a mapping from contigs to indexes (index for contig order in the current current_min_path)
        current_path_instance, current_contigs_to_indexes, current_indexes_to_contigs, current_index_observations = calculate_path_LP(current_path, Scaffolds, small_scaffolds, observations, param, moved_forward=ctg_to_move, shifted=ctg_after)
        if not current_path_instance:
            continue
    ## 3 Check of current current_min_path is better than previous
        #print 'Current objective: {0}, best objective: {1}'.format(current_path_instance.objective, current_min_path_instance.objective)
        # if  -1 < current_path_instance.objective - current_min_path_instance.objective < 1:
        #     print 'Very close:', current_path_instance.objective, current_min_path_instance.objective
        #     print current_path_instance.objective < current_min_path_instance.objective
        #print 'obj:', current_path_instance.objective,current_min_path_instance.objective
        if current_path_instance.objective < current_min_path_instance.objective - 1:
            current_min_path = copy.deepcopy(current_path)
            current_min_path_instance = copy.deepcopy(current_path_instance)
            current_min_contigs_to_indexes = current_contigs_to_indexes
            current_min_indexes_to_contigs = current_indexes_to_contigs
            current_min_index_observations = current_index_observations
            ctgs_perm_chosen = (edge_end1, edge_end2)
    #print "IN FCN:", ctgs_perm_chosen, current_min_path_instance.objective
    return current_min_path, current_min_path_instance, current_min_contigs_to_indexes, current_min_indexes_to_contigs, current_min_index_observations, ctgs_perm_chosen



def estimate_path_gaps(Contigs, small_contigs, path,Scaffolds,small_scaffolds, G_prime, param):
    #print "New path!"
    ## ACCURATE GAP EST HERE

    #print G_prime.subgraph(path).edges(data=True)
    sub_graph = G_prime.subgraph(path)


    #print G_prime.subgraph(path).edges(data=True)
    sub_graph_reduced = filter(lambda x: sub_graph[x[0]][x[1]]['nr_links'] != None and x[0][0] in sub_graph[x[0]][x[1]] , sub_graph.edges())

    observations = dict(map(lambda x: (x, [i+j for i,j in zip(sub_graph[x[0]][x[1]][x[0][0]], sub_graph[x[0]][x[1]][x[1][0]] )]), sub_graph_reduced))
    sub_graph_small_to_large_ctgs = filter(lambda x: sub_graph[x[0]][x[1]]['nr_links'] != None and x[0][0] not in sub_graph[x[0]][x[1]] , sub_graph.edges())
    for c1,c2 in sub_graph_small_to_large_ctgs:
        #observations[(c1,c2)] = (sub_graph[c1][c2]['obs']/ sub_graph[c1][c2]['nr_links'], sub_graph[c1][c2]['nr_links']) #[sub_graph[c1][c2]['obs']/ sub_graph[c1][c2]['nr_links']]*sub_graph[c1][c2]['nr_links']
        n = sub_graph[c1][c2]['nr_links']
        mean_obs = sub_graph[c1][c2]['obs'] / n
        if n > 1:
            std_dev_obs = ((sub_graph[c1][c2]['obs_sq'] - n * mean_obs ** 2) / (n - 1)) ** 0.5
        else:
            std_dev_obs = 0
        obs_list = sub_graph[c1][c2]['observations']
        observations[(c1, c2)] = (mean_obs, n, std_dev_obs, obs_list)

    # for c1,c2 in observations:
    #     if (other_end(c2),other_end(c1)) in observations:
    #         print observations
    # print param.contamination_ratio 



    #########################
    ##########################
    ## Now we have all info that we need to send to module for calculating optimal path


    # only one contig, nothing to permute
    if len(path) <= 4 or not param.contamination_ratio or param.NO_ILP:

        # constraint_pairs = kmer_overlaps(path, Scaffolds, small_scaffolds, Contigs, small_contigs)
        # print 'small', constraint_pairs

        final_path_instance, final_contigs_to_indexes, final_indexes_to_contigs, final_index_observations = calculate_path_LP(path,Scaffolds,small_scaffolds,observations,param)
        final_path = path
        #print final_path

    ## algm here
    else:


        # new algorithm
        # final_path_instance, final_contigs_to_indexes, final_indexes_to_contigs, final_index_observations = calculate_path_LP(path,Scaffolds,small_scaffolds,observations,param)
        # final_path = copy.deepcopy(path)
        # original_path = copy.deepcopy(path)

        # for p in path_permutations_iterator(path, 3, original_path):
        #     print p
        #      ## 1 Get a mapping from contigs to indexes (index for contig order in the current path)
        #     current_path_instance, current_contigs_to_indexes, current_indexes_to_contigs, current_index_observations = calculate_path_LP(p,Scaffolds,small_scaffolds,observations,param)
        #     if not current_path_instance:
        #         continue
        #     ## 3 Check of current path is better than previous            
        #     if current_path_instance.objective < final_path_instance.objective:
        #         final_path = copy.deepcopy(p)
        #         final_path_instance = copy.deepcopy(current_path_instance)
        #         final_contigs_to_indexes = current_contigs_to_indexes 
        #         final_indexes_to_contigs = current_indexes_to_contigs
        #         final_index_observations = current_index_observations
        ############################



        ## overlap idea ##

        # original_path = copy.deepcopy(path)
        # constraint_dict = kmer_overlaps(path, Scaffolds, small_scaffolds, Contigs, small_contigs)
        # print 'big', constraint_dict
        # original_path_with_constraints = path_permutations_with_overlap_constraints(original_path, constraint_dict)

        # final_path_instance, final_contigs_to_indexes, final_indexes_to_contigs, final_index_observations = calculate_path_LP(original_path_with_constraints, Scaffolds, small_scaffolds, observations,param)
        # final_path = copy.deepcopy(original_path_with_constraints)

        # end_constraints = map(lambda x: x[0], constraint_dict.values())
        # # inital valid permutations
        # for i in range(3, len(path) - 1, 2):
        #     # switch positions of two contigs
        #     contig_beginning = original_path[i]
        #     contig_end = original_path[i-2]
        #     if contig_beginning in constraint_dict or contig_end in end_constraints:
        #         continue
        #     else:
        #         # switch positions of two contigs
        #         current_path = copy.deepcopy(final_path)
        #         print current_path
        #         ctg_to_move = original_path[i][0]
        #         contig_after = original_path[i-2][0]
        #         current_path  = permute_path(current_path, ctg_to_move, contig_after)

        #         #print 'Current path:',current_path
        #     ## 1 Get a mapping from contigs to indexes (index for contig order in the current path)
        #         current_path_instance, current_contigs_to_indexes, current_indexes_to_contigs, current_index_observations = calculate_path_LP(current_path,Scaffolds,small_scaffolds,observations,param)
        #         if not current_path_instance:
        #             continue
        #     ## 3 Check of current path is better than previous
        #         #print 'Current objective: {0}, best objective: {1}'.format(current_path_instance.objective, final_path_instance.objective)

        #         if current_path_instance.objective < final_path_instance.objective:
        #             final_path = copy.deepcopy(current_path)
        #             final_path_instance = copy.deepcopy(current_path_instance)
        #             final_contigs_to_indexes = current_contigs_to_indexes
        #             final_indexes_to_contigs = current_indexes_to_contigs
        #             final_index_observations = current_index_observations



        #########################



        if param.FASTER_ILP:

            #### ORIGINAL ILP ##############

            final_path_instance, final_contigs_to_indexes, final_indexes_to_contigs, final_index_observations = calculate_path_LP(path,Scaffolds,small_scaffolds,observations,param)
            final_path = copy.deepcopy(path)
            original_path = copy.deepcopy(path)
            #print 'WORK IS DONE'

            for i in range(3, len(path) - 1, 2):
                # switch positions of two contigs
                current_path = copy.deepcopy(final_path)
                ctg_to_move = original_path[i][0]
                contig_after = original_path[i-2][0]
                current_path  = permute_path(current_path, ctg_to_move, contig_after)

                #print 'Current path:',current_path
            ## 1 Get a mapping from contigs to indexes (index for contig order in the current path)
                current_path_instance, current_contigs_to_indexes, current_indexes_to_contigs, current_index_observations = calculate_path_LP(current_path, Scaffolds, small_scaffolds, observations, param, moved_forward=ctg_to_move, shifted=contig_after)
                if not current_path_instance:
                    continue
            ## 3 Check of current path is better than previous
                #print 'Current objective: {0}, best objective: {1}'.format(current_path_instance.objective, final_path_instance.objective)
                if  -1 < current_path_instance.objective - final_path_instance.objective < 1:
                    print 'Very close:', current_path_instance.objective, final_path_instance.objective
                    print current_path_instance.objective < final_path_instance.objective

                if current_path_instance.objective < final_path_instance.objective - 1:
                    final_path = copy.deepcopy(current_path)
                    final_path_instance = copy.deepcopy(current_path_instance)
                    final_contigs_to_indexes = current_contigs_to_indexes 
                    final_indexes_to_contigs = current_indexes_to_contigs
                    final_index_observations = current_index_observations

        #########################################################################


        ### O(n*(n-1)/2) solution idea ######

        else:
            max_nr_rounds = len(range(3, len(path) - 1, 2))
            current_min_path = copy.deepcopy(path)

            # calculate initial round with all edges
            valid_permutations = set([(current_min_path[i], current_min_path[i-2]) for i in range(3, len(path) - 1, 2)])
            #print "PERM set:", valid_permutations

            # store what are valid adjacent contigs, we need it as info when we permute contigs in a path
            valid_connections = set([(current_min_path[i][0], current_min_path[i-2][0]) for i in range(3, len(path) - 1, 2)])

            final_path, final_path_instance, final_contigs_to_indexes, final_indexes_to_contigs, final_index_observations, ctgs_perm_chosen = calculate_min_objective(current_min_path, Scaffolds, small_scaffolds, observations, param, valid_permutations, valid_connections)

            # In initial step, at least one solution was better than original path
            if ctgs_perm_chosen:
                valid_permutations.remove(ctgs_perm_chosen)
                #print "PERM set:", valid_permutations
                i = 1
                while i < max_nr_rounds:
                    #print 'WE are here!!!!! round:', i, max_nr_rounds
                    min_path, min_solution_path_instance, min_solution_contigs_to_indexes, min_solution_indexes_to_contigs, min_solution_index_observations, ctgs_perm_chosen = calculate_min_objective(final_path,Scaffolds,small_scaffolds,observations,param, valid_permutations, valid_connections)
                    #print "IN WHILE:", ctgs_perm_chosen
                    if ctgs_perm_chosen:
                        valid_permutations.remove(ctgs_perm_chosen)
                    else:
                        break
                    #print "PERM set LOLLL:", valid_permutations
                    # allow smal approximation error of LP solver, hence at least -1 lower obj of new instance is accepted
                    if min_solution_path_instance.objective < final_path_instance.objective - 1:
                        final_path = copy.deepcopy(min_path)
                        final_path_instance = copy.deepcopy(min_solution_path_instance)
                        final_contigs_to_indexes = min_solution_contigs_to_indexes
                        final_indexes_to_contigs = min_solution_indexes_to_contigs
                        final_index_observations = min_solution_index_observations
                    else:
                        break
                    i += 1


        ######################################





    ## 5 Calculate some stats on the path that we have chosen

    param.path_gaps_estimated += len(final_path_instance.gaps)
    path_dict_index = final_path_instance.make_path_dict_for_besst()
    #print 'FINAL:', final_path
    #print 'OBJECTIVE: {0}'.format(final_path_instance.objective)

    # for ctg in final_path:
    #     if ctg[0] in Scaffolds:
    #         for c in Scaffolds[ctg[0]].contigs:
    #             print 'ctg name, scaffold{0}: {1}'.format(ctg[0],c.name)
    #     else:
    #         print 'small ctg name:', small_scaffolds[ctg[0]].contigs[0].name


    #print final_path_instance
    #print path_dict_index
    path_dict = map(lambda x: (final_indexes_to_contigs[x[0].index],final_indexes_to_contigs[x[1].index], path_dict_index[x]), path_dict_index)
    #print path_dict
    # for ctg in final_path_instance.ctgs:
    #     contig = indexes_to_contigs[ctg.index]
    #     ctg.length
    #     ctg.position

    G_ = nx.Graph()
    final_path.insert(0, (final_path[0][0], 'R')) if final_path[0][1] == 'L' else final_path.insert(0, (final_path[0][0], 'L'))
    final_path.insert(len(final_path), (final_path[-1][0], 'R'))  if final_path[-1][1] == 'L' else final_path.insert(len(final_path), (final_path[-1][0], 'L'))
    G_.add_edges_from(zip(final_path[::1], final_path[1::]))

    for c1,c2,gap in path_dict:
        if (c2,'L') in G_[(c1,'L')]:
            G_[(c1,'L')][(c2,'L')]['avg_gap'] = gap
        elif (c2,'R') in G_[(c1,'L')]:
            G_[(c1,'L')][(c2,'R')]['avg_gap'] = gap
        elif (c2,'L') in G_[(c1,'R')]:
            G_[(c1,'R')][(c2,'L')]['avg_gap'] = gap
        elif (c2,'R') in G_[(c1,'R')]:
            G_[(c1,'R')][(c2,'R')]['avg_gap'] = gap
        else:
            print 'Could not find edge!'
            sys.exit()
    #print G_.edges(data=True)

    return(G_, final_path)


def PROBetweenScaf(G_prime, Contigs, small_contigs, Scaffolds, small_scaffolds, param, dValuesTable, Information):
    start_scaf_index = param.scaffold_indexer
    G = nx.Graph()
    for node in G_prime:
        if node[0] in Scaffolds: # meets the length criteria
            G.add_node(node)

    # Filtering and heuristic here to reduce computation if needed O(n^2) in contigs on pathfinder

    #remove all solated contigs
    for node in G.nodes():
        if node in G:
            nbr = G_prime.neighbors(node)[0]
            if len(G_prime.neighbors(node)) == 1 and len(G_prime.neighbors(nbr)) == 1:
                G.remove_nodes_from([node, nbr])


    if len(G.nodes()) / 2.0 > 10000:
        # Too few short contigs compared to long (ratio set to 0.1) or lib ins size + 2*std_dev - 2*read_len < 200 ) and too many large contigs (> 10 000) do not enter path extension algm since to low payoff:

        if len(small_scaffolds) / float(len(Scaffolds)) < 0.1 and not param.no_score:
            print >> Information, "Did not enter path seartching algorithm between scaffolds due to too small fraction of small scaffolds, fraction were: ", len(small_scaffolds) / float(len(Scaffolds))
            return(start_scaf_index)

    ########### Find paths between scaffolds here ###############

    # Multi Processing (if available), check nr of available cores
    num_cores = multiprocessing.cpu_count()
    #TODO: If we get too many paths back and run into memory issues we could change so that only paths with score over 0 are stored in ELS module
    if param.multiprocess and num_cores > 1:
        import workerprocess
        import heapq
        print >> Information, 'Entering ELS.BetweenScaffolds parallelized with ', num_cores, ' cores.'
        start = time.time()
        # load up work queue
        work_queue = multiprocessing.Queue()
        end = set()
        for node in G:
            end.add(node)
        nodes = G.nodes()
        nr_jobs = len(nodes)
        chunk = nr_jobs / (num_cores)
        counter = 0
        nr_processes = 0
        # partition equally many nodes in G to each core
        while counter < nr_jobs:
            work_queue.put((set(nodes[counter:counter + chunk]), G_prime, end, param))
            nr_processes += 1
            print >> Information, 'node nr', counter, 'to', counter + chunk - 1, 'added'
            #print work_queue.get()
            counter += chunk

        # create a queue to pass to workers to store the results
        result_queue = multiprocessing.Queue()

        # spawn workers
        while not work_queue.empty():
            worker = workerprocess.Worker(work_queue.get(), result_queue)
            worker.start()

        # collect the results off the queue
        results = []
        for i in range(nr_processes):
            res = result_queue.get()
            results.append(res)

        def wrapper(func, args):
            return(func(*args))
        all_paths_sorted_wrt_score_itr = wrapper(heapq.merge, results) #tot_result
        all_paths_sorted_wrt_score = [i for i in all_paths_sorted_wrt_score_itr]
        elapsed = time.time() - start
        print >> Information, "Elapsed time multiprocessing: ", elapsed

    else:
        start = time.time()
        end = set()
        for node in G:
            end.add(node)
        iter_nodes = end.copy()
        print >> Information, 'Entering ELS.BetweenScaffolds single core'
        all_paths_sorted_wrt_score = ELS.BetweenScaffolds(G_prime, end, iter_nodes, param)
        elapsed = time.time() - start
        print >> Information, "Elapsed time single core pathfinder: ", elapsed
        print "Elapsed time single core pathfinder: ", elapsed

    ################################################################

    nr_paths = len(all_paths_sorted_wrt_score)
    start_end_node_update_storage = {}
    print >> Information, '{0} paths detected are with score greater or equal to {1} '.format(nr_paths, param.score_cutoff)
    for sublist in reversed(all_paths_sorted_wrt_score):
        path = sublist[2]
        bad_links = sublist[1]
        score = sublist[0]
        path_len = sublist[3]
        if nr_paths <= 100000:
            print >> Information, 'Path: path length: {0}, nr bad links: {1}, score: {2} '.format((path_len - 2) / 2.0, bad_links, score)

        ## Need something here that keeps track on which contigs that are added to Scaffolds so that a
        ## contig is only present once in each path

        #print start_end_node_update_storage
        # Either a small contig/scaffold has been included in a path earlier and thus has moved it's object to Scaffolds (and changed index) 
        small_scaf_is_already_in = 0
        for scaf_ in path[1:-1]:
            if scaf_[0] not in small_scaffolds:
                small_scaf_is_already_in = 1
                #print 'At least one of the contigs is already in another scaffold'
                break
        if small_scaf_is_already_in:
            continue

        # A very special corner case (circular paths)
        if path[0][0] not in Scaffolds and path[-1][0] not in Scaffolds:
            try:
                strt = start_end_node_update_storage[path[0]][0]
                nd = start_end_node_update_storage[path[-1]][0]
                if strt[0] == nd[0]:
                    print >> Information, 'Rare case (circular paths) detected and treated. '
                    continue
            except KeyError:
                pass

        # Or a large scaffold/contig has changed scaffold index due to one of it's sides is present in another path (we still want to allow for paths from the other side)
        case1 = 0
        case2 = 0
        if path[0][0] not in Scaffolds:
            if path[0] in start_end_node_update_storage:
                case1 = 1
            else:
                print >> Information, 'Beginning is already in path'
                continue

        if path[-1][0] not in Scaffolds:
            if path[-1] in start_end_node_update_storage:
                case2 = 1
            else:
                print >> Information, 'End is already in path'
                continue


        original_start_node = path[0]

        if path[0][0] not in Scaffolds:
            #large scaffold has changed index before. This suggested path is however from it's other side
            node_to_remove1 = path[0]
            path[0] = start_end_node_update_storage[node_to_remove1][0]
            #update the node on the other end of the end scaffold to point at the newest index
            node_to_refresh1 = start_end_node_update_storage[node_to_remove1][1]
            #print 'Enter 1'
            try:
                node_ptr = start_end_node_update_storage[ path[-1] ][1]
                #print '1.1', node_ptr,start_end_node_update_storage[ path[-1] ]
            except KeyError:
                other_side = 'L' if path[-1][1] == 'R' else 'R'
                node_ptr = (path[-1][0], other_side)
                #print '1.2', node_ptr, path[-1]
            start_end_node_update_storage[node_to_refresh1] = [(param.scaffold_indexer + 1, 'L'), node_ptr  ]
            #path pointer can be accesed only once needs to be destroyed after
            del start_end_node_update_storage[node_to_remove1]


        if path[-1][0] not in Scaffolds:
            #large scaffold has changed index before. This suggested path is however from it's other side
            #print 'case2.2'
            node_to_remove2 = path[-1]
            path[-1] = start_end_node_update_storage[node_to_remove2][0]
            #update the node on the other end of the end scaffold to point at the newest index
            node_to_refresh2 = start_end_node_update_storage[node_to_remove2][1]
            #print 'Enter 2'
            try:
                node_ptr = start_end_node_update_storage[ original_start_node ][1]
                #print '2.1', node_ptr, start_end_node_update_storage[ original_start_node ]
            except KeyError:
                other_side = 'L' if original_start_node[1] == 'R' else 'R'
                node_ptr = (original_start_node[0], other_side)
                #print '2.2', node_ptr,original_start_node          
            start_end_node_update_storage[node_to_refresh2] = [(param.scaffold_indexer + 1, 'R'), node_ptr ]
            #path pointer can be accesed only once needs to be destroyed after
            del start_end_node_update_storage[node_to_remove2]


        # Here we update the contigs that lies in small_contigs to Contigs. We need to do this here because
        # we update the scaffold index below

        # move all contig and scaffold objects from "small" structure to large structure to fit with UpdateInfo structure

        small_scafs = map(lambda i: path[i], filter(lambda i: i % 2 == 1, range(len(path) - 1)))
        for item in small_scafs:
            scaf_obj = small_scaffolds[item[0]]
            Scaffolds[item[0]] = scaf_obj
            cont_objects = scaf_obj.contigs
            for obj_ in cont_objects:
                ctg_name = obj_.name
                Contigs[ctg_name] = obj_
                del small_contigs[ctg_name]
            del small_scaffolds[item[0]]
        ## Here we do the "joining of two scaffolds with the new path if no contig/scaffold is present
        ## in another path, we need to update "Scaffolds" structure here along as we go in order for
        ## the above dublette checking function to work


        G_, path = estimate_path_gaps(Contigs, small_contigs, path,Scaffolds,small_scaffolds, G_prime,param)


        # #make the path a small linear graph
        # G_ = nx.Graph()
        # path.insert(0, (path[0][0], 'R')) if path[0][1] == 'L' else path.insert(0, (path[0][0], 'L'))
        # path.insert(len(path), (path[-1][0], 'R'))  if path[-1][1] == 'L' else path.insert(len(path), (path[-1][0], 'L'))
        # G_.add_edges_from(zip(path[::1], path[1::]))

        start_end_node_update_storage[path[0]] = 0
        start_end_node_update_storage[path[-1]] = 0


        start = path[0]
        end = path[-1]
        prev_node = ('', '')
        pos = 0
        scaffold_length = 0
        contig_list = []
        param.scaffold_indexer += 1

        values = UpdateInfo(G_, Contigs, small_contigs, Scaffolds, small_scaffolds, start, prev_node, pos, contig_list, scaffold_length, dValuesTable, param)
        #print values
        while len(values) !=2:
            values = UpdateInfo(G_,*values)
    
        (contig_list, scaffold_length) = values[0],values[1]
        S = Scaffold.scaffold(param.scaffold_indexer, contig_list, scaffold_length)  #Create the new scaffold object 

        if nr_paths <= 100000:
            print >> Information, 'Path taken! path length: {0}, nr bad links: {1}, score: {2} '.format((path_len - 2) / 2.0, bad_links, score)

        Scaffolds[S.name] = S        #include in scaffold library
        #add the new scaffold object to G_prime

        G_prime.add_node((S.name, 'L'))  #start node
        G_prime.add_node((S.name, 'R'))  # end node
        G_prime.add_edge((S.name, 'L'), (S.name, 'R'), nr_links=None)
        for nbr in G_prime.neighbors(start):
            nr_links_ = G_prime[start][nbr]['nr_links']
            if nr_links_:
                obs_ = G_prime[start][nbr]['obs']
                observations_ = G_prime[start][nbr]['observations']
                #G_prime.add_edge((S.name, 'L'), nbr, nr_links=nr_links_, obs=obs_)
                obs_sq_ = G_prime[start][nbr]['obs_sq']
                G_prime.add_edge((S.name, 'L'), nbr, nr_links=nr_links_, obs=obs_, obs_sq=obs_sq_, observations=observations_)

        for nbr in G_prime.neighbors(end):
            nr_links_ = G_prime[end][nbr]['nr_links']
            if nr_links_:
                obs_ = G_prime[end][nbr]['obs']
                observations_ = G_prime[end][nbr]['observations']
                # G_prime.add_edge((S.name, 'R'), nbr, nr_links=nr_links_, obs=obs_)
                obs_sq_ = G_prime[end][nbr]['obs_sq']
                G_prime.add_edge((S.name, 'R'), nbr, nr_links=nr_links_, obs=obs_, obs_sq=obs_sq_, observations=observations_)

        #remove the old scaffold objects from G_prime
        G_prime.remove_nodes_from(path)

        #updated beginning
        if case1 and not case2:
            start_end_node_update_storage[node_to_refresh1] = [(S.name, 'L'), path[-1] ]
            start_end_node_update_storage[path[-1]] = [(S.name, 'R'), node_to_refresh1 ]
        elif case2 and not case1:
            start_end_node_update_storage[path[0]] = [(S.name, 'L'), node_to_refresh2 ]
            start_end_node_update_storage[node_to_refresh2] = [(S.name, 'R'), path[0] ]
        elif case1 and case2:
            start_end_node_update_storage[node_to_refresh1] = [(S.name, 'L'), node_to_refresh2 ]
            start_end_node_update_storage[node_to_refresh2] = [(S.name, 'R'), node_to_refresh1 ]
        else:
            start_end_node_update_storage[path[0]] = [(S.name, 'L'), path[-1] ]
            start_end_node_update_storage[path[-1]] = [(S.name, 'R'), path[0] ]

    return(start_scaf_index)

def GiveLinkConnection(Contigs, contig_objects1, contig_objects2, side1, side2, nr_links):
    if side1 == 'R' and side2 == 'L':
        max_pos = 0
        for contig in contig_objects1:
            if contig.position >= max_pos:
                linking_contig1 = contig
                max_pos = contig.position
        min_pos = sys.maxint
        for contig in contig_objects2:
            if contig.position <= min_pos:
                linking_contig2 = contig
                min_pos = contig.position

        linking_contig1.links[linking_contig2.name] = nr_links
        linking_contig2.links[linking_contig1.name] = nr_links


    elif side1 == 'L' and side2 == 'R':
        max_pos = 0
        for contig in contig_objects2:
            if contig.position >= max_pos:
                linking_contig2 = contig
                max_pos = contig.position

        min_pos = sys.maxint
        for contig in contig_objects1:
            if contig.position <= min_pos:
                linking_contig1 = contig
                min_pos = contig.position
        linking_contig1.links[linking_contig2.name] = nr_links
        linking_contig2.links[linking_contig1.name] = nr_links


    elif side1 == 'R' and side2 == 'R':
        max_pos = 0
        for contig in contig_objects1:
            if contig.position >= max_pos:
                linking_contig1 = contig
                max_pos = contig.position

        max_pos = 0
        for contig in contig_objects2:
            if contig.position >= max_pos:
                linking_contig2 = contig
                max_pos = contig.position
        linking_contig1.links[linking_contig2.name] = nr_links
        linking_contig2.links[linking_contig1.name] = nr_links


    elif side1 == 'L' and side2 == 'L':
        min_pos = sys.maxint
        for contig in contig_objects1:
            if contig.position <= min_pos:
                linking_contig1 = contig
                min_pos = contig.position
        min_pos = sys.maxint
        for contig in contig_objects2:
            if contig.position <= min_pos:
                linking_contig2 = contig
                min_pos = contig.position

        linking_contig1.links[linking_contig2.name] = nr_links
        linking_contig2.links[linking_contig1.name] = nr_links

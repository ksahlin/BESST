'''
Created on Jun 17, 2013

@author: ksahlin
'''
'''
Created on Mar 15, 2013

@author: ksahlin
'''

import os
import sys
import time
from collections import defaultdict
import itertools
from operator import itemgetter
import gc

import networkx as nx

import smith_waterman
import wrapper_sw
import output_contigs



# generate the reverse complement of a sequence
def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def create_kmers(k, m, sequence, acc, k_mer_hash):
    #get length of contig
    read_length = len(sequence)
    # k_mer_hash =(contig accession, beginning/end, forw/rev)
    # beginning = 0
    # end = 1
    # forward = 0
    # reverse complement = 1

    contig_begin = sequence[:k ]
    contig_end = sequence[read_length - k:]

    # iterate over all m-mers in the k subsequence of the contig
    for i in range(len(contig_begin) - m + 1):
        # fill up k-mer hash table with k-mer as key and 
        # tuple (contig accession, beginning/end, forw/rev) as value
        k_mer_hash[contig_begin[i : m + i]].append((acc, 0 , 0))
        k_mer_hash[reverse_complement(contig_begin[i : m + i])].append((acc, 0 , 1))
        k_mer_hash[contig_end[i : m + i]].append((acc, 1 , 0))
        k_mer_hash[reverse_complement(contig_end[i : m + i])].append((acc, 1 , 1))
        # fill up contig hash table with contig accession as key and 
        # k-mer hash table  (contig accession, beginning/end, forw/rev) as value

    return()

def get_kmer_index(contigs, k, m):
    k_mer_hash = defaultdict(list)
    cntr = 0
    for accession, sequence in contigs.iteritems():
        cntr += 1
        create_kmers(k, m, sequence, accession, k_mer_hash)
        if cntr % 100000 == 0:
            print 'Creating k-mers of ', cntr, ' contigs.'
    return(k_mer_hash)

# given an iterable of pairs return the key corresponding to the greatest value
def argmax(pairs):
    return max(pairs, key=itemgetter(1))[0]

# given an iterable of values return the index of the greatest value
def argmax_index(values):
    return argmax(enumerate(values))

def is_haplotype_region(contigs, haplotypes, similarity_threshold):

    # Check if only one contig on this side, it should then be merged if
    # the other region is haplotypic. We therefore return True

    if len(haplotypes) == 1:
        return(1)

    # Check if no contig on this side, this could be the case if e.g. low coverage on one s
    # side. We want to merge the other side anyway if it is similar. We therefore return True
    if len(haplotypes) == 0:
        print 'No contig on one side!'
        return(1)
    seqs = itertools.combinations(haplotypes, 2)
    score = 0
    nr_comparisons = 0
    nr_matching = 0
    for pair in seqs:
        nr_comparisons += 1
        len_ctg1 = len(contigs[pair[0][0]])
        len_ctg2 = len(contigs[pair[1][0]])
        if pair[0][2] == 1:
            contig1 = reverse_complement(contigs[pair[0][0]])
        else:
            contig1 = contigs[pair[0][0]]

        if pair[1][2] == 1:
            contig2 = reverse_complement(contigs[pair[1][0]])
        else:
            contig2 = contigs[pair[1][0]]

        score, max_i, max_j = wrapper_sw.wrap_sw(contig1, contig2, 0, 0)

        #print score, max_i, max_j

        ## TODO:
        # Wrong in formula!! need too make sure that we start aligning from the similar k-mer 
        # and then take min(max_i, max_j) or max(....)?? Figure out
        procent_similarity = score / float(max(min(max_i, max_j), min(len_ctg1, len_ctg2)))  #float(min(max_i, max_j, len_ctg1, len_ctg2))
        print 'Score', procent_similarity, 'length alignment:', float(min(max_i, max_j, len_ctg1, len_ctg2)), 'length ctgs:', len_ctg1, len_ctg2
        if procent_similarity > similarity_threshold:
            nr_matching += 1

    if nr_matching == nr_comparisons:
        return(1)
    else:
        return()

def graph_updater(connectedness_graph, contigs, nodes, remove_nodes):
    # One side might be empty this is a special case. (see comments in "is_haplotype_region" function)
    if len(nodes) == 0:
        return False


    # Remove all reverse complement nodes (they are redundant) after we
    # have classified the region (they are removed after the loop in search_regions)
    nodes_rev_comp = map(lambda x: (x[0], x[1], 1 - x[2]) , nodes)
    remove_nodes.update(nodes_rev_comp)
    longest_contig_index = argmax_index([len(contigs[ctg[0]]) for ctg in nodes])
    hapl_lengths = map(lambda x: len(contigs[x[0]]), nodes)

    # Need to pick the longest haplotype. This is a corner case where
    # there are at least two longest contigs. We then need to stick to one consistently
    if hapl_lengths.count(max(hapl_lengths)) > 1:
        #print 'OMGGGGG Same length of haplotypes (need to choose one):'
        #print nodes
        pass

        index_ = 0
        for hapl in nodes:
            if len(contigs[hapl[0]]) == max(hapl_lengths):
                try:
                    connectedness_graph.node[hapl]['r']
                except KeyError:
                    longest_contig_index = index_
            index_ += 1


    #print 'INDEX:', longest_contig_index
    for index, ctg in enumerate(nodes):
        if index != longest_contig_index:
            connectedness_graph.node[ctg]['r'] = 1

            # need to update other end of contig as well
            # if it was added to the k-mer connection graph
            # note: not all k-mers of contig ends need to have
            # multiple hits. If only one hit, it is not added.
            # therefore, not all contig ends might be represented in 
            # connectedness_graph
            try:
                connectedness_graph.node[(ctg[0], 1 - ctg[1], 1)]['r'] = 1
            except KeyError:
                pass
            try:
                connectedness_graph.node[(ctg[0], 1 - ctg[1], 0)]['r'] = 1
            except KeyError:
                pass

        else:
            ctg_to_join = ctg




    return(ctg_to_join)


def haplotype_detect(connectedness_graph, contigs, nodes, similarity_threshold, remove_nodes):

    ## check to see that node has not been merged previously
    #if not all([ctg[0] in contigs for ctg in nodes]):
    #    return(0)

    ctg_left, ctg_right = None, None
    before_kmer = filter(lambda x: (x[1] + x[2]) % 2 == 1, nodes)
    if is_haplotype_region(contigs, before_kmer, similarity_threshold):
        ctg_left = graph_updater(connectedness_graph, contigs, before_kmer, remove_nodes)



    after_kmer = filter(lambda x: (x[1] + x[2]) % 2 == 0, nodes)
    if is_haplotype_region(contigs, after_kmer, similarity_threshold):
        ctg_right = graph_updater(connectedness_graph, contigs, after_kmer, remove_nodes)

    if ctg_right:
        for nbr in connectedness_graph.neighbors(ctg_right):
            if nbr in nodes:
                connectedness_graph[ctg_right][nbr]['a'] = 1

    if ctg_left:
        for nbr in connectedness_graph.neighbors(ctg_left):
            if nbr in nodes:
                connectedness_graph[ctg_left][nbr]['a'] = 1

    if ctg_left and ctg_right:
        if connectedness_graph.has_edge(ctg_left, ctg_right):
            connectedness_graph[ctg_left][ctg_right]['m'] = 1
        else:
            connectedness_graph.add_edge(ctg_left, ctg_right, m=1)




    # Save information of what contigs that should be saved, 
    # what contigs that shoud be removed and what contigs that 
    # should be "joined" with SW. Also store the coordinates for merging


    return()

def search_regions(k_mer_hash, contigs, similarity_threshold):

    # create graph with connections
    cntr = 0
    cntr2 = 0
    #cntr3 = 0
    connectedness_graph = nx.Graph()
    print 'Number of m-mers in m-mer hash table: ', len(k_mer_hash)

    # We will not look at nodes with more than 10 triangles = clique < 6 nodes = 
    filtered_list = filter(lambda x: 1 < len(x) < 6, k_mer_hash.values())
    #print filtered_list
    for ctgs in filtered_list:
        cntr2 += 1
        #cntr3 += 1
        #if cntr3 >= 7000000:
        #    break
        if cntr2 % 100000 == 0:
            print cntr2, ' items have been treated.'
        #print k_mer
        if len(ctgs) > 1:
            cntr += 1
        #    print 'here'
            connections = itertools.combinations(ctgs, 2)
            connectedness_graph.add_edges_from(connections)

            if cntr % 10000 == 0:
                print 'Constructed graph from ', cntr, ' k-mers.', 'nr ctgs for k-mer: ', len(ctgs)
                #print 'Total nr of nodes: ', len(connectedness_graph.nodes())
                #print 'Total nr of edges: ', len(connectedness_graph.edges())

    # Look for interesting (looking like split haplotypes) regions in connectedness_graph
    interesting_regions = set([1, 3, 6])
    visited_nodes = set()
    remove_nodes = set()
    cntr = 0
    print 'len', len(contigs)
    for node in connectedness_graph:
        # We modify (remove nodes) as we go so we need to check if the node is still there
        if node not in visited_nodes:
            # visited nodes are both the actual visited triangle and its reverse complement
            if nx.algorithms.triangles(connectedness_graph, node) in interesting_regions:
                nodes = connectedness_graph.neighbors(node) + [node]
                # add speed up here, mark visited nodes / triangles and there reverse complements
                nodes_rev_comp = map(lambda x: (x[0], x[1], 1 - x[2]) , nodes)
                #print nodes , nodes_rev_comp
                visited_nodes.update(nodes)
                visited_nodes.update(nodes_rev_comp)
                cntr += 1
                if cntr % 10 == 0:
                    print cntr
                #if cntr > 20:
                #    break
                haplotype_detect(connectedness_graph, contigs, nodes, similarity_threshold, remove_nodes)



    print 'REMOVED:', len(remove_nodes)

    # Remove all reverse complement nodes to regions that have been classified as haplotypic regions
    connectedness_graph.remove_nodes_from(remove_nodes)


    # Remove all connections that did not pass as haplotyic regions 
    # They should have no key attributes
    rmv_edges = []
    for edge in connectedness_graph.edges_iter():
        if not connectedness_graph[edge[0]][edge[1]]:
            rmv_edges.append(edge)
    connectedness_graph.remove_edges_from(rmv_edges)

    rmv_nodes = []
    for node in connectedness_graph:
        if connectedness_graph.neighbors(node) == 0:
            rmv_nodes.append(node)
    connectedness_graph.remove_nodes_from(rmv_nodes)

    for node in connectedness_graph:
        try:
            connectedness_graph.node[node]['r']
        except KeyError:
            connectedness_graph.node[node]['r'] = 0


    print 'nr nodes left:', len(connectedness_graph.nodes())
    print 'nr edges left:', len(connectedness_graph.edges())

    return(connectedness_graph)


def main(output_dir, contigs, k, m, similarity_threshold):
    # start timing
    start_time = time.time()
    k_mer_hash = get_kmer_index(contigs, k, m)
    print 'Done with m-mer graph.'
    connectedness_graph = search_regions(k_mer_hash, contigs, similarity_threshold)
    output_contigs.generate_fasta(connectedness_graph, contigs, output_dir, k)
    # stop timing
    end_time = time.time()
    # print timing
    print "time:", end_time - start_time
    return()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Outputs k-mers for a given input sequence')
    parser.add_argument('k', type=int, help='assembly k-mer size.')
    parser.add_argument('m', type=int, help='Length of matching seeds (m-mer size) <= k-mer.')
    parser.add_argument('t', type=float, default=0.7, help='Procent similarity threshold (number between 0 and 1)')
    #parser.add_argument('sequence', type=str, help='Sequence to split into k-mers')
    args = parser.parse_args()
    contigs = {'ctg1' : 'TTTGCGCGTCGTGCTGCAAGGTTGTTAC', 'ctg2' : 'AAGGTTGTTACTTTTTTTGGGGGGTTTTGAAACAACCAA',
                'ctg3' : 'AAGGTTGTTACTCGTAAATAAAACAACCAA', 'ctg4' : 'AAACAACCAAGGTGTGGTGTAGTCGTCGTGTCGTC'}
    # start timing
    start_time = time.time()

    k_mer_hash = get_kmer_index(contigs, args.k, args.m)
    search_regions(k_mer_hash, contigs, args.t)


    # stop timing
    end_time = time.time()


    # print timing
    print "time:", end_time - start_time


'''
Created on Jun 17, 2013

@author: ksahlin
'''
'''
Created on Mar 15, 2013

@author: ksahlin
'''

import sys
import time
from collections import defaultdict
import itertools
from operator import itemgetter

import networkx as nx

import smith_waterman
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
    for accession, sequence in contigs.iteritems():
        create_kmers(k, m, sequence, accession, k_mer_hash)
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

        score, max_i, max_j, traceback_matrix = smith_waterman.SW(contig1, contig2, -1, -1)
        #print score, max_i, max_j
        procent_similarity = score / float(min(max_i, max_j, len_ctg1, len_ctg2))
        #print 'Score', procent_similarity
        if procent_similarity > similarity_threshold:
            nr_matching += 1

    if nr_matching == nr_comparisons:
        return(1)
    else:
        return()

def graph_updater(connectedness_graph, contigs, nodes, remove_nodes):

    # Remove all reverse complement nodes (they are redundant) after we
    # have classified the region (they are removed after the loop in search_regions)
    nodes_rev_comp = map(lambda x: (x[0], x[1], 1 - x[2]) , nodes)
    remove_nodes.update(nodes_rev_comp)

    longest_contig_index = argmax_index([len(contigs[ctg[0]]) for ctg in nodes])
    for index, ctg in enumerate(nodes):
        if index != longest_contig_index:
            connectedness_graph.node[ctg]['r'] = 1
        else:
            ctg_to_join = ctg


    return(ctg_to_join)


def haplotype_detect(connectedness_graph, contigs, nodes, similarity_threshold, remove_nodes):

    ## check to see that node has not been merged previously
    if not all([ctg[0] in contigs for ctg in nodes]):
        return(0)

    ctg_left, ctg_right = None, None
    before_kmer = filter(lambda x: (x[1] + x[2]) % 2 == 1, nodes)
    if is_haplotype_region(contigs, before_kmer, similarity_threshold):
        ctg_left = graph_updater(connectedness_graph, contigs, before_kmer, remove_nodes)



    after_kmer = filter(lambda x: (x[1] + x[2]) % 2 == 0, nodes)
    if is_haplotype_region(contigs, after_kmer, similarity_threshold):
        ctg_right = graph_updater(connectedness_graph, contigs, after_kmer, remove_nodes)

    if ctg_left and ctg_right:
        #print connectedness_graph[ctg_left]
        #print connectedness_graph[ctg_right]
        #print ctg_left, ctg_right
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
    connectedness_graph = nx.Graph()
    for k_mer in k_mer_hash:
        #print k_mer
        if len(k_mer_hash[k_mer]) > 1:
        #    print 'here'
            #ctgs = map(lambda x : x[0], k_mer_hash[k_mer])
            #for ctg in contig_hash[k_mer]:

            connections = itertools.combinations(k_mer_hash[k_mer], 2)
            connectedness_graph.add_edges_from(connections)

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
                print cntr
                if cntr > 20:
                    break
                haplotype_detect(connectedness_graph, contigs, nodes, similarity_threshold, remove_nodes)


    connectedness_graph.remove_nodes_from(remove_nodes)

    for node in connectedness_graph:
        try:
            connectedness_graph.node[node]['r']
        except KeyError:
            connectedness_graph.node[node]['r'] = 0


    print 'nr nodes left:', len(connectedness_graph.nodes())
    print 'nr edges left:', len(connectedness_graph.edges())


    #triangle_connections = nx.algorithms.triangles(connectedness_graph)
#    for node, nr_triangles in triangle_connections.iteritems():
#        if nr_triangles == 1 or nr_triangles == 10:
#            print node, connectedness_graph.neighbors(node)
#            max_score, max_i, max_j, traceback_matrix = smith_waterman.SW(contigs[node[0]], contigs[connectedness_graph.neighbors(node)[0][0]], -1, -1)



    return(connectedness_graph)


def main(output_dir, contigs, k, m, similarity_threshold):
    k_mer_hash = get_kmer_index(contigs, k, m)
    connectedness_graph = search_regions(k_mer_hash, contigs, similarity_threshold)
    output_contigs.generate_fasta(connectedness_graph, contigs, output_dir)
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


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


class HaplotypeGraph():
    def __init__(self):
        self.graph
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
        print score, max_i, max_j
        procent_similarity = score / float(min(max_i, max_j, len_ctg1, len_ctg2))
        print 'Score', procent_similarity
        if procent_similarity > similarity_threshold:
            nr_matching += 1

    if nr_matching == nr_comparisons:
        return(1)
    else:
        return()



def decide_if_haplotype_region(contigs, nodes, similarity_threshold):

    ## check to see that node has not been merged previously
    if not all([ctg[0] in contigs for ctg in nodes]):
        return(0)

    before_kmer = filter(lambda x: (x[1] + x[2]) % 2 == 1, nodes)
    if is_haplotype_region(contigs, before_kmer, similarity_threshold):
        longest_contig_index = argmax_index([len(contigs[ctg[0]]) for ctg in before_kmer])

        for index, ctg in enumerate(before_kmer):
            if index != longest_contig_index:
                pass
                #del contigs[ctg[0]]


    after_kmer = filter(lambda x: (x[1] + x[2]) % 2 == 0, nodes)
    if is_haplotype_region(contigs, after_kmer, similarity_threshold):
        longest_contig_index = argmax_index([len(contigs[ctg[0]]) for ctg in after_kmer])

        for index, ctg in enumerate(after_kmer):
            if index != longest_contig_index:
                pass
                #del contigs[ctg[0]]




    return()

def search_regions(k_mer_hash, contigs, similarity_threshold):
    connectedness_graph = nx.Graph()
    for k_mer in k_mer_hash:
        #print k_mer
        if len(k_mer_hash[k_mer]) > 1:
        #    print 'here'
            #ctgs = map(lambda x : x[0], k_mer_hash[k_mer])
            #for ctg in contig_hash[k_mer]:

            connections = itertools.combinations(k_mer_hash[k_mer], 2)
            connectedness_graph.add_edges_from(connections)
    #print connectedness_graph.edges()

    interesting_regions = set([1, 3, 6])
    visited_nodes = set()
    for node in connectedness_graph:
        if node not in visited_nodes:
            #print nx.algorithms.triangles(connectedness_graph, node)
            visited_nodes.add(node)
            if nx.algorithms.triangles(connectedness_graph, node) in interesting_regions:
                nodes = connectedness_graph.neighbors(node) + [node]
                if decide_if_haplotype_region(contigs, nodes, similarity_threshold):
                    pass



    #triangle_connections = nx.algorithms.triangles(connectedness_graph)
#    for node, nr_triangles in triangle_connections.iteritems():
#        if nr_triangles == 1 or nr_triangles == 10:
#            print node, connectedness_graph.neighbors(node)
#            max_score, max_i, max_j, traceback_matrix = smith_waterman.SW(contigs[node[0]], contigs[connectedness_graph.neighbors(node)[0][0]], -1, -1)



    return()


def main(contigs, k, m, similarity_threshold):
    k_mer_hash = get_kmer_index(contigs, k, m)
    search_regions(k_mer_hash, contigs, similarity_threshold)
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


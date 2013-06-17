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
import smith_waterman

import networkx as nx

# generate the reverse complement of a sequence
def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def create_kmers(k, m, sequence, acc, k_mer_hash):
    #get length of contig
    read_length = len(sequence)
    # k_mer_hash =(contig accession, beginning/end, forw/rev)
    # forward = 0
    # reverse complement = 1
    # beginning = 1
    # end = 0
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

def search_regions(k_mer_hash, contigs):
    connectedness_graph = nx.Graph()
    for k_mer in k_mer_hash:
        print k_mer
        if len(k_mer_hash[k_mer]) > 1:
            print 'here'
            #ctgs = map(lambda x : x[0], k_mer_hash[k_mer])
            #for ctg in contig_hash[k_mer]:

            connections = itertools.combinations(k_mer_hash[k_mer], 2)
            connectedness_graph.add_edges_from(connections)
    #print connectedness_graph.edges()

    for node in connectedness_graph:
        print nx.algorithms.triangles(connectedness_graph, node)

    #triangle_connections = nx.algorithms.triangles(connectedness_graph)
#    for node, nr_triangles in triangle_connections.iteritems():
#        if nr_triangles == 1 or nr_triangles == 10:
#            print node, connectedness_graph.neighbors(node)
#            max_score, max_i, max_j, traceback_matrix = smith_waterman.SW(contigs[node[0]], contigs[connectedness_graph.neighbors(node)[0][0]], -1, -1)



    return()


def main(contigs, k, m):
    k_mer_hash = get_kmer_index(contigs, k, m)
    search_regions(k_mer_hash, contigs)
    return()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Outputs k-mers for a given input sequence')
    parser.add_argument('k', type=int, help='assembly k-mer size.')
    parser.add_argument('m', type=int, help='Length of matching seeds (m-mer size) <= k-mer.')
    #parser.add_argument('sequence', type=str, help='Sequence to split into k-mers')
    args = parser.parse_args()
    contigs = {'ctg1' : 'TTTGCGCGTCGTGCTGCAAGGTTGTTAC', 'ctg2' : 'AAGGTTGTTACTTTTTTTGGGGGGTTTTGAAACAACCAA',
                'ctg3' : 'AAGGTTGTTACTCGTAAATAAAACAACCAA', 'ctg4' : 'AAACAACCAAGGTGTGGTGTAGTCGTCGTGTCGTC'}
    # start timing
    start_time = time.time()

    k_mer_hash = get_kmer_index(contigs, args.k, args.m)
    search_regions(k_mer_hash, contigs)


    # stop timing
    end_time = time.time()


    # print timing
    print "time:", end_time - start_time


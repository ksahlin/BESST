'''
Created on Jun 19, 2013

@author: ksahlin
'''

import os

import networkx as nx

def generate_fasta(graph, contigs, output_dir):

    haplotypes = open(os.path.join(output_dir, 'haplotypes.fa'), 'w')
    nodes_to_remove = set()
    for node in graph:
        if graph.node[node]['r']:
            nodes_to_remove.add(node)
            print 'remove'

    graph.remove_nodes_from(nodes_to_remove)
    for node in nodes_to_remove:
        print >> haplotypes, '>' + node[0].strip() + '\n' + contigs[node[0]]


    print len(nx.connected_components(graph))

    for component in nx.connected_components(graph):
        pass

#        try:
#            graph.edge[edge[0]][edge[1]]['m']
#            print 'merge'
#        except KeyError:
#            pass

    return()

'''
Created on Mar 11, 2013

@author: ksahlin
'''


import networkx as nx


def HaplotypicRegions(G, G_prime, Contigs, Scaffolds, param, Information):
    CB = nx.cycle_basis(G)
    cb_6 = 0
    cb_6_true = 0
    cb_else_true = 0
    potentially_merged = 0
    tot_nr_contigs_merged = 0
    strange_cases = 0

    for cycle in CB:
        str_case_abort = False
        contigs = [node[0] for node in cycle]
        d = {}

####### Very temporary implementation of dealing with haplotypes!!! #########            
        if len(cycle) >= 6:
            for i in contigs: d[i] = d.has_key(i)
            singles = [k for k in d.keys() if not d[k]]
            if len(singles) == 2:
                cb_else_true += 1
                #find length of the two paths between the source and sink
                haplotype_region = True
                #first path
                #get sink and source
                try:
                    cycle.index((singles[0], 'L'))
                    source_node = (singles[0], 'L')
                except ValueError:
                    source_node = (singles[0], 'R')
                try:
                    cycle.index((singles[1], 'L'))
                    sink_node = (singles[1], 'L')
                except ValueError:
                    sink_node = (singles[1], 'R')

                #check if region has been removed in some previuos step: this suggests some strange region
                for scaf in cycle:
                    if not scaf[0] in Scaffolds:
                        strange_cases += 1
                        #print 'STRANGE!'
                        str_case_abort = True
                if str_case_abort:
                    continue

                sub_G = nx.subgraph(G, cycle)
                path = nx.algorithms.shortest_path(sub_G, source=source_node, target=sink_node)


                ## Get length of path (OBS: gaps not implemented yet!!) ##
                #print 'First path' ,path
                length_path1 = 0
                nr_contigs_path1 = 0
                for scaffold_ in path:
                    scaffold = scaffold_[0]
                    if sink_node != scaffold_ and source_node != scaffold_:
                        length_path1 += Scaffolds[scaffold].s_length
                        nr_contigs = len(Scaffolds[scaffold].contigs)
                        #if nr_contigs > 1:
                        #    print 'More than one contig in scaffold, contigs are:'
                        for cont_obj in Scaffolds[scaffold].contigs:
                            nr_contigs_path1 += 1
                            #print 'Haplotype: ', cont_obj.is_haplotype, cont_obj.coverage
                            if not cont_obj.is_haplotype:
                        #        print 'Not haplotype'
                                haplotype_region = False



                #print 'Total length of path 1: ', length_path1/2.0

                set_of_nodes = set(path)
                start_end = set([source_node, sink_node])
                tot_set = set(cycle)
                nodes_to_remove = set_of_nodes.symmetric_difference(start_end)
                remaining_path = tot_set.symmetric_difference(nodes_to_remove)
                #print 'Secont path', remaining_path
                length_path2 = 0
                nr_contigs_path2 = 0
                for scaffold_ in remaining_path:
                    scaffold = scaffold_[0]
                    if sink_node != scaffold_ and source_node != scaffold_:
                        length_path2 += Scaffolds[scaffold].s_length
                        nr_contigs = len(Scaffolds[scaffold].contigs)
                        #if nr_contigs > 1:
                        #    print 'More than one contig in scaffold'
                        for cont_obj in Scaffolds[scaffold].contigs:
                            nr_contigs_path2 += 1
                            #print 'Haplotype: ', cont_obj.is_haplotype, cont_obj.coverage
                            if not cont_obj.is_haplotype:
                                haplotype_region = False

                try:
                    if length_path2 / float(length_path1) < param.hapl_ratio or length_path1 / float(length_path2) > 1 / param.hapl_ratio and haplotype_region:
                        potentially_merged += 1
                        tot_nr_contigs_merged += nr_contigs_path2 / 2.0
                        #Remove all contigs from path 2
                        to_remove = remaining_path.symmetric_difference(start_end)
                        G.remove_nodes_from(to_remove)
                        G_prime.remove_nodes_from(to_remove)
                        for node in to_remove:
                            try: #remove scaffold with all contigs
                                for contig_obj in Scaffolds[node[0]].contigs:
                                    del Contigs[contig_obj.name]
                                del Scaffolds[node[0]]
                            except KeyError: #scaffold and all contigs has already been removed since to_revove-path contains two nodes for each scaffold
                                pass
                except ZeroDivisionError:
                    pass

    print >> Information, 'NR of other interesting cycles: ', cb_else_true
    print >> Information, 'Potential hapl regions treated: ', potentially_merged
    print >> Information, 'Potential hapl contigs "removed": ', tot_nr_contigs_merged
    print >> Information, 'Nr of strange cases (contigs occurring in multiple regions): ', strange_cases
    return()

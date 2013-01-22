'''
Created on Sep 23, 2011

@author: ksahlin
'''
import Contig,Scaffold,Parameter

import GenerateOutput as GO
import networkx as nx
from networkx import algorithms
import sys
import GapCalculator as GC
from Norm import normpdf,normcdf
import ExtendLargeScaffolds as ELS


    

def Algorithm(G,G_prime,Contigs,small_contigs,Scaffolds,small_scaffolds,Information,C_dict,param):
    #search for linear streches in the graph, remove all cliques >2 and all contigs having more than three neighbors
    nr_edges=0
    for edge in G.edges_iter():
        if G[edge[0]][edge[1]]['nr_links']:
            nr_edges+=1

    print >> Information, str(nr_edges)+' link edges created.'
    print 'Perform inference on scaffold graph...'          
    #VizualizeGraph(G_prime,param)

    #TODO: should this come after removal of spurious edges?
    if param.detect_haplotype: 
        HaplotypicRegions(G,G_prime,Contigs,Scaffolds,param,C_dict)
    #save graph in dot format to file here
            
    ##If sigma specified. Pre calculate a look up table for every possible gap estimate in the common case
    ##where we have two long contigs that are linked. We do this one time per library to save computation time.
    if param.std_dev_ins_size:
        dValuesTable=GC.PreCalcMLvaluesOfdLongContigs(param.mean_ins_size,param.std_dev_ins_size,param.read_len)
    else :
        dValuesTable=None
    already_visited = set()
    if param.extend_paths:
        for node in G:
            already_visited.add(node)
        
        
    ##### Here is the scaffolding algorithm #######

    G = RemoveIsolatedContigs(G,Information)     #step1
    G = RemoveSpuriousEdges(G,G_prime,Scaffolds,Contigs,Information,param)    #step2 Remove edges that has very few edges compared to what the coverage and gap estimate says (should I also remove those with too many compared to what GapEst says? -> implement later)
    #G,Contigs,Scaffolds = DetectHaplotypes(G,Contigs,Scaffolds,param,C_dict) #remove haplotypes from scaffolds (to get more linear regions), one haplotype is left in graph and one in outputted to a file haplotypes.fa (with reference to which contig it is haplotypic to) 
    G = RemoveAmbiguousRegions(G,G_prime,Information,param) #step2
    G = RemoveIsolatedContigs(G,Information) #there are probably new isolated nodes created from step 2
    G,Contigs,Scaffolds = RemoveLoops(G,G_prime,Scaffolds,Contigs,Information,param)    #step4
    
    #The contigs that made it to proper scaffolds
    (Contigs,Scaffolds,param)=NewContigsScaffolds(G,G_prime,Contigs,small_contigs,Scaffolds,small_scaffolds,Information,C_dict,dValuesTable,param,already_visited)   #step5

    ##Here PathExtension algorithm between created scaffolds is called if PRO is activated
    
    if param.extend_paths:
        print '\n\n\n Searching for paths BETWEEN scaffolds\n\n\n'
        #TODO: Whenever we have removed an edge in G in algm above (deduced it as spurious), we should remove the same edge in G_prime. 
        #Also: we must update G_prime and G with the new scaffold objects created in this step. The isolated nodes removed here should not
        #be removed from G_prime
        PROBetweenScaf(G_prime,Contigs,small_contigs,Scaffolds,small_scaffolds,param,dValuesTable)
#        #TODO: Error here because we need to update G_prime and small_scaffolds with new scaffolds created in G
#        # in NewContigsScaffolds
        #TODO: G is reduced a lot when we reach here.. We must have all initial contigs in G.
        #already_visited contains all these nodes that need to be included here
        # Use it to "fill up" G again. Or do it in some other way??
        # remove isolated nodes and remove cycles step.
#        all_paths_sorted_wrt_score = ELS.BetweenScaffolds(G, G_prime,small_scaffolds)
#        for path in all_paths_sorted_wrt_score:
#            for node_ in path:
#                print node_
#                try:
#                    cont = Scaffolds[node_[0]].contigs[0].name
#                except KeyError:
#                    cont = small_scaffolds[node_[0]].contigs[0].name
#                    
#                print 'Contig:', cont
    
    
    ####### End of algorithm #####################

    
    return()

def VizualizeGraph(G,param):
    import os
    try: 
        import matplotlib
        matplotlib.use('Agg')
        
        try:
            os.mkdir(param.output_directory +'/graph_regions'+str(int(param.mean_ins_size)))
        except OSError:
            #directory is already created
            pass
        counter = 1
        #CB = nx.connected_component_subgraphs(G)
        nx.draw(G)
        matplotlib.pyplot.savefig(param.output_directory +'graph_regions'+str(int(param.mean_ins_size))+'/'+str(counter)+'.png')
#        CB = nx.cycle_basis(G)
#        print 'NR of SubG: ' ,len(CB)
#
#        for cycle in CB:
#            if len(cycle) >= 6: # at leats 6 nodes if proper cycle (haplotypic region)
#                subgraph = nx.Graph()
#                subgraph.add_cycle(cycle)
#                nx.draw(subgraph)
#                matplotlib.pyplot.savefig(param.output_directory +'graph_regions'+str(int(param.mean_ins_size))+'/'+str(counter)+'.png')
#                matplotlib.pyplot.clf()
#                counter += 1
    except ImportError:
        pass
    return()
    
def RemoveIsolatedContigs(G,Information):
    print 'Remove isolated nodes.'
    counter=0
    for node in G.nodes():
        if node in G:
            nbr=G.neighbors(node)[0]
            if len(G.neighbors(node)) == 1 and len(G.neighbors(nbr)) == 1:
                counter+=1
                G.remove_nodes_from([node,nbr])
    print >> Information, str(counter)+' isolated contigs removed from graph.'
    return(G)

def RemoveSpuriousEdges(G,G_prime,Scaffolds,Contigs,Information,param):
    ## ToDo: IMPROVE THE MODEL THIS FUNCTION USES ##
    # step 1 : infer gap distance for an edge
    # step 2 : Calculate what should be the expected nr of links between the contigs/super contigs given the coverage and the estimated distance (done a very naively at the moment, not considering haplotypes)
    # step 3 : Remove edge if the total number of links is less than E[nr of links]/3 (maybe also remove edges with more than 3*E[nr of links] links?)
    
    #remove links with less than user set edge_support
#    for node in G:
#        nbrs=G.neighbors(node)
#        for nbr in nbrs:
#            if G[node][nbr]['nr_links'] and G[node][nbr]['nr_links'] < param.edgesupport:
#                G.remove_edge(node,nbr)

    #n=0.0
    #x_links=0
    #x_links2 =0
    #true_links=0
    
    ####### Remove edges that has way too few links given what the gap observation is #######
    fishy_edges=[]
    for edge in G.edges_iter():
        if G[edge[0]][edge[1]]['nr_links']:
            len1 = Scaffolds[ edge[0][0] ].s_length 
            len2 = Scaffolds[ edge[1][0] ].s_length
            contig_coverages1 = [ Contigs[cont_obj.name].coverage for cont_obj in Scaffolds[ edge[0][0] ].contigs ]
            mean_cov1 = sum(contig_coverages1)/float(len(contig_coverages1))
            
            contig_coverages2 = [ Contigs[cont_obj.name].coverage for cont_obj in Scaffolds[ edge[1][0] ].contigs ]
            mean_cov2 = sum(contig_coverages2)/float(len(contig_coverages2))
            #CALCULATE GAP
            ## have precalc table here as well for time issues!!
            
            #here I merge the individual gap observations o_1 and o_2 into one o.
            #obs_list = map(lambda (x, y): x + y, G[edge[0]][edge[1]]['gap_dist'])
            #G[edge[0]][edge[1]]['gap_dist'] = obs_list
            tot_obs=G[edge[0]][edge[1]]['gap_dist']
            nr_links=G[edge[0]][edge[1]]['nr_links']
            data_observation=(nr_links*param.mean_ins_size -tot_obs)/float(nr_links)

            if nr_links >= 5:  
                gap=GC.GapEstimator(param.mean_ins_size,param.std_dev_ins_size,param.read_len,data_observation,len1,len2)  
                #FIND THE EXPECTED NR OF LINKS BASED ON THE MEAN COVERAGE AND THE INFERRED DISTANCE 
                #n +=1
                #k is the average distance between two successive reads on the genome for the current library
                try:
                    k = 2*param.read_len/float(min(mean_cov1,mean_cov2))
                except ZeroDivisionError:
                    print "Contig with link had 0 coverage, this should not happen and is due to BWAs buggy BAM file, removing this link obviously.."
                    fishy_edges.append((edge[0], edge[1]))
                    #print 'scaffold 1: ', contig_coverages1
                    #print 'scaffold 2: ', contig_coverages2
                    continue
#                    print 'len1', len1
#                    print 'len2', len2
#                    print 'Contigs in first scaffold: '
#                    for cont_obj in Scaffolds[ edge[0][0] ].contigs:
#                        print cont_obj.name
#                    print 'Contigs in second scaffold: '
#                    for cont_obj in Scaffolds[ edge[1][0] ].contigs:
#                        print cont_obj.name
#                    sys.exit(0)
                #E_links2 = (param.mean_ins_size-gap-2*param.read_len)/(k)
                E_links = ExpectedLinks(len1,len2,gap,k,param)
                #print 'NEW WAY: ', E_links, 'True: ', nr_links,'gap: ',gap                
                #x_links += E_links                    
                #x_links2 += E_links2
                #true_links += nr_links
                
                # the threshold_decider varies with the estimated gap, it is more loose for larger gaps 
                # since for large gaps, both the gap estimate and the E_ links is more uncertain. 
                try: 
                    relation = max((param.mean_ins_size - max(0,gap)) , 0) # if gap is larger than mean set to 0 , next line then becomes division error
                    threshold_decider = 8/( relation / float(param.mean_ins_size) )
                except ZeroDivisionError:
                    #gap is the same as mean ins size so set threshold to large nr (no discarding of edge)
                    threshold_decider = 2**32                
                if nr_links < E_links/threshold_decider:
                    fishy_edges.append((edge[0], edge[1]))
                    print 'FISHY EDGE: Nr_links: ', nr_links, ' E_links: ', E_links,' Gap: ', gap, 'threshold decider: ',threshold_decider
#                    print 'scaf1: '
#                    for cont in Scaffolds[edge[0][0]].contigs:
#                        print cont.name,cont.coverage
#                    print 'scaf2: '
#                    for cont in Scaffolds[edge[1][0]].contigs:
#                        print cont.name,cont.coverage                   
    #print 'Naive mean: ', x_links2/n, 'Pro mean: ', x_links/n,'true: ', true_links/n
    #print 'Average position to next read-pair based on coverage:', k
    counter=0
    for edge in fishy_edges:
        G.remove_edge(edge[0], edge[1])
        counter +=1   
        if param.extend_paths:
            try:
                G_prime.remove_edge(edge[0], edge[1])
            except nx.exception.NetworkXError: #this edge might have been removed in the earlier filtering step of G_prime (the very heuristic < 5 edges)
                pass  
    print 'Removed ', counter, ' spurious edges.'
    
#        for contig in Scaffolds[edge[0][0]].contigs:
#            print 'Scaf 1: ', contig.name
#        for contig in Scaffolds[edge[1][0]].contigs:
#            print 'Scaf 2: ', contig.name
#        print 'removed: ', (edge[0], edge[1])     
    
    return(G)

# SLOOOW expected links
#def ExpectedLinks(len1,len2,gap,k,param):
#    from scipy.stats import norm
#    mean = param.mean_ins_size
#    std_dev = param.std_dev_ins_size 
#    read_len = param.read_len  
#    nr_iterations = int( min((min(len1,len2) - param.read_len)/k,mean+3*std_dev) )
#    E_links=0
#    for i in range(0,nr_iterations):
#        E_links += norm.cdf( (len1 + len2 + gap + -i*k-mean)/float(std_dev) ) - norm.cdf( (min(len1, len2) + gap + read_len -i*k-mean)/float(std_dev) )
#        
#    return(E_links)

def ExpectedLinks(len1,len2,gap_est,k,param):
    #Working good for contigs at the moment (not working good for scaffolds with lots of N's) 
    std_dev = float(param.std_dev_ins_size)
    #We should not expect more links for a negative gap since the negative gap is most likely:
    # 1) a bad gap estimate 2) a strange haplotypic region 
    # Even if the overlap is correct we should not expect too many links in the overlapping region since part of
    # the PE/MP will map to the same contig so this will create a problem anyway, therefor we count with gap = max(gap_est,0)
    gap = max(gap_est,0) 
    #Specifying input arguments
    b1 = (len1 + len2 + gap - param.mean_ins_size)/std_dev
    a1 = ( max(len1,len2) + gap - param.read_len - param.mean_ins_size  ) / std_dev
    b2 = ( min(len1,len2) + gap + param.read_len - param.mean_ins_size ) / std_dev
    a2 = ( gap  + param.read_len - param.read_len - param.mean_ins_size  ) / std_dev
    def Part(a,b): 
        expr1 = (min(len1,len2)-param.read_len)/k * normcdf(a, 0, 1)
        expr2 = (b*std_dev)/k*( normcdf(b, 0, 1) - normcdf(a, 0, 1) ) 
        expr3 = (std_dev/k)*(normpdf(b, 0, 1) - normpdf(a, 0, 1))
        value = expr1 + expr2 + expr3
        return value 
    E_links = Part(a1,b1) - Part(a2,b2)    
    return(E_links)


def DetectHaplotypes(G,Contigs,Scaffolds,param,C_dict):
    ##FIX THIS FCN
    ## We detect haplotypes only if it is two single contigs that is suspicious (that is, if we detect a fishy region between two scaffolds, we are not trying to infer that but leave it as it is)
    Haplotypes={}
    
    #VERY NASTY IMPLEMENTATION FOR FINDING SIMPLE BUBBLES!!
    for node in G:
        nbrs=G.neighbors(node)
        common_contig=[]
        hapl_cont=[]
        if len(nbrs) == 3: #exactly two neighbors:
            for nbr in nbrs: #for both nbrs
                if G[node][nbr]['nr_links'] and len(Scaffolds[ nbr[0] ].contigs) == 1: #they should both be single contigs
                    if Scaffolds[ nbr[0] ].contigs[0].is_haplotype:
                        nbrs2=G.neighbors(nbr)
                        if len(nbrs2) == 2: #exactly one neighbor
                            for nbr2 in nbrs2: #for both nbrs
                                if not G[nbr][nbr2]['nr_links']: #the other node of the contig
                                    nbrs3=G.neighbors(nbr2)
                                    if len(nbrs3) == 2: #exactly one neighbor
                                        for nbr3 in nbrs3: #for both nbrs
                                            if G[nbr2][nbr3]['nr_links']:
                                                common_contig.append(nbr3)
                                                hapl_cont.append([Scaffolds[ nbr[0] ].contigs[0].name, nbr]) #[contig name, node label in G]
                                                #print [Scaffolds[ nbr[0] ].contigs[0].name, nbr]
        if len(common_contig) == 2 and common_contig[0] == common_contig[1]:
            if hapl_cont[0][0] not in Haplotypes and hapl_cont[1][0] not in Haplotypes:
                Haplotypes[hapl_cont[0][0]] = (hapl_cont[1][0],hapl_cont[1][1])
                #print (hapl_cont[1][0],hapl_cont[1][1])
           
    for key in Haplotypes:
        print 'Haplotype region: ', key,Haplotypes[key][0] ,
        if ( Contigs[key].length / float( Contigs[Haplotypes[key][0]].length ) ) < 1.3 and ( Contigs[key].length / float( Contigs[Haplotypes[key][0]].length )) >0.77 : #about the same length 
            print 'OK', Contigs[key].length, Contigs[Haplotypes[key][0]].length
        else: 
            print 'More than 20% length difference' ,Contigs[key].length, Contigs[Haplotypes[key][0]].length
        G.remove_node( Haplotypes[key][1] )
        if Haplotypes[key][1][1] == 'L': 
            G.remove_node( (Haplotypes[key][1][0],'R') )
        else:
            G.remove_node( (Haplotypes[key][1][0],'L') )
        del Scaffolds[Haplotypes[key][1][0]]
    nr_haplotypes = len(Haplotypes)
    print 'Resolved ', nr_haplotypes, ' simple haplotypic regions (simple bubble in graph)' 
    Contigs = GO.PrintOutHaplotypes(Haplotypes,Contigs,param.output_directory,C_dict)

                
        #    for contig in Scaffolds[node[0]].contigs:
        #        print 'REGION: ',contig.name, ' nbrs: ', 
        #        for nbr in nbrs:
        #            for contig_nbr in Scaffolds[nbr[0]].contigs:
        #                print contig_nbr.name
                
        #for contig in Scaffolds[node[0]].contigs:
        #    pass
            #print Contigs[contig.name].coverage, Contigs[contig.name].name
    return(G,Contigs,Scaffolds)

def HaplotypicRegions(G,G_prime,Contigs,Scaffolds,param,C_dict):
    CB = nx.cycle_basis(G)
    print 'NR of cycles bases: ' ,len(CB)
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
                haplotype_region =True
                #first path
                #get sink and source
                try:
                    cycle.index((singles[0],'L'))
                    source_node = (singles[0],'L')
                except ValueError:
                    source_node = (singles[0],'R')
                try:
                    cycle.index((singles[1],'L'))
                    sink_node = (singles[1],'L')
                except ValueError:
                    sink_node = (singles[1],'R')                   
            
                #check if region has been removed in some previuos step: this suggests some strange region
                for scaf in cycle:
                    if not scaf[0] in Scaffolds:
                        strange_cases +=1
                        #print 'STRANGE!'
                        str_case_abort = True
                if str_case_abort:
                    continue
                
                sub_G = nx.subgraph(G, cycle)                
                path = nx.algorithms.shortest_path(sub_G, source = source_node, target = sink_node)
                
                        
                ## Get length of path (OBS: gaps not implemented yet!!) ##
                #print 'First path' ,path
                length_path1 = 0
                nr_contigs_path1 =0
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
                                haplotype_region =False
                                
                                

                #print 'Total length of path 1: ', length_path1/2.0
                        
                set_of_nodes =set(path)
                start_end= set([source_node,sink_node])
                tot_set =set(cycle)
                nodes_to_remove = set_of_nodes.symmetric_difference(start_end)
                remaining_path = tot_set.symmetric_difference(nodes_to_remove)
                #print 'Secont path', remaining_path
                length_path2 = 0
                nr_contigs_path2 =0
                for scaffold_ in remaining_path:
                    scaffold = scaffold_[0]
                    if sink_node != scaffold_ and source_node != scaffold_:
                        length_path2 += Scaffolds[scaffold].s_length
                        nr_contigs = len(Scaffolds[scaffold].contigs)
                        #if nr_contigs > 1:
                        #    print 'More than one contig in scaffold'
                        for cont_obj in Scaffolds[scaffold].contigs:
                            nr_contigs_path2 +=1
                            #print 'Haplotype: ', cont_obj.is_haplotype, cont_obj.coverage
                            if not cont_obj.is_haplotype:
                        #        print 'Not haplotype'
                                haplotype_region =False
                #print 'Total length of path 2: ', length_path2/2.0 
                #print 'Cycle involving ', (len(contigs)-2)/2.0, 'contigs.'
                #print 'Nr contigs path one: ', nr_contigs_path1/2.0
                #print 'Nr contigs path two: ', nr_contigs_path2/2.0
                #print 'Length 1st path: ',  length_path1/2.0
                #print 'Length 2nd path: ',  length_path2/2.0
                #print 'All haplotypes in region: ', haplotype_region 
                    
                try: 
                    if length_path2/float(length_path1) < param.hapl_ratio or length_path1/float(length_path2) > 1/param.hapl_ratio and haplotype_region:
                        potentially_merged +=1
                        tot_nr_contigs_merged += nr_contigs_path2/2.0  
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


                                
                              
    #print 'NR of 6 cycles: ', cb_6, 'NR of interesting 6 cycles: ', cb_6_true
    print 'NR of other interesting cycles: ', cb_else_true
    print 'Potential hapl regions treated: ', potentially_merged
    print 'Potential hapl contigs "removed": ',tot_nr_contigs_merged
    print 'Nr of strange cases (contigs occurring in multiple regions): ', strange_cases
    return()                


def RemoveAmbiguousRegions(G,G_prime,Information,param):    
### remove edges from node if more than two edges 
### Keep the edge with max support e* if all the other edges of that node satisfies links(e*)>l*links(e_i) (l input threshold) all i
### and links(e*)>edge_support (input value)
    edge_support = param.edgesupport
    ratio = param.rel_weight
    print 'Remove edges from node if more than two edges'
    counter1=0
    counter2=0
    counter3 = 0
    for node in G:
        nbrs=G.neighbors(node)
        #Remove ambiguous edges
        if len(nbrs) > 2:
            nr_link_list=[]
            for nbr in nbrs:
                if G[node][nbr]['nr_links']:
                    nr_link_list.append((G[node][nbr]['nr_links'],nbr)) 
            nr_link_list.sort()
            
        ###### this if statement is used when scaffolding with fosmidpools, i.e. -p <int> is specified ###
            if param.fosmidpool != None and nr_link_list[-2][0] > param.fosmidpool:
            ### remove all edges on this side of the contig
                nr_nbrs=len(nr_link_list)
                for i in xrange(0,nr_nbrs):
                    G.remove_edge(node,nr_link_list[i][1]) 
                    if param.extend_paths:
                        try: #we might have been removed this edge from G_prime when we did individual filtering of G_prime in CreateGraph module
                            G_prime.remove_edge(node,nr_link_list[i][1])
                        except nx.exception.NetworkXError:
                            pass 
                counter3+=1
                continue
            ##################################################################
            
            if nr_link_list[-1][0] >= (ratio*nr_link_list[-2][0]) and nr_link_list[-1][0] >= edge_support:
            ### save the dominating link edge on this side of the contig
                #print 'Linklist: ', nr_link_list
                nr_nbrs=len(nr_link_list)
                for i in xrange(0,nr_nbrs-1):
                    G.remove_edge(node,nr_link_list[i][1]) 
                    if param.extend_paths:
                        try: #we might have been removed this edge from G_prime when we did individual filtering of G_prime in CreateGraph module
                            G_prime.remove_edge(node,nr_link_list[i][1])
                        except nx.exception.NetworkXError:
                            pass
                counter1+=1      
            else:
            ### remove all edges on this side of the contig
                nr_nbrs=len(nr_link_list)
                for i in xrange(0,nr_nbrs):
                    G.remove_edge(node,nr_link_list[i][1]) 
                    if param.extend_paths:
                        try: #we might have been removed this edge from G_prime when we did individual filtering of G_prime in CreateGraph module
                            G_prime.remove_edge(node,nr_link_list[i][1])
                        except nx.exception.NetworkXError:
                            pass
                counter2+=1
                
        #Remove low support edges e < k in remaining linear regions             
        else:
            for nbr in nbrs:
                if G[node][nbr]['nr_links'] and G[node][nbr]['nr_links'] < edge_support:
                    G.remove_edge(node,nbr)
                    if param.extend_paths:
                        try: #we might have been removed this edge from G_prime when we did individual filtering of G_prime in CreateGraph module
                            G_prime.remove_edge(node,nbr)
                        except nx.exception.NetworkXError:
                            pass
                            
                
    print >> Information,str(counter1+counter2)+ ' ambiguous regions in graph ( a contig with more than 2 neighbors).'
    print >> Information,str(counter1)+ ' of these regions solved (one dominating edge)'
    print >> Information,str(counter2)+ ' of these regions unsolved (removed all edges)'
    if param.fosmidpool != None:
        print str(counter3), 'of these regions were unsolved due to fosmid pool specific settings (removed all edges)'
    return(G)
 


    

def RemoveLoops(G,G_prime,Scaffolds,Contigs,Information,param):    
#### After the proceure above, we hope that the graph is almost perfectly linear but we can still be encountering cycles (because of repeats or haplotypic contigs that has slipped through our conditions). Thus we finally search for loops
    print 'Contigs/scaffolds left:', len(G.nodes())/2
    print 'Remove remaining cycles...'
    graphs=nx.connected_component_subgraphs(G)
    #print 'Nr connected components',len(graphs)
    counter=0
    for graph in graphs:
        list_of_cycles=algorithms.cycles.cycle_basis(graph)
        for cycle in list_of_cycles:
            print >> Information, 'A cycle in the scaffold graph: '+str(cycle)+'\n'
            print 'A cycle in the scaffold graph: '+str(cycle), graph.edges()
            counter+=1
            for node in cycle:
                if node in G:
                    #we split up the whole cycle into separate contigs and send them to F
                    scaffold_=node[0]
                    G.remove_nodes_from([(scaffold_,'L'),(scaffold_,'R')])
                    if param.extend_paths:
                        G_prime.remove_nodes_from([(scaffold_,'L'),(scaffold_,'R')])
#                    S_obj=Scaffolds[scaffold_]
#                    list_of_contigs=S_obj.contigs   #list of contig objects contained in scaffold object
#                    Contigs, F = GO.WriteToF(F,Contigs,list_of_contigs)
#                    del Scaffolds[scaffold_]
    print >> Information, str(counter)+' cycles removed from graph.'
    return(G,Contigs, Scaffolds)

def NewContigsScaffolds(G,G_prime,Contigs,small_contigs,Scaffolds,small_scaffolds,Information,C_dict,dValuesTable,param,already_visited):    
### Remaining scaffolds are true sensible scaffolds, we must now update both the library of scaffold objects and the library of contig objects
    new_scaffolds_=nx.connected_component_subgraphs(G)
    print 'Nr of new scaffolds created: '+str(len(new_scaffolds_))
    print>>Information, 'Nr of new scaffolds created in this step: '+str(len(new_scaffolds_))
    for new_scaffold_ in new_scaffolds_:
        param.scaffold_indexer+=1
        #scaf_size=len(new_scaffold_)
        scaffold_length=0
        contig_list=[]
        
        #Store nr_of links between contigs before "destroying" the graph
        for edge in new_scaffold_.edges_iter():
            nr_links=G[edge[0]][edge[1]]['nr_links']
            side1=edge[0][1]
            side2=edge[1][1]
            if nr_links:
                contig_objects1=Scaffolds[edge[0][0]].contigs
                contig_objects2=Scaffolds[edge[1][0]].contigs
                GiveLinkConnection(Contigs,contig_objects1,contig_objects2,side1,side2,nr_links) 
                             
        ##### Here PathExtension algorithm is called if PRO is activated #####
        if param.extend_paths:
            PROWithinScaf(G,G_prime,Contigs,small_contigs,Scaffolds,small_scaffolds,param,new_scaffold_,dValuesTable,already_visited)
                                    
                                            
        for node in new_scaffold_:   
            if len(G.neighbors(node))==1:
                start = node
                break
        for node in new_scaffold_:   
            if len(G.neighbors(node)) == 1 and node != start:
                end = node        
        #Create info to new scaffold object such as total length and the contig objects included
        
        prev_node=('','')
        pos=0
        (G,contig_list,scaffold_length)=UpdateInfo(G,Contigs,small_contigs,Scaffolds,small_scaffolds,start,prev_node,pos,contig_list,scaffold_length,dValuesTable,param)
        S=Scaffold.scaffold(param.scaffold_indexer,contig_list,scaffold_length,{},{})  #Create the new scaffold object 

        Scaffolds[S.name]=S        #include in scaffold library

        if param.extend_paths:
            # Find the ends of the old subgraph new_scaffold_. We want them to be able to relabel these end nodes as the new sides on the new scaffold object created
            #only these ends are allowed to have links because they are of size mean+ 4*sigma so nothing is supposed to span over.
            
            #add the new scaffold object to G_prime

            G_prime.add_node((S.name,'L'))  #start node
            G_prime.add_node((S.name,'R'))  # end node
            G_prime.add_edge((S.name,'L'),(S.name,'R'),nr_links = None )
            for nbr in G_prime.neighbors(start):
                nr_links_ = G_prime[start][nbr]['nr_links']
                if nr_links_:
                    gap_dist_ = G_prime[start][nbr]['gap_dist']
                    G_prime.add_edge((S.name,'L'),nbr,nr_links = nr_links_,gap_dist = gap_dist_ )
                    
            for nbr in G_prime.neighbors(end):
                nr_links_ = G_prime[end][nbr]['nr_links']
                if nr_links_:
                    gap_dist_ = G_prime[end][nbr]['gap_dist']
                    G_prime.add_edge((S.name,'R'),nbr,nr_links = nr_links_,gap_dist = gap_dist_ )                    
                
            #remove the old scaffold objects from G_prime
            G_prime.remove_nodes_from(new_scaffold_)      

    return(Contigs,Scaffolds,param)


def UpdateInfo(G,Contigs,small_contigs,Scaffolds,small_scaffolds,node,prev_node,pos,contig_list,scaffold_length,dValuesTable,param):
    scaf=node[0]
    side=node[1]
    prev_scaf=prev_node[0]
    if len(G.neighbors((scaf,side))) ==0:  #reached end of scaffol
        #find the contig with the largest position
        object_with_largest_pos_in_scaffold=max(contig_list, key=lambda object: object.position+object.length)
        scaffold_length=object_with_largest_pos_in_scaffold.position+object_with_largest_pos_in_scaffold.length
        del Scaffolds[scaf] #finally, delete the old scaffold object
        return(G,contig_list,scaffold_length)
    else:
        nbr_node=G.neighbors((scaf,side)) 
        nbr_scaf = nbr_node[0][0]
        nbr_side = nbr_node[0][1]
        if scaf != prev_scaf:
            if side =='L':    #Contig/scaffold still has same orientation as in previous iteration, just update position in scaffold                                           
                #want to assign nr of links to contig object, note that in case of a "multiple contigs"-scaffold object, only the outermost contig holds the information of the total nr of links between the two scaffold objects
                contig_objects=Scaffolds[scaf].contigs #list of contig objects
                
                #Update just update position in scaffold 
                for contig in contig_objects:
                    contig.scaffold =param.scaffold_indexer
                    contig.position +=pos
                    #direction unchanged
                    contig_list.append(contig)
                G.remove_node((scaf,side))   
                prev_node=node
                node=(nbr_scaf,nbr_side)                   
                pos+=Scaffolds[scaf].s_length  #update position before sending it to next scaffold
                G,contig_list,scaffold_length=UpdateInfo(G,Contigs,small_contigs,Scaffolds,small_scaffolds,node,prev_node,pos,contig_list,scaffold_length,dValuesTable,param)               
          
            else:  #Contig/scaffold need to change orientation as well as modify orientation in this case
                contig_objects=Scaffolds[scaf].contigs #list of contig objects
                for contig in contig_objects:
                    contig.scaffold =param.scaffold_indexer
                    curr_scaf_length=Scaffolds[scaf].s_length
                    curr_pos_within_scaf=contig.position
                    contig_length = contig.length
                    contig.position = pos +  (curr_scaf_length - curr_pos_within_scaf) -contig_length #updates the position within scaf
                    contig.direction = bool(True - contig.direction) #changes the direction
                    contig_list.append(contig)

                G.remove_node((scaf,side))   
                prev_node=node
                node= (nbr_scaf,nbr_side)                    
                pos+=Scaffolds[scaf].s_length  #update position before sending it to next scaffold
                G,contig_list,scaffold_length=UpdateInfo(G,Contigs,small_contigs,Scaffolds,small_scaffolds,node,prev_node,pos,contig_list,scaffold_length,dValuesTable,param)  
        else: 
            #calculate gap to next scaffold
            sum_obs=G[(scaf,side)][(nbr_scaf,nbr_side)]['gap_dist']
            nr_links=G[(scaf,side)][(nbr_scaf,nbr_side)]['nr_links']
            data_observation=(nr_links*param.mean_ins_size -sum_obs)/float(nr_links)
            c1_len=Scaffolds[scaf].s_length
            c2_len=Scaffolds[nbr_scaf].s_length
            #do fancy gap estimation by the bias estimator formula
            if param.std_dev_ins_size and nr_links >= 5:  
                #pre calculated value in lookup table 
                if c1_len > param.mean_ins_size+ 4*param.std_dev_ins_size and c2_len > param.mean_ins_size+ 4*param.std_dev_ins_size:
                    #(heuristic scale down of table to gaps of at most 2 stddevs away from mean)
                    try:
                        avg_gap=dValuesTable[int(round(data_observation,0))]
                    except KeyError:
                        avg_gap=GC.GapEstimator(param.mean_ins_size,param.std_dev_ins_size,param.read_len,data_observation,c1_len,c2_len)
                        print 'Gap estimate was outside the boundary of the precalculated table, obs were: ',data_observation,'binary search gave: ', avg_gap
                        #print 'Gap estimate was outside the boundary of the precalculated table'
                        #print 'Boundaries were [' ,-2*param.std_dev_ins_size, ' , ',param.mean_ins_size+2*param.std_dev_ins_size, ' ]. Observation were: ', data_observation
                        #if data_observation < -2*param.std_dev_ins_size:                            
                        #    avg_gap=int(-2*param.std_dev_ins_size)
                        #    print 'Setting gap estimate to min value:', avg_gap,'number of links: ' , nr_links
                        #else:
                        #    avg_gap=int(param.mean_ins_size+2*param.std_dev_ins_size-2*param.read_len)
                        #    print 'Setting gap estimate to max value:', avg_gap,'number of links: ' , nr_links
                #Do binary search for ML estimate of gap
                else:
                    avg_gap=GC.GapEstimator(param.mean_ins_size,param.std_dev_ins_size,param.read_len,data_observation,c1_len,c2_len)
                    #print avg_gap
            #do naive gap estimation
            else: 
                avg_gap=int(data_observation)
                
            #print avg_gap, nr_links, 'OBS: ',data_observation #, obs_list
            
            ##Here PathExtension algorithm is called if PRO is activated
#            if param.extend_paths:
#                start = (scaf,side)
#                #print G_prime.edges()
#                end = (nbr_scaf,nbr_side)
#                #if Scaffolds[scaf].contigs[0].name == 'contig95,159692,0' or Scaffolds[nbr_scaf].contigs[0].name == 'contig95,159692,0':
#                    #print 'NEIGHBORS',G_prime.neighbors(((51,'L'))), G_prime[(51,'L')], G.neighbors(((51,'L'))), G[(51,'L')]
#                high_score_path,bad_links,score,path_len = ELS.WithinScaffolds(G, G_prime,small_scaffolds, start, end, already_visited,param.ins_size_threshold )
#                if len(high_score_path) > 1:
#                    print avg_gap,high_score_path, 'nr bad links path:',bad_links, 'Score:',score,'path length:', path_len
#                    if score >= 2.0: 
#                        for small_scaf in high_score_path:
                            
                        

#                    for node_ in high_score_path:
#                        try:
#                            cont = Scaffolds[node_[0]].contigs[0].name
#                        except KeyError:
#                            cont = small_scaffolds[node_[0]].contigs[0].name
#                            
#                        print 'Contig:', cont
#                if Scaffolds[scaf].contigs[0].name == 'contig95,159692,0' or Scaffolds[nbr_scaf].contigs[0].name == 'contig95,159692,0':
#                    print 'NEIGHBORS',G_prime.neighbors(((30,'L'))) #,G_prime[(30,'L')][(76,'R')]
#                    high_score_paths = ELS.WithinScaffolds(G, G_prime, start, end, already_visited )
#                    for path_ in high_score_path:
#                        print 'Score:',path_[0]
#                        for node_ in path_[1]:
#                            try:
#                                cont = Scaffolds[node_[0]].contigs[0].name
#                            except KeyError:
#                                cont = small_scaffolds[node_[0]].contigs[0].name
#                                
#                            print 'Contig:', cont                
                
            if avg_gap < 0:
                #print avg_gap,'number of links: ' , nr_links
                #TODO: Eventually implement SW algm to find ML overlap
                avg_gap=0 
            pos += int(avg_gap)
            G.remove_node((scaf,side))   
            prev_node=node
            node= (nbr_scaf,nbr_side) 
            del Scaffolds[scaf] #finally, delete the old scaffold object
            G,contig_list,scaffold_length=UpdateInfo(G,Contigs,small_contigs,Scaffolds,small_scaffolds,node,prev_node,pos,contig_list,scaffold_length,dValuesTable,param) 
    return(G,contig_list,scaffold_length)

def PROWithinScaf(G,G_prime,Contigs,small_contigs,Scaffolds,small_scaffolds,param,new_scaffold_,dValuesTable,already_visited):
    #loc_count = 0
    for edge in new_scaffold_.edges_iter():
        nr_links_=G[edge[0]][edge[1]]['nr_links']
        if nr_links_:
            start = edge[0]
            end = edge[1]
            c1_len = Scaffolds[start[0]].s_length
            c2_len = Scaffolds[end[0]].s_length
            sum_obs=G[start][end]['gap_dist']
            data_observation=(nr_links_*param.mean_ins_size -sum_obs)/float(nr_links_)
            avg_gap=GC.GapEstimator(param.mean_ins_size,param.std_dev_ins_size,param.read_len,data_observation,c1_len,c2_len)
#            try:
#                avg_gap=dValuesTable[int(round(data_observation,0))]
#            except KeyError:
#                print 'Gap estimate was outside the boundary of the precalculated table'
#                print 'Boundaries were [' ,-2*param.std_dev_ins_size, ' , ',param.mean_ins_size+2*param.std_dev_ins_size, ' ]. Observation were: ', data_observation
#                if data_observation < -2*param.std_dev_ins_size:                            
#                    avg_gap=int(-2*param.std_dev_ins_size)
#                    print 'Setting gap estimate to min value:', avg_gap,'number of links: ' , nr_links_
#                else:
#                    avg_gap=int(param.mean_ins_size+2*param.std_dev_ins_size-2*param.read_len)
#                    print 'Setting gap estimate to max value:', avg_gap,'number of links: ' , nr_links_
                    
            high_score_path,bad_links,score,path_len = ELS.WithinScaffolds(G, G_prime,small_scaffolds, start, end, already_visited,param.ins_size_threshold )
            if len(high_score_path) > 1:
                #print 'Start scaf:',start, 'End scaf:', end, 'Avg gap between big:', avg_gap,'Path:',high_score_path, 'nr bad links path:',bad_links, 'Score:',score,'path length:', path_len
                if score >= 0.0: 
                    #loc_count += 1
                    #remove edge in G to fill in the small scaffolds
                    G.remove_edge(start,end)
                    #add small scaffolds to G
                    for i in range(0,len(high_score_path)-1):
                        nr_lin= G_prime[high_score_path[i]][high_score_path[i+1]]['nr_links']
                        try:
                            total_dist= G_prime[high_score_path[i]][high_score_path[i+1]]['gap_dist']
                            G.add_edge(high_score_path[i],high_score_path[i+1],nr_links=nr_lin,gap_dist=total_dist)
                        except KeyError:
                            G.add_edge(high_score_path[i],high_score_path[i+1],nr_links=nr_lin)
                    #remove the small contigs from G_prime
                    G_prime.remove_nodes_from(high_score_path[1:-1])
                    # move all contig and scaffold objects from "small" structure to large structure to fit with UpdateInfo structure
                    small_scafs = map(lambda i: high_score_path[i],filter(lambda i: i%2 == 1,range(len(high_score_path)-1)))
                    for item in small_scafs:
                        scaf_obj = small_scaffolds[item[0]]
                        Scaffolds[item[0]] = scaf_obj
                        cont_objects = scaf_obj.contigs
                        for obj_ in cont_objects:
                            ctg_name = obj_.name
                            Contigs[ctg_name] = obj_
                            del small_contigs[ctg_name]
                        del small_scaffolds[item[0]]  
                        
    ####################################################################
    return()

def PROBetweenScaf(G_prime,Contigs,small_contigs,Scaffolds,small_scaffolds,param,dValuesTable):
    # Must fill up G with all the new scaffold objects here, they will act as starts and ends for paths in G_prime
#    for edge in G_prime.edges():
#        if edge[0][0] == 118 or edge[1][0] == 118:
#            print edge
    G = nx.Graph()
    for node in G_prime:
        if node[0] in Scaffolds: #meets the length criteria
            G.add_node(node)
            
    #Find paths between scaffolds here
    print 'Entering ELS.BetweenScaffolds'
    all_paths_sorted_wrt_score = ELS.BetweenScaffolds(G, G_prime,small_scaffolds)
    #TODO: If we get too many paths back and run into memory issues we could change so that only paths with score over 2.0 are stored in ELS module
    #print 'Total number of paths between scaffolds detected:', len(all_paths_sorted_wrt_score)
    for sublist in reversed(all_paths_sorted_wrt_score):
        path = sublist[2]
        bad_links = sublist[1]
        score = sublist[0]
        path_len = sublist[3]
        #if score < 0:
        #    break
        print 'Path:',path, 'nr bad links path:',bad_links, 'Score:',score,'path length:', path_len
        
        ## Need something here that keeps track on which contigs that are added to Scaffolds so that a
        ## contig is only present once in each path
        ## !!!! this does not work!! need to keep tack of paths added since scaffolds change name below
        for scaf_ in path[1:-1]:
            if scaf_[0] in Scaffolds:
                continue
        
        
        ## Here we do the "joining of two scaffolds with the new path if no contig/scaffold is present
        ## in another path, we need to update "Scaffolds" structure here along as we go in order for
        ## the above dublette checking function to work
        
        #make the path a small linear graph
        G_ = nx.Graph()
        path.insert(0,(path[0][0],'R')) if path[0][1] == 'L' else path.insert(0,(path[0][0],'L'))
        path.insert(len(path),(path[0][0],'R')) if path[len(path)-1][1] == 'L' else path.insert(len(path),(path[len(path)-1][0],'L'))
        #print '\n', path,'\n'
        G_.add_edges_from(zip(path[::1], path[1::]))
        for edge in G_.edges():
            #print edge[0],edge[1]
            #print G_prime[edge[0]][edge[1]] #,(G_prime.node[edge[0]],G_prime.node[edge[1]]),G_prime[edge[0]][edge[1]]['nr_links'] ,G_prime[edge[0]][edge[1]]['gap_dist']
            G_[edge[0]][edge[1]]['nr_links'] = G_prime[edge[0]][edge[1]]['nr_links'] 
            
            try:
                G_[edge[0]][edge[1]]['gap_dist'] = G_prime[edge[0]][edge[1]]['gap_dist']
            except KeyError:
                #may be the two different sides of a cotig (has no gap dist)
                pass
            #G_.add_edge([(G_prime[edge[0]],G_prime[edge[1]])])
            #print G_[edge[0]][edge[1]],G_prime[edge[0]][edge[1]]['nr_links'] ,G_prime[edge[0]][edge[1]]['gap_dist']
        #print G_.edges()
        for node in G_:   
            if len(G_.neighbors(node))==1:
                start = node
                break
        #print '\n',G_.edges(),'\n' 
        prev_node = ('','') 
        pos = 0
        scaffold_length=0
        contig_list=[]
        #(G,contig_list,scaffold_length)=UpdateInfo(G_,Contigs,small_contigs,Scaffolds,small_scaffolds,start,prev_node,pos,contig_list,scaffold_length,dValuesTable,param)
        #S=Scaffold.scaffold(param.scaffold_indexer,contig_list,scaffold_length,{},{})  #Create the new scaffold object 

        #Scaffolds[S.name]=S        #include in scaffold library
        
#        for node_ in path:
#            try:
#                contigs = Scaffolds[node_[0]].contigs
#                for cont in contigs:
#                    print cont.name
#            except KeyError:
#                contigs = small_scaffolds[node_[0]].contigs
#                for cont in contigs:
#                    print cont.name
        

        

    #high_score_path,bad_links,score,path_len = ELS.BetweenScaffolds(G, G_prime,small_scaffolds, start, end, already_visited,param.ins_size_threshold )
    
    return()

def GiveLinkConnection(Contigs,contig_objects1, contig_objects2, side1, side2, nr_links):
    if side1 == 'R' and side2 == 'L':        
        max_pos=0
        for contig in contig_objects1:
            if contig.position >= max_pos:            
                linking_contig1=contig
                max_pos=contig.position
        min_pos=sys.maxint
        for contig in contig_objects2:
            if contig.position <= min_pos:            
                linking_contig2=contig  
                min_pos=contig.position
        
        linking_contig1.links[linking_contig2.name]=nr_links
        linking_contig2.links[linking_contig1.name]=nr_links 
        
        #print linking_contig1.name, linking_contig2.name, nr_links
        
                        
    elif side1 == 'L' and side2 == 'R':        
        max_pos=0
        for contig in contig_objects2:
            if contig.position >= max_pos:            
                linking_contig2=contig
                max_pos=contig.position
                
        min_pos=sys.maxint
        for contig in contig_objects1:
            if contig.position <= min_pos:            
                linking_contig1=contig  
                min_pos=contig.position       
        linking_contig1.links[linking_contig2.name]=nr_links 
        linking_contig2.links[linking_contig1.name]=nr_links 
        
        #print linking_contig1.name, linking_contig2.name,nr_links
        
    elif side1 == 'R' and side2 == 'R':        
        max_pos=0
        for contig in contig_objects1:
            if contig.position >= max_pos:            
                linking_contig1=contig
                max_pos=contig.position
                
        max_pos=0
        for contig in contig_objects2:
            if contig.position >= max_pos:            
                linking_contig2=contig
                max_pos=contig.position     
        linking_contig1.links[linking_contig2.name]=nr_links 
        linking_contig2.links[linking_contig1.name]=nr_links
        
        #print linking_contig1.name, linking_contig2.name , nr_links
  
    elif side1 == 'L' and side2 == 'L':        
        min_pos=sys.maxint
        for contig in contig_objects1:
            if contig.position <= min_pos:            
                linking_contig1=contig  
                min_pos=contig.position              
        min_pos=sys.maxint
        for contig in contig_objects2:
            if contig.position <= min_pos:            
                linking_contig2=contig
                min_pos=contig.position  
        
        linking_contig1.links[linking_contig2.name]=nr_links 
        linking_contig2.links[linking_contig1.name]=nr_links 

             
        #print linking_contig1.name, linking_contig2.name , nr_links

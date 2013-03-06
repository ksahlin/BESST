'''
Created on Mar 5, 2013

@author: ksahlin
'''

from Norm import normpdf,normcdf
import networkx as nx
#G = RemoveSpuriousEdges(G,G_prime,Scaffolds,Contigs,Information,param)    #step2 Remove edges that has very few edges compared to what the coverage and gap estimate says (should I also remove those with too many compared to what GapEst says? -> implement later)
#G = RemoveAmbiguousRegions(G,G_prime,Information,param) #step2

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
            #obs_list = map(lambda (x, y): x + y, G[edge[0]][edge[1]]['obs'])
            #G[edge[0]][edge[1]]['obs'] = obs_list
            tot_obs=G[edge[0]][edge[1]]['obs']
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
 



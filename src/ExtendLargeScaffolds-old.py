'''
Created on Jun 4, 2012

@author: ksahlin
'''

from time import time
import networkx as nx


def Insert_path(a, x, path, path_len, lo=0, hi=None): 
    if len(a) == 0:
        a.append([x,path_len, path])
        return()        
    if x > a[-1][0]:
        a.append([x,path_len, path])
        return()
    if x < a[0][0]:
        a.insert(0, [x,path_len, path])
        return()
    if hi is None: 
        hi = len(a) 
    while lo < hi: 
        mid = (lo + hi) // 2 
        midval = a[mid][0] 
        if midval < x: 
            lo = mid + 1 
        elif midval > x:  
            hi = mid 
        else: 
            #exact same score as some other path, look what path has the most nr of contigs and prefer that one
            if len(a[mid][2]) > len(path):
                #path to be inserted is shorter
                hi = mid
                #a.insert(mid, [x, path])
            else:
                lo = mid + 1
                #a.insert(mid + 1, [x, path])
            #return () 
    if lo == len(a):
        a.append([x,path_len, path])
    elif x > a[lo]:
        a.insert(lo + 1, [x,path_len, path])
        #print mid
    else:
        a.insert(lo, [x,path_len, path])
    return () 

def ScorePaths(G, nodes_present_in_path, paths, all_paths_sorted_wrt_score):
    if len(paths) == 0:
        return 0
    def CalculateSpanScore(path, G): 
        contig_end_prev = path[-1] #(-1, 'O')
        #contig_start_prev = path[0]  #(-1, 'O') 
        spanned_contigs_in_path = set([])
        for j, contig_end in reversed(list(enumerate(path))):
            if contig_end_prev[0] != contig_end[0]:
                contig_end_prev = contig_end
                continue
            contig_start_prev = path[0]
            for i, contig_start in list(enumerate(path)):
                if i < j:
                    #print contig_start                
                    if contig_start_prev[0] != contig_start[0]:
                        contig_start_prev = contig_start
                        continue
                    #print i, contig_start
                    if contig_end in G[contig_start]:
                        for contig in path[i + 1:j-1]:
                            spanned_contigs_in_path.add(contig[0])
                else:
                    break
                contig_start_prev = contig_start                           
            contig_end_prev = contig_end

            
        #print spanned_contigs_in_path 
        score = len(spanned_contigs_in_path)        
        return(score)
    
    def CalculateNbrScore(path, nodes_present_in_path, path_start, path_end, G):
        tot_nbr_ratio = 0
        #print path
        tot_nbrs = 0
        edges_not_in_path = set()
        
        # need to traverse the path once forward and once reverse to check for inconsistent edges
        
        #forward traversal
        edges_already_considered = set()
        prev_contig_end = path[0]
        for contig_end in path: #[1:-1]:
            for nbr in G.neighbors(contig_end):
                if nbr[0] != contig_end[0] and (nbr, contig_end) not in edges_already_considered: 
                    #print (contig_end, nbr)
                    tot_nbrs += 1
                    edges_already_considered.add((contig_end, nbr))
                    edges_already_considered.add((nbr, contig_end))
                    if nbr not in nodes_present_in_path[(path_start, path_end)] or (prev_contig_end[0] != contig_end[0] and contig_end != path[-1] and nbr in path):
                        # end criteria in this if-statement is to get rid of edges that goes from a contig start
                        # over to some contig later in the path considered (but the and statement is to include the last contig in path). Example: path (1R->2R->2L->3L->3R->4L) and there is an
                        #edge from 2L to any of 3L,3R or 4L, or an edge from 3L to 4L. (suggests inversions or other strange regions?)
                        #print contig_end, nbr
                        #tot_nbrs_not_in_path += 1 
                        edges_not_in_path.add((contig_end, nbr))  
                        edges_not_in_path.add((nbr, contig_end))
            prev_contig_end = contig_end    
            
        
        # Reverse traversal
        edges_already_considered = set()
        prev_contig_end = path[-1]
        for contig_end in reversed(path):
            for nbr in G.neighbors(contig_end):
                if nbr[0] != contig_end[0] and (nbr, contig_end) not in edges_already_considered: 
                    #print (contig_end, nbr)
                    #tot_nbrs += 1
                    edges_already_considered.add((contig_end, nbr))
                    edges_already_considered.add((nbr, contig_end))
                    if nbr not in nodes_present_in_path[(path_start, path_end)] or (prev_contig_end[0] != contig_end[0] and contig_end != path[0] and nbr in path):
                        # end criteria in this if-statement is to get rid of edges that goes from a contig start
                        # over to some contig later in the path considered (but the and statement is to include the last contig in path). Example: path (1R->2R->2L->3L->3R->4L) and there is an
                        #edge from 2L to any of 3L,3R or 4L, or an edge from 3L to 4L. (suggests inversions or other strange regions?)
                        #print contig_end, nbr
                        #tot_nbrs_not_in_path += 1 
                        edges_not_in_path.add((contig_end, nbr))  
                        edges_not_in_path.add((nbr, contig_end))
            prev_contig_end = contig_end           
            
        tot_nbrs_not_in_path =  len(edges_not_in_path)/2.0             
        if tot_nbrs != 0:  
            #print tot_nbrs, tot_nbrs_not_in_path
            tot_nbr_ratio += (tot_nbrs - tot_nbrs_not_in_path) / float(tot_nbrs)    
        else:
            #direct path between two scaffolds with larger contigs
            print 'These contigs/scaffolds should have been joined earlier.. why did we get here?'
            
        return tot_nbr_ratio
    
    for path_ in paths:
        path = path_[0]
        path_start = path[0]
        path_end = path[-1]
        path_len = path_[1]
        #calculate spanning score s_ci
        span_score = CalculateSpanScore(path, G)
        #calculate neighbour explaining path score r_ci
        nbr_score = CalculateNbrScore(path, nodes_present_in_path, path_start, path_end, G)
        #print nbr_score, span_score
        if len(path) > 2: #startnode and end node are not directly connected
            tot_score = (span_score) /  ((len(path) - 2) / 2.0) + nbr_score
        else:
            tot_score = 0 #if they are directly connected, they should by definition already be in scaffold together
            
        #insert path in sorted list using binary search
        Insert_path(all_paths_sorted_wrt_score, tot_score, path,path_len)

    return ()

def find_all_paths_for_start_node(graph, small_scaffolds, start, end, nodes_present_in_path, already_visited,is_withing_scaf,max_path_length_allowed): 
    path = [] 
    paths = [] 
    if start[1] == 'L':
        forbidden = set()
        forbidden.add( (start[0], 'R') )
    else:
        forbidden = set()
        forbidden.add( (start[0], 'L') )

    #Joining within scaffolds
    if is_withing_scaf:
        element = end.pop()
        end.add(element)
        if element[1] == 'L':
            forbidden.add( (element[0], 'R') )
        else:
            forbidden.add( (element[0], 'L') ) 
    
    start_node = start
    
    #TODO: Have length criteria that limits the path lenght due to complecity reasons. Can also identify strange
    #links by looking how many neighbors a contig has and how mych the library actually can span
    path_len = 0
    queue = [(start, end, path, path_len)]#, sum_path)]
    #prev_node = start
    counter = 0
    while queue: 
        #prev_node = start
        counter += 1
        if counter % 100000 == 0:
            print 'Paths:', counter
        start, end, path, path_len = queue.pop() #start, end, path, sum_path = queue.pop()  
        try:
            prev_node = path[-1]
        except IndexError:
            prev_node = start
        path = path + [start]          
        #print 'PATH', path ,'end', end 
        if path_len > max_path_length_allowed: #All possible paths can be exponential!! need something to stop algorithm in time
            continue
        #if score < score_best_path: # need something to stop a bad path
        #    continue
        if start in already_visited or start in forbidden:
            continue
        
        if start in end:
            if (start_node, start) in nodes_present_in_path:
                nodes_present_in_path[(start_node, start)] = nodes_present_in_path[(start_node, start)].union(path)
            else:
                nodes_present_in_path[(start_node, start)] = set(path)
                
            #print 'FOUND', path , path_len#, sum_path
            #score = ScorePath(G, path)
            paths.append((path,path_len)) 
            #nodes_present_in_path[(start_node, start)] = nodes_present_in_path[(start_node, start)].union(path)
            continue          
        if  prev_node[0] != start[0]:
            if start[1] == 'L' and (start[0], 'R') not in forbidden: 
                queue.append(((start[0], 'R'), end, path,path_len)) #, sum_path + graph[start][(start[0], 'R')]['nr_links']))
            elif start[1] == 'R' and (start[0], 'L') not in forbidden: 
                queue.append(((start[0], 'L'), end, path,path_len))#, sum_path + graph[start][(start[0], 'L')]['nr_links']))                
        else:
            for node in set(graph[start]).difference(path): 
                if node not in forbidden: # and node not in already_visited: 
                    try: # if last node (i.e. "end") it is not present in small_scaffolds and it should not be included in the length
                        queue.append((node, end, path, path_len + small_scaffolds[node[0]].s_length))   
                    except KeyError:
                        queue.append((node, end, path, path_len ))
                          
    return paths

def ExtendScaffolds(all_paths_sorted_wrt_score):
#    for score_and_path in reversed(all_paths_sorted_wrt_score):
#        print 'Score: ',score_and_path[0] #,'Path: ' , score_and_path[1]
        #MakeScaffolds()?
    n =len(all_paths_sorted_wrt_score)
    if n > 0:
        print'Total nr of paths found: ', n,
        for i in reversed(xrange(0,n)):
            if all_paths_sorted_wrt_score[i][0] < 1.5:
                print 'Paths with score equal or over 1.5: ', n-(i+1)
                break
        else:
            print 'All paths had score equal or over 1.5: '
            
    return(all_paths_sorted_wrt_score)


def BetweenScaffolds(G,G_prime,small_scaffolds):
    end = set()
    for node in G:
        end.add(node)
    
    # here we should have a for loop looping over all start nodes. Start nodes already examined should be removed in a nice way to skip over counting
    iter_nodes = end.copy()
    already_visited = set()
    
    all_paths_sorted_wrt_score = []
    while len(iter_nodes) > 1: 
        start_node = iter_nodes.pop()
        #print 'START NODE: ', start_node 
        nodes_present_in_path = {}
        paths = find_all_paths_for_start_node(G_prime,small_scaffolds, start_node, end.difference(set([start_node])), nodes_present_in_path, already_visited,0,2**32)
        already_visited.add(start_node) 
        #print 'START NODE: ', start_node, 'Tot nr of paths for this start node: ', len(paths)
        ScorePaths(G_prime, nodes_present_in_path, paths, all_paths_sorted_wrt_score)
        #print  all_paths_sorted_wrt_score
    all_paths_sorted_wrt_score = ExtendScaffolds(all_paths_sorted_wrt_score)
    all_paths_found = [all_paths_sorted_wrt_score[i][2] for i in range(0,len(all_paths_sorted_wrt_score))]
    for path in all_paths_found:
        print path
    return(all_paths_found)

def WithinScaffolds(G,G_prime,small_scaffolds,start,end_node,already_visited,max_path_length):
    end = set()
    end.add(end_node)
    nodes_present_in_path = {}
    all_paths_sorted_wrt_score = []
    already_visited.difference_update(set([start,end_node]))
    paths = find_all_paths_for_start_node(G_prime,small_scaffolds, start, end, nodes_present_in_path, already_visited,1,max_path_length)
    already_visited.add(start)
    already_visited.add(end_node)
    #print paths
    if len(paths) > 1:
        ScorePaths(G_prime, nodes_present_in_path, paths, all_paths_sorted_wrt_score)
        all_paths_sorted_wrt_score = ExtendScaffolds(all_paths_sorted_wrt_score)  
#        for path in all_paths_sorted_wrt_score:
#            print path
        return(all_paths_sorted_wrt_score[-1][2],all_paths_sorted_wrt_score[-1][1],all_paths_sorted_wrt_score[-1][0]) #return(all_paths_sorted_wrt_score) #
    return([],0,0)

if __name__ == '__main__':
    import Scaffold 
    small_scaffolds_test ={} 
    for i in range(1,7):
        S = Scaffold.scaffold(i,0,0,{},{})
        small_scaffolds_test[S.name] = S 
    start = time()   
    G_prime = nx.Graph() 
    #G.add_nodes_from([(1, 'L'), (1, 'R'), (2, 'L'), (2, 'R'), (3, 'L'), (3, 'R'), (4, 'L'), (4, 'R'), (5, 'L'), (5, 'R')]) 
    for i in range(1, 7):
        G_prime.add_edge((i, 'L'), (i, 'R'), {'nr_links':0})
    G_prime.add_edges_from([((1, 'R'), (2, 'R'), {'nr_links':1}), ((3, 'L'), (4, 'L'), {'nr_links':1}), ((2, 'L'), (3, 'R'), {'nr_links':1}), ((1, 'R'), (5, 'L'), {'nr_links':2}),
                       ((5, 'R'), (4, 'L'), {'nr_links':3}), ((2, 'L'), (5, 'L'), {'nr_links':2}), ((1, 'R'), (4, 'L'), {'nr_links':8}), ((2, 'L'), (6, 'L'), {'nr_links':3}),
                       ((1, 'L'), (4, 'R'), {'nr_links':1}), ((1, 'L'), (4, 'L'), {'nr_links':1}), ((3, 'L'), (4, 'R'), {'nr_links':1}),
                        ((1, 'R'), (2, 'L'), {'nr_links':1}), ((1, 'R'), (5, 'R'), {'nr_links':1}),((2, 'L'), (5, 'R'), {'nr_links':1})]) 
    G = nx.Graph()
    G.add_nodes_from([(1, 'L'), (1, 'R'), (4, 'L'), (4, 'R'), (6, 'L'), (6, 'R')])
    contigs = [1, 2, 3, 4, 5, 6]
    
    print 'Between'
    BetweenScaffolds(G,G_prime,small_scaffolds_test) 
    start_node = (1, 'R')
    end_node = (4, 'L')
    print 'Within'
    already_visited = set(G.nodes())
    print already_visited
    WithinScaffolds(G,G_prime,small_scaffolds_test,start_node,end_node,already_visited,0)
    elapsed = time() - start
    print 'time all paths: ', elapsed



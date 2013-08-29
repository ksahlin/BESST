
'''
    Created on Jun 21, 2012
    
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

import time
import networkx as nx
from collections import defaultdict
import multiprocessing as mp

def Insert_path(a, x, path, bad_link_weight, path_len, lo=0, hi=None):
    if len(a) == 0:
        a.append([x, bad_link_weight, path, path_len])
        return()
    if x > a[-1][0]:
        a.append([x, bad_link_weight, path, path_len])
        return()
    if x < a[0][0]:
        a.insert(0, [x, bad_link_weight, path, path_len])
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
            #exact same score as some other path, look what path has the  nr of bad links and prefer that one
            if a[mid][1] < bad_link_weight:
                #path to be inserted has more bad links
                hi = mid
                #a.insert(mid, [x, path])
            else:
                lo = mid + 1
                #a.insert(mid + 1, [x, path])
            #return () 

    #print 'insert on pos:', lo, len(a)
    if lo == len(a):
        a.append([x, bad_link_weight, path, path_len])
    elif x > a[lo]:
        a.insert(lo + 1, [x, bad_link_weight, path, path_len])
        #print mid
    else:
        a.insert(lo, [x, bad_link_weight, path, path_len])
    return ()

def ScorePaths(G, nodes_present_in_path, paths, all_paths_sorted_wrt_score):
    if len(paths) == 0:
        return ()

    def CalculateConnectivity(path, G):
        good_edges = set(map(lambda i: path[i], filter(lambda i: i % 2 == 1, range(len(path)))))
        #print good_edges, 'path len:',len(path)- 2, 'nr of good edges: ', len(good_edges)*(len(good_edges)+1)/2.0
        good_edges_count = 0
        bad_edges_count = 0
        prev_node = path[0]
        good_edges_already_considered = set()
        bad_edges_already_considered = set()
        bad_link_weight = 0
        good_link_weight = 0
        link_weights = defaultdict(lambda : defaultdict(float))

        for node in path:
            if node[0] == prev_node[0]:
                #TODO: eventually calculate if global optima here. We then need a hash table good_e ={node:[tot_edge = set(nbr1,nbr2..), bad_edge = set(nbr3,..)]}
                for nbr in G.neighbors(node):
                    if nbr in good_edges and not nbr[0] == node[0]:
                        good_edges_count += 1
                        good_link_weight += G[node][nbr]['nr_links']
                        link_weights[node[0]]['good'] += G[node][nbr]['nr_links']
                        good_edges_already_considered.add((nbr, node))
                        good_edges_already_considered.add((node, nbr))
                    elif not nbr[0] == node[0] and (nbr, node) not in bad_edges_already_considered:
                        bad_edges_count += 1
                        bad_link_weight += G[node][nbr]['nr_links']
                        link_weights[node[0]]['bad'] += G[node][nbr]['nr_links']
                        bad_edges_already_considered.add((nbr, node))
                        bad_edges_already_considered.add((node, nbr))
            else:
                good_edges.remove(node)
                for nbr in G.neighbors(node):
                    if (nbr, node) not in good_edges_already_considered and not nbr[0] == node[0] and (nbr, node) not in bad_edges_already_considered:
                        bad_edges_count += 1
                        bad_link_weight += G[node][nbr]['nr_links']
                        link_weights[node[0]]['bad'] += G[node][nbr]['nr_links']
                        bad_edges_already_considered.add((nbr, node))
                        bad_edges_already_considered.add((node, nbr))
                    elif (nbr, node) in good_edges_already_considered:
                        link_weights[node[0]]['good'] += G[node][nbr]['nr_links']
            prev_node = node

        #print 'Edges good/bad:',good_edges_count,bad_edges_count,'links good/bad:', good_link_weight, bad_link_weight
        #print 'contig weights:'
        for key in link_weights:
            try:
                weight = link_weights[key]['good'] / link_weights[key]['bad']
                #print 'contig_weight:', weight
            except ZeroDivisionError:
                weight = 2 ** 16
                #print 'contig_weight: inf (no bad edges)'
            if weight < 2.0:
                return 0, 0

        try:
            score = good_link_weight / float(bad_link_weight)
        except ZeroDivisionError:
            score = good_link_weight

        #return score, bad_link_weight
        return good_link_weight - bad_link_weight, bad_link_weight
        #return good_edges_count - bad_edges_count, bad_link_weight


    #print '\nSTARTING scoring paths:'
    for path_ in paths:
        path = path_[0]
        path_len = path_[1]
        #calculate spanning score s_ci
        score, bad_link_weight = CalculateConnectivity(path, G)
        if len(path) > 2 and score > 0: #startnode and end node are not directly connected
            Insert_path(all_paths_sorted_wrt_score, score, path , bad_link_weight, path_len)

#        if len(path) > 2: #startnode and end node are not directly connected
#            tot_score = (score)*2 /  ((len(path) / 2.0))  # formula for score
#        else:
#            tot_score = 0 #if they are directly connected, they should by definition already be in scaffold together
#            
#        #insert path in sorted list using binary search
#        Insert_path(all_paths_sorted_wrt_score, tot_score, path,bad_link_weight,path_len)
    #print 'This is the stats for the best path: ',all_paths_sorted_wrt_score[-1]
    return ()

def find_all_paths_for_start_node(graph, start, end, nodes_present_in_path, already_visited, is_withing_scaf, max_path_length_allowed, param):
    path = []
    paths = []
    if start[1] == 'L':
        forbidden = set()
        forbidden.add((start[0], 'R'))
    else:
        forbidden = set()
        forbidden.add((start[0], 'L'))

    #Joining within scaffolds
    if is_withing_scaf:
        element = end.pop()
        end.add(element)
        if element[1] == 'L':
            forbidden.add((element[0], 'R'))
        else:
            forbidden.add((element[0], 'L'))

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
        #if counter % 100 == 0:
        #    print 'Potential paths:', counter, 'paths found: ', len(paths)
        if counter > param.path_threshold or len(path) > 100:
            longest_path = 0
            if len(paths) > 0:
                longest_path = len(paths[-1][0])
            #print 'Stopping at', counter, 'iterations..', 'longest path added so far (nr_of_contigs/2 incl in path):', longest_path
            break
        start, end, path, path_len = queue.pop() #start, end, path, sum_path = queue.pop()  
        try:
            prev_node = path[-1]
        except IndexError:
            prev_node = start
        path = path + [start]
        path_len = len(path)
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
            paths.append((path, path_len))
            #nodes_present_in_path[(start_node, start)] = nodes_present_in_path[(start_node, start)].union(path)
            continue
        if  prev_node[0] != start[0]:
            if start[1] == 'L' and (start[0], 'R') not in forbidden:
                queue.append(((start[0], 'R'), end, path, path_len)) #, sum_path + graph[start][(start[0], 'R')]['nr_links']))
            elif start[1] == 'R' and (start[0], 'L') not in forbidden:
                queue.append(((start[0], 'L'), end, path, path_len))#, sum_path + graph[start][(start[0], 'L')]['nr_links']))                
        else:
            for node in set(graph[start]).difference(path):
                if node not in forbidden: # and node not in already_visited: 
                    try: # if last node (i.e. "end") it is not present in small_scaffolds and it should not be included in the length
                        queue.append((node, end, path, path_len + graph[node[0]]['length'])) #  small_scaffolds[node[0]].s_length))   #
                    except KeyError:
                        queue.append((node, end, path, path_len))

    return paths

def ExtendScaffolds(all_paths_sorted_wrt_score):
#    for score_and_path in reversed(all_paths_sorted_wrt_score):
#        print 'Score: ',score_and_path[0] #,'Path: ' , score_and_path[1]
        #MakeScaffolds()?
    n = len(all_paths_sorted_wrt_score)
    if n > 0:
        print'Total nr of paths found: ', n
#        for i in reversed(xrange(0,n)):
#            if all_paths_sorted_wrt_score[i][0] < 0:
#                print 'Paths with score equal or over 0: ', n-(i+1)
#                break
#        else:
#            print 'All paths had score equal or over 0: '

    return(all_paths_sorted_wrt_score)




def BetweenScaffolds(G_prime, end, iter_nodes, param):
    # here we should have a for loop looping over all start nodes. Start nodes already examined should be removed in a nice way to skip over counting
    already_visited = set()
    all_paths_sorted_wrt_score = []
    print 'Entering "find_all_paths_for_start_node" '
    iter_threshold = 0
    cnter = 0
    while len(iter_nodes) > 0 and iter_threshold <= 100000:
        iter_threshold += 1
        start_node = iter_nodes.pop()
        #print 'START NODE: ', start_node 
        nodes_present_in_path = {}
        cnter += 1
        if cnter % 100 == 0:
            print 'enter Betwween scaf node: ', cnter
        paths = find_all_paths_for_start_node(G_prime, start_node, end.difference(set([start_node])), nodes_present_in_path, already_visited, 0, 2 ** 32, param)
        already_visited.add(start_node)
        #print 'START NODE: ', start_node, 'Tot nr of paths for this start node: ', len(paths)
        ScorePaths(G_prime, nodes_present_in_path, paths, all_paths_sorted_wrt_score)
        #print  all_paths_sorted_wrt_score
    all_paths_sorted_wrt_score = ExtendScaffolds(all_paths_sorted_wrt_score)
#    all_paths_found = [all_paths_sorted_wrt_score[i][2] for i in range(0,len(all_paths_sorted_wrt_score))]
#    for path in all_paths_found:
#        print path

    return(all_paths_sorted_wrt_score)

def WithinScaffolds(G, G_prime, start, end_node, already_visited, max_path_length, param):
    end = set()
    end.add(end_node)
    nodes_present_in_path = {}
    all_paths_sorted_wrt_score = []
    already_visited.difference_update(set([start, end_node]))
    paths = find_all_paths_for_start_node(G_prime, start, end, nodes_present_in_path, already_visited, 1, max_path_length, param)
    already_visited.add(start)
    already_visited.add(end_node)
    #print paths
    if len(paths) > 1:
        ScorePaths(G_prime, nodes_present_in_path, paths, all_paths_sorted_wrt_score)
#        for path in all_paths_sorted_wrt_score:
#            print path
        if len(all_paths_sorted_wrt_score) > 0:
            #all_paths_sorted_wrt_score = ExtendScaffolds(all_paths_sorted_wrt_score)
            return(all_paths_sorted_wrt_score[-1][2], all_paths_sorted_wrt_score[-1][1], all_paths_sorted_wrt_score[-1][0], all_paths_sorted_wrt_score[-1][3]) #return(all_paths_sorted_wrt_score) #
    return([], 0, 0, 0)

if __name__ == '__main__':
    import Scaffold
    small_scaffolds_test = {}
    for i in range(1, 7):
        S = Scaffold.scaffold(i, 0, 0, {}, {})
        small_scaffolds_test[S.name] = S
    start = time()
    G_prime = nx.Graph()
    #G.add_nodes_from([(1, 'L'), (1, 'R'), (2, 'L'), (2, 'R'), (3, 'L'), (3, 'R'), (4, 'L'), (4, 'R'), (5, 'L'), (5, 'R')]) 
    for i in range(1, 7):
        G_prime.add_edge((i, 'L'), (i, 'R'), {'nr_links':0})
    G_prime.add_edges_from([((1, 'R'), (2, 'R'), {'nr_links':1}), ((3, 'L'), (4, 'L'), {'nr_links':1}), ((2, 'L'), (3, 'R'), {'nr_links':1}), ((1, 'R'), (5, 'L'), {'nr_links':2}),
                       ((5, 'R'), (4, 'L'), {'nr_links':3}), ((2, 'L'), (5, 'L'), {'nr_links':2}), ((1, 'R'), (4, 'L'), {'nr_links':8}), ((2, 'L'), (6, 'L'), {'nr_links':3}),
                       ((1, 'L'), (4, 'R'), {'nr_links':1}), ((1, 'L'), (4, 'L'), {'nr_links':1}), ((3, 'L'), (4, 'R'), {'nr_links':1}),
                        ((1, 'R'), (2, 'L'), {'nr_links':1}), ((1, 'R'), (5, 'R'), {'nr_links':1}), ((2, 'L'), (5, 'R'), {'nr_links':1})])
    G = nx.Graph()
    G.add_nodes_from([(1, 'L'), (1, 'R'), (4, 'L'), (4, 'R'), (6, 'L'), (6, 'R')])
    contigs = [1, 2, 3, 4, 5, 6]

    print 'Between'
    BetweenScaffolds(G, G_prime, small_scaffolds_test)
    start_node = (1, 'R')
    end_node = (4, 'L')
    print 'Within'
    already_visited = set(G.nodes())
    print already_visited
    WithinScaffolds(G, G_prime, small_scaffolds_test, start_node, end_node, already_visited, 0)
    elapsed = time() - start
    print 'time all paths: ', elapsed



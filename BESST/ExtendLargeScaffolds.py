
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


def ScorePaths(G, paths, all_paths, param):
    if len(paths) == 0:
        return ()

    def calculate_connectivity(path, G):
        """
            Contig ends has even or odd numbered placements (index) in the path list. In case
            of contamination, good links are between an even to an odd indexed contig-end, or vice versa.

        """
        bad_link_weight = 0
        good_link_weight = 0
        links_not_in_path = 0
        links_wrong_orientation_in_path = 0
        even, odd = path[::2], path[1::2]
        nodes_even = set(even)
        nodes_odd = set(odd)
        visited = set()
        for i, node in enumerate(path):
            for nbr in G.neighbors(node):
                if i % 2 == 0 and node[0] != nbr[0]:
                    if nbr in nodes_odd:
                        if nbr not in visited:
                            good_link_weight += G[node][nbr]['nr_links']
                        else:
                            pass
                    else:
                        bad_link_weight += G[node][nbr]['nr_links']
                elif i % 2 == 1 and node[0] != nbr[0]:
                    if nbr not in nodes_even:
                        bad_link_weight += G[node][nbr]['nr_links']
                    elif nbr not in visited:
                        bad_link_weight += G[node][nbr]['nr_links']
            visited.add(node)
        good_link_weight = good_link_weight
        bad_link_weight = bad_link_weight
        try:
            score = good_link_weight / float(bad_link_weight)
        except ZeroDivisionError:
            score = good_link_weight

        return score, bad_link_weight


    def calculate_connectivity_contamination(path, G):
        """
            Contig ends has even or odd numbered placements (index) in the path list. In case
            of contamination, good links are between an even to an odd indexed contig-end, or vice versa.

        """
        bad_link_weight = 0
        good_link_weight = 0
        links_not_in_path = 0
        links_wrong_orientation_in_path = 0
        even, odd = path[::2], path[1::2]
        nodes_even = set(even)
        nodes_odd = set(odd)
        for i, node in enumerate(path):
            for nbr in G.neighbors(node):
                if i % 2 == 0 and node[0] != nbr[0]:
                    if nbr in nodes_odd:
                        good_link_weight += G[node][nbr]['nr_links']
                    else:
                        bad_link_weight += G[node][nbr]['nr_links']
                elif i % 2 == 1 and node[0] != nbr[0]:
                    if nbr in nodes_even:
                        good_link_weight += G[node][nbr]['nr_links']
                    else:
                        bad_link_weight += G[node][nbr]['nr_links']
        good_link_weight = good_link_weight/2
        bad_link_weight = bad_link_weight

        # ###############


        try:
            score = good_link_weight / float(bad_link_weight)
        except ZeroDivisionError:
            score = good_link_weight

        return score, bad_link_weight


    #print '\nSTARTING scoring paths:'
    for path_ in paths:
        path = path_[0]
        path_len = path_[1]
        #calculate spanning score s_ci
        if param.contamination_ratio:
            score, bad_link_weight = calculate_connectivity_contamination(path, G)
        else:
            score, bad_link_weight = calculate_connectivity(path, G)


        if param.no_score and score >= param.score_cutoff:
            all_paths.append([score, bad_link_weight, path, path_len])
            #Insert_path(all_paths, score, path , bad_link_weight, path_len)
        elif len(path) > 2 and score >= param.score_cutoff: #startnode and end node are not directly connected
            all_paths.append([score, bad_link_weight, path, path_len])
            #Insert_path(all_paths, score, path , bad_link_weight, path_len)

    return ()

def find_all_paths_for_start_node(graph, start, end, already_visited, is_withing_scaf, max_path_length_allowed, param):
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


    #TODO: Have length criteria that limits the path lenght due to complecity reasons. Can also identify strange
    #links by looking how many neighbors a contig has and how mych the library actually can span
    path_len = 0
    queue = [(start, path, path_len)]#, sum_path)]
    #prev_node = start
    counter = 0
    while queue:
        #prev_node = start
        counter += 1
        #if counter % 100 == 0:
        #    print 'Potential paths:', counter, 'paths found: ', len(paths)
        if counter > param.path_threshold or len(path) > 100:
            print 'Hit path_threshold of {0} iterations! consider increase --iter <int> parameter to over {0} if speed of BESST is not a problem. Standard increase is, e.g., 2-10x of current value'.format(param.path_threshold)
            break
            
        start, path, path_len = queue.pop() #start, end, path, sum_path = queue.pop()  
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
            # if (start_node, start) in nodes_present_in_path:
            #     nodes_present_in_path[(start_node, start)] = nodes_present_in_path[(start_node, start)].union(path)
            # else:
            #     nodes_present_in_path[(start_node, start)] = set(path)
            paths.append((path, path_len))
            continue


        if  prev_node[0] != start[0]:
            if start[1] == 'L' and (start[0], 'R') not in forbidden:
                queue.append(((start[0], 'R'), path, path_len)) #, sum_path + graph[start][(start[0], 'R')]['nr_links']))
            elif start[1] == 'R' and (start[0], 'L') not in forbidden:
                queue.append(((start[0], 'L'), path, path_len))#, sum_path + graph[start][(start[0], 'L')]['nr_links']))                
        else:
            for node in set(graph[start]).difference(path):
                if node not in forbidden: # and node not in already_visited: 
                    try: # if last node (i.e. "end") it is not present in small_scaffolds and it should not be included in the length
                        queue.append((node, path, path_len + graph[node[0]]['length'])) #  small_scaffolds[node[0]].s_length))   #
                    except KeyError:
                        queue.append((node, path, path_len))

    return paths



def BetweenScaffolds(G_prime, end, iter_nodes, param):
    # here we should have a for loop looping over all start nodes. Start nodes already examined should be removed in a nice way to skip over counting
    already_visited = set()
    all_paths = []
    print 'Entering "find_all_paths_for_start_node" '
    iter_count = 0
    cnter = 0
    if param.max_extensions:
        iter_threshold = param.max_extensions
    else: 
        iter_threshold = len(end)

    print 'iterating until maximum of {0} extensions.'.format(iter_threshold) 
    print 'nodes:{0}, edges: {1}'.format(len(G_prime.nodes()), len(G_prime.edges()))
    while len(iter_nodes) > 0 and iter_count <= iter_threshold:
        iter_count += 1
        start_node = iter_nodes.pop()
        if cnter % 100 == 0:
            print 'enter Betwween scaf node: ', cnter
        end.difference_update(set([start_node]))
        paths = find_all_paths_for_start_node(G_prime, start_node, end, already_visited, 0, 2 ** 32, param)
        already_visited.add(start_node)
        ScorePaths(G_prime, paths, all_paths, param)
        cnter += 1
    #all_paths = ExtendScaffolds(all_paths)
    #print all_paths
    print 'Total nr of paths found: {0} with score larger than: {1}'.format(len(all_paths), param.score_cutoff)
    all_paths.sort(key=lambda list_: list_[0]) 
    #print all_paths
    return(all_paths)

def WithinScaffolds(G, G_prime, start, end_node, already_visited, max_path_length, param):
    end = set()
    end.add(end_node)
    all_paths = []
    already_visited.difference_update(set([start, end_node]))
    paths = find_all_paths_for_start_node(G_prime, start, end, already_visited, 1, max_path_length, param)
    already_visited.add(start)
    already_visited.add(end_node)
    #print paths
    if len(paths) > 1:
        ScorePaths(G_prime, paths, all_paths,param)
        all_paths.sort(key=lambda list_: list_[0]) 

        if len(all_paths) > 0:
            return all_paths

    return []

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



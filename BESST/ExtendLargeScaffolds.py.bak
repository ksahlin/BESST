
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
from collections import defaultdict, deque
import multiprocessing as mp
import heapq

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
        # print path
        # print even, odd
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
    #print 'scores:',
    for path in paths:
        #path = path_[0]
        #path_len = path_[1]
        #calculate spanning score s_ci
        if param.contamination_ratio:
            score, bad_link_weight = calculate_connectivity_contamination(path, G)
        else:
            score, bad_link_weight = calculate_connectivity(path, G)

        #print int(score),' ',
        if param.no_score and score >= param.score_cutoff:
            all_paths.append([score, bad_link_weight, path, len(path)])
            #Insert_path(all_paths, score, path , bad_link_weight, path_len)
        elif len(path) > 2 and score >= param.score_cutoff: #startnode and end node are not directly connected
            all_paths.append([score, bad_link_weight, path, len(path)])
            #Insert_path(all_paths, score, path , bad_link_weight, path_len)

    return ()

# def find_all_paths_for_start_node_BFS(graph, start, end, already_visited, is_withing_scaf, max_path_length_allowed, param):
#     path = []
#     paths = []
#     if start[1] == 'L':
#         forbidden = set()
#         forbidden.add((start[0], 'R'))
#     else:
#         forbidden = set()
#         forbidden.add((start[0], 'L'))

#     #Joining within scaffolds
#     if is_withing_scaf:
#         element = end.pop()
#         end.add(element)
#         if element[1] == 'L':
#             forbidden.add((element[0], 'R'))
#         else:
#             forbidden.add((element[0], 'L'))


#     #TODO: Have length criteria that limits the path lenght due to complecity reasons. Can also identify strange
#     #links by looking how many neighbors a contig has and how mych the library actually can span
#     path_len = 0
#     queue = [(start, path, path_len)]#, sum_path)]
#     #prev_node = start
#     counter = 0
#     while queue:
#         #prev_node = start
#         counter += 1
#         #if counter % 100 == 0:
#         #    print 'Potential paths:', counter, 'paths found: ', len(paths)
#         if counter > param.path_threshold or len(path) > 100:
#             #print 'Hit path_threshold of {0} iterations! consider increase --iter <int> parameter to over {0} if speed of BESST is not a problem. Standard increase is, e.g., 2-10x of current value'.format(param.path_threshold)
#             param.hit_path_threshold = True
#             break
            
#         start, path, path_len = queue.pop() #start, end, path, sum_path = queue.pop()  
#         try:
#             prev_node = path[-1]
#         except IndexError:
#             prev_node = start
#         path = path + [start]
#         path_len = len(path)
#         #print 'PATH', path ,'end', end 
#         if path_len > max_path_length_allowed: #All possible paths can be exponential!! need something to stop algorithm in time
#             continue
#         #if score < score_best_path: # need something to stop a bad path
#         #    continue
#         if start in already_visited or start in forbidden:
#             continue

#         if start in end:
#             # if (start_node, start) in nodes_present_in_path:
#             #     nodes_present_in_path[(start_node, start)] = nodes_present_in_path[(start_node, start)].union(path)
#             # else:
#             #     nodes_present_in_path[(start_node, start)] = set(path)
#             paths.append((path, path_len))
#             continue


#         if  prev_node[0] != start[0]:
#             if start[1] == 'L' and (start[0], 'R') not in forbidden:
#                 queue.append(((start[0], 'R'), path, path_len)) #, sum_path + graph[start][(start[0], 'R')]['nr_links']))
#             elif start[1] == 'R' and (start[0], 'L') not in forbidden:
#                 queue.append(((start[0], 'L'), path, path_len))#, sum_path + graph[start][(start[0], 'L')]['nr_links']))                
#         else:
#             for node in set(graph[start]).difference(path):
#                 if node not in forbidden: # and node not in already_visited: 
#                     try: # if last node (i.e. "end") it is not present in small_scaffolds and it should not be included in the length
#                         queue.append((node, path, path_len + graph[node[0]]['length'])) #  small_scaffolds[node[0]].s_length))   #
#                     except KeyError:
#                         queue.append((node, path, path_len))

#     return paths

def find_all_paths_for_start_node_BFS_Dynamic_Programming_ish(graph, start, end, already_visited, is_withing_scaf, max_path_length_allowed, param):
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
    queue = deque(((start, path, path_len, (set([start]), set(), 0, 0)),))#, sum_path)]
    #prev_node = start
    counter = 0
    head_dict = {} # contains current best paths for 3_prime contigs

    while queue:
        #prev_node = start
        counter += 1
        #if counter % 100 == 0:
        #    print 'Potential paths:', counter, 'paths found: ', len(paths)
        if counter > param.path_threshold or len(path) > 100:
            #print 'Hit path_threshold of {0} iterations! consider increase --iter <int> parameter to over {0} if speed of BESST is not a problem. Standard increase is, e.g., 2-10x of current value'.format(param.path_threshold)
            param.hit_path_threshold = True
            break
            
        start, path, path_len, (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs, bad_link_count) = queue.popleft() #start, end, path, sum_path = queue.pop()

        if start in head_dict:
            # strictly worse path than previous seen one
            if nr_bad_nbrs > head_dict[start][0] and bad_link_count > head_dict[start][1]:
                continue

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
            paths.append(path)
            continue


        if  prev_node[0] != start[0]:


            if start[1] == 'L' and (start[0], 'R') not in forbidden:
                additional_bad_nbrs = 0 
                additional_bad_link_count = 0
                for nbr in graph.neighbors(start):
                    if nbr[0] != start[0] and nbr not in ctg_ends_in_path and nbr not in bad_ctgs:
                        additional_bad_nbrs += 1
                        additional_bad_link_count += graph[start][nbr]['nr_links']
                        bad_ctgs.add(nbr)

                if (start[0], 'R') not in head_dict:
                    head_dict[(start[0], 'R')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                    head_dict[(start[0], 'L')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                    # print "Here", (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count), start[0]
                    ctg_ends_in_path.add((start[0], 'R'))
                    queue.append(((start[0], 'R'), path, path_len, (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)))

                elif (start[0], 'R') in head_dict:
                    # print 'Were in!!'
                    # print "OK", (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count), start[0]
                    if nr_bad_nbrs + additional_bad_nbrs > head_dict[(start[0], 'R')][0] and bad_link_count + additional_bad_link_count > head_dict[(start[0], 'R')][1]:
                        #print "Skipping path here"
                        continue
                    # stictly better than any previous path, continue exploring  and add new strict minimum     
                    elif nr_bad_nbrs + additional_bad_nbrs < head_dict[(start[0], 'R')][0] and bad_link_count + additional_bad_link_count < head_dict[(start[0], 'R')][1]:
                        head_dict[(start[0], 'R')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                        head_dict[(start[0], 'L')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                        ctg_ends_in_path.add((start[0], 'R'))
                        queue.append(((start[0], 'R'), path, path_len, (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)))
                    else:  
                        # print "FINAL LEVEL"
                        ctg_ends_in_path.add((start[0], 'R'))
                        queue.append(((start[0], 'R'), path, path_len, (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)))


            elif start[1] == 'R' and (start[0], 'L') not in forbidden:
                additional_bad_nbrs = 0 
                additional_bad_link_count = 0
                for nbr in graph.neighbors(start):
                    if nbr[0] != start[0] and nbr not in ctg_ends_in_path and nbr not in bad_ctgs:
                        additional_bad_nbrs += 1 
                        additional_bad_link_count += graph[start][nbr]['nr_links']
                        bad_ctgs.add(nbr)
                        
                if (start[0], 'L') not in head_dict:
                    head_dict[(start[0], 'L')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                    head_dict[(start[0], 'R')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                    # print "Here", (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count), start[0]
                    ctg_ends_in_path.add((start[0], 'L'))
                    queue.append(((start[0], 'L'), path, path_len, (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)))

                elif (start[0], 'L') in head_dict:
                    # print 'Were in!!'
                    # stictly worse than any previous path 
                    if nr_bad_nbrs + additional_bad_nbrs > head_dict[(start[0], 'L')][0] and bad_link_count + additional_bad_link_count > head_dict[(start[0], 'L')][1]:
                        #print "Skipping path here"
                        continue

                    # stictly better than any previous path, continue exploring  and add new strict minimum 
                    elif nr_bad_nbrs + additional_bad_nbrs < head_dict[(start[0], 'L')][0] and bad_link_count + additional_bad_link_count < head_dict[(start[0], 'L')][1]:
                        head_dict[(start[0], 'L')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                        head_dict[(start[0], 'R')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                        ctg_ends_in_path.add((start[0], 'L'))
                        queue.append(((start[0], 'L'), path, path_len, (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)) ) 

                    # ambiguous continue exploring traversal
                    else:
                        # print "FINAL LEVEL"
                        ctg_ends_in_path.add((start[0], 'L'))
                        queue.append(((start[0], 'L'), path, path_len, (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)) ) 



            # if start[1] == 'L' and (start[0], 'R') not in forbidden:
            #     queue.append(((start[0], 'R'), path, path_len)) 
            # elif start[1] == 'R' and (start[0], 'L') not in forbidden:
            #     queue.append(((start[0], 'L'), path, path_len))          
        else:
            # nbrs = set(graph[start]).difference(path)
            # nr_nbrs = len(nbrs)
            # if nr_nbrs > 1:
            #     # print 'Complex area'
            #     least_nr_of_links = param.edgesupport
            # else:
            #     least_nr_of_links = 0
            for node in set(graph[start]).difference(path): #nbrs:
                if node not in forbidden: # and node not in already_visited:
                    # nr_links = graph[start][node]['nr_links']
                    # if  nr_links >= least_nr_of_links:
                    path_new = path + [node]
                    ctg_ends_in_path_new = set(path_new)
                    queue.append((node, path, path_len, (ctg_ends_in_path_new, bad_ctgs, nr_bad_nbrs, bad_link_count) ))

            # for node in set(graph[start]).difference(path):
            #     if node not in forbidden: # and node not in already_visited: 
            #         try: # if last node (i.e. "end") it is not present in small_scaffolds and it should not be included in the length
            #             queue.append((node, path, path_len + graph[node[0]]['length'])) #  small_scaffolds[node[0]].s_length))   #
            #         except KeyError:
            #             queue.append((node, path, path_len))
    # print "paths:", len(paths)
    return paths


# def find_all_paths_for_start_node_BFS_improved(graph, start, end, already_visited, is_withing_scaf, max_path_length_allowed, param):
#     paths = []
#     if start[1] == 'L':
#         forbidden = set()
#         forbidden.add((start[0], 'R'))
#     else:
#         forbidden = set()
#         forbidden.add((start[0], 'L'))

#     path = [start]
#     #Joining within scaffolds
#     if is_withing_scaf:
#         element = end.pop()
#         end.add(element)
#         if element[1] == 'L':
#             forbidden.add((element[0], 'R'))
#         else:
#             forbidden.add((element[0], 'L'))


#     #TODO: Have length criteria that limits the path lenght due to complecity reasons. Can also identify strange
#     #links by looking how many neighbors a contig has and how mych the library actually can span
#     queue = [path]#, sum_path)]
#     #prev_node = start
#     counter = 0
#     while queue:
#         counter += 1
#         if counter > param.path_threshold or len(path) > 100:
#             #print 'Hit path_threshold of {0} iterations! consider increase --iter <int> parameter to over {0} if speed of BESST is not a problem. Standard increase is, e.g., 2-10x of current value'.format(param.path_threshold)
#             param.hit_path_threshold = True
#             break
            
#         path = queue.pop() 
#         prev_node = path[-1]

#         if len(path) > 1 and prev_node in already_visited or prev_node in forbidden:
#             continue

#         if len(path) > 1  and path[-2] in end:
#             path.pop()
#             paths.append(path)
#             continue


#         # if prev_node[1] == 'L' and (prev_node[0], 'R') not in forbidden:
#         #     path.append( (prev_node[0], 'R') )
#         #     #queue.append(path)
#         # elif prev_node[1] == 'R' and (prev_node[0], 'L') not in forbidden:
#         #     path.append( (prev_node[0], 'L') )
#         #     #queue.append(path)                

#         for node in set(graph[prev_node]).difference(path):
#             if node not in forbidden: # and node not in already_visited: 
#                 new_path = path + [node]
#                 # path.append(node)
#                 if node[1] == 'L':
#                      new_path.append( (node[0], 'R') )
#                 else:
#                      new_path.append( (node[0], 'L') )

#                 queue.append(new_path)

#     return paths

# def find_all_paths_for_start_node_DFS(graph, start, end, already_visited, is_withing_scaf, max_path_length_allowed, param):
#     path = []
#     paths = []
#     max_size_heap = 0
#     if start[1] == 'L':
#         forbidden = set()
#         forbidden.add((start[0], 'R'))
#     else:
#         forbidden = set()
#         forbidden.add((start[0], 'L'))

#     #Joining within scaffolds
#     if is_withing_scaf:
#         element = end.pop()
#         end.add(element)
#         if element[1] == 'L':
#             forbidden.add((element[0], 'R'))
#         else:
#             forbidden.add((element[0], 'L'))


#     #TODO: Have length criteria that limits the path lenght due to complecity reasons. Can also identify strange
#     #links by looking how many neighbors a contig has and how mych the library actually can span
#     heap = [(0,(start, path))]#, sum_path)]
#     #prev_node = start
#     counter = 0
#     while heap:
#         #prev_node = start
#         counter += 1
#         #if counter % 100 == 0:
#         #    print 'Potential paths:', counter, 'paths found: ', len(paths)
#         if counter > param.path_threshold or len(path) > 100:
#             #print 'Hit path_threshold of {0} iterations! consider increase --iter <int> parameter to over {0} if speed of BESST is not a problem. Standard increase is, e.g., 2-10x of current value'.format(param.path_threshold)
#             param.hit_path_threshold = True
#             break
            
#         nr_links, (start, path) = heapq.heappop(heap) #start, end, path, sum_path = heapq.pop()  
#         try:
#             prev_node = path[-1]
#         except IndexError:
#             prev_node = start
#         path = path + [start]
#         path_len = len(path)
#         #print 'PATH', path ,'end', end 
#         if path_len > max_path_length_allowed: #All possible paths can be exponential!! need something to stop algorithm in time
#             continue
#         #if score < score_best_path: # need something to stop a bad path
#         #    continue
#         if start in already_visited or start in forbidden:
#             continue

#         if start in end:
#             # if (start_node, start) in nodes_present_in_path:
#             #     nodes_present_in_path[(start_node, start)] = nodes_present_in_path[(start_node, start)].union(path)
#             # else:
#             #     nodes_present_in_path[(start_node, start)] = set(path)
#             paths.append(path)
#             continue


#         if  prev_node[0] != start[0]:
#             nr_links = 2**16 # large number to give high priority to this intra-contig edge in this two node representation
#             if start[1] == 'L' and (start[0], 'R') not in forbidden:
#                 heapq.heappush(heap, (nr_links, ((start[0], 'R'), path))) #, sum_path + graph[start][(start[0], 'R')]['nr_links']))
#             elif start[1] == 'R' and (start[0], 'L') not in forbidden:
#                 heapq.heappush(heap, (nr_links, ((start[0], 'L'), path)))#, sum_path + graph[start][(start[0], 'L')]['nr_links']))                
#         else:
#             for node in set(graph[start]).difference(path):
#                 if node not in forbidden: # and node not in already_visited: 
#                     # try: # if last node (i.e. "end") it is not present in small_scaffolds and it should not be included in the length
#                     #     heapq.heappush(heap, (nr_links, node, path, path_len + graph[node[0]]['length'])) #  small_scaffolds[node[0]].s_length))   #
#                     # except KeyError:
#                     nr_links = graph[start][node]['nr_links']
#                     heapq.heappush(heap, (nr_links, (node, path)))
#             if len(heap) > max_size_heap:
#                 max_size_heap = len(heap)
#     #print "max lenght heap:", max_size_heap, "paths:", len(paths)
#     return paths

def find_all_paths_for_start_node_DFS_dynamic_programming_ish(graph, start, end, already_visited, is_withing_scaf, max_path_length_allowed, param):
    path = []
    paths = []
    max_size_heap = 0
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
    heap = [(0,(start, path), (set([start]), set(), 0, 0))] #last tuple is: current contigs in path, and how many bad neighbours (skipped contigs hanging from contigs included in path on 3_prime end)
    #prev_node = start
    counter = 0
    head_dict = {} # contains current best paths for 3_prime contigs
    while heap:
        #prev_node = start
        # print head_dict
        counter += 1
        #if counter % 100 == 0:
        #    print 'Potential paths:', counter, 'paths found: ', len(paths)
        if counter > param.path_threshold or len(path) > 100:
            #print 'Hit path_threshold of {0} iterations! consider increase --iter <int> parameter to over {0} if speed of BESST is not a problem. Standard increase is, e.g., 2-10x of current value'.format(param.path_threshold)
            param.hit_path_threshold = True
            break
            
        nr_links, (start, path), (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs, bad_link_count)  = heapq.heappop(heap) #start, end, path, sum_path = heapq.pop()

        if start in head_dict:
            # strictly worse path than previous seen one
            if nr_bad_nbrs > head_dict[start][0] and bad_link_count > head_dict[start][1]:
                continue

        # print path
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
            paths.append(path)
            continue


        if  prev_node[0] != start[0]:
            nr_links = 2**16 # large number to give high priority to this intra-contig edge in this two node representation
            if start[1] == 'L' and (start[0], 'R') not in forbidden:
                additional_bad_nbrs = 0 
                additional_bad_link_count = 0
                for nbr in graph.neighbors(start):
                    if nbr[0] != start[0] and nbr not in ctg_ends_in_path and nbr not in bad_ctgs:
                        additional_bad_nbrs += 1
                        additional_bad_link_count += graph[start][nbr]['nr_links']
                        bad_ctgs.add(nbr)

                if (start[0], 'R') not in head_dict:
                    head_dict[(start[0], 'R')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                    head_dict[(start[0], 'L')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                    # print "Here", (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count), start[0]
                    ctg_ends_in_path.add((start[0], 'R'))
                    heapq.heappush(heap, (nr_links, ((start[0], 'R'), path), (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)) ) 

                elif (start[0], 'R') in head_dict:
                    # print 'Were in!!'
                    # print "OK", (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count), start[0]
                    if nr_bad_nbrs + additional_bad_nbrs > head_dict[(start[0], 'R')][0] and bad_link_count + additional_bad_link_count > head_dict[(start[0], 'R')][1]:
                        #print "Skipping path here"
                        continue
                    else:  
                        # print "FINAL LEVEL"
                        ctg_ends_in_path.add((start[0], 'R'))
                        heapq.heappush(heap, (nr_links, ((start[0], 'R'), path), (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)) ) 


            elif start[1] == 'R' and (start[0], 'L') not in forbidden:
                additional_bad_nbrs = 0 
                additional_bad_link_count = 0
                for nbr in graph.neighbors(start):
                    if nbr[0] != start[0] and nbr not in ctg_ends_in_path and nbr not in bad_ctgs:
                        additional_bad_nbrs += 1 
                        additional_bad_link_count += graph[start][nbr]['nr_links']
                        bad_ctgs.add(nbr)
                        
                if (start[0], 'L') not in head_dict:
                    head_dict[(start[0], 'L')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                    head_dict[(start[0], 'R')] = (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)
                    # print "Here", (nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count), start[0]
                    ctg_ends_in_path.add((start[0], 'L'))
                    heapq.heappush(heap, (nr_links, ((start[0], 'L'), path), (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)) ) 

                elif (start[0], 'L') in head_dict:
                    # print 'Were in!!'
                    if nr_bad_nbrs + additional_bad_nbrs > head_dict[(start[0], 'L')][0] and bad_link_count + additional_bad_link_count > head_dict[(start[0], 'L')][1]:
                        #print "Skipping path here"
                        continue
                    else:
                        # print "FINAL LEVEL"
                        ctg_ends_in_path.add((start[0], 'L'))
                        heapq.heappush(heap, (nr_links, ((start[0], 'L'), path), (ctg_ends_in_path, bad_ctgs, nr_bad_nbrs + additional_bad_nbrs, bad_link_count + additional_bad_link_count)) ) 

        else:
            for node in set(graph[start]).difference(path):
                path_new = path + [node]
                ctg_ends_in_path_new = set(path_new)
                if node not in forbidden: # and node not in already_visited: 
                    nr_links = graph[start][node]['nr_links']
                    heapq.heappush(heap, (nr_links, (node, path), (ctg_ends_in_path_new, bad_ctgs, nr_bad_nbrs, bad_link_count) ))
                    # print 'Added', ctg_ends_in_path_new
            # if len(heap) > max_size_heap:
                # max_size_heap = len(heap)
    # print "max lenght heap:", max_size_heap, "paths:", len(paths)
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

    param.hit_path_threshold = False
    print 'iterating until maximum of {0} extensions.'.format(iter_threshold) 
    print 'Number of nodes:{0}, Number of edges: {1}'.format(len(G_prime.nodes()), len(G_prime.edges()))
    print >> param.information_file,  'iterating until maximum of {0} extensions.'.format(iter_threshold)
    print >> param.information_file, 'Number of nodes:{0}, Number of edges: {1}'.format(len(G_prime.nodes()), len(G_prime.edges()))
    while len(iter_nodes) > 0 and iter_count <= iter_threshold:
        iter_count += 1
        start_node = iter_nodes.pop()
        if cnter % 100 == 0:
            print 'enter Between scaf node:{0}, scaffold progression {1}%. '.format(cnter, round( cnter / float(iter_threshold)*100, 1 ))
        end.difference_update(set([start_node]))
        if param.dfs_traversal:
            paths = find_all_paths_for_start_node_DFS_dynamic_programming_ish(G_prime, start_node, end, already_visited, 0, 2 ** 32, param)
            #paths = find_all_paths_for_start_node_BFS_improved(G_prime, start_node, end, already_visited, 0, 2 ** 32, param)
            #paths = map(lambda x: x[0], paths)
            # print 'NEW PATHS', paths2
            # print
            # print
            # print "OLD PATHS", p
            #print len(paths)
            #print len(p)
            #assert p == paths2
        else:
            paths = find_all_paths_for_start_node_BFS_Dynamic_Programming_ish(G_prime, start_node, end, already_visited, 0, 2 ** 32, param)

        already_visited.add(start_node)
        ScorePaths(G_prime, paths, all_paths, param)
        cnter += 1
    #all_paths = ExtendScaffolds(all_paths)
    #print all_paths
    print 'Total nr of paths found: {0} with score larger than: {1}'.format(len(all_paths), param.score_cutoff)
    all_paths.sort(key=lambda list_: list_[0])
    if param.hit_path_threshold:
        print 'Hit path_threshold of {0} iterations! consider increase --iter <int> parameter to over {0} if speed of BESST is not a problem. Standard increase is, e.g., 2-10x of current value'.format(param.path_threshold)
        print >> param.information_file, 'Hit path_threshold of {0} iterations! consider increase --iter <int> parameter to over {0} if speed of BESST is not a problem. Standard increase is, e.g., 2-10x of current value'.format(param.path_threshold)
    #print all_paths
    return(all_paths)

def WithinScaffolds(G, G_prime, start, end_node, already_visited, max_path_length, param):
    end = set()
    end.add(end_node)
    all_paths = []
    already_visited.difference_update(set([start, end_node]))
    if param.dfs_traversal:
        paths = find_all_paths_for_start_node_DFS_dynamic_programming_ish(G_prime, start, end, already_visited, 1, max_path_length, param)
    else:
        paths = find_all_paths_for_start_node_BFS_Dynamic_Programming_ish(G_prime, start, end, already_visited, 1, max_path_length, param)


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



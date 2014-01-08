'''
    Created on Jun 17, 2012
    
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


import networkx as nx
from itertools import combinations
import random
import multiprocessing , Queue
import time
import ExtendLargeScaffolds as ELS


class Worker(multiprocessing.Process):

    def __init__(self, work_queue, result_queue):

        # base class initialization
        multiprocessing.Process.__init__(self)

        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
    def run(self):
        #while not self.kill_received:

            # get a task
            #job = self.work_queue.get_nowait()
        #try:
        job = self.work_queue #.get_nowait()
        start_nodes = job[0]
        G_prime = job[1]
        #small_scaffolds = job[2]
        end = job[2]
        param = job[3]
            #all_paths_sorted_wrt_score = [(1,[]),(2,[])]

            # except Queue.Empty:
            #    print 'Noeee! empty!!', multiprocessing.current_process().name
            #    break

        print 'Enter ', multiprocessing.current_process().name
            # the actual processing
        all_paths_sorted_wrt_score = ELS.BetweenScaffolds(G_prime, end, start_nodes, param)
        #all_paths_sorted_wrt_score = [(1,[]),(2,[])]
            #nr_paths = find_all_paths(graph, start_node, end_node)
            # store the result
        self.result_queue.put(all_paths_sorted_wrt_score)
        print 'Exit', multiprocessing.current_process().name




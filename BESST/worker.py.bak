'''
    Created on Oct 30, 2012
    
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


import multiprocessing, Queue
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
        while not self.kill_received:

                # get a task
                #job = self.work_queue.get_nowait()
            worker_name = multiprocessing.current_process().name
            try:
                jobs, G, G_prime, small_scaffolds, end = self.work_queue.get_nowait()
            except Queue.Empty:
                break

            # the actual processing
            all_paths_sorted_wrt_score = ELS.BetweenScaffolds(G, G_prime, small_scaffolds, end, jobs)
            self.result_queue.put(all_paths_sorted_wrt_score)
            print 'Exited from ', worker_name

        #return()

#            for job in jobs:
#                start_node = job[0]
#                end_node = job[1]
#                #nr_paths = find_all_paths(graph, start_node, end_node)
#                    #print nr_paths, multiprocessing.current_process().name
#                # store the result
#                self.result_queue.put(nr_paths)
#import os

    #result_queue.put(all_paths_sorted_wrt_score)

    #os.kill()

result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    print 'Going in here'
    result_list.append(result)
def worker(jobs, G, G_prime, small_scaffolds, end):
    worker_name = multiprocessing.current_process().name
    # the actual processing
    all_paths_sorted_wrt_score = ELS.BetweenScaffolds(G, G_prime, small_scaffolds, end, jobs)
    print 'Exited from ', worker_name
    return(all_paths_sorted_wrt_score)

def BetweenScaffoldsParallellized(G, G_prime, small_scaffolds, num_processes, end):
    import multiprocessing , Queue
    #import worker #(multi processing)
    import heapq

    #import random
#    end = set()
#    for node in G:
#        end.add(node)

    # load up work queue
    work_queue = multiprocessing.Queue()
    nodes = G.nodes()
    #random.shuffle(nodes)
    nr_jobs = len(nodes)
    chunk = nr_jobs / num_processes
    #print jobs, nr_jobs, chunk
    counter = 0

    # partition equally many nodes in G to each core
    while counter < nr_jobs:
        work_queue.put((set(nodes[counter:counter + chunk]), G, G_prime, small_scaffolds, end))
        print 'node nr', counter, 'to', counter + chunk - 1, 'added'
        counter += chunk

    # create a queue to pass to workers to store the results
    #result_queue = multiprocessing.Queue()    

    pool = multiprocessing.Pool()
    #num_processes = multiprocessing.cpu_count()
    #for i in range(0,num_processes):
    while not work_queue.empty():
        nodes, G, G_prime, small_scaffolds, end = work_queue.get_nowait()
        pool.apply_async(worker, args=(nodes, G, G_prime, small_scaffolds, end), callback=log_result)
    pool.close()
    pool.join()
    #print(result_list)

#
#    jobs = []
#    for i in range(num_processes):
#        try:
#            nodes, G, G_prime, small_scaffolds, end = work_queue.get_nowait()
#            print 'Time to work'
#        except Queue.Empty:
#            print 'lol no work'
#        p = multiprocessing.Process(target=worker.worker, args=(nodes, G, G_prime, small_scaffolds, end,result_queue))
#        jobs.append(p)
#        print 'job nr ',i, 'start'
#        p.start()

#    
#    # spawn workers (equally many as available cores)
#    workers =[]
#    for i in range(0,num_processes):
#        workers.append(MP.Worker(work_queue, result_queue))
#        workers[i].start()
    #print 'hejhej!'
    work_queue.close() # Indicate that no more data will be put on this queue by the current process
    # Wait for the workers to finish
#    for i in range(0,num_processes):
#        jobs[i].join()
#        print "I'm waiting..."
#        #workers[i].join()

    # collect the results off the queue
    def wrapper(func, args):
        return(func(*args))
    #tot_result = []

    #while not result_queue.empty():
    #    tot_result.append(result_queue.get())
    #print  'tot result', tot_result
    #Merge the scored paths from each of the subprocesses into one consensus path
    itr = wrapper(heapq.merge, result_list) #tot_result)

    return(itr)

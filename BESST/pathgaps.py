'''
    Created on May 30, 2014

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

import copy
import math
import random

from mathstats.normaldist.normal import normpdf


class Contig(object):
    """docstring for Contig"""
    def __init__(self, index, length):
        super(Contig, self).__init__()
        self.index = index
        self.length = length
        self.position = None
        

class Path(object):
    """docstring for Path"""
    def __init__(self, ctg_lengths, observations):
        super(Path, self).__init__()
        self.ctgs = []
        for i,length in enumerate(ctg_lengths):
            self.ctgs.append(Contig(i, length))
        self.ctgs = tuple(self.ctgs)
        self.gaps = [0]*(len(ctg_lengths)-1) # n contigs has n-1 gaps between them, start with gap size 0  
        
        # get positions for when all gaps are 0
        self.update_positions()

        self.observations = observations


    def get_distance(self,start_index,stop_index):
        total_contig_length = sum(map(lambda x: x.length ,filter(lambda x: start_index <= x.index < stop_index, self.ctgs) ))
        total_gap_length = sum(self.gaps[start_index:stop_index])
        index_adjusting = len(self.gaps[start_index:stop_index]) # one extra bp shifted each time
        return (total_contig_length, total_gap_length, index_adjusting)

    def update_positions(self):
        for ctg in self.ctgs:
            index = ctg.index
            ctg.position = sum(self.get_distance(0,index))

    def get_inferred_isizes(self):
        self.isizes = {}
        for (c1,c2) in self.observations:
            gap = self.ctgs[c2].position - (self.ctgs[c1].position + self.ctgs[c1].length) - (c2-c1) # last thing is an index thing
            x = map(lambda obs: obs + gap , self.observations[(c1,c2)]) # inferr isizes
            self.isizes[(c1,c2)] = x

    def propose_new_state(self,mean,stddev):
        path_proposed = copy.deepcopy(self) # create a new state
        #(index,gap) = random.choice([(i,gap) for i,gap in enumerate(self.gaps)]) # choose a gap to change
        (c1,c2) = random.choice(self.observations.keys())
        print 'CHOSEN:', (c1,c2)
        obs = self.observations[(c1,c2)] # take out observations and

        #obs = self.observations[(index,index+1)] # take out observations and
        exp_mean_over_bp = mean + stddev**2/float(mean+1)
        mean_obs = sum(obs)/float(len(obs)) # get the mean
        proposed_distance = exp_mean_over_bp - mean_obs # choose what value to set between c1 and c2 
        (total_contig_length, total_gap_length, index_adjusting) = self.get_distance(c1+1,c2)
        avg_suggested_gap = (proposed_distance - total_contig_length)/ c2-c1
        for index in range(c1,c2):
            path_proposed.gaps[index] = avg_suggested_gap
        #path_proposed.gaps[index] = proposed_distance
        return path_proposed

    def calc_log_likelihood(self,mean,stddev):
        log_likelihood_value = 0
        for (c1,c2) in self.isizes:
            for isize in self.isizes[(c1,c2)]:
                log_likelihood_value += math.log( normpdf(isize,mean,stddev) )
            
        return log_likelihood_value

    def make_path_dict_for_besst(self):
        path_dict = {}
        for ctg1,ctg2 in zip(self.ctgs[:-1],self.ctgs[1:]):
            path_dict[(ctg1,ctg2)] = ctg2.position - (ctg1.position + ctg1.length) - 1 
        print path_dict
        return path_dict


    def __str__(self):
        string= ''
        for ctg in self.ctgs:
            string += 'c'+str(ctg.index) +',startpos:'+ str(ctg.position)+',endpos:'+str(ctg.position + ctg.length)+'\n'

        return string


def position_maximum_likelihood(path, mean,stddev):
    
    path.get_inferred_isizes()
    iteration = 0 
    curr_state_count = 1 # how many iterations we have remained in state
    while True:
        #print 'CURRENT PATH:'
        #print path
        #print 'curr highest likelihood: ', path.calc_log_likelihood(mean,stddev)

        # change gaps
            # MCMC based sampling:
            # take the most deviant mean insert size observation (x_hat = o_hat + current_gap)
            # given the current positioning. 
            # Adjust the current_gap to new_gap so that o_hat + new_gap = mean

        suggested_path = path.propose_new_state(mean,stddev)
        suggested_path.update_positions()
        suggested_path.get_inferred_isizes()
        suggested_path.calc_log_likelihood(mean,stddev)

        #print 'SUGGESTED PATH:'
        #print suggested_path
        #print 'Suggested path likelihood: ', suggested_path.calc_log_likelihood(mean,stddev)
        # compare likelihood values ML value for path
        if suggested_path.calc_log_likelihood(mean,stddev) > path.calc_log_likelihood(mean,stddev):
            path = suggested_path
            curr_state_count = 1
        else:
            print 'PATH not taken!'
            curr_state_count += 1

            


        # see if "converged" or not
        if curr_state_count >= 5:
            break

        
        iteration += 1
        if iteration >= 5:
            break
        

    return path
        



def main(contig_lenghts, observations, mean, stddev):
    """
    contig_lenghts: Ordered list of integers which is contig lengths (ordered as contigs comes in the path)
    observations:  dictionary oflist of observations, eg for contigs c1,c2,c3
                    we can have [(c1,c2):[23, 33, 21],(c1,c3):[12,14,11],(c2,c3):[11,34,32]]
    """

    path = Path(contig_lenghts,observations)

    ML_path = position_maximum_likelihood(path,mean,stddev)
    return ML_path


if __name__ == '__main__':
    contig_lenghts = [3000,500,500,500,3000]
    observations = {(0,1):[1800,2000], (0,2):[1500,1800,1400,1700], (0,3):[1200,800,1000], 
                    (1,2):[750,800],(1,3):[600,700], (1,4):[300,600,700],
                    (2,3):[700,750], (2,4):[1400,1570],
                    (3,4):[2000,1750], }
    mean = 1500
    stddev = 500
    main(contig_lenghts,observations,mean,stddev)



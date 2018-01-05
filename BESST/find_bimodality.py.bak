'''
Created on Sep 29, 2011

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
import sys

import numpy as np
import math



def params(x):
    n = float(len(x)) 
    mu = sum(x)/ float(len(x))
    x_sq = sum(map(lambda x: x ** 2, x))
    var = (x_sq - n * mu ** 2) / (n - 1)
    return mu, var


def checkEqualIvo(lst):
    return not lst or lst.count(lst[0]) == len(lst)

def split_distribution(all_isizes):
    """ Split a population into sub-populations

        Based on binning the data into n_bins and finding contigous groups of non-empty bins.

        Returns [lowest, ..., highest] all of which are sorted sequences
        
        EXAMPLE: 
        all_isizes = [2, 1, 2, 1, 2, 3, 3, 3, 2, 1, 2,
                     3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 11,
                    12, 12, 11, 13, 11, 12, 13, 11, 12, 
                    13, 1, 2, 3, 2, 1, 3]

        find_bimodality.split_population2(all_isizes,20)

        bins = [  1.    1.6   2.2   2.8   3.4   4.    4.6   5.2   5.8   6.4   7.    7.6
                8.2   8.8   9.4  10.   10.6  11.2  11.8  12.4  13. ]

        # n bins (set by n_bins parameter) will always create n+2 sublists as seen below
        binned = [[], [1, 1, 1, 1, 1], [2, 2, 2, 2, 2, 2, 2], [], [3, 3, 3, 3, 3, 3],
                 [], [4], [5], [], [6, 6], [], [7], [8], [], [9], [], [10],
                [11, 11, 11, 11, 11], [], [12, 12, 12, 12], [], [13, 13, 13]]

        # we join them to clusters
        joined_bins = [[1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2], [3, 3, 3, 3, 3, 3],
                        [4, 5], [6, 6], [7, 8], [9], [10, 11, 11, 11, 11, 11],
                        [12, 12, 12, 12], [13, 13, 13]]

        # joined_bins contains 9 sublists. Now we split joined_bins into two sublists
        # by splitting joined_bins at each position 1,...,8. The variance of the two sublists
        # are thereafter calculated and summed up. We return the two sublists (i.e., clusters/populations)
        # that gives the lowest sum. For this example we return:

        cluster1 = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 5, 6, 6, 7, 8]
        cluster2 = [9, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13]
        mean1 = 3.0416666666666665
        stddev1 = 1.9886453039512946
        mean2 = 11.5
        stddev2 = 1.1602387022306428

        Values will slightly change with value on n_bins. It is recommended to choose n_bins
        too high rather than too low. For isize data, choose e.g., n_bins = max_isize

        ###########################################################################
        Some other tests

        We don't want to split l5 (approx normal dist) and l2,l4 since theu have one mode and fat tailed.
        l1 and l3 have two modes.
        >>> l
        [2, 1, 2, 1, 2, 3, 3, 3, 2, 1, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 11, 12, 12, 11, 13, 11, 12, 13, 11, 12, 13, 1, 2, 3, 2, 1, 3]
        >>> l2
        [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 11]
        >>> l3
        [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 7, 7, 7, 8, 8, 9, 9, 10, 10, 11]
        >>> l4
        [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 10, 11]
        >>> l5
        [1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9]
        >>> find_bimodality.split_distribution(l)
        ([1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 5, 6, 6, 7], [8, 9, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13], 2.8260869565217392, 1.7228985595908761, 11.266666666666667, 1.4375905768565118)
        >>> find_bimodality.split_distribution(l2)
        ([1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 9, 10, 10, 11], [], 4.621621621621622, 2.8901302484235636, 0, 0)
        >>> find_bimodality.split_distribution(l3)
        ([1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3], [7, 7, 7, 8, 8, 9, 9, 10, 10, 11], 2.1176470588235294, 0.781213234429025, 8.6, 1.4298407059684846)
        >>> find_bimodality.split_distribution(l4)
        ([1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 10, 11], [], 4.264705882352941, 2.711466376496939, 0, 0)
        >>> find_bimodality.split_distribution(l5)
        ([1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 9], [], 5.0, 1.8754309849499513, 0, 0)
    """

    sorted_isizes = sorted(all_isizes)
    n_bins = sorted_isizes[-1]*2
    # bin the data into n_bins in a 2d structure, one sequence for each bin:
    _ , bins = np.histogram(sorted_isizes, bins=n_bins)
    # print bins
    bin_indices = np.digitize(sorted_isizes, bins)
    # print bin_indices
    binned = []
    for i in range(len(bins)+1):
        binned.append([])
    for ix_bin, v in zip(bin_indices, sorted_isizes):
        binned[ix_bin].append(v)
    # print binned

    # now join-up non-empty bins
    joined_bins = [[]] # so 2D, with initially 1 sub-list
    len_last_bin = 0
    for bin in binned:
        len_bin = len(bin)
        if len_bin == 0 and len_last_bin != 0: # will correctly handle the case where bin 0 is empty
            joined_bins.append([])
        if len_bin != 0:
            joined_bins[-1].extend(bin)
        len_last_bin = len_bin


    # special case: only one bin was detected, probably due to a tight cluster
    # i.e. no split distributions
    if len(joined_bins) < 2:
        observations = [item for sublist in joined_bins for item in sublist]
        mu, var = params(observations)
        stddev = math.sqrt(var)
        mean = mu
        return cluster1, [], mean1, stddev1, 0, 0



    # if we are here we passed the special case of only one cluster
    # we calculate the mean and stddev for the base case now.
    # the base case is that there is in fact only one cluster of isizes
    # and there is no need to separete observations. We will later on compare this
    # base case varaince to variances when splitting the distribution.

    observations = [item for sublist in joined_bins for item in sublist]
    filter_observations = filter(lambda x: x>100, observations)
    if len(filter_observations) < 2:
        return [],[],0,0,0,0
    base_mu, base_var = params(filter_observations)
    base_stddev = math.sqrt(base_var)
    mean = base_mu
    # print "base", base_mu, base_stddev, base_var
    
    #calculate sum of two varinaces from the two sublists when splitting
    # this list of list into two. We ar guaranteed to calculate the varance for 
    # at most n+1 pairs of lists. In practice often fewer times.
    # print 'JOINED BINS:', joined_bins


    # initialze lowest var to the variance given that there is only one cluster of isizes
    lowest_var = base_var
    lowest_stddev = base_stddev
    # print lowest_var, lowest_stddev
    split_index = 0

    for i in range(1,len(joined_bins)):
        observations1 = [item for sublist in joined_bins[:i] for item in sublist]
        observations2 = [item for sublist in joined_bins[i:] for item in sublist]
        # print observations1
        # print observations2
        if len(observations1) < 2 or len(observations2) < 2 or checkEqualIvo(observations1) or checkEqualIvo(observations2):
            continue
        mu1, var1 = params(observations1)
        mu2, var2 = params(observations2)
        # print 'mu1', mu1
        # print 'mu2', mu2

        sum_var = var1 + var2
        sum_stddev = math.sqrt(var1) + math.sqrt(var2)
        # print sum_var, sum_stddev
        if sum_var < lowest_var and sum_stddev < lowest_stddev:
            lowest_var = sum_var
            lowest_stddev = sum_stddev
            stddev1 = math.sqrt(var1)
            stddev2 = math.sqrt(var2)
            mean1 = mu1
            mean2 = mu2
            split_index = i

    if lowest_var < base_var and lowest_stddev < base_stddev:
        cluster1 = [item for sublist in joined_bins[:split_index] for item in sublist]
        cluster2 = [item for sublist in joined_bins[split_index:] for item in sublist]
        return cluster1, cluster2, mean1, stddev1, mean2, stddev2

    else:
        cluster1 = [item for sublist in joined_bins for item in sublist]
        base_stddev = math.sqrt(base_var)
        return cluster1, [], base_mu, base_stddev, 0, 0

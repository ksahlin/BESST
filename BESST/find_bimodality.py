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

def split_list(a_list, position):
    # """
    #     loc: a rational number between 0 and 1 
    #     (preferrably specified as a multiple of1/len(a_list) parts)
    #     to ensure split on a specific item. where the list should be split.
    #     0 is the first item and 1 is the last. 
    # """
    # if not 0<= loc <= 1:
    #     sys.error("Invalid range on loc. Exiting.")
    #     sys.exit()

    # position = len(a_list)*loc
    return a_list[:position], a_list[position:]



def params(x):
    n = float(len(x)) 
    mu = sum(x)/ float(len(x))
    x_sq = map(lambda x: x ** 2, x)
    var = (x_sq - n * mu ** 2) / (n - 1)
    return mu,var

def find_bimodal_fit(sorted_observations, separation_point=len(observations)/2):
    """
        Simpler algorithm to find roughly the two normal-like distributions
        in a bimodal distribution. We don't care about exact fit, just a rough separation
        value k, that is a x-coordinate separating the two distributions.
        This only works if the two distributions are far from eachother
    """


    # # initializations
    # observations.sort()
    # separation_coordinate = observations[len(observations)/2]
    # x1, x2 = split_list(observations)
    # mu1, var1 = params(x1)
    # mu2, var2 = params(x2)

    # sum_var = var1 + var2

    if sum_var < get_sum_var(sorted_observations, separation_point)
    return min(sum_var ) observations[separation_point], )

def get_sum_var(sorted_observations, separation_point):
    x1, x2 = split_list(sorted_observations, separation_point)
    mu1, var1 = params(x1)
    mu2, var2 = params(x2)
    sum_var = var1 + var2
    return sum_var


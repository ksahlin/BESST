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


class scaffold(object):
    '''
    classdocs
    '''

    __slots__ = ('name', 'contigs', 's_length', 'lower_left_nbrs_obs',
                 'lower_right_nbrs_obs', 'upper_left_nbrs_obs', 'upper_right_nbrs_obs')
    def __init__(self, scaffold_name, scaffold_contigs, scaffold_length,
                 scaffold_lower_left_nbrs, scaffold_lower_right_nbrs,
                 scaffold_upper_left_nbrs, scaffold_upper_right_nbrs):
        '''
        Constructor
        '''
        self.name = scaffold_name           # String name
        self.contigs = scaffold_contigs     # list of contig objects that are ordered ass they should be placed in the scaffold
        self.s_length = scaffold_length     # integer of total length of the scaffold   
        self.lower_left_nbrs_obs = scaffold_lower_left_nbrs
        self.lower_right_nbrs_obs = scaffold_lower_right_nbrs
        self.upper_left_nbrs_obs = scaffold_upper_left_nbrs
        self.upper_right_nbrs_obs = scaffold_upper_right_nbrs



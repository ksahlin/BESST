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

    __slots__ = ('name', 'contigs', 's_length')
    def __init__(self, scaffold_name, scaffold_contigs, scaffold_length):
        '''
        Constructor
        '''
        self.name = scaffold_name           # String name
        self.contigs = scaffold_contigs     # list of contig objects that are ordered as they should be placed in the scaffold
        self.s_length = scaffold_length     # integer of total length of the scaffold   




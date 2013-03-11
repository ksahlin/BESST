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


class contig(object):
    __slots__ = ('name', 'scaffold', 'direction', 'position', 'length',
                  'coverage', 'repeat', 'is_haplotype', 'sequence')

    def __init__(self, contig_name, contig_scaffold=None , contig_direction=None,
                  contig_position=None, contig_length=None, contig_coverage=None,
                  contig_repeat=False, contig_haplotype=False, contig_sequence=None):
        self.name = contig_name
        self.scaffold = contig_scaffold
        self.direction = contig_direction
        self.position = contig_position
        self.length = contig_length
        self.sequence = contig_sequence
        self.coverage = contig_coverage
        self.repeat = contig_repeat
        self.is_haplotype = contig_haplotype

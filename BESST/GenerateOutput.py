'''
    Created on Sep 29, 2011
    @author: ksahlin

    Updates on march, 2017
    Updates on august, 2018
    @contributor: sletort

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
from __future__ import print_function
import os
import time

rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}


def longest_kmer_overlap(s1, s2):
    i = len(s1)
    while i > 0:
        if s1[-i:] == s2[:i]:
            return i
        i -= 1
    return i

def PrintOutHaplotypes(Haplotypes, Contigs, output_dest):
    haplotype_file = open(output_dest + '/haplotypes.fa', 'a')
    for cont_hapl in Haplotypes:
        print('>' + Haplotypes[cont_hapl][0] + ' (variant of ' + cont_hapl + ')', file=haplotype_file)
        hapl_contig = Contigs[Haplotypes[cont_hapl][0]].sequence
        for i in range(0, len(hapl_contig), 60):
            print(hapl_contig[i:i + 60], file=haplotype_file)
        #remove the output haplotype, the other one is left in graph
        del Contigs[Haplotypes[cont_hapl][0]]
    return(Contigs)

def PrintOutRepeats(Repeats, Contigs, output_dest, small_contigs):
    repeat_file = open(output_dest + '/repeats.fa', 'w')
    for cont_obj in Repeats:
        print('>' + cont_obj.name, file=repeat_file)
        repeat_contig = cont_obj.sequence
        for i in range(0, len(repeat_contig), 60):
            print(repeat_contig[i:i + 60], file=repeat_file)
        try:
            del Contigs[cont_obj.name]
        except KeyError:
            del small_contigs[cont_obj.name]
    return()

def repeat_contigs_logger(Repeats, Contigs, output_dest, small_contigs, param):
    repeat_logger_file = open(output_dest + '/repeats_log.tsv', 'w')
    print("contig_accession\tlength\tcoverage\tcov/mean_cov(exp number of placements)\tlib_mean\tplacable", file=repeat_logger_file)
    repeats_sorted = sorted(Repeats, key=lambda x: x.coverage, reverse=True)
    for cont_obj in repeats_sorted:
        placable = "Yes" if param.mean_ins_size > cont_obj.length else 'No'
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(cont_obj.name, cont_obj.length, round(cont_obj.coverage,1), round(cont_obj.coverage/param.mean_coverage,0), round(param.mean_ins_size,0), placable), file=repeat_logger_file)

def PrintOut_low_cowerage_contigs(low_coverage_contigs, Contigs, output_dest, small_contigs):
    low_coverage_contigs_file = open(output_dest + '/low_coverage_contigs.fa', 'w')
    for cont_obj in low_coverage_contigs:
        print('>' + cont_obj.name, file=low_coverage_contigs_file)
        low_coverage_contig = cont_obj.sequence
        for i in range(0, len(low_coverage_contig), 60):
            print(low_coverage_contig[i:i + 60], file=low_coverage_contigs_file)
        try:
            del Contigs[cont_obj.name]
        except KeyError:
            del small_contigs[cont_obj.name]
    return()

def ChangeToSmallContigs(Contigs, list_of_contigs, small_contigs):
    #change from Contigs to small_contige (contigs are included in scaffolds that does not meet the lenght criteria, they will be used in path extension algorithm)
    for cont_obj in list_of_contigs:
        del Contigs[cont_obj.name]
        small_contigs[cont_obj.name] = cont_obj
    return()

def WriteToF(F, Contigs, list_of_contigs):
    info_list = []
    for cont_obj in list_of_contigs:
        info_list.append((cont_obj.name, cont_obj.direction, cont_obj.position, cont_obj.length, cont_obj.sequence)) #,cont_obj.links
        if cont_obj.position < 0:
            print('Write to F: Position is negative!', cont_obj.position, cont_obj.name, cont_obj.direction)
        #del Contigs[cont_obj.name]
    F.append(info_list)
    return(F)

# ==================================
def compute_kmer_overlap( end1,end2 ):
	# should be a class method of Scaffold ?
    i = len(end1)
    while i > 0:
        if end1[-i:] == end2[:i]:
            return i
        i -= 1
    return i

# thoses classes are based on AGP concept
# Scaffold is 'Object' in AGP concept
class Component(object):
    """'Component' is for what compose the 'Object', 'Contig' and 'Gap'."""
    i = 0 # an instance counter, a way to have id for component.

    def __init__( self, id, len, seq ):
        super( Component, self ).__init__()
        Component.i += 1

        self._len = len
        self._id  = id if( id != None ) else Component.i
        self._seq = seq

        # position of the first/last contig base that belong to the scaffold
        # contig coordinate
        self._scaff_start = 1
        self._scaff_end   = len

    def __del__( self ):
        Component.i -= 1

    @property
    def seq( self ):
        return self._seq

    @property
    def len( self ):
        return self._len

    @property
    def id( self ):
        return self._id

    @property
    def scaff_start( self ):
        return self._scaff_start
    @scaff_start.setter
    def scaff_start( self, value ):
        self._scaff_start = value
        #~ self._scaff_end  -= value

    @property
    def scaff_end( self ):
        return self._scaff_end

    def get_scaff_ctg_coords( self ):
        return [ self.scaff_start, self.scaff_end ]

    def __str__( self ):
        s  = "\t".join([ str(x) for x in [self._id, self._len] ])
        s += "\t".join([ str(x) for x in [ '\ncoords' ]+self.get_scaff_ctg_coords() ])
        s += "\n" + self._seq
        return s

class Contig( Component ):
    '''a Contig is a Component with a direction(=strand).

    This class provide a builder and sequence is returned regarding direction.
    '''
    def __init__( self, id, len, direction, seq ):
        super( Contig, self ).__init__( id, len, seq )
        self.__direction = direction

    @classmethod
    def buildFromTuple( cls, t_ctg ):
        return cls( t_ctg[0], t_ctg[3], t_ctg[1], t_ctg[4] )

    @property
    def seq( self ):
        if self.direction:
            return self._seq
        else:
            return RevComp( self._seq,rev_nuc )

    @property
    def direction( self ):
        return self.__direction

class Gap( Component ):
    '''a Gap is a component whose sequence is a repetition of N or n.'''

    def __init__( self, char, len ):
        super( Gap, self ).__init__( None, len, char*len )

# ==================================

class Scaffold(object):
    """docstring for Scaffold"""
    def __init__(self, name, param, info_tuple):
        super( Scaffold, self ).__init__()
        self.name  = name
        self.param = param

        self.l_components = []
        self.__scaff_len  = 0
        self.__buildComponents( info_tuple )

    @property
    def scaff_len( self ):
        return self.__scaff_len

    def __incorporateComponent( self, o_compo, compo_pos=0 ):
        ctg_offset = 0

        if( 0 == len( self.l_components ) ):
            # the first contig align all along the scaffold
            s_start = compo_pos + 1  # if the first contig start on base 1 of the scaffold, use 1 make it easier to understand
                # compo_pos + 1 to turn to 1-based.
        else:
            s_start = self.scaff_len + 1 # the base following the last one

        s_end = s_start + o_compo.len - o_compo.scaff_start
        self.__scaff_len = s_end
        self.l_components.append([ o_compo, s_start, s_end ])

    def __buildComponents( self, l_info_tuples ):
        o_ctg = Contig.buildFromTuple( l_info_tuples[0] )

        self.__incorporateComponent( o_ctg, l_info_tuples[0][2] )

        max_gap_size = 2*self.param.std_dev_ins_size # not good var name, max_small_gap_size ?
        max_overlap  = self.param.max_contig_overlap

        for i in range( 1, len( l_info_tuples ) ):
            o_next   = Contig.buildFromTuple( l_info_tuples[i] )
            next_pos = l_info_tuples[i][2] +1 # 1-based

            gap_estimate = next_pos - ( l_info_tuples[i-1][2]+1 + o_ctg.len )
            if( max_gap_size < gap_estimate ):
                #~ print "***Big Gap"
                o_gap = Gap( 'N', gap_estimate )
            else:
                # compute the real overlap
                prev_seq = o_ctg.seq[-max_overlap:]
                next_seq = o_next.seq[:max_overlap]
                overlap  = compute_kmer_overlap( prev_seq, next_seq )

                #~ print "Overlap = " + str( overlap )
                if( 20 <= overlap ): # big overlap
                    #~ print "***Big Overlap"
                    print( 'merging {0} bp here'.format(overlap), file=self.param.information_file )
                    o_gap = Gap( 'n', 1 )
                    o_next.scaff_start = overlap + 1 # +1 to turn 1-based
                elif( gap_estimate <= 1 ): # supposed overlap, but not found
                    #~ print "***Supposed overlap"
                    o_gap = Gap( 'n', 1 )
                else:
                    o_gap = Gap( 'N', gap_estimate )

            self.__incorporateComponent( o_gap )
            self.__incorporateComponent( o_next )
            o_ctg = o_next

    def make_fasta_string( self,fasta_file ):
        '''write a fasta sequence from the object.

        Builds a header and concatenates Component sequences.
        When no overlapping, offset will be 0.
        '''
        l_fasta = [ '>'+str(self.name)+'\n' ]
        for o_compo in [ l_[0] for l_ in self.l_components ]:
                offset = o_compo.scaff_start - 1 # scaff_start is 1-based
                l_fasta.append( o_compo.seq[offset:] )

        print( ''.join( l_fasta ), file=fasta_file )

    def make_AGP_string( self, AGP_file ):
        '''write AGP rows from the object.

        AGP row starts with scaff_id, scaff_start, scaff_end.
        Follows the component number, type.
        Then depending on type we have
            for type=contig id, scaff_start, scaff_end, orientation
            for type=gap    len, gap_type, linkage, links_evidence
        '''
        component_count = 0

        for l_compo in self.l_components:
            component_count += 1
            o_compo = l_compo[0]
            l_elts = [ self.name, l_compo[1], l_compo[2], component_count ]

            if( isinstance( o_compo, Gap ) ):
                l_component  = [ 'N', o_compo.len, 'scaffold', 'yes', 'paired-ends' ]
            else:
                strand = '+' if o_compo.direction else '-'
                l_component  = [ 'W', o_compo.id ]
                l_component += o_compo.get_scaff_ctg_coords()
                l_component.append( strand )

            print( '\t'.join([ str(x) for x in l_elts + l_component ]), file=AGP_file )

    def make_GFF_string( self, gff_file ):
        '''write GFF rows from the object.

        GFF row starts with scaff_id, source, component type, scaff_start, scaff_end.
        Follows score, strand, phase and ends with attributes.
        Note: we only provide ID and Name as attributes.
        '''
        source = 'besst_assembly'

        for l_compo in self.l_components:
            o_compo = l_compo[0]

            ( score, strand, phase ) = ( '.', '.', '.' )
            if( isinstance( o_compo, Gap ) ):
                type_   = 'gap'
                l_attrs = []
            else:
                type_   = 'contig'
                strand  = '+' if o_compo.direction else '-'
                name    = "_".join( o_compo.id.split('_',2)[:2] ) # only kept NODE_XXXX from ID
                l_attrs = [ 'ID='+ o_compo.id, 'Name='+name ]

            l_elts  = [ self.name, source, type_, l_compo[1], l_compo[2] ]
            l_elts += [ score, strand, phase, ';'.join( l_attrs ) ]

            print( '\t'.join([ str(x) for x in l_elts ]), file=gff_file )


def PrintOutput(F, param, pass_nr):
    try:
        os.mkdir(param.output_directory + '/pass' + str(pass_nr))
    except OSError:
        #directory is already created
        pass
    #contigs_before=len(C_dict)
    contigs_after = len(F)
    print('(super)Contigs after scaffolding: ' + str(contigs_after) + '\n', file=param.information_file)
    gff_file = open(param.output_directory + '/pass' + str(pass_nr) + '/info-pass' + str(pass_nr) + '.gff', 'w')
    AGP_file = open(param.output_directory + '/pass' + str(pass_nr) + '/info-pass' + str(pass_nr) + '.agp', 'w')
    print('##gff-version 3', file=gff_file)
    print('##agp-version 2.0\n#lw-scaffolder output', file=AGP_file)
    fasta_file = open(param.output_directory + '/pass' + str(pass_nr) + '/Scaffolds-pass' + str(pass_nr) + '.fa', 'w')
    header_index = 0

    unique_id = str(int(time.time()))

    for scaf_ in reversed(F):
        #sort contigs in scaf w.r.t position here
        scaf = sorted(scaf_, key=lambda tuple: tuple[2])
        header_index += 1


        s = Scaffold('scaffold_' + str(header_index) + "_uid_" + unique_id , param, scaf)
        s.make_fasta_string( fasta_file )
        s.make_AGP_string( AGP_file )
        s.make_GFF_string( gff_file )

    return()


def RevComp(string, rev_nuc):
    #rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

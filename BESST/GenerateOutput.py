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

class Scaffold(object):
    """docstring for Scaffold"""
    def __init__(self, name, param, info_tuple):
        super(Scaffold, self).__init__()
        self.name = name
        self.param = param

        self.seqs = [x[4] for x in info_tuple]
        self.gaps = list(map(lambda x,y : y[2] - (x[2] + x[3]), info_tuple[:-1],info_tuple[1:]))
        self.directions = [x[1] for x in info_tuple]
        self.positions = [(x[2],x[2] + x[3] - 1) for x in info_tuple]
        self.contigs = [x[0] for x in info_tuple]

    def check_kmer_overlap(self,end1,end2):
        i = len(end1)
        while i > 0:
            if end1[-i:] == end2[:i]:
                return i
            i -= 1
        return i

    def get_sequence(self, string, direction):
        if direction:
            return string
        else:
            return RevComp(string,rev_nuc)


    def make_fasta_string(self,fasta_file):
        fasta = []
        #fasta.append('>{0}\n'.format(self.name))
        fasta.append('>'+str(self.name)+'\n')
        # first contig

        fasta.append( self.get_sequence(self.seqs[0], self.directions[0]))

        for i in range(len(self.seqs)-1):
            gap = self.gaps[i]
            if gap <= 2*self.param.std_dev_ins_size:
                overlap = self.check_kmer_overlap( self.get_sequence(self.seqs[i], self.directions[i])[-self.param.max_contig_overlap:], self.get_sequence(self.seqs[i+1], self.directions[i+1])[:self.param.max_contig_overlap])
                if overlap >= 20:
                    fasta.append('n' + self.get_sequence(self.seqs[i+1], self.directions[i+1])[overlap:])
                    print('merging {0} bp here'.format(overlap), file=self.param.information_file)
                else:
                    #print gap
                    if gap <= 1:
                        fasta.append('n' + self.get_sequence(self.seqs[i+1], self.directions[i+1]))
                    else:
                        fasta.append('N'*int(gap) + self.get_sequence(self.seqs[i+1], self.directions[i+1]))
            else:
                if gap <= 1:
                    fasta.append('n' + self.get_sequence(self.seqs[i+1], self.directions[i+1]))
                else:
                    fasta.append('N'*int(gap) + self.get_sequence(self.seqs[i+1], self.directions[i+1]))


        print(''.join([ x for x in fasta]), file=fasta_file)

    def make_AGP_string(self, AGP_file):
        component_count = 0
        for i in range( len(self.seqs) ):
            sign = '+' if self.directions[i] else '-'
            if i > 0 and self.gaps[i-1] > 0:
                component_count += 1
                obj_start = self.positions[i-1][1] + 2
                obj_end   = self.positions[i][0]
                compo_len = self.gaps[i-1]
                l_elts = [ self.name, obj_start, obj_end, component_count ]
                l_elts +=[ 'N', compo_len, 'scaffold', 'yes', 'paired-ends' ]
                print('\t'.join([ str(x) for x in l_elts ]), file=AGP_file)

            component_count += 1
            obj_start = self.positions[i][0] + 1
            obj_end   = self.positions[i][1] + 1
            compo_len = self.positions[i][1] - self.positions[i][0] + 1
            l_elts = [ self.name, obj_start, obj_end, component_count ]
            l_elts +=[ 'W', self.contigs[i],'1', compo_len, sign ]
            print('\t'.join([ str(x) for x in l_elts ]), file=AGP_file)

    def make_GFF_string( self, gff_file ):
        source = 'besst_assembly'
        for i in range( len(self.seqs) ):
            sign = '+' if self.directions[i] else '-'
            if i > 0 and self.gaps[i-1] > 0:
                obj_start = self.positions[i-1][1] + 2
                obj_end   = self.positions[i][0]
                l_attrs   = []
                l_elts  = [ self.name, source, 'gap', obj_start, obj_end ]
                l_elts += [ '.', '.', '.', ';'.join( l_attrs ) ]
                print('\t'.join([ str(x) for x in l_elts ]), file=gff_file)

            obj_start = self.positions[i][0] + 1
            obj_end   = self.positions[i][1] + 1
            name      = "_".join( self.contigs[i].split('_',2)[:2] ) # only kept NODE_XXXX from ID
            l_attrs   = [ 'ID='+self.contigs[i], 'Name='+name ]
            l_elts  = [ self.name, source, 'contig', obj_start, obj_end ]
            l_elts += [ '.', sign, '.', ';'.join( l_attrs ) ]
            print('\t'.join([ str(x) for x in l_elts ]), file=gff_file)

def PrintOutput(F, Information, output_dest, param, pass_nr):
    try:
        os.mkdir(param.output_directory + '/pass' + str(pass_nr))
    except OSError:
        #directory is already created
        pass
    #contigs_before=len(C_dict)
    contigs_after = len(F)
    print('(super)Contigs after scaffolding: ' + str(contigs_after) + '\n', file=Information)
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
        s.make_fasta_string(fasta_file)
        s.make_AGP_string(AGP_file)
        s.make_GFF_string( gff_file )

    return()


def RevComp(string, rev_nuc):
    #rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)



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


def PrintOutHaplotypes(Haplotypes, Contigs, output_dest):
    haplotype_file = open(output_dest + '/haplotypes.fa', 'a')
    for cont_hapl in Haplotypes:
        print >> haplotype_file, '>' + Haplotypes[cont_hapl][0] + ' (variant of ' + cont_hapl + ')'
        hapl_contig = Contigs[Haplotypes[cont_hapl][0]].sequence
        for i in range(0, len(hapl_contig), 60):
            print >> haplotype_file, hapl_contig[i:i + 60]
        #remove the output haplotype, the other one is left in graph
        del Contigs[Haplotypes[cont_hapl][0]]
    return(Contigs)

def PrintOutRepeats(Repeats, Contigs, output_dest, small_contigs):
    repeat_file = open(output_dest + '/repeats.fa', 'w')
    for cont_obj in Repeats:
        print >> repeat_file, '>' + cont_obj.name
        repeat_contig = cont_obj.sequence
        for i in range(0, len(repeat_contig), 60):
            print >> repeat_file, repeat_contig[i:i + 60]
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
            print 'Write to F: Position is negative!', cont_obj.position, cont_obj.name, cont_obj.direction
        #del Contigs[cont_obj.name]        
    F.append(info_list)
    return(F)

def PrintOutput(F, Information, output_dest, param, pass_nr):
    import os
    try:
        os.mkdir(param.output_directory + '/pass' + str(pass_nr))
    except OSError:
        #directory is already created
        pass
    #contigs_before=len(C_dict)
    contigs_after = len(F)
    ## this is the madness reverse complement table that needs to be specified when working with AbySS...
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    out_scaf_file = open(param.output_directory + '/pass' + str(pass_nr) + '/scaffold_info-pass' + str(pass_nr), 'w')
    #print >>Information, 'Contigs before scaffolding: '+ str(contigs_before)+'\n'
    print >> Information, '(super)Contigs after scaffolding: ' + str(contigs_after) + '\n'
    gff_file = open(param.output_directory + '/pass' + str(pass_nr) + '/info-pass' + str(pass_nr) + '.gff', 'w')
    AGP_file = open(param.output_directory + '/pass' + str(pass_nr) + '/info-pass' + str(pass_nr) + '.agp', 'w')
    print >> gff_file, '#gff-version 3'
    print >> AGP_file, '#APG file\n#lw-scaffolder output'
    output = open(param.output_directory + '/pass' + str(pass_nr) + '/Scaffolds-pass' + str(pass_nr) + '.fa', 'w')
    header_index = 0
    prev_pos = 0
    print_list = []
    print >> out_scaf_file, 'Scaffold name/length/nr contigs:\n '
    component_count = 0
    for scaf_ in reversed(F):
        #sort contigs in scaf w.r.t position here
        scaf = sorted(scaf_, key=lambda tuple: tuple[2])
        header_index += 1
        print_list = []
        print >> output, '>scaffold_' + str(header_index)
        #print_list.append('>scaffold_'+str(header_index)+'\n')
        scaf_len = len(scaf)
        prev_pos = 0
        component_count = 0
        for i in range (0, scaf_len):
            name = scaf[i][0]
            direction = scaf[i][1]
            pos = scaf[i][2]
            length = scaf[i][3]
            sequence = scaf[i][4]

#            #Get the number of left and right links            
#            if i-1>=0:
#                left_name=scaf[i-1][0]
#                try:
#                    nr_links_left=scaf[i][4][left_name]
#                except KeyError:
#                    nr_links_left='!'                    
#            else:
#                nr_links_left=0
#                
#            if i+1<=scaf_len-1:
#                right_name=scaf[i+1][0]
#                try:
#                    nr_links_right=scaf[i][4][right_name]
#                except KeyError:
#                    nr_links_right='!'                    
#            else:
#                nr_links_right=0

            gap = pos - prev_pos
            N_seq = 'N' * gap
            if gap < 0: #remove a part of the last contig so we can append the next one (overwrite it's sequence)
                last_contig = print_list.pop()
                last_contig = last_contig[0:gap]
                print_list.append(last_contig)

            print_list.append(N_seq)
            if direction:
                print_list.append(sequence)
            else:  #needs to be reverse complemented before outputed
                rev_comp = RevComp(sequence, rev_nuc)
                print_list.append(rev_comp)

            sign = '+' if direction else '-'
            #FOR GFF file
            if scaf_len > 1 and  i > 0 and gap > 0:
                print >> gff_file, 'scaffold_' + str(header_index) + '\tBESST V0.7 \tfragment\t' + str(prev_pos + 1) + '\t' + str(pos) + '\t.\t' + sign + '\t.\t' + 'Parent=scaffold_' + str(header_index)
            #FOR AGP file            
            if i > 0 and gap > 0:
                component_count += 1
                print >> AGP_file, 'scaffold_' + str(header_index) + '\t' + str(prev_pos + 1) + '\t' + str(pos) + '\t' + str(component_count) + '\t' + 'N\t' + str(gap) + '\tfragment\tyes\t'
            component_count += 1
            prev_pos = pos + length
            print >> AGP_file, 'scaffold_' + str(header_index) + '\t' + str(pos + 1) + '\t' + str(prev_pos) + '\t' + str(component_count) + '\t' + 'W\t' + name + '\t1\t' + str(length) + '\t' + sign

            #FOR GFF file
            if scaf_len > 1:
                print >> gff_file, 'scaffold_' + str(header_index) + '\tBESST V0.7\tcontig\t' + str(pos + 1) + '\t' + str(prev_pos) + '\t.\t' + sign + '\t.\t' + 'Parent=scaffold_' + str(header_index) + ';ID=' + name # +';Left-links='+str(nr_links_left) +';Right-links='+str(nr_links_right) 

        if scaf_len > 1:
            print >> out_scaf_file, 'scaffold_' + str(header_index) + '\t' + str(prev_pos) + '\t' + str(scaf_len)
            print >> gff_file, 'scaffold_' + str(header_index) + '\tBESST V0.7\tscaffold\t' + '1' + '\t' + str(prev_pos) + '\t.\t' + '+' + '\t.\t' + 'ID=scaffold_' + str(header_index)
        items = len(print_list)
        output_scaf = ''.join([print_list[i] for i in range(0, items)])
        for i in range(0, len(output_scaf), 60):
            print >> output, output_scaf[i:i + 60]

        #print >>output, ''.join([print_list[i] for i in range(0,items)])
    return()


def RevComp(string, rev_nuc):
    #rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)



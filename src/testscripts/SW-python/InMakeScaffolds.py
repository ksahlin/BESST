

def CalculateOverlap(C_dict,scaf,side,nbr_scaf,nbr_side,gap,Scaffolds):
    #print 'Neg'
    ########## If negative gap, do Smith Waterman ###################
    ## check if negative gap, if gap is negative: do SW for the gap +300bp
    ## and retrieve a score of the best overlap
    def RevComp(string):
        rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}
        rev_comp=''.join([rev_nuc[nucl] for nucl in reversed(string)])
        return(rev_comp)
    ### GET SEQUENCES ####
    
    #First sequence
    #Need to obtain the outermost sequence within a scaffold object
    contigs_in_scaf=Scaffolds[scaf].contigs
    if side == 'R':        
        max_pos=0
        for contig in contigs_in_scaf:
            if contig.position >= max_pos:            
                overlapping_contig=contig
                max_pos=contig.position
        #if contig direction == -1 do rev comp
        #else seq1=C_dict[overlapping_contig.name][-(abs(gap)+200):]
        seq1=C_dict[overlapping_contig.name][-(abs(gap)+200):] #last part of current sequence
    else:
        min_pos=sys.maxint
        for contig in contigs_in_scaf:
            if contig.position <= min_pos:            
                overlapping_contig=contig  
                min_pos=contig.position
        seq1=RevComp(C_dict[overlapping_contig.name])
        seq1=seq1[-(abs(gap)+200):]  # last part of current sequence
    #if contig direction == -1 do rev comp
    #else seq1=C_dict[overlapping_contig.name][-(abs(gap)+200):]
    
    
    #piece of second sequence
    
    contigs_in_scaf_nbr=Scaffolds[nbr_scaf].contigs
    if nbr_side == 'R':        
        max_pos=0
        for contig in contigs_in_scaf_nbr:
            if contig.position >= max_pos:            
                overlapping_contig=contig
                max_pos=contig.position
        seq2=C_dict[overlapping_contig.name][0:abs(gap)+200] #first part of next sequence
    else:
        min_pos=sys.maxint
        for contig in contigs_in_scaf_nbr:
            if contig.position <= min_pos:            
                overlapping_contig=contig  
                min_pos=contig.position
        seq2=RevComp(C_dict[overlapping_contig.name])
        seq2=seq2[0:abs(gap)+200]  #first part of next sequence        
    
    ### START ALIGNMENT ####   
    
    score,row,column,traceback_matrix=SW.SW(seq1,seq2,-4,-4)
    #print score,row,column
    ## if score is significant and overlap is over (say) 16 bp, merge the contigs on the positions 
    ## suggested. Otherwise: set gap to 0
    if score > 0 and ( column > 16 and column < len(seq2) ):  ## no use to have: and column < len(seq1)-16)
        gap=SW.traceback_parser(score,row,column,traceback_matrix,seq1,seq2)
    elif score > 0 and (row > 16 and column == len(seq2)):
        ##error, we have chosen too short sequences or repeat structure within contigs, or by chanse got a good alignment (unlikely) 
        ## just append the contigs too eachother (0 gap distance)
        gap=0
    else:
        ## If no significant alignment was found, the most likely gap estimate between the contigs is zero
        ## since we got an negative distance estimate from the start. Set the gap to 0, abutting contigs 
        gap=0
    return(gap)
#################################################################### 


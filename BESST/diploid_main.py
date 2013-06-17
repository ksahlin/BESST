'''
Created on Jun 17, 2013

@author: ksahlin
'''

import diploid

def main(contigs):
    kmer_dict, contig_dict = diploid.get_kmers(contigs)
    potential_haplotypes = diploid.get_haplotype_regions(kmer_dict, contig_dict)
    diploid.smith_waterman(potential_haplotypes)

    return


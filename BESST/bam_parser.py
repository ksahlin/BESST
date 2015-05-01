import os
import sys
import pysam

##
# Opens a .bam or .sam file and returns the file
# object.
#
# @param bam_file_path Path to the .bam or .sam.
# 
# @return File object for .bam or .sam.
#
def open_bam_file(bam_file_path):
    bam_file_name, bam_file_ext = os.path.splitext(bam_file_path)
    if bam_file_ext == ".bam":
        return pysam.Samfile(bam_file_path, 'rb')
    elif bam_file_ext == ".sam":
        return pysam.Samfile(bam_file_path, 'r')
    else:
        return IOError("open_bam_file: File must be either .bam or .sam.")

def is_proper_aligned_unique_innie(read):
    return (read.is_reverse and not read.mate_is_reverse and read.is_read2 and read.tlen < 0 and read.rname == read.mrnm) or \
                (not read.is_reverse and read.mate_is_reverse and read.is_read2 and read.tlen > 0 and read.rname == read.mrnm ) \
                and not read.mate_is_unmapped and read.mapq > 10 and not read.is_secondary
def is_proper_aligned_unique_outie(read):
    return (read.is_reverse and not read.mate_is_reverse and read.is_read2 and read.tlen > 0 and read.rname == read.mrnm) or \
                (not read.is_reverse and read.mate_is_reverse and read.is_read2 and read.tlen < 0 and read.rname == read.mrnm ) \
                and not read.mate_is_unmapped and read.mapq > 10 and not read.is_secondary

def is_unique_read_link(read):
    # if  not read.is_unmapped and not read.mate_is_unmapped and read.rname != read.mrnm \
    # and read.opt('XT')=='U' and not read.is_secondary and read.rlen != read.alen:
    #     print read
    return not read.is_unmapped and not read.mate_is_unmapped and read.rname != read.mrnm \
    and read.mapq > 10 and not read.is_secondary


### BOWTIE ####
def proper_unique_alignment_innie_bowtie(read1,read2):
    unique1 = not read1.is_unmapped and not read1.is_secondary 
    unique2 = not read2.is_unmapped and not read2.is_secondary 
    is_innie = not read1.is_reverse and read2.is_reverse and read1.pos < read2.pos \
                or read1.is_reverse and not read2.is_reverse and read1.pos > read2.pos
    return unique1 and unique2 and is_innie 

def proper_unique_alignment_outie_bowtie(read1,read2):
    unique1 = not read1.is_unmapped and not read1.is_secondary 
    unique2 = not read2.is_unmapped and not read2.is_secondary  
    is_outie = read1.is_reverse and not read2.is_reverse and read1.pos < read2.pos \
                or not read1.is_reverse and read2.is_reverse and read1.pos > read2.pos
    return unique1 and unique2 and is_outie 

def unique_link_alignment_bowtie(read1,read2):
    unique1 = not read1.is_unmapped and not read1.is_secondary 
    unique2 = not read2.is_unmapped and not read2.is_secondary  
    return unique1 and unique2 

##############

def get_mp_observation(read1, read2, ctg_len1, ctg_len2):
    if not read1.is_reverse:
        obs1 = read1.pos + read1.qlen
    else:
        obs1 = ctg_len1 - read1.pos

    if not read2.is_reverse:
        obs2 = read2.pos + read2.qlen
    else:
        obs2 = ctg_len2 - read2.pos
    return obs1, obs2

def get_pe_observation(read1, read2, ctg_len1, ctg_len2):
    if read1.is_reverse:
        obs1 = read1.pos + read1.qlen
    else:
        obs1 = ctg_len1 - read1.pos

    if read2.is_reverse:
        obs2 = read2.pos + read2.qlen
    else:
        obs2 = ctg_len2 - read2.pos
    return obs1, obs2

class BamParser(object):
    """docstring for BamParser"""
    def __init__(self, bam_file,bam_path2 = None):
        super(BamParser, self).__init__()
        self.bam_file = open_bam_file(bam_file)
        if bam_path2:
            self.bam_file2 = open_bam_file(bam_path2) 
        self.contig_lengths = dict(zip(self.bam_file.references,self.bam_file.lengths))

    def proper_aligned_unique_pairs(self,aligner, samples=2**32):
        nr_samples = 0
        if aligner == 'bwa' or aligner == 'bwa_mem':
            for ref, length in sorted(zip(self.bam_file.references,self.bam_file.lengths),key = lambda x: x[1], reverse = True):
                try:
                    iter_ = self.bam_file.fetch(ref)
                except ValueError:
                    sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
                    sys.exit(0)
                for read in iter_:
                    if read.is_read2 and is_proper_aligned_unique_innie(read):
                        nr_samples += 1
                        yield 'innie',read
                    elif read.is_read2 and is_proper_aligned_unique_outie(read): 
                        nr_samples += 1
                        yield 'outie',read

                    if nr_samples >= samples:
                        break
            self.bam_file.seek(0)


        elif aligner == 'bowtie':
            ctgs_largest_first = [ref for ref, length in  sorted(zip(self.bam_file.references,self.bam_file.lengths), key = lambda x: x[1], reverse = True)]
            long_ctgs = set(ctgs_largest_first[:10])
            for read1,read2 in zip(self.bam_file,self.bam_file2):
                assert read1.qname == read2.qname
                same_ref = read1.rname == read2.rname
                if same_ref and proper_unique_alignment_innie_bowtie(read1,read2): 
                    if not self.bam_file.getrname( read1.tid) in long_ctgs and  not self.bam_file2.getrname(read2.tid) in long_ctgs:
                        continue
                    nr_samples += 1
                    read1.tlen = abs(read1.pos - read2.pos) + read1.rlen
                    yield 'innie',read1
                elif same_ref and proper_unique_alignment_outie_bowtie(read1,read2): 
                    if not self.bam_file.getrname( read1.tid) in long_ctgs and not self.bam_file2.getrname(read2.tid) in long_ctgs:
                        continue
                    nr_samples += 1
                    read1.tlen = abs(read1.pos - read2.pos) + 2*read1.rlen
                    yield 'outie',read1

                if nr_samples >= samples:
                    break
            self.bam_file.seek(0)      

    def aligned_reads(self,aligner):
        if aligner == 'bwa' or aligner == 'bwa_mem':
            for read in self.bam_file:
                if not read.is_unmapped:
                    yield read 
            self.bam_file.seek(0)
        elif aligner == 'bowtie':
            for read1,read2 in zip(self.bam_file,self.bam_file2):
                if not read1.is_unmapped and not read2.is_unmapped:
                    yield read1,read2           
                elif not read1.is_unmapped:
                    yield read1,False
                elif not read2.is_unmapped:
                    yield False, read2

    def unique_reads_on_different_references(self, aligner):
        read_pairs = {}
        if aligner == 'bwa_mem':
            for read in self.bam_file:
                if is_unique_read_link(read):
                    #tmp for tests:
                    #print read.qname[:-1]
                    #print read_pairs
                    if read.qname in read_pairs:
                        #print 'lol'
                        #read2 = read_pairs[read.qname]
                        #if read.tid == read2.tid:    
                        yield read, read_pairs[read.qname]
                        #else: 
                        #    pass
                        del read_pairs[read.qname]
                    else:
                        read_pairs[read.qname] = read
                    


        elif aligner == 'bwa':
            pass

        elif aligner == 'bowtie':
            for read1,read2 in zip(self.bam_file,self.bam_file2):
                assert read1.qname == read2.qname
                same_ref = read1.rname == read2.rname
                if not same_ref and unique_link_alignment_bowtie(read1,read2): 
                    read1.tlen = abs(read1.pos - read2.pos) + read1.rlen
                    yield read1, read2
                elif not same_ref and unique_link_alignment_bowtie(read1,read2): 
                    read1.tlen = abs(read1.pos - read2.pos) + 2*read1.rlen
                    yield read1, read2

        self.bam_file.seek(0)

    def long_reads_for_coverage(self):
        pass
    def long_reads_scaffold_links(self):
        pass

def main(path):
    for read in BamParser(path).reads_for_coverage():
        if read.qlen != 100:
        #print read.tags
            print read.qlen
        #print read.opt('XT')


if __name__ == '__main__':
    main(sys.argv[1])

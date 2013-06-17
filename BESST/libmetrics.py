'''
Created on Mar 7, 2013

@author: ksahlin
'''

import sys

from heapq import nlargest

from mathstats.normaldist.normal import MaxObsDistr

#with pysam.Samfile(param.bamfile, 'rb') as bam_file:

def get_metrics(bam_file, param, Information):
    def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
        k = MaxObsDistr(len(insert_list), 0.95)
        filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
        if len(insert_list) > len(filtered_list):
            return(True, filtered_list)
        else:
            return(False, filtered_list)
    informative_pair = set([147, 163]) #161,145,129,177,
    cont_names = bam_file.references
    cont_lengths = bam_file.lengths
    #cont_lengths=[int(nr) for nr in cont_lengths]  #convert long to int object
    cont_lengths_list = list(cont_lengths)
    indexes = [i for i in range(0, len(cont_lengths_list))]
    largest_contigs_indexes = nlargest(1000, indexes, key=lambda i: cont_lengths_list[i]) #get indexes of the 1000 longest contigs


    if not param.read_len: # user has not specified read len  
        #get read length
        try:
            iter_ = bam_file.fetch(cont_names[largest_contigs_indexes[0]])
        except ValueError:
            sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
            sys.exit(0)
        nr_reads = 0
        tot_read_len = 0
        for read in iter_:
            if read.rlen != 0:
                tot_read_len += read.rlen
                nr_reads += 1
            else:
                tot_read_len += read.alen
                nr_reads += 1
        param.read_len = tot_read_len / float(nr_reads)


    if param.mean_ins_size and param.std_dev_ins_size and not param.ins_size_threshold: # user has specified mean and std dev but no thresholds
        param.ins_size_threshold = param.mean_ins_size + 4 * param.std_dev_ins_size
        if param.extend_paths:
            param.contig_threshold = param.ins_size_threshold
        else:
            param.contig_threshold = param.mean_ins_size + (param.std_dev_ins_size / float(param.mean_ins_size)) * param.std_dev_ins_size
        print >> Information, '-T', param.ins_size_threshold, '-t', param.contig_threshold

    if not param.mean_ins_size: # user has not specified mean and std dev. (and no thresholds)
        counter = 1
        ins_size_reads = []
        for index in largest_contigs_indexes:
            try:
                iter_ = bam_file.fetch(cont_names[index])
            except ValueError:
                sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
                sys.exit(0)
            for read in iter_:
                if read.flag in informative_pair and read.rname == read.mrnm  and not read.mate_is_unmapped and not read.is_unmapped and read.mapq > 10: ##read.tid == read.rnext and not read.mate_is_unmapped and not read.is_unmapped: #
                    ins_size_reads.append(abs(read.tlen))
                    counter += 1
                if counter > 1000000:
                    break
            if counter > 1000000:
                break

        #get mean and std dev here. 
        #Assure that there were enough reads  for computation of mean and variance
        if len(ins_size_reads) <= 1000:
            sys.stderr.write('To few valid read alignments exists to compute mean and variance of library (need at least 1000 observations). Got only ' + str(len(ins_size_reads)) + ' valid alignments. Please specify -m and -s to the program. \nPrinting out scaffolds produced in earlier steps...')
            sys.stderr.write('\nterminating...\n')
            sys.exit(0)

        ## SMOOTH OUT THE MEAN HERE by removing extreme observations## 
        n = float(len(ins_size_reads))
        mean_isize = sum(ins_size_reads) / n
        std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5
        print >> Information, 'Mean before filtering :', mean_isize
        print >> Information, 'Std_est  before filtering: ', std_dev_isize
        extreme_obs_occur = True
        while extreme_obs_occur:
            extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, ins_size_reads)
            n = float(len(filtered_list))
            mean_isize = sum(filtered_list) / n
            std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n - 1)) ** 0.5
            ins_size_reads = filtered_list

        n = float(len(ins_size_reads))
        mean_isize = sum(ins_size_reads) / n
        std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5

        print >> Information, 'Mean converged:', mean_isize
        print >> Information, 'Std_est converged: ', std_dev_isize

        param.mean_ins_size = mean_isize
        param.std_dev_ins_size = std_dev_isize
    if not param.ins_size_threshold:
        param.ins_size_threshold = param.mean_ins_size + 4 * param.std_dev_ins_size
        if param.extend_paths:
            param.contig_threshold = param.ins_size_threshold
        else:
            param.contig_threshold = param.mean_ins_size + (param.std_dev_ins_size / float(param.mean_ins_size)) * param.std_dev_ins_size

    print >> Information, ''
    print >> Information, 'Mean of library set to:', param.mean_ins_size
    print >> Information, 'Standard deviation of library set to: ', param.std_dev_ins_size
    print >> Information, '-T (library insert size threshold) set to: ', param.ins_size_threshold
    print >> Information, '-k set to (Scaffolding with contigs larger than): ', param.contig_threshold
    print >> Information, 'Number of links required to create an edge: ', param.edgesupport
    print >> Information, 'Read length set to: ', param.read_len
    print >> Information, 'Relative weight of dominating link set to (default=3): ', param.rel_weight
    print >> Information, ''
    return()


    #### Get parameters -r, -m, -s, -T, -t for library ####
    #print >> Information, 'Computing parameters not set by user...'
    #GetParams(bam_file, param, Information)

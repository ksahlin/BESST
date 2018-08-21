'''
Created on Mar 7, 2013

@author: ksahlin
'''

from __future__ import print_function
import sys
import math

from heapq import nlargest
from collections import Counter

from mathstats.normaldist.normal import MaxObsDistr
from BESST import bam_parser
from BESST import find_bimodality
from BESST import plots




def AdjustInsertsizeDist(param, mean_insert, std_dev_insert, insert_list): 
    k = 1.5*MaxObsDistr(len(insert_list), 0.95) # allow 1.5x for thicker tails and skewed distributions
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)

def remove_outliers(ins_size_reads):
    ## SMOOTH OUT THE MEAN HERE by removing extreme observations
    if not ins_size_reads:
        return 0,0
    n = float(len(ins_size_reads))
    mean_isize = sum(ins_size_reads) / n
    std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5
    extreme_obs_occur = True
    while extreme_obs_occur:
        extreme_obs_occur, filtered_list = AdjustInsertsizeDist(param, mean_isize, std_dev_isize, ins_size_reads)
        n = float(len(filtered_list))
        mean_isize = sum(filtered_list) / n
        std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n - 1)) ** 0.5
        ins_size_reads = filtered_list
    n = float(len(ins_size_reads))
    mean_isize = sum(ins_size_reads) / n
    std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5
    return mean_isize, std_dev_isize

def get_contamination_metrics(largest_contigs_indexes, bam_file, cont_names, param, Information):
    iter_threshold = 1000000
    counter_total= 0
    count_contamine = 0

    contamination_reads = []
    sample_counter = 0

    # for index in largest_contigs_indexes:
    #     try:
    #         iter_ = bam_file.fetch(cont_names[index])
    #     except ValueError:
    #         sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
    #         sys.exit(0)
    for read in bam_file:
        # all reads mapping
        if read.rname in largest_contigs_indexes:
            sample_counter += 1
            if not read.is_unmapped: ##read.tid == read.rnext and not read.mate_is_unmapped and not read.is_unmapped: #
                counter_total += 1

            # contamination reads (mapped in reverse complemented orientation)
            if  param.orientation == 'fr' and bam_parser.is_proper_aligned_unique_outie(read, param.min_mapq):
                frag_size = abs(read.tlen)+2*param.read_len
                if param.read_len < frag_size:
                    contamination_reads.append(frag_size)
                    count_contamine += 2

            if  param.orientation == 'rf' and bam_parser.is_proper_aligned_unique_innie(read, param.min_mapq):
                frag_size = abs(read.tlen)
                if param.read_len < frag_size:
                    contamination_reads.append(frag_size)
                    count_contamine += 2

            if sample_counter >= iter_threshold:
                    break

    ## SMOOTH OUT contamine distribution here by removing extreme observations##

    n_contamine = float(len(contamination_reads))
    mean_isize = 0
    std_dev_isize = 0
    if n_contamine > 2:
        mean_isize = sum(contamination_reads) / n_contamine
        std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), contamination_reads))) / (n_contamine - 1)) ** 0.5
        print('Contamine mean before filtering :', mean_isize, file=Information)
        print('Contamine stddev before filtering: ', std_dev_isize, file=Information)
        # #new method to find distribution
        # dist1, dist2, mean_isize, std_dev_isize, outliers_mean, outliers_stddev = find_bimodality.split_distribution(contamination_reads)
        # n_contamine = len(dist1)
        extreme_obs_occur = True
        while extreme_obs_occur:
            extreme_obs_occur, filtered_list = AdjustInsertsizeDist(param, mean_isize, std_dev_isize, contamination_reads)
            n_contamine = float(len(filtered_list))
            if n_contamine > 2:
                mean_isize = sum(filtered_list) / n_contamine
                std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n_contamine - 1)) ** 0.5
                contamination_reads = filtered_list
            else:
                break
        print('Contamine mean converged:', mean_isize, file=Information)
        print('Contamine std_est converged: ', std_dev_isize, file=Information)

    if counter_total > 0:
        contamination_ratio = 2*n_contamine / float(counter_total)
    else:
        contamination_ratio = 0

    if mean_isize >= param.mean_ins_size or std_dev_isize >= param.std_dev_ins_size or contamination_ratio <= 0.05:
        # either contamine mean or stddev is higher than MP lb mean which means it's spurious alignments or
        # other wierd thing -> no true PE-contamine as artifact discribed in the illumina MP construction protocol.
        # or we have less than 5% of contamine reads. Then we skip dealing with them since they constitute such a small
        # fraction of the total reads and introduce more complexity when orienting scaffolds in pathfinder module
        param.contamination_ratio = False
        param.contamination_mean = 0
        param.contamination_stddev = 0
    else:
        param.contamination_mean = mean_isize
        param.contamination_stddev = std_dev_isize
        param.contamination_ratio = contamination_ratio

    bam_file.reset()
    return n_contamine

def argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

def sum_chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield sum(l[i:i+n])

def getdistr(ins_size_reads, cont_lengths_list, param, Information):
    largest_contigs = [int(x) for x in sorted(nlargest(1000, cont_lengths_list))]
    #print largest_contigs
    #sorted_lengths = sorted(cont_lengths_list)

    min_ctg_length = min(largest_contigs)
    max_isize = int(max(ins_size_reads))
    adjusted_distribution = [0]*int(max_isize+1)
    ## Make hasmap of the number of contigs and their sum 
    ## that is larger than a given isize
    current_sum_ctgs = sum(largest_contigs)
    current_nr_ctgs = len(largest_contigs)
    ctgs_larger_than = [] # array indexed by isize, item is tuple is a tuple (n_ctgs, their nucleotide sum)

    ctgs_larger_than.append( (current_nr_ctgs, current_sum_ctgs) )
    current_smallest_contig = largest_contigs[0]
    current_smallest_contig_index = 0
    upper_isize = min(max_isize+1, largest_contigs[-1])
    for isize in range(upper_isize):
        if isize <= current_smallest_contig:
            ctgs_larger_than.append( (current_nr_ctgs, current_sum_ctgs) )
        else:
            #current_smallest_contig_index += 1
            while isize > largest_contigs[current_smallest_contig_index]:
                current_smallest_contig_index += 1
                current_nr_ctgs -= 1
                current_sum_ctgs -= current_smallest_contig

            ctgs_larger_than.append( (current_nr_ctgs, current_sum_ctgs) )
            current_smallest_contig = largest_contigs[current_smallest_contig_index]

    # print min_ctg_length, max_isize
    # print ctgs_larger_than
    for o in ins_size_reads:
        obs = int(o)
        if obs > upper_isize:
            continue
        nr_ctgs = ctgs_larger_than[obs][0]
        sum_ctgs = ctgs_larger_than[obs][1]
        w = float(max(sum_ctgs - (obs-1)*nr_ctgs, 10000))
        adjusted_distribution[obs] += 1/w

    # if min_ctg_length >= max_isize:
    #     for o in ins_size_reads:
    #         obs = int(o)
    #         w = max(float(sum_ctgs - (obs-1)*nr_ctgs), 1)
    #         adjusted_distribution[obs] += 1/w

    tot_density = float(sum(adjusted_distribution))
    cum_sum = 0
    curr_isize = 0
    median_density = tot_density/2.0
    while cum_sum <= median_density:
         cum_sum += adjusted_distribution[curr_isize]
         curr_isize +=1

    median_adj = curr_isize

    if param.plots:
        plots.histogram(ins_size_reads, param, bins=100, x_label='fragment length', y_label='frequency', title='Frag_length_distribuion' + '.' + param.bamfile.split('/')[-1])
        plots.dot_plot(list(range(len(adjusted_distribution))), adjusted_distribution, param, x_label='fragment length', y_label='frequency', title='Frag_length_distribuion_adjusted' + '.' + param.bamfile.split('/')[-1], set_marker= '.')

    ## find a stable mode of the fitted distribution since our distribution is a sample.
    ## It is not very stable just to chose the isize with the highest count
    mode_for_different_windows = []
    for chunk_size in range(1,102, 5):
        adj_distr_chunked = list(sum_chunks(adjusted_distribution, chunk_size))
        mode_adj = (argmax(adj_distr_chunked) + 0.5)*chunk_size
        mode_for_different_windows.append(int(mode_adj))
        print("mode for chunk size ", chunk_size, " : ", mode_adj, file=Information)
    mode_adj = sorted(mode_for_different_windows)[int(len(mode_for_different_windows)/2)]
    print("Choosing mode:", mode_adj)

    #mode_adj = argmax(adjusted_distribution)
    mu_adj = sum([i_f_x1[0]*i_f_x1[1] for i_f_x1 in enumerate(adjusted_distribution)])/tot_density
    sigma_adj = math.sqrt(sum([(i_f_x[0]-mu_adj)**2 * i_f_x[1] for i_f_x in enumerate(adjusted_distribution)]) / tot_density)
    m_3 = sum([(i_f_x2[0]-mu_adj)**3 * i_f_x2[1] for i_f_x2 in enumerate(adjusted_distribution)]) / tot_density

    # m_3 = sum(map(lambda x: (x - mean_isize) ** 3, ins_size_reads))/n
    skew_adj = m_3 / sigma_adj**3

    print('mu_adjusted:{0}, sigma_adjusted:{1}, skewness_adjusted:{2}'.format(mu_adj, sigma_adj, skew_adj))
    return adjusted_distribution, mu_adj, sigma_adj, skew_adj, median_adj, mode_adj
#with pysam.Samfile(param.bamfile, 'rb') as bam_file:

def get_metrics(bam_file, param, Information):
    #informative_pair = set([147, 163]) #161,145,129,177,
    cont_names = bam_file.references
    cont_lengths = bam_file.lengths
    #cont_lengths=[int(nr) for nr in cont_lengths]  #convert long to int object
    cont_lengths_list = list(cont_lengths)
    indexes = [i for i in range(0, len(cont_lengths_list))]
    largest_contigs_indexes = set(nlargest(1000, indexes, key=lambda i: cont_lengths_list[i])) #get indexes of the 1000 longest contigs

    param.lognormal = False #default as False, but cna change below dependant on skew

    try:
        bam_file.fetch(cont_names[0])
    except ValueError:
        sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
        sys.exit(0)
    #largest_contigs_indexes = nlargest(1000, indexes, key=lambda i: cont_lengths_list[i]) #get indexes of the 1000 longest contigs
    
    #print largest_contigs_indexes

    if not param.read_len: # user has not specified read len  
        #get read length
        nr_reads = 0
        tot_read_len = 0
        # for index in largest_contigs_indexes:
        #     try:
        #         iter_ = bam_file.fetch(cont_names[index])
        #     except ValueError:
        #         sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
        #         sys.exit(0)

        for read in bam_file:
            if read.rlen != 0:
                tot_read_len += read.rlen
                nr_reads += 1
            else:
                tot_read_len += read.alen
                nr_reads += 1
            if nr_reads >= 1000:
                param.read_len = tot_read_len / float(nr_reads)
                break        
        else:
            sys.stderr.write('Did not get sufficient readmappings to calculate\
             read_length from mappings. Got {0} mappings. Please provide this parameter or more importantly\
             check why almost no reads are mapping to the contigs.\nterminating..\n'.format(nr_reads))
            sys.exit(0)

        bam_file.reset()

    if param.mean_ins_size and param.std_dev_ins_size and not param.ins_size_threshold: # user has specified mean and std dev but no thresholds
        param.ins_size_threshold = param.mean_ins_size + 6 * param.std_dev_ins_size
        if param.extend_paths:
            param.contig_threshold = param.mean_ins_size + 4 * param.std_dev_ins_size
        else:
            param.contig_threshold = param.mean_ins_size + (param.std_dev_ins_size / float(param.mean_ins_size)) * param.std_dev_ins_size
        print('-T', param.ins_size_threshold, '-t', param.contig_threshold, file=Information)

    if not param.mean_ins_size: # user has not specified mean and std dev. (and no thresholds)
        #total_reads_iterated_through = 0
        counter = 1
        ins_size_reads = []
        # for index in largest_contigs_indexes:
        #     try:
        #         iter_ = bam_file.fetch(cont_names[index])
        #     except ValueError:
        #         sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
        #         sys.exit(0)
        for read in bam_file:
            if param.orientation == 'fr' and bam_parser.is_proper_aligned_unique_innie(read, param.min_mapq):
                if read.rname in largest_contigs_indexes:
                    ins_size_reads.append(abs(read.tlen))
                    counter += 1
            if param.orientation == 'rf' and bam_parser.is_proper_aligned_unique_outie(read, param.min_mapq):
                if read.rname in largest_contigs_indexes:
                    ins_size_reads.append(abs(read.tlen) + 2*param.read_len)
                    counter += 1
            if counter > 1000000:
                break
        bam_file.reset()
        # if counter > 1000000:
        #     break
        print("Estimating insert size from {0} mappings with quality over --min_mapq {1}.".format(counter, param.min_mapq))
        print("Estimating insert size from {0} mappings with quality over --min_mapq {1}.".format(counter, param.min_mapq), file=Information)
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
        print('Mean before filtering :', mean_isize, file=Information)
        print('Std_est  before filtering: ', std_dev_isize, file=Information)
        extreme_obs_occur = True
        while extreme_obs_occur:
            extreme_obs_occur, filtered_list = AdjustInsertsizeDist(param, mean_isize, std_dev_isize, ins_size_reads)
            n = float(len(filtered_list))
            mean_isize = sum(filtered_list) / n
            std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n - 1)) ** 0.5
            ins_size_reads = filtered_list

        n = float(len(ins_size_reads))
        mean_isize = sum(ins_size_reads) / n
        std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5

        print('Mean converged:', mean_isize, file=Information)
        print('Std_est converged: ', std_dev_isize, file=Information)

        param.mean_ins_size = mean_isize
        param.std_dev_ins_size = std_dev_isize

        m_3 = sum([(x - mean_isize) ** 3 for x in ins_size_reads])/n
        skewness = m_3 / std_dev_isize**3
        param.skewness = skewness
        print('Skewness of distribution: ', param.skewness, file=Information)

        # weight each observation with how likely it is to see it
        adj_distr, mu_adj, sigma_adj, skew_adj, median_adj, mode_adj = getdistr(ins_size_reads, cont_lengths_list, param, Information)
        param.skew_adj = skew_adj
        param.empirical_distribution = dict(list(zip(list(range(len(adj_distr))), adj_distr))) #Counter(ins_size_reads)
        print('Mean of getdistr adjusted distribution: ', mu_adj, file=Information)
        print('Sigma of getdistr adjusted distribution: ', sigma_adj, file=Information)
        print('Skewness of getdistr adjusted distribution: ', skew_adj, file=Information)
        print('Median of getdistr adjusted distribution: ', median_adj, file=Information)
        print('Mode of getdistr adjusted distribution: ', mode_adj, file=Information)
        print('Using mean and stddev of getdistr adjusted distribution from here: ', mu_adj, sigma_adj, file=Information)
        param.mean_ins_size = mu_adj
        param.std_dev_ins_size = sigma_adj

        #### If skewness (of original - not the getdistr)is positive and larger than 0.5 
        #### (big enough skew to have impact), we fit to the lognormal distribution 
        if param.skew_adj > 0.5 and math.log(median_adj) > math.log(mode_adj):

            ###################################################
            #### NOTE: Fitting lognormal of original sample, not getdistr adjusted for 
            #### smaller isizes observation bias
            #### because I don't know how yet. If the two distributions are not too unsimilar
            #### it should be a good approximation in practice
            # median = sorted(ins_size_reads)[len(ins_size_reads)/2]
            # mode = mean_isize - 3*(mean_isize - median)
            # print >> Information, 'Mode on initial sample (not getdistr adjusted): ', mode
            # print >> Information, "Median on initial sample (not getdistr adjusted)", median
            # print "mode:", mode
            # print "median", median
            # param.lognormal_mean = math.log(median)
            # param.lognormal_sigma = math.sqrt(param.lognormal_mean - math.log(mode))
            # print >> Information, 'Lognormal mean (not getdistr adjusted): ', param.lognormal_mean
            # print >> Information, "Lognormal stddev (not getdistr adjusted)", param.lognormal_sigma
            #################################################

            ## the statistics for the getdistr adjusted distribution
            #mode_adj = mu_adj - 3*(mu_adj - median_adj)
            print('Mode on getdistr adjusted: ', mode_adj, file=Information)
            print("Median on getdistr adjusted:", median_adj, file=Information)
            print("mode adj:", mode_adj)
            print("median adj", median_adj)
            param.lognormal_mean = math.log(median_adj)
            param.lognormal_sigma = math.sqrt(param.lognormal_mean - math.log(mode_adj))
            print('Lognormal mean getdistr adjusted: ', param.lognormal_mean, file=Information)
            print("Lognormal stddev getdistr adjusted", param.lognormal_sigma, file=Information)

            param.lognormal = True
            #sys.exit()

        # TODO: calculate skew of contamination distribution

        # stddev_fit = (sum(map(lambda x: (x - mode)**2 , ins_size_reads))/n)**0.5
        # print stddev_fit
        
        # import matplotlib.pyplot as plt
        # plt.hist(ins_size_reads,100)
        # plt.savefig("/Users/ksahlin/_tmp/BESST_ILP/ARABI_27_statistical_score_no_paths/isize_plot")
    else:
        total_reads_iterated_through = None

    if not param.ins_size_threshold:
        param.ins_size_threshold = param.mean_ins_size + 6 * param.std_dev_ins_size
        if param.extend_paths:
            param.contig_threshold = param.mean_ins_size + 4 * param.std_dev_ins_size
        else:
            param.contig_threshold = param.mean_ins_size + (param.std_dev_ins_size / float(param.mean_ins_size)) * param.std_dev_ins_size

    ## finally, get a reverse complemented read contamination distribution from the MP library if it exists
    n_contamine = get_contamination_metrics(largest_contigs_indexes,bam_file, cont_names, param, Information)


    print('', file=Information)
    print('LIBRARY STATISTICS', file=Information)
    print('Mean of library set to:', param.mean_ins_size, file=Information)
    print('Standard deviation of library set to: ', param.std_dev_ins_size, file=Information)
    print('MP library PE contamination:', file=Information)
    print('Contamine rate (rev comp oriented) estimated to: ',  param.contamination_ratio, file=Information)
    print('lib contamine mean (avg fragmentation size): ', param.contamination_mean, file=Information)
    print('lib contamine stddev: ', param.contamination_stddev, file=Information) 
    print('Number of contamined reads used for this calculation: ', n_contamine, file=Information)


    print('-T (library insert size threshold) set to: ', param.ins_size_threshold, file=Information)
    print('-k set to (Scaffolding with contigs larger than): ', param.contig_threshold, file=Information)
    print('Number of links required to create an edge: ', param.edgesupport, file=Information)
    print('Maximum identical contig-end overlap-length to merge of contigs that are adjacent in a scaffold: ', param.max_contig_overlap, file=Information)
    print('Read length set to: ', param.read_len, file=Information)
    print('', file=Information)
    return()


    #### Get parameters -r, -m, -s, -T, -t for library ####
    #print >> Information, 'Computing parameters not set by user...'
    #GetParams(bam_file, param, Information)

'''
    Created on Mar 31, 2012

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


class parameter(object):
    '''
    classdocs
    '''
    def __init__(self, parameter_mean_coverage=None, parameter_std_dev_coverage=None,
                 parameter_mean_ins_size=None , parameter_std_dev_ins_size=None,
                 parameter_output_directory=None, parameter_bamfile=None,
                 parameter_read_len=None,
                 parameter_ins_size_threshold=None, parameter_contigfile=None,
                 parameter_edgesupport=None, parameter_contig_threshold=None,
                 parameter_scaffold_indexer=0, parameter_first_lib=None,
                 parameter_cov_cutoff=None, parameter_tot_assembly_length=None,
                 parameter_current_N50=None, parameter_current_L50=None,
                 parameter_hapl_ratio=None, parameter_hapl_threshold=None,
                 parameter_detect_haplotype=None, parameter_detect_duplicate=None,
                 parameter_gff_file=None, parameter_information_file=None,
                 parameter_extend_paths=None, FASTER_ILP=None,
                 parameter_development=None, parameter_plots=None, parameter_path_threshold=None,
                 path_gaps_estimated =0, parameter_gap_estimations = [],
                 contamination_mean=None, contamination_stddev = None, contamination_ratio=0,
                 no_score=None, orientation = None, contig_index= None,
                 score_cutoff = None, max_extensions = None, NO_ILP=None, print_scores=None,
                 pass_number=None, max_contig_overlap=None, lower_cov_cutoff=None):



        # Library information
        self.mean_ins_size = parameter_mean_ins_size
        self.std_dev_ins_size = parameter_std_dev_ins_size
        self.read_len = parameter_read_len
        self.mean_coverage = parameter_mean_coverage
        self.lower_cov_cutoff = lower_cov_cutoff
        self.cov_cutoff = parameter_cov_cutoff

        # Contig information
        self.contig_index = contig_index
        # Algorithm parameters
        self.score_cutoff = score_cutoff
        self.max_extensions = max_extensions
        self.max_contig_overlap = max_contig_overlap
        # Assembly information
        # Algorithm information
        # Output information
        self.std_dev_coverage = parameter_std_dev_coverage
        self.output_directory = parameter_output_directory
        self.bamfile = parameter_bamfile
        self.ins_size_threshold = parameter_ins_size_threshold
        self.contigfile = parameter_contigfile
        self.edgesupport = parameter_edgesupport
        self.contig_threshold = parameter_contig_threshold
        self.scaffold_indexer = parameter_scaffold_indexer
        self.first_lib = parameter_first_lib
        self.tot_assembly_length = parameter_tot_assembly_length
        self.current_N50 = parameter_current_N50
        self.current_L50 = parameter_current_L50
        self.hapl_ratio = parameter_hapl_ratio
        self.hapl_threshold = parameter_hapl_threshold
        self.detect_haplotype = parameter_detect_haplotype
        self.detect_duplicate = parameter_detect_duplicate
        self.gff_file = parameter_gff_file
        self.information_file = parameter_information_file
        self.extend_paths = parameter_extend_paths
        self.development = parameter_development
        self.plots = parameter_plots
        self.path_threshold = parameter_path_threshold
        self.no_score = no_score
        self.orientation = orientation
        self.pass_number = pass_number
        # debugging and logging
        self.print_scores = print_scores


        # related to gap distances
        self.path_gaps_estimated = path_gaps_estimated
        self.gap_estimations = parameter_gap_estimations

        # If MP library (PE contamination)
        self.contamination_ratio = contamination_ratio
        self.contamination_mean = contamination_mean
        self.contamination_stddev = contamination_stddev
        self.NO_ILP = NO_ILP
        self.FASTER_ILP = FASTER_ILP

    def get_params(self):
        output = "param\tvalue\n"
        values = "".join([ "{0}\t{1}\n".format(attr, value) if value == None or type(value) in [bool, int,float, file] or len(value) < 5 else "" for attr, value in self.__dict__.items() if not callable(value)])
        output += values
        return output

class counters(object):
    def __init__(self, param_count=None, param_non_unique=None, param_non_unique_for_scaf=None,
                 param_nr_of_duplicates=None, param_prev_obs1=None, param_prev_obs2=None,
                 param_reads_with_too_long_insert=None):

        self.count = param_count
        self.non_unique = param_non_unique
        self.non_unique_for_scaf = param_non_unique_for_scaf
        self.nr_of_duplicates = param_nr_of_duplicates
        self.prev_obs1 = param_prev_obs1
        self.prev_obs2 = param_prev_obs2
        self.reads_with_too_long_insert = param_reads_with_too_long_insert

#
#def MemoryUsage():
#    percent_usage = subprocess.Popen("ps -v | grep 'python Main.py' | awk '{sum+=$12} END {print sum}'",
#                                     shell=True,
#                                     stdout=subprocess.PIPE,
#                                     )
#    stdout_list = percent_usage.communicate()[0]
#    print 'Percentage of memory occupied by BESST from total memory: ', stdout_list
#    return()




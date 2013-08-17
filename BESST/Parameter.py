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
import subprocess


class parameter(object):
    '''
    classdocs
    '''
    def __init__(self, parameter_mean_coverage=None, parameter_std_dev_coverage=None,
                 parameter_mean_ins_size=None , parameter_std_dev_ins_size=None,
                 parameter_output_directory=None, parameter_bamfile=None,
                 parameter_read_len=None, parameter_rel_weight=None,
                 parameter_ins_size_threshold=None, parameter_contigfile=None,
                 parameter_edgesupport=None, parameter_contig_threshold=None,
                 parameter_scaffold_indexer=0, parameter_first_lib=None,
                 parameter_cov_cutoff=None, parameter_tot_assembly_length=None,
                 parameter_current_NG50=None, parameter_current_LG50=None,
                 parameter_hapl_ratio=None, parameter_hapl_threshold=None,
                 parameter_detect_haplotype=None, parameter_detect_duplicate=None,
                 parameter_gff_file=None, parameter_information_file=None,
                 parameter_fosmidpool=None, parameter_extend_paths=None,
                 parameter_development=None, parameter_plots=None, parameter_path_threshold=None):

        # Library information
        # Contig information
        # Algorithm parameters
        # Assembly information
        # Algorithm information
        # Output information
        self.mean_coverage = parameter_mean_coverage
        self.std_dev_coverage = parameter_std_dev_coverage
        self.mean_ins_size = parameter_mean_ins_size
        self.std_dev_ins_size = parameter_std_dev_ins_size
        self.output_directory = parameter_output_directory
        self.bamfile = parameter_bamfile
        self.read_len = parameter_read_len
        self.rel_weight = parameter_rel_weight
        self.ins_size_threshold = parameter_ins_size_threshold
        self.contigfile = parameter_contigfile
        self.edgesupport = parameter_edgesupport
        self.contig_threshold = parameter_contig_threshold
        self.scaffold_indexer = parameter_scaffold_indexer
        self.first_lib = parameter_first_lib
        self.cov_cutoff = parameter_cov_cutoff
        self.tot_assembly_length = parameter_tot_assembly_length
        self.current_NG50 = parameter_current_NG50
        self.current_LG50 = parameter_current_LG50
        self.hapl_ratio = parameter_hapl_ratio
        self.hapl_threshold = parameter_hapl_threshold
        self.detect_haplotype = parameter_detect_haplotype
        self.detect_duplicate = parameter_detect_duplicate
        self.gff_file = parameter_gff_file
        self.information_file = parameter_information_file
        self.fosmidpool = parameter_fosmidpool
        self.extend_paths = parameter_extend_paths
        self.development = parameter_development
        self.plots = parameter_plots
        self.path_threshold = parameter_path_threshold


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




#!/usr/bin/env python3

# Copyright (C) 2017, Weizhi Song, Torsten Thomas.
# modified from https://github.com/songweizhi/Binning_refiner/blob/master/bin/Binning_refiner
# with a copy of the GNU Affero General Public License along with this program.

import shutil
import os
from Bio import SeqIO
import logging
import numpy as np
import copy

# package logger
logger = logging.getLogger(__name__)


class FindShare:
    def __init__(self , hic_folder , shotgun_folder , ref_len , min_binsize , min_frac , min_complete_size):
        self.min_binsize = min_binsize
        self.ref_len = ref_len
        self.hic_bin = []
        self.shotgun_bin = []
        self.shared_bin = []
        self.shared_bin_plus = None ##add relatively complete bins
        self.frac = min_frac
        self.complete = min_complete_size

        input_bin_subfolders = [hic_folder , shotgun_folder]

        for bin_subfolder in input_bin_subfolders:
            # get bin file list in each input bin subfolder
            pwd_bin_subfolder = '%s' % bin_subfolder
            bin_file_list = self.get_no_hidden_folder_list(pwd_bin_subfolder)

            bin_ext_list = set()
            for each_bin in bin_file_list:
                each_bin_ext = each_bin.split('.')[-1]
                bin_ext_list.add(each_bin_ext)

            if len(bin_ext_list) > 1:
                logger.info('Program exited, please make sure all bins within %s folder have the same extension' % bin_subfolder)
                raise IOError('Extensions are not consistent')


        ctg_length_dict = {}
        ctg_to_bin_dict = {}

        for bin_subfolder in input_bin_subfolders:
        # bin_subfolder is just one folder of input files
        # get bin file list in each input bin subfolder
            bin_file_list = self.get_no_hidden_folder_list(bin_subfolder)

            logger.info('Read in %s %s bins' % (len(bin_file_list), bin_subfolder))
            
            # each bin is one bin in the bin_subfolder
            for each_bin in bin_file_list:
                pwd_each_bin = '%s/%s' % (bin_subfolder, each_bin)
                gp_temp = []

                for each_seq in SeqIO.parse(pwd_each_bin, 'fasta'):
                    each_seq_id = str(each_seq.id)
                    each_seq_len = len(each_seq.seq)

                    # store contig length into dict
                    if each_seq_id not in ctg_length_dict:
                        ctg_length_dict[each_seq_id] = each_seq_len

                    # store contig to bin info into dict
                    if each_seq_id not in ctg_to_bin_dict:
                        ctg_to_bin_dict[each_seq_id] = [each_bin]
                    else:
                        ctg_to_bin_dict[each_seq_id].append(each_bin)

                    gp_temp.append(each_seq_id)

                gp_temp = list(np.unique(gp_temp)) 
                if bin_subfolder == hic_folder:
                    self.hic_bin.append(gp_temp)
                else:
                    self.shotgun_bin.append(gp_temp)

        ##########find shared contigs####################
        ctg_to_bin_dict_shared = {}
        for each_ctg in ctg_to_bin_dict:
            if len(ctg_to_bin_dict[each_ctg]) == len(input_bin_subfolders):
                ctg_to_bin_dict_shared[each_ctg] = '___'.join(ctg_to_bin_dict[each_ctg])


        concatenated_bins_to_ctg_dict = {}
        concatenated_bins_length_dict = {}

        for each_shared_ctg in ctg_to_bin_dict_shared:
            concatenated_bins = ctg_to_bin_dict_shared[each_shared_ctg]

            if concatenated_bins not in concatenated_bins_to_ctg_dict:
                concatenated_bins_to_ctg_dict[concatenated_bins] = [each_shared_ctg]
                concatenated_bins_length_dict[concatenated_bins] = ctg_length_dict[each_shared_ctg]
            else:
                concatenated_bins_to_ctg_dict[concatenated_bins].append(each_shared_ctg)
                concatenated_bins_length_dict[concatenated_bins] += ctg_length_dict[each_shared_ctg]

        # remove short length; maximum_length is the largest shared bins

        for concatenated_bins in concatenated_bins_length_dict:
            concatenated_bins_length = concatenated_bins_length_dict[concatenated_bins]
            if concatenated_bins_length >= self.min_binsize:
                self.shared_bin.append(concatenated_bins_to_ctg_dict[concatenated_bins])


        logger.info('Got %s shared bins with size larger than %s Kbp' % (len(self.shared_bin), self.min_binsize))

        conINshare = []
        for bin in self.shared_bin:
            for i in bin:
                conINshare.append(i)

        #######################################################
        addbin = []
        addbin_size = []
        #######################################################
        #######################################################
        for i in range(len(self.hic_bin)):
            uninclude_size = 0
            total_size = 0
            temp_bin = []
            for temp_name in self.hic_bin[i]:
                temp_len = self.ref_len[temp_name]
                total_size += temp_len
                if temp_name not in conINshare:
                    uninclude_size += temp_len
                    temp_bin.append(temp_name)
            if uninclude_size/total_size > self.frac and uninclude_size > self.complete:
                addbin.append(temp_bin)
                addbin_size.append(uninclude_size)
    
        #######################################################
        #######################################################
        for bin in addbin:
            for contig in bin:
                conINshare.append(contig)
        #######################################################
        #######################################################
        for i in range(len(self.shotgun_bin)):
            uninclude_size = 0
            total_size = 0
            temp_bin = []
            for temp_name in self.shotgun_bin[i]:
                temp_len = self.ref_len[temp_name]
                total_size += temp_len
                if temp_name not in conINshare:
                    uninclude_size += temp_len
                    temp_bin.append(temp_name)
            if uninclude_size/total_size > self.frac and uninclude_size > self.complete:
                addbin.append(temp_bin)
                addbin_size.append(uninclude_size)

        #######################################################
        #######################################################
        ##########insert the unselected bin into shared bins##############
        binsize_shared_bin = []
        for i in range(len(self.shared_bin)):
            total_size = 0
            for temp_name in self.shared_bin[i]:
                temp_len = self.ref_len[temp_name]
                total_size += temp_len
            binsize_shared_bin.append(total_size)
            
    
        self.shared_bin_plus = copy.deepcopy(self.shared_bin)
        binsize_shared_bin_plus = copy.deepcopy(binsize_shared_bin)


        for i in range(len(addbin_size)):
            temp_bin = addbin[i]
            temp_binsize = addbin_size[i]
            temp_index = self.search(binsize_shared_bin_plus , temp_binsize)
            self.shared_bin_plus.insert(temp_index , temp_bin)
            binsize_shared_bin_plus.insert(temp_index , temp_binsize)


    def get_no_hidden_folder_list(self , wd):
        folder_list = []
        for each_folder in os.listdir(wd):
            if not each_folder.startswith('.'):
                folder_list.append(each_folder)

        folder_list_sorte = sorted(folder_list)
        return folder_list_sorte

    def search(self , size , item):
        for i in range(len(size)):
            if item >= size[i]:
                break
        return(i)




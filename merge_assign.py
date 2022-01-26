import numpy as np
import copy
import logging

# package logger
logger = logging.getLogger(__name__)

def merge(gp_share_plus , ref_ind , map_matrix , alpha):
########step 2: merge fragmented bins############
    num_bin = len(gp_share_plus)
    binned_contig = []
    set_contig = {}
    label_contig = []
    ref_local = {}
    temp_ind = 0
    for i in range(num_bin):
        for contig in gp_share_plus[i]:
            set_contig[contig] = i
            binned_contig.append(contig)
            label_contig.append(ref_ind[contig])
            ref_local[contig] = temp_ind
            temp_ind += 1
            
    map_filter = map_matrix.tocsr()
    map_filter = map_filter[label_contig , :]
    map_filter = map_filter.tocsc()
    map_filter = map_filter[: , label_contig]

    map_filter = map_filter.tocoo()
    iter_dist = np.zeros((num_bin,num_bin))

    map_row = map_filter.row
    map_col = map_filter.col
    map_data = map_filter.data

    index = map_row < map_col
    map_row = map_row[index]
    map_col = map_col[index]
    map_data = map_data[index]

    ref_ind_rev = {}
    for key , value in ref_local.items():
        ref_ind_rev[value] = key
        

    for temp_row,temp_col,temp_data in zip(map_row, map_col , map_data):
        contig1 = ref_ind_rev[temp_row]
        contig2 = ref_ind_rev[temp_col]
        if (contig1 in binned_contig) and (contig2 in binned_contig):
            id1con = set_contig[contig1]
            id2con = set_contig[contig2]
            if id1con < id2con:
                iter_dist[id1con , id2con] = iter_dist[id1con , id2con] + temp_data
            else:
                iter_dist[id2con , id1con] = iter_dist[id2con , id1con] + temp_data
        
    for i in range(num_bin):
        for j in range(i , num_bin):
            if i == j:
                iter_dist[i , j] = 2*iter_dist[i , j]/(len(gp_share_plus[i]) * len(gp_share_plus[j]))
            else:
                iter_dist[i , j] = iter_dist[i , j]/(len(gp_share_plus[i]) * len(gp_share_plus[j]))
                iter_dist[j , i] = iter_dist[i , j]

    combine = {}
    for i in range(num_bin):
        temp_dist = iter_dist[i , :]
        temp_combine = []
        for j in range((i+1),num_bin):
            if np.bitwise_and(temp_dist[j] > alpha * temp_dist[i] , i != j):
                temp_combine.append(j)
        combine[i] = temp_combine


    assign = np.ones(num_bin)*num_bin
    for key , value in combine.items():
        if assign[key] == num_bin:
            temp = []
            for i in value:
                if assign[i] != num_bin:
                    temp.append(assign[i])
            if len(temp) == 0:
                assign[key] = key
                for i in value:
                    assign[i] =key
            elif len(set(temp)) == 1:
                assign[key] = temp[0]
        else:
            temp = []
            uninc = []
            for i in value:
                if assign[i] != num_bin:
                    temp.append(assign[i])
                else:
                    uninc.append(i)
            if len(temp) == 0:
                for i in value:
                    assign[i] = assign[key]
            elif len(set(temp)) == 1:
                assign[key] = temp[0]
                for i in uninc:
                    assign[i] = temp[0]

    gp_merge = [[] for i in range(len(np.unique(assign)))]
    order_bin = list(np.unique(assign))
    for i in range(len(assign)):
        idbin = assign[i]
        pos = order_bin.index(idbin)
        for contig in gp_share_plus[i]:
            gp_merge[pos].append(contig)
    if gp_merge[-1] == num_bin:
        gp_merge = gp_merge[0:(len(gp_merge)-1)]
    return gp_merge



def assign(gp_merge , map_matrix , ref_ind , beta):
########step 3:  assigned unbinned contigs############
    map_filter = map_matrix.tolil()
    bin2bin_dist = np.zeros(len(gp_merge))
    for i in range(len(gp_merge)):
        info_avg = 0
        for contig1 in gp_merge[i]:
            for contig2 in gp_merge[i]:
                index1 = ref_ind[contig1]
                index2 = ref_ind[contig2]
                temp = map_filter[index1 , index2]
                info_avg += temp
        bin2bin_dist[i] = info_avg/(len(gp_merge[i]) * len(gp_merge[i]))

    ######Select unbinned contigs###########
    binned_name = []
    for bin in gp_merge:
        for i in bin:
            binned_name.append(i)

    unbinned = []
    for i in ref_ind.keys():
        if i not in binned_name:
            unbinned.append(i)

    gp_final_refine = copy.deepcopy(gp_merge)
    count = 0
    drop = 0
    for conName in unbinned:
        index1 = ref_ind[conName]
        con2bin_dist = np.zeros(len(gp_merge))
        for i in range(len(gp_merge)):
            info = 0
            for contig1 in gp_merge[i]:
                index2 = ref_ind[contig1]
                info += map_filter[index1 , index2]
            con2bin_dist[i] = info/len(gp_merge[i])
        if np.max(con2bin_dist) < beta * bin2bin_dist[np.argmax(con2bin_dist)]:
            drop += 1
            continue
        count += 1
        gp_final_refine[np.argmax(con2bin_dist)].append(conName)
    logger.info('{} contigs are dropped and {} contigs are added into the bins'.format(drop,count))
    return gp_final_refine

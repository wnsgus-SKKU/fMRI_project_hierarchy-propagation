import numpy as np
from scipy import io
from scipy.io import savemat
from scipy.signal import find_peaks
from function import *principal_delay

if __name__ == '__main__':

    path ="/home/cjh/skku/2021_second_semester/FMRI/project/Data"
   

    epi = loadMat(path+"/epi1msk","epi1msk")
    global_mean_sig = epi.mean(axis=0) # calculate the global mean of input data - Figure 2A upper graph

    # calculate delay profile across time segments: idx_dly

    locs_seg_edge = find_peaks(-1 * global_mean_sig)[0] # cut the global mean signal into time segments - Figure 2A upper graph
    locs_global_peaks = find_peaks(global_mean_sig)[0] # find global mean signal peak
    
    delay_mat = np.zeros((epi.shape[0], locs_seg_edge.shape[0]-1),object) # delay matrix: vertices * segments

    for li in np.arange(locs_seg_edge.shape[0]-1):
        tmp_seg = np.array([]) 
        tmp_seg = epi[:,locs_seg_edge[li]:locs_seg_edge[li+1]] # get one time segment

        # locate max peaks in each vertex
        for lj in np.arange(epi.shape[0]):
            tmp_vertex = tmp_seg[lj,:] # get one vertex in one time segment
            value_vertex_peak, locs_vertex_peak = findPeaks(tmp_vertex) # and find peak of that vertex - Figure 2A Delay profile #1
            # value_vertex_peak: the value of each peak
            # locs_vertex_peak: index of the peak in each 
            
            threshold = 0

            # find location of largest peak in each vertex
            if locs_vertex_peak.size == 0: # if no peak
                delay_mat[lj,li] = np.NaN; # nan

            elif locs_vertex_peak.shape[0] > 1: # if more than one peaks exist
                value_max = np.max(value_vertex_peak) # pick just one peak that has the highest value
                idx_max = np.argmax(value_vertex_peak) 

                if value_max <= threshold: # if the peak value is lower than the threshold value,
                    delay_mat[lj,li] = np.NaN; # nan
                else: # if the peak value is above the threshold value,
                    delay_mat[lj,li] = locs_vertex_peak[idx_max] - (locs_global_peaks[li] - locs_seg_edge[li]) # save the index of the peak (as a relative index to the global peak)

            elif value_vertex_peak > threshold: # one peak, and above the threshold
                delay_mat[lj,li] = locs_vertex_peak[0] - (locs_global_peaks[li] - locs_seg_edge[li]) # save the index of the peak (as a relative index to the global peak)
            else :
                delay_mat[lj,li] = np.NaN

    ## applying SVD on delay profiles to calculate the principal delay profile: pd
      
    # fill nan
    delay_mat_filtered = np.array([])
    for li in np.arange(delay_mat.shape[1]): # for the number of segments
        tmp_seg_delay = delay_mat[:,li] # all vertices of one segment
        if np.sum(np.isnan(tmp_seg_delay.astype(np.float))) < delay_mat.shape[0] * 0.2: # if the number of nan vertices in the segment is lower than 20%,
            if delay_mat_filtered.shape[0] == 0:
                delay_mat_filtered = inpaintNans(tmp_seg_delay)
            else:
                delay_mat_filtered = np.concatenate((delay_mat_filtered, inpaintNans(tmp_seg_delay)),axis=1) # fill nan with neighbor average and save it to the delay matrix

    # apply svd
    U,S,V = np.linalg.svd(delay_mat_filtered, full_matrices = True)

    PD1 = U[:, 1] # PD profile: the first column of U, which explains the largest variance of X

    # variance explained for each components
    temp = np.diag(S)
    var_exp = temp**2/np.sum(temp**2)
    print(var_exp)
    
    pg1 = loadMat(path+"/pg1","pg1") # principle gradient 1 (vertices * 1)
    pg1 = pg1.flatten()
    r_btw_pd1_pg1 = np.corrcoef(PD1, pg1)[0,1]
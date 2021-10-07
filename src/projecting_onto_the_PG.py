import numpy as np
from scipy import io
from scipy.io import savemat
from scipy.signal import find_peaks
import pandas as pd
import matplotlib.pyplot as plt
from function import *

if __name__ == '__main__':

    path ="/home/cjh/skku/2021_second_semester/FMRI/project/Data"
   
    pg1 = loadMat(path+"/pg1","pg1") # principle gradient 1 (vertices * 1)
    epi = loadMat(path+"/epi1msk","epi1msk") # spatially & temporally smoothed -> zscored data (vertices * time points)s

    ##

    global_mean_sig = epi.mean(axis=0) # calculate the global mean of input data - Figure 1B middle graph
    
    pg1_70bins = np.digitize(pg1, prctile(pg1,0,100/70,100)) # divide pg1 vertices into 70 bins (vertices * 1 - the value in each row indicates which bin the vertex belongs to)
    pg1_70bins = pg1_70bins-np.min(pg1_70bins)+1 # make sure that pg1_70bins values are in 1~70
    time_pos = grpstats(epi, pg1_70bins) # % the time-position plot (time_pos) along the brain map (postion(70) * time(1200)) - Figure 1B upper graph

    ## find local peak at each position in each segment: idx_peak_seg

    locs_seg_edge = find_peaks(-1 * global_mean_sig)[0] # cut the global mean signal into time segments - Figure 1B middle graph
    # locs_seg_edge: location of negative peaks - indicate edge time points of time segments

    ## find delay of local peaks relative to global peak at each position: idx_peak_seg

    idx_peak_seg = np.zeros((time_pos.shape[0], locs_seg_edge.shape[0]-1),object)  # position * segment

    for li in np.arange(locs_seg_edge.shape[0]-1): 
        tmp_seg = np.array([])
        tmp_seg = epi[:,locs_seg_edge[li]:locs_seg_edge[li+1]] # get one time segment

        for lj in np.arange(time_pos.shape[0]): 
            tmp_pos = tmp_seg[lj,:] # get one position
            value_pos_peak, locs_pos_peak = findPeaks(tmp_pos) # and find peak of that position - Figure 1C red line
                                                    # value_pos_peak: the value of each peak
                                                    # locs_pos_peak: index of the peak
            threshold = 0       

            # find a peak in each segment - Figure 1C red line

            if locs_pos_peak.size == 0: # no peak
                idx_peak_seg[lj,li] = np.NaN; 

            elif locs_pos_peak.shape[0] > 1: # more than one peaks exist
                value_max = np.max(value_pos_peak) # pick just one peak that has the highest value
                idx_max = np.argmin(value_pos_peak)

                if value_max <= threshold: # if the peak value is lower than the threshold value, make it nan
                    idx_peak_seg[lj,li] = np.NaN; 

                else: # save the index of the peak only if the value is above the threshold
                    idx_peak_seg[lj,li] = locs_pos_peak[idx_max]

            elif value_pos_peak > threshold: # one peak, and above the threshold -> save the index
                idx_peak_seg[lj,li] = locs_pos_peak[0]

            else:
                idx_peak_seg[lj,li] = np.NaN; 

    # calculate time-position correlation at each segment: r_time_pos (Figure 1C r value, Figure 1D Principal gradient(PG) graph)

    sz = np.diff(locs_seg_edge) # width of each segment
    r_time_pos = np.zeros((1,locs_seg_edge.shape[0]-1),object) # 1 * segments
    
    for ln in np.arange(locs_seg_edge.shape[0]-1):

        #pos = np.arange(70) # 70 positions
        x_ax = np.array([])
        x_ax = np.arange(0,(sz[ln]-1),(sz[ln]-1)/idx_peak_seg.shape[0]) # make an array of 70 bins
        
        if np.sum(np.isnan(idx_peak_seg[:,ln].astype(np.float))) > 14: # if more than 14 peak indices in each time point are nan,
            r_time_pos[:,ln] = np.NaN # correlation of that time point is nan

        else: # if not, calculate time-position correlation at each segment - Figure 1C r value   
            r_time_pos[:,ln] = np.corrcoef(x_ax, idx_peak_seg[:,ln].astype(np.float))[0,1]

    r_time_pos = r_time_pos.astype(np.float)
    r_time_pos = r_time_pos[~np.isnan(r_time_pos)]

    # histogram (Figure 1D Principal gradient(PG) graph)
    plt.hist(r_time_pos,bins=30)
    plt.show()
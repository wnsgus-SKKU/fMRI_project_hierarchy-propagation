import numpy as np
from scipy import io
from scipy.io import savemat
from scipy.signal import find_peaks
from function import *

if __name__ == '__main__':

    path ="/home/cjh/skku/2021_second_semester/FMRI/project/Data"
   

    epi = loadMat(path+"/epi1msk","epi1msk")

    # calculate delay profile across time segments: idx_dly

    gls_neg_pk,locs = findPeaks(epi.mean(axis=0)) # cut the global mean signal into time segments - Figure 2A upper graph

    idx_dly = np.zeros((size(epi,0),size(locs,1)-1),object) # vertices * segments

    for li in range(0,size(locs,1)-1):
        tmp_prin = np.array([]) 
        tmp_prin = epi[:,locs[li]:locs[li+1]] # get one time segment

        # locate max peaks in each bin/electrode
        for lj in (0,size(epi,1)):
            tmp2_prin = tmp_prin[lj,:] # get one vertex in one time segment
            pks_prin, a_prin = findPeaks(tmp2_prin) # and find peak of that vertex - Figure 2A Delay profile #1
            # pks_prin: the value of each peak
            # a_prin: index of the peak
            thre1 = 0

            # find location of largest peak in each vertex
            if a_prin.size == 0: # no peak
                idx_dly[lj,li] = np.NaN; 

            elif size(a_prin,1) > 1: # more than one peaks exist
                valmax1= np.max(pks_prin) # pick just one peak that has the highest value
                id_prin = np.argmin(pks_prin) 

                if valmax1 <= thre1:
                    idx_dly[lj,li] = np.NaN; # if the peak value is lower than the threshold value, make it nan
                    
                else:
                    idx_dly[lj,li] = a_prin[id_prin] # save the index of the peak only if the value is above the threshold

            elif pks_prin > thre1: # one peak, and above the threshold -> save the index
                idx_dly[lj,li] = a_prin[0]
            else :
                idx_dly[lj,li] = np.NaN

    ## applying SVD on delay profiles to calculate the principal delay profile: pd
        
    # fill nan
    seg_all = idx_dly
    seg_all2 = np.array([])
    for l in range(0,size(seg_all,1)):
        temp1 = seg_all[:,li]
        if np.sum(np.isnan(temp1.astype(np.float))) < size(seg_all,0) * 0.2:
            if seg_all2.size == 0:
                seg_all2 = inpaintNans(temp1)
            else:
                seg_all2 = np.concatenate((seg_all2,inpaintNans(temp1)),axis=0)

    delay_test3 = seg_all2 - np.tile(np.mean(seg_all2),(size(seg_all,1),1))
    U,S,V = np.linalg.svd(delay_test3, full_matrices = True)

    pd = U[:, 1]

    # variance explained for each components
    temp2 = np.diag(S)
    var_exp = temp2**2/np.sum(temp2**2)
    print(var_exp)
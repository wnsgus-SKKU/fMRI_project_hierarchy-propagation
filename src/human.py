import numpy as np
from scipy import io
from scipy.io import savemat
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from function import *

if __name__ == '__main__':

    path ="/home/cjh/skku/2021_second_semester/FMRI/project/Data"
   
    pg1 = loadMat(path+"/pg1","pg1") #principle gradient 1 (vertices * 1)
    epi1msk = loadMat(path+"/epi1msk","epi1msk") # spatially & temporally smoothed -> zscored data (vertices * time points)s

    ##

    gs_LR1 = epi1msk.mean(axis=0) # calculate the global mean of input data - Figure 1B middle graph
    
    pgd1 = np.digitize(pg1, prctile(pg1,0,100/70,100)) # divide pg1 vertices into 70 bins (vertices * 1 - the value in each row indicates which bin the vertex belongs to)
    pgd1 = pgd1-np.min(pgd1)+1 # make sure that pgd1 values are in 1~70
    tw1 = grpstats(epi1msk,pgd1) # % the time-position plot (tw1) along the brain map (postion(70) * time(1200)) - Figure 1B upper graph

    ## find delay of local peaks relative to global peak at each position: idx_tem_prin

    gls_neg_pk,locs = findPeaks(-1 * gs_LR1) # cut the global mean signal into time segments - Figure 1B middle graph
    # locs: location of negative peaks - indicate edge time points of time segments

    idx_tem_prin = np.zeros((size(tw1),size(locs,1)-1),object)  # position * segment

    for li in range(0,size(locs,1)-1): 
        tmp_prin = np.array([])
        tmp_prin = tw1[:,locs[li]:locs[li+1]]

        for lj in range(1,size(tw1)): 
            tmp2_prin = tmp_prin[lj-1,:] # get one position
            pks_prin, a_prin = findPeaks(tmp2_prin) # and find peak of that position - Figure 1C red line
                                                    # pks_prin: the value of each peak
                                                    # a_prin: index of the peak
            thre1 = 0       

            # find a peak in each segment - Figure 1C red line

            if a_prin.size == 0: # no peak
                idx_tem_prin[lj,li] = np.NaN; 

            elif size(a_prin,1) > 1: # more than one peaks exist
                valmax1= np.max(pks_prin) # pick just one peak that has the highest value
                id_prin = np.argmin(pks_prin)

                if valmax1 <= thre1: # if the peak value is lower than the threshold value, make it nan
                    idx_tem_prin[lj,li] = np.NaN; 

                else: # save the index of the peak only if the value is above the threshold
                    idx_tem_prin[lj,li] = a_prin[id_prin]

            elif pks_prin > thre1: # one peak, and above the threshold -> save the index
                idx_tem_prin[lj,li] = a_prin[0]

            else:
                idx_tem_prin[lj,li] = np.NaN; 

    # calculate time-position correlation at each segment: rval_prin_2 (Figure 1D Principal gradient(PG) graph)

    sz = np.diff(locs) # width of each segment
    rval_prin_2 = np.zeros((1,size(locs,1)-1),object)
    for ln in range(0,size(locs,1)-1):

        x_ax = np.array([])
        x_ax = np.arange(0,(sz[ln]-1),(sz[ln]-1)/size(idx_tem_prin,0)) # divide each segment into 70 bins(-1 means not using signals on the segment edges)

        if np.sum(np.isnan(idx_tem_prin[:,ln].astype(np.float))) > 14:
            rval_prin_2[:,ln] = np.NaN

        else: # if not, calculate time-position correlation at each segment - Figure 1D Principal gradient(PG) graph   
            rval_prin_2[:,ln] = np.corrcoef(x_ax.astype(np.float),idx_tem_prin[:,ln].astype(np.float))[0,1]

    rval_prin_2_tmp = rval_prin_2.astype(np.float)
    rval_prin_2_tmp = rval_prin_2_tmp[~np.isnan(rval_prin_2_tmp)]

    # ploting 
    plt.hist(rval_prin_2_tmp,bins=10)
    plt.show()

    # save generated matrices

    mat = {'locs' : locs,'tw1':tw1,'idx_tem_prin':idx_tem_prin,'rval_prin_2':rval_prin_2}
    io.savemat(path+"/results.mat",mat)

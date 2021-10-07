import numpy as np
from scipy import io
from scipy.io import savemat
from scipy.signal import find_peaks
from function import *

if __name__ == '__main__':

    path ="/home/cjh/skku/2021_second_semester/FMRI/project/Data"
   
    pg1 = loadMat(path+"/pg1","pg1")
    epi1msk = loadMat(path+"/epi1msk","epi1msk")
    gs_LR1 = epi1msk.mean(axis=0)
    
    pgd1 = np.digitize(pg1, prctile(pg1,0,100/70,100))
    pgd1 = pgd1-np.min(pgd1)+1
    tw1 = grpstats(epi1msk,pgd1)

    gls_neg_pk,locs = findPeaks(-1 * gs_LR1)
    idx_tem_prin = np.zeros((size(tw1),size(locs,1)-1),object) # need to check range 

    for li in range(0,size(locs,1)-1): # need to check range 
        tmp_prin = np.array([])
        tmp_prin = tw1[:,locs[li]:locs[li+1]]

        for lj in range(1,size(tw1)):
            tmp2_prin = tmp_prin[lj-1,:]
            pks_prin, a_prin = findPeaks(tmp2_prin)
            thre1 = 0

            if a_prin.size == 0:
                idx_tem_prin[lj,li] = np.NaN; 

            elif size(a_prin,1) > 1:
                valmax1= np.max(pks_prin)
                id_prin = np.argmin(pks_prin)

                if valmax1 <= thre1:
                    idx_tem_prin[lj,li] = np.NaN; 
                else:
                    idx_tem_prin[lj,li] = a_prin[id_prin]

            elif pks_prin > thre1:
                idx_tem_prin[lj,li] = a_prin[0]

            else:
                idx_tem_prin[lj,li] = np.NaN; 

    sz = np.diff(locs)
    rval_prin_2 = np.zeros((1,size(locs,1)-1),object)
    for ln in range(0,size(locs,1)-1): # need to check range 

        x_ax = np.array([])
        x_ax = np.arange(0,(sz[ln]-1),(sz[ln]-1)/size(idx_tem_prin,0))

        if np.sum(np.isnan(idx_tem_prin[:,ln].astype(np.float))) > 14:
            rval_prin_2[:,ln] = np.NaN
        else:   
            rval_prin_2[:,ln] = np.corrcoef(x_ax.astype(np.float),idx_tem_prin[:,ln].astype(np.float))[0,1]
    print(rval_prin_2)
    # mat = {'locs' : locs,'tw1':tw1,'idx_tem_prin':idx_tem_prin,'rval_prin_2':rval_prin_2}
    mat = {'locs' : locs,'tw1':tw1}
    io.savemat(path+"/results.mat",mat) # Error in saving because some arrays are empty

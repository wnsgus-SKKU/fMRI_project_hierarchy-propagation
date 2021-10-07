import numpy as np
from scipy import io
from scipy.io import savemat
from scipy.signal import find_peaks
from function import *

if __name__ == '__main__':

    path ="/home/cjh/skku/2021_second_semester/FMRI/project/Data"
   

    epi = loadMat(path+"/epi1msk","epi1msk")
    gls_neg_pk,locs = findPeaks(epi.mean(axis=0))
    print(size(locs,1))
    idx_dly = np.zeros(size(epi,1),size(locs,2)-1)
    print(idx_dly.size)
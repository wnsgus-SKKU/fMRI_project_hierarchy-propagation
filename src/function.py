import numpy as np
from scipy import io
from scipy.io import savemat
from scipy.signal import find_peaks
import pandas as pd

def loadMat(path,name):
    mat_file1 = io.loadmat(path)
    data = mat_file1[name]
    return data


def prctile(data,start,interval,end):
    prctile = np.array([])
    percentage = start
    while(percentage <= 100):
        prctile = np.append(prctile,np.percentile(data,percentage))
        percentage += interval
    return prctile


def grpstats(epi1msk,pgd1):
    max = np.max(pgd1)
    min = np.min(pgd1)
    result = np.empty((max-min+1 ,epi1msk.shape[1]),object)
    nVtxInBin = np.zeros(max-min+1) # the number of vertices in each bin (1->70)
    for i in range(pgd1.size):
        nVtxInBin[pgd1[i][0]-1] +=1
        if(result[pgd1[i][0]-1][0] == None):
            temp = epi1msk[i,:]
        else:
            temp = result[pgd1[i][0]-1,:] + epi1msk[i,:]
        result[pgd1[i][0]-1,:] = temp
    for j in range(max-min+1):
        result[j,:] = result[j,:]/ nVtxInBin[j]
    return result


def size(data,axis=0):
    if data.ndim == 1:
        if axis == 0:   
            return 1
        elif axis == 1:
            return data.size
    else :
        return data.shape[axis]


def findPeaks(data):
    locs= find_peaks(data)[0]
    peak = np.array([])
    for i in range(locs.size):
        peak = np.append(peak,data[locs[i]]) 
    return peak,locs


def inpaintNans(array):
    if np.isnan(array[0]):
        array[0] = array[1]
    if np.isnan(array[-1]):
        array[-1] = array[-2]
    df = pd.DataFrame(array.astype(np.float)).interpolate()
    return df.to_numpy()
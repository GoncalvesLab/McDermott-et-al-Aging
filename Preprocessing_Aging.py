# -*- coding: utf-8 -*-
"""
Created on Tue May 31 11:11:16 2022

@author: tiagolab
"""


import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
import os
import numpy as np
import glob
#from scipy import signal
#from scipy.signal import lfilter
#from scipy import stats
#import peakutils
#from scipy.signal import find_peaks
import h5py 
from skimage.transform import resize
#import xml.etree.ElementTree as ET
#import seaborn as sns
#import cmath
#import csv


def extract_zones(textures_trace, t_voltages = [0.6, 1.2, 1.8, 2.4]):
        
        delay = 15000
        transitions = [0]
        textures = []
        mean_values = []
        thre_values = []
        t_voltages = np.array(t_voltages)
        
        while True:
        
            last_transition = transitions[-1]
            
            ref = textures_trace[last_transition+delay:last_transition+delay+1000, 0]
            if ref.size == 0:
                #print("Last transition at the end, assuming next texture")
                textures.append((textures[-1] + 1) % len(t_voltages))
                #print(f"There are {len(transitions)} zones.")
                break
            current_mean = np.mean(ref)
            textures.append(np.argmin(np.abs(t_voltages - current_mean)))
            mean_values.append(current_mean)
            current_threshold = 0.8 * (np.max(ref) - np.min(ref))
            thre_values.append(current_threshold)
            outside_this_texture = np.where(np.abs(textures_trace[last_transition+delay:, 0] - current_mean) > current_threshold)[0]
            if outside_this_texture.size == 0:
                #print(f"There are {len(transitions)} zones.")
                break
            transitions.append(outside_this_texture[0] + last_transition+delay)
        # Return transitions and textures as a tuple of two lists:
        return np.array(transitions), np.array(textures)

filelist = glob.glob(r"filelist") #"\\data.einsteinmed.org\users\Maria Frechou\Ex7\*EE*\mossyconc")
brain_region = "DG"

average_response_by_cell = []
average_amplitud_by_cell = []

for filename in filelist:

    '''
    LOAD FILES
    '''
    mouse_dir = filename
    cells = np.load(os.path.join(mouse_dir, 'suite2p\plane0\iscell.npy'))
    ftraces = np.load(os.path.join(mouse_dir, 'suite2p\plane0\F.npy'))
    fneu = np.load(os.path.join(mouse_dir, 'suite2p\plane0\Fneu.npy'))
    deconvtraces = np.load(os.path.join(mouse_dir, 'suite2p\plane0\spks.npy'))
    
    # Frames per Second:
    fps = 15.253
    
    #Total frames of video
    tf = len(deconvtraces[1])
    
    '''
    SELECT USEFUL ROIs
    '''
    #select good ROI from iscell and get only the good traces (transpose values because if not it plots them the other way around):
    ftraces_good = ftraces[cells[:,0]==1]
    ftraces_noneu = ftraces - fneu 
    ftraces_good_noneu = ftraces_noneu[cells[:,0]==1]
    dtraces_good = deconvtraces[cells[:,0]==1]
    
    dataf_mean = np.mean(ftraces_good, axis=1)[:, None]
    dataf2 = ftraces_good - dataf_mean
    # Threshold of "dataf2"
    threshold2 = (np.std(ftraces_good, axis=1) * 2)[:, None]
    threshold = threshold2 + dataf_mean
    # Fluorescence signal above mean+2*sigma:
    above_threshold = (dataf2 > threshold2) & (dtraces_good > 0) #binary data == (data2 > 0). True are cells with activity above the threshold
    
    data2 = dtraces_good.copy()
    # Use calculated filter on "Suite2P" cleaned signal:
    # set all values to 0 where NOT above_theshold:
    data2[~above_threshold] = 0
    
    # Plot each cleaned&re-thresholded trace in one figure:
    
    dataf_mean = np.mean(ftraces_good, axis=1)[:, None]
    threshold2 = (np.std(ftraces_good, axis=1) * 2)[:, None]
    threshold = threshold2 + dataf_mean
    
    np.save(os.path.join(mouse_dir, 'cleaned_traces_binary.npy'), above_threshold) #number of active cells (yes/no type of data)
    np.save(os.path.join(mouse_dir, 'cleaned_traces_data.npy'), data2) #actual values of active cleaned cells 
       
#%%
        
    '''
    spatial information
    '''
    
    '''
    load files
    '''
    
    #using the length of the video folder
    a = os.path.split(mouse_dir)[0]
    b = glob.glob(os.path.join(a, brain_region))
    lenexp = []
    for file in b:
        c = glob.glob(os.path.join(file, "*.tif"))
        lens = len(c)
        lenexp.append(lens)
    
    print(lenexp)
    
    #spatial_dir = glob.glob(os.path.join(mouse_dir, 'SyncData*'))
    spatial_dir = glob.glob(os.path.join(a, brain_region, 'SyncData*')) 
    
    '''
    clean files
    '''
    normedlist_all = []
    trans1_all = []
    indexes_final = []
    count = 0
    
    mset_all = []
    tset_all = []
    
    for index in range(len(spatial_dir)):
        
        f = h5py.File(os.path.join(spatial_dir[index], 'Episode001.h5'), 'r+')
        list(f.keys())
        a_group_key = list(f.keys())[0]
    
        mset = np.array(f['AI/movement']) #converts the h5 file to numpy array
        tset = np.array(f['AI/texture'])
        frame_out = np.array(f['DI/FrameOut'])
        
        
        non_zero_indexes = np.nonzero(frame_out) #find non zero indexes of frame_out. The first and last index are the relevant ones
        t = np.arange(mset.size) / 30000 # 
    
        mset = mset[non_zero_indexes[0][0]:non_zero_indexes[0][-1]] #cuts mset and leaves only the part that we need
        tset = tset[non_zero_indexes[0][0]:non_zero_indexes[0][-1]] #cuts tset and leaves only the part that we need 
        frame_out = frame_out[non_zero_indexes[0][0]:non_zero_indexes[0][-1]] #cuts frame_out and leaves only the part that we need 
        t = t[non_zero_indexes[0][0]:non_zero_indexes[0][-1]] #cuts t

        # Unpack tuple into two lists: get transitions and zones
        transitions, zones = extract_zones(tset)
        
        mset_all.extend(mset)
        tset_all.extend(tset)
        

#%%
        '''
        new approach
        '''
        cum = []
        # resample mset 
        res1 = resize(mset, (lenexp[index],1), order=1, preserve_range=True)
        res11 = res1 - 1.372462080290498
        res11[res11<0] = 0
        #resize transitions
        trans = (transitions/len(mset))*lenexp[index]
        #convert floats in trans to integers so you can use them to slice
        trans1 = [int(i) for i in trans]
        #do the cumsum for every zone
        for k in range(1, len(trans1)-2): # this selects all values from the transitions except the first and last
            cum.append(np.cumsum(res11[trans1[k]:trans1[k+1]]))
        
        [l.tolist() for l in cum]
        #normalize
        normed = []
        for i in range(len(cum)):
            if zones[i] == 3:
                normed.append((cum[i]/np.max(cum[i]))*25 + 75)
            if zones[i] == 2:
                normed.append((cum[i]/np.max(cum[i]))*25 + 50)
            if zones[i] == 1:
                normed.append((cum[i]/np.max(cum[i]))*25 + 25)
            if zones[i] == 0:
                normed.append((cum[i]/np.max(cum[i]))*25)

        normed = [l.tolist() for l in normed]
        normed = np.array(normed)

        flatten = lambda l: [item for sublist in l for item in sublist]
        normedlist = flatten(normed)
        normedlist = np.array(normedlist)

        normedlist_all.extend(normedlist)
        trans1_all.append(trans1)

        lens = list(range(lenexp[index] + 1))
        vid = lens[trans1[1]:trans1[-2]]
        video = [x + count for x in vid]
        count += lenexp[index] + 1
        indexes_final.extend(video)


    final_fsubneu = ftraces_good_noneu[:,indexes_final]  #all traces cut without including the 1st and last trans.
    
    final_f = ftraces_good[:,indexes_final]
    
    fneu_good = fneu[:,indexes_final]
    
    final_spks = dtraces_good[:,indexes_final]

    np.save(os.path.join(mouse_dir, 'positions.npy'), normedlist_all)

    np.save(os.path.join(mouse_dir, 'fluorescencesubneuropil.npy'), final_fsubneu)  #with neuropil subtracted

    np.save(os.path.join(mouse_dir, 'spks_final.npy'), final_spks)
    
    np.save(os.path.join(mouse_dir, 'fluorescence.npy'), final_f)
    
    np.save(os.path.join(mouse_dir, 'fneuropil.npy'), fneu_good)
    
    np.save(os.path.join(mouse_dir, 'mset'), mset_all)
    np.save(os.path.join(mouse_dir, 'tset'), tset_all)
    np.save(os.path.join(mouse_dir, 'frame_out'), frame_out)
    
    print('Processed file:' + filename)
    
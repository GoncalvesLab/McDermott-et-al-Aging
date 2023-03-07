# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 13:51:51 2022

@author: M. Agustina Frechou
"""

'''

Bootstrapping code and preprocessing to analyze tuning data
 
'''

import os
import glob
import scipy.io
import numpy as np
import random
from random import choices
import csv

#%%

'''
Functions
'''


def get_fi(filelist):
    
    activity_all = []
    for filename in filelist:
        
        load_file = open(os.path.join(filename, 'FI per cell.npy'), 'r')
        activity = load_file.read()
        activity1= activity.splitlines()
        activity2 = [eval(i) for i in activity1]
        activity_all.append(activity2)
    return activity_all



def get_activity(filelist):
    
    activity_all = []
    for filename in filelist:
        
        load_file = open(os.path.join(filename, 'activity_area_under_curve_normed.npy'), 'r')
        activity = load_file.read()
        activity1= activity.splitlines()
        activity2 = [eval(i) for i in activity1]
        activity_all.append(activity2)
    return activity_all

def get_tuning(filelist):
    
    tuning_final = []
    for filename in filelist:
        
        file = scipy.io.loadmat(os.path.join(filename, 'CellTuning.mat'))
        tuning = file['CellTuning']
        
        tuningg = []
        for ii in range(len(tuning)):
            tuningg.append(tuning[ii][2][0][0])
        
        #remove NaNs
        
        tuning_cleaned = [x for x in tuningg if str(x) != 'nan']
        
        tuning_final.append(tuning_cleaned)
    
    return tuning_final

def bootstrap(group1, group2, groupALL, bootstraps=100000):

    #you need list containing the lists of cells for each mice for each experimental group and one list containing all (both groups in one list)
    
    #find total number of cells per group
    total_cells_group1 = 0
    for i in group1:
        total_cells_group1 +=len(i)
        
    total_cells_group2 = 0
    for i in group2:
        total_cells_group2 +=len(i)
        
    #find total amount of mice per group    
    
    mice_group1 = len(group1)
    mice_group2 = len(group2)
    
    #to determine GROUP1NULL and GROUP2NULL
    
    #sample with replacement animal number a total equal to the total amount of cells in group1 and group2.
    
    DIFFNULL_mean_all = []
    DIFFNULL_median_all = []
    
    for ii in range(bootstraps):
    
        y_group1 = choices(range(mice_group1), k=total_cells_group1)
        
        y_group2 = choices(range(mice_group2), k=total_cells_group2)
        
        #Sample a cell value at random with replacement from the animals chosen in y_group1 and y_group2
        
        GROUP1NULL = []
        for i in y_group1:
            GROUP1NULL.extend(random.choices(groupALL[i], k=1))
         
        GROUP2NULL = []
        for i in y_group2:
            GROUP2NULL.extend(random.choices(groupALL[i], k=1))
        
        #calculate mean and median for each NULL group and compute the difference (should be zero because of the way in which it was constructed)
        
        DIFFNULL_mean = np.mean(GROUP2NULL) - np.mean(GROUP1NULL)
        
        DIFFNULL_median = np.median(GROUP2NULL) - np.median(GROUP1NULL)
    
        DIFFNULL_mean_all.append(DIFFNULL_mean)
        DIFFNULL_median_all.append(DIFFNULL_median)
    
    #calculate the measured mean of each group 
    
    mean_group1 = np.mean([item for sublist in group1 for item in sublist])
    mean_group2 = np.mean([item for sublist in group2 for item in sublist])
    
    #calculate difference
    
    meanDIFF_group1group2 = mean_group2 - mean_group1
    
    #p value
    
    if meanDIFF_group1group2 < 0:
    
        count = 0
        for i in DIFFNULL_mean_all:
            if i < meanDIFF_group1group2:
                count += 1
    
        p_value = count/len(DIFFNULL_mean_all)
        return [p_value, mean_group1, mean_group2]
    else:
        count = 0
        for i in DIFFNULL_mean_all:
            if i > meanDIFF_group1group2:
                count += 1
    
        p_value = count/len(DIFFNULL_mean_all)
        return [p_value, mean_group1, mean_group2]
    
def load_csvs(file):
    #Gets all values for Amp, Width and R2 > 50%
    with open(file, "r") as read_obj:
      
        # Return a reader object which will iterate over lines in the given csv file
        csv_reader = csv.reader(read_obj)
      
        # convert string to list
        list_of_csv = list(csv_reader)
        
    listcsv = [[ele for ele in sub if ele != ''] for sub in list_of_csv]
    
    for i in listcsv:
        del i[0]
        
    listcsv2 = [[eval(ele) for ele in sub] for sub in listcsv]
    
    return listcsv2

#%%


'''
Load files
'''

filelist1 = glob.glob(r"directory_group1")
filelist2 = glob.glob(r"directory_group2")


'''
Tuning Bootstrap
'''

''' Preprocessing '''

# 1) Save first experimental group in variable "group 1".

group1 = get_tuning(filelist1)

# 2) Save second experimental group in variable "group 2".

group2 = get_tuning(filelist2)

# 3) Save both experimental groups combined in variable groupALL

groupALL = group1 + group2


''' BOOTSTRAP '''
# 4) Run bootstrap and print p_value

print("Tuning p value: ", bootstrap(group1, group2, groupALL))


'''
Activity Bootstrap
'''

''' Preprocessing '''

# 1) Save first experimental group in variable "group 1".

group1 = get_activity(filelist1)

# 2) Save second experimental group in variable "group 2".

group2 = get_activity(filelist2)

# 3) Save both experimental groups combined in variable groupALL

groupALL = group1 + group2


'''
BOOTSTRAP
'''

# 4) Run bootstrap and print p_value

print("Activity p value: ", bootstrap(group1, group2, groupALL))

    
'''
FI Bootstrap
'''

''' Preprocessing '''

# 1) Save first experimental group in variable "group 1".

group1 = get_fi(filelist1)

# 2) Save second experimental group in variable "group 2".

group2 = get_fi(filelist2)

# 3) Save both experimental groups combined in variable groupALL

groupALL = group1 + group2


'''
BOOTSTRAP
'''
# 4) Run bootstrap and print p_value

print("FI p value: ", bootstrap(group1, group2, groupALL))
    

'''
subsets of data e.g. matched cells
'''

file1 = r"directory\group1.csv"
file2 = r"directory\group2.csv"

 
group1 = load_csvs(file1)
group2 = load_csvs(file2)
groupALL = group1 + group2

print("p value, mean group 1 , mean group 2: ", bootstrap(group1, group2, groupALL))
    
    

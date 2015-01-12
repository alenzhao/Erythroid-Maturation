#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon 12 Jan 2015 11:13:18 AM PST 
Measures the experimental differences betten spade_gated cases
@author: ngdavid
"""

from SPADE_Gating import SpadeGating
from scipy.spatial.distance import cdist
import numpy as np
from itertools import product
import pandas as pd
import matplotlib.pyplot as plt
from munkres import Munkres
from pyemd import emd
from Analysis_Variables import *

output = {}
for (Normal,Base) in product(NormNames,AbnNames):
    print("Comparing {} with {}".format(Normal,Base))
    t = pd.read_table(table_of_contents)
    fl = t.loc[t['ShortName'] == Normal,'Folder'].values[0]
    A_Filename = '/home/ngdavid/FCS_Data/MDS_Plates/{}/SynData_v3.hdf5'.format(fl)
    fl = t.loc[t['ShortName'] == Base,'Folder'].values[0]
    B_Filename = '/home/ngdavid/FCS_Data/MDS_Plates/{}/SynData_v3.hdf5'.format(fl)
    try:
        A=SpadeGating(filename=A_Filename,
                         gates=Gates,
                         parameters=Parameters,
                         source=Source,
                         target=Target,
                         removal_list=Removal_list,
                         isotypes=Isotypes)
        
        B=SpadeGating(filename=B_Filename,
                         gates=Gates,
                         parameters=Parameters,
                         source=Source,
                         target=Target,
                         removal_list=Removal_list,
                         isotypes=Isotypes)
        
        plate_plate=pd.Series()
        
        for antigen in A.data.columns:
            A_dat= A.data[A.erythroid_mask]
            idx = np.random.choice(A_dat.index,size=1000,replace=False)
            A_comp = A_dat[antigen][idx].values.astype(np.float)
            B_dat = B.data[B.erythroid_mask]
            idx = np.random.choice(B_dat.index,size=1000,replace=False)
            B_comp = B_dat[antigen][idx].values.astype(np.float)
            
            print antigen
            sample = 1000
            dmtx=np.zeros((sample,sample),dtype=np.float)
            for p in product(range(sample),range(sample)):
                dmtx[p[0]][p[1]]=abs(A_comp[p[0]]-B_comp[p[1]])
            
            cost = emd(A_comp,B_comp,dmtx)
            plate_plate[antigen]=cost
            
        output["{}vs{}".format(Normal,Base)] = plate_plate
    except:
        output["{}vs{}_failed".format(Normal,Base)] = plate_plate
    writing_out = pd.DataFrame(output)
    writing_out.to_csv("/home/ngdavid/Desktop/Ubuntu_Dropbox/Interplate Comparison/EMD_Analysis/complete.csv")

c = [i for i in writing_out.columns if i not in ['Untitled: 0']]
writing_out['mean'] = writing_out.mean(axis=1)
writing_out['rms'] = sqrt(sum(output[c]*output[c],axis=1))


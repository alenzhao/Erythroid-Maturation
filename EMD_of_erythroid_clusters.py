#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon 12 Jan 2015 11:13:18 AM PST 
Measures the experimental differences betten spade_gated cases
@author: ngdavid
"""

from ERY_subroutines.SPADE_Gating import SpadeGating
from ERY_subroutines.clustering import kmeans_clustering
from scipy.spatial.distance import cdist
import numpy as np
from itertools import product
import pandas as pd
import matplotlib.pyplot as plt
import os

from pyemd import emd
from Analysis_Variables import *

output = {}

for (Normal,Base) in product(NormNames,NormNames+AbnNames):
    print("Comparing {} with {}".format(Normal,Base))
    t = pd.read_table(table_of_contents)
    fl = t.loc[t['ShortName'] == Normal,'Folder'].values[0]
    A_Filename = os.path.join(root_dir,"{}/SynData_v3.hdf5".format(fl))
    fl = t.loc[t['ShortName'] == Base,'Folder'].values[0]
    B_Filename = os.path.join(root_dir,"{}/SynData_v3.hdf5".format(fl))

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
    antigens_to_test = [i for i in A.data.columns if i not in backbone_antigens]
     
    for antigen in antigens_to_test:
        atg_list = backbone_antigens+[antigen]
        A_dat= A.data[atg_list][A.erythroid_mask]
        B_dat = B.data[atg_list][B.erythroid_mask]
        A_cluster = kmeans_clustering(A_dat,n_clusters=40,n_jobs=4)
        
        B_cluster = kmeans_clustering(B_dat,n_clusters=40,n_jobs=4)
        
        distance_matrix = cdist(A_cluster.output[atg_list],
                                B_cluster.output[atg_list],
                                'cityblock')
        
        XA = A_cluster.output["cluster_size"].values.astype(np.double)
        XA = XA/np.sum(XA)
        XB = B_cluster.output["cluster_size"].values.astype(np.double)
        XB = XB/np.sum(XB)
        cost = emd(XA,XB,distance_matrix)
        plate_plate[antigen]=cost
        print("{} has emd of {}".format(antigen,cost))
    output["{}vs{}".format(Normal,Base)] = plate_plate

    writing_out = pd.DataFrame(output)
    writing_out.to_csv("/home/ngdavid/Desktop/Ubuntu_Dropbox/Interplate Comparison/EMD_Analysis/BackboneSimulation.csv")

c = [i for i in writing_out.columns if i not in ['Untitled: 0']]
writing_out['mean'] = writing_out.mean(axis=1)
writing_out['rms'] = np.sqrt(sum(output[c]*output[c],axis=1))

writing_out.to_csv("/home/ngdavid/Desktop/Ubuntu_Dropbox/Interplate Comparison/EMD_Analysis/BackboneSimulation.csv")


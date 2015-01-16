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
from multiprocessing import pool
from multiprocessing.dummy import Pool as ThreadPool 

from pyemd import emd
from Analysis_Variables import *

output = {}

kmeans_jobs = 1
num_workers = 52
clusters = 40

def worker(antigen):
    print antigen
    atg_list = backbone_antigens+[antigen]
    A_dat= A.data[atg_list][A.erythroid_mask]
    B_dat = B.data[atg_list][B.erythroid_mask]
    A_cluster = kmeans_clustering(A_dat,n_clusters=clusters,n_jobs=kmeans_jobs)
    
    B_cluster = kmeans_clustering(B_dat,n_clusters=clusters,n_jobs=kmeans_jobs)
    
    distance_matrix = cdist(A_cluster.output[atg_list],
                            B_cluster.output[atg_list],
                            'cityblock')
    
    XA = A_cluster.output["cluster_size"].values.astype(np.double)
    XA = XA/np.sum(XA)
    XB = B_cluster.output["cluster_size"].values.astype(np.double)
    XB = XB/np.sum(XB)
    cost = emd(XA,XB,distance_matrix)

    print("{} has emd of {}".format(antigen,cost))
    return cost
    
        

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
    
    pool = ThreadPool(num_workers)
    costs = pool.map(worker,antigens_to_test)
    pool.close()
    pool.join()        
    output["{}vs{}".format(Normal,Base)] = dict(zip(antigens_to_test,costs))

    writing_out = pd.DataFrame(output)
    writing_out.to_csv(output_file)

c = [i for i in writing_out.columns if i not in ['Untitled: 0']]
#writing_out['mean'] = writing_out.mean(axis=1)
#writing_out['rms'] = np.sqrt(sum(output[c]*output[c],axis=1))

#writing_out.to_csv(output_file)


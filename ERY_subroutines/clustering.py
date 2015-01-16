"""
This file describes wrappers for various clustering functions
"""
from sklearn.cluster import KMeans
import pandas as pd
import numpy as np

class kmeans_clustering(object):

    def __init__(self,data,**kwargs):
        """
        Wrapper for kmeans clustering
        """
        clusterer = KMeans(**kwargs)
        clusterer.fit_predict(data)
        self.centroid = clusterer.cluster_centers_
        self.labels = clusterer.labels_
        self.output = pd.DataFrame(columns=data.columns.tolist()+['cluster_size'])
        
        for k in np.unique(self.labels):
            self.output.loc[k]=data[self.labels==k].median(axis=0)
            self.output.loc[k]['cluster_size'] = np.sum(self.labels==k)
            

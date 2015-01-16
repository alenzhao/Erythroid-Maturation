# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 11:03:16 2014
Defines a class for spade gating
@author: ngdavid
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
import networkx as nx

from sklearn.cluster import KMeans, MiniBatchKMeans
from scipy.interpolate import interp1d
from scipy.spatial.distance import pdist, squareform,cdist
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.signal import wiener
from matplotlib.path import Path

class SpadeGating:
    def __init__(self, filename, source, target, gates, isotypes, \
                       removal_list, parameters):
        """
        filename = string of names
        source = source immunophenotype coords
        targe = end stage immunophenotype coords
        isotypes = list out isotype controls
        removal list = list of antigens that are worthless
        parameters = dictionary to set various parameters
        """
        self.data, self.Antigen_list = self._reconstitute_data(filename,isotypes)    # load an antigen list and data files
        self.data = self._apply_gates(gates)   # apply gates
        self.centroids, self.labels_mask, self.Antigen_list_filtered = self._full_clustering(parameters,isotypes,removal_list)    #apply clustering
        self.centroids, self.kept_indicies = self._filter_centroids(parameters) #remove clusters

        self.source_graph_idx,self.target_graph_idx = self._find_ends(parameters,source,target) # generates the indicies of the source and sinks

        self.erythroid_mask = self._path_gating(parameters) # generates a mask based on the selected path
        self.erythroid_cluster_path = self.erythroid_clusters(parameters)
        
    def _reconstitute_data(self,Filename,isotypes):
        print "Loading Data"
        hdf_handle = h5py.File(Filename, "r")
        if 'v3' in Filename:
            print "Version 3"
            antigen_list = set(hdf_handle['Antigens'])
            data_header = list(antigen_list-set(isotypes))
            data = [hdf_handle['Data'][a][:][np.newaxis] for a in data_header]
            data = np.concatenate(data,axis=0).T
            data = pd.DataFrame(data = data,columns = data_header)
        elif 'v4' in Filename:
            print "Version 4"
            data = pd.DataFrame(data = hdf_handle['data\values'],
                                index = hdf_handle['data\index'],
                                columns = hdf_handle['data\columns'])
            data_header = pd.DataFrame(data = hdf_handle['parameters\values'],
                                       index = hdf_handle['parameters\index'],
                                       columns = hdf_handle['parameters\columns'])
        hdf_handle.close()
        print "Data Frame loaded"
        return data,data_header 

    def erythroid_clusters(self,parameters):
        dist_metric=parameters['distance_metric']
        erythroid_clusters=parameters['erythroid_clusters']
        ery_cluster_antigens=self.Antigen_list

        ery_model=MiniBatchKMeans(n_clusters=erythroid_clusters)
        erythroid_DF=self.data[self.erythroid_mask]
        labels=ery_model.fit_predict(erythroid_DF[ery_cluster_antigens])
        columns=[ery_cluster_antigens+['cluster_size']+['centroid_disp']]
        centroids=pd.DataFrame(columns=columns)
        for k in np.unique(labels):
            centroids.loc[k]=erythroid_DF[labels==k].median(axis=0)
            centroids.loc[k]['cluster_size'] = np.sum(labels==k)
            centroids.loc[k]['centroid_disp'] = np.sum( cdist(centroids.loc[k][self.Antigen_list].values[np.newaxis],\
                                                        erythroid_DF[labels==k], \
                                                        dist_metric) \
                                                        /centroids.loc[k]['cluster_size'])

        dist_mtx=squareform(pdist(centroids[ery_cluster_antigens],dist_metric))
        G_ery=nx.Graph(dist_mtx)
        G_ery_mst=nx.minimum_spanning_tree(G_ery)

        all_paths=nx.all_pairs_dijkstra_path_length(G_ery_mst)
        max_path=0

        for l in range(erythroid_clusters):
            for m in range(erythroid_clusters):
                if all_paths[l][m] > max_path:
                    max_path=all_paths[l][m]
                    max_l=l
                    max_m=m

        if centroids['CD34'][max_l] > centroids['CD34'][max_m]:
            sub_path=nx.dijkstra_path(G_ery_mst,source=max_l,target=max_m)
        else:
            sub_path=nx.dijkstra_path(G_ery_mst,source=max_m,target=max_l)
        return centroids.iloc[sub_path]
        
       
    def stage_plot(self,stages,figsize=(20,14),title='default',**kwargs):
        plt.figure(figsize=figsize)
        plt.title(title)
        output_array=self._stage_averages(self.disp_array,stages,axes=('CD34','CD71 FITC'))
        plt.imshow(output_array.values,cmap=plt.cm.bwr,interpolation='none')
        plt.xticks(range(output_array.shape[1]),output_array.columns,rotation=90)

    def _delinate_stages(self,maturation_DF,stages,axes):
        path=np.array([maturation_DF[axes[0]],maturation_DF[axes[1]]]).T
        stage_array=np.array(stages)
        return np.argmin(cdist(path,stage_array),axis=1)
    
    def _stage_averages(self,maturation_DF,stages,axes):
        mask=self._delinate_stages(maturation_DF,stages,axes=axes)
        output=[]
        for i in np.unique(mask):
            output.append(np.median(maturation_DF[mask == i],axis=0))
        return pd.DataFrame(data=output,columns=maturation_DF.columns)

    def waterfall_plot(self,figsize=(20,14),title='default',plot=True,normalize=False,smooth=True,ptp=0.1,var=0.05,**kwargs):
        """
        Creates waterfall plot of antigen levels with
        If Antigen_list is provided, no correlative clustering will be performed
        """
        fig=plt.figure(figsize=figsize)
        fig.suptitle(title)
               
        if kwargs.has_key("Antigen_list"):
            Antigen_list = kwargs.get("Antigen_list")
            if set(Antigen_list).issubset(set(self.Antigen_list)) == False:
                print Antigen_list
                print self.Antigen_list
                print "The Difference is: "
                print set(Antigen_list)-set(self.Antigen_list)
                raise RuntimeError('Provided Antigen_list does not match')
        else:
            Antigen_list = self.Antigen_list_filtered

        ery_antigen_path=self.erythroid_cluster_path[Antigen_list]
        self.test=ery_antigen_path
        if smooth == True:
            ery_antigen_path = self._smooth_columns(ery_antigen_path)
        
        if kwargs.has_key("Antigen_list"):
            antigen_dist=pdist(ery_antigen_path[Antigen_list].T,'correlation')
            ax2 = fig.add_axes([0.025,0.1,0.9,0.8])
            a_list=Antigen_list
        else:
            ptp_mask = ery_antigen_path.max(axis=0) - ery_antigen_path.min(axis=0) > ptp #select peak-to-peak height
            var_mask = ery_antigen_path.var(axis=0) > var #select variance greater than...
            selection_mask = ptp_mask | var_mask
            selection_mask_index = [key for key,value in selection_mask.iteritems() if value==True]
            antigen_dist=pdist(ery_antigen_path[selection_mask_index].T,'correlation')
            Y_a = linkage(antigen_dist)
            ax1 = fig.add_axes([0.025,0.6,0.9,0.3])
            Z_a = dendrogram(Y_a,no_plot=False,orientation='top')
            link_sort = Z_a['leaves']
            a_list=[selection_mask_index[i] for i in link_sort]
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax2 = fig.add_axes([0.025,0.1,0.9,0.5])

            
        if kwargs.has_key("num"):
            self.disp_array=self._upsample(ery_antigen_path[a_list],num=kwargs.get('num'))
        else:
            self.disp_array=ery_antigen_path[a_list]
        print a_list
        print self.disp_array.shape
        if normalize == True:
            for col in self.disp_array:
                self.disp_array[col]=self._rescale_array(self.disp_array[col].values)

        ax2.imshow(self.disp_array.values.astype(np.float64),cmap=plt.cm.bwr,interpolation='none')
        plt.xticks(np.arange(len(a_list)),a_list,rotation=90)
        if kwargs.has_key("savefile"):
            fig.savefig(kwargs.get('savefile'),bbox_inches='tight')

    def _antigen_filtering(self,ptp,var):
        ery_antigen_path = self.erythroid_cluster_path[Antigen_list]
        ptp_mask = ery_antigen_path.max(axis=0) - ery_antigen_path.min(axis=0) > ptp #select peak-to-peak height
        var_mask = ery_antigen_path.var(axis=0) > var #select variance greater than...
        selection_mask = ptp_mask | var_mask
        selection_mask_index = [key for key,value in selection_mask.iteritems() if value==True]
        return selection_mask_index
        
    def _antigen_sorting(self):
        pass
    
    def _upsample(self,input_df,num=40,kind='linear'):
        pd.options.mode.chained_assignment = None  # default='warn'
        length = input_df.shape[0]
        output_df=[]
        columns = []
        self.test2=input_df
        for col in input_df:
                self.test1=col
                f=interp1d(range(length),input_df[col].values,kind=kind)
                
                temp = f(np.linspace(0,length-1,num=num))
                output_df.append(temp)
                columns.append(col)
        return pd.DataFrame(data=np.array(output_df).T,columns=columns)
    
    def dotplot(self,axes,display_pts,**kwargs):
        k=0
        if display_pts < 0 | display_pts > self.data.shape[0]:
            display_pts = self.data.shape[0]

        if kwargs.has_key('figsize'):
            figsize = kwargs.get('figsize')
        else:
            figsize = (5*len(axes),6)
            
        fig,ax_mat=plt.subplots(1,len(axes),sharex=True,sharey=True,figsize=figsize)

        if kwargs.has_key('title'):
            title = kwargs.get('title')        
            fig.suptitle(title)

        for ax in axes:
            x_ax=ax[0]
            y_ax=ax[1]
            ax_h=ax_mat[k]
            x=self.data[x_ax][:display_pts]
            y=self.data[y_ax][:display_pts]
            ax_h.scatter(x,y,lw=0,s=1,c='r',marker=',',alpha=0.15) # red background
            x=self.data[x_ax][self.erythroid_mask]
            y=self.data[y_ax][self.erythroid_mask]
            ax_h.scatter(x,y,lw=0,s=1,c='b',marker=',')
            x=self.centroids[x_ax][:]
            y=self.centroids[y_ax][:]
            ax_h.scatter(x,y,lw=2,s=10,c='k',marker='o')
            ax_h.set_xlabel(x_ax)
            ax_h.set_ylabel(y_ax)
            ax_h.set_xlim(0,1)
            ax_h.set_ylim(0,1)

            k+=1
        fig.tight_layout()
        if kwargs.has_key("savefile"):
            fig.savefig(kwargs.get('savefile'),bbox_inches='tight')
            
    def _smooth_columns(self,maturation_paths):
        pd.options.mode.chained_assignment = None  # default='warn'
        smoothed_maturation_paths=maturation_paths.copy()
        for col in maturation_paths:
                temp = wiener(maturation_paths[col])
                smoothed_maturation_paths.loc[:,col]=temp
        return smoothed_maturation_paths
        
    def _rescale_array(self,input_array,inwardcount=0):
        """
        This function takes an array and rescales all values to the range [0,1]
        inwardcount will rescale to I-th less than max and I-th greater than min.
        """
        sort_array=np.sort(input_array)
        maximum=sort_array[-1-inwardcount]
        minimum=sort_array[0+inwardcount]
        return (input_array-minimum)/float(maximum-minimum)

    def _path_gating(self,parameters):
        """
        Internal function for SpadeGating
        Generates a mask of events between the source and target
        """
        dist_metric = parameters['distance_metric']
        dist_mtx=squareform(pdist(self.centroids[self.Antigen_list],dist_metric))
        G=nx.Graph(dist_mtx)
        G_mst=nx.minimum_spanning_tree(G)
        maturation_path=nx.dijkstra_path(G_mst,self.source_graph_idx,self.target_graph_idx)

        maturation_mask=self.labels_mask<0
        for l in maturation_path:
            maturation_mask+= self.labels_mask==self.kept_indicies[l]
        return maturation_mask

    def _find_ends(self,parameters,source,target):
        dist_metric = parameters['distance_metric']

        XA = np.array([source.values()])
        XB = self.centroids[source.keys()]
        source_graph_idx = np.argmin(cdist(XA,XB,dist_metric))
        XA = np.array([target.values()])
        XB = self.centroids[target.keys()]
        target_graph_idx = np.argmin(cdist(XA,XB,dist_metric))
        return source_graph_idx,target_graph_idx
        
    def _filter_centroids(self,parameters):
        elements2remove=parameters['centroid_filter']
        sorted_indicies=self.centroids.sort('centroid_disp').index
        kept_indicies=sorted_indicies[:-elements2remove]
        #print kept_indicies
        #print self.centroids.shape
        return self.centroids.iloc[kept_indicies],kept_indicies

    def _full_clustering(self,parameters,isotypes,removal_list):
        if parameters['clustering'].lower() == 'MiniBatchKmeans'.lower():
            model=MiniBatchKMeans(n_clusters=parameters['clusters_per_file'])
        elif parameters['clustering'].lower() == 'Kmeans'.lower():
            model=KMeans(n_clusters=parameters['clusters_per_file'])
        else:
            raise RuntimeError('clustering parameter is undefined')
        Antigen_to_cluster=[str(key) for key in self.data.keys() if key not in removal_list+isotypes]

        labels=model.fit_predict(self.data[Antigen_to_cluster])

        columns=[self.Antigen_list+['cluster_size']+['centroid_disp']]
        centroids=pd.DataFrame(columns=columns)
        for k in np.unique(labels):
            centroids.loc[k]=self.data[labels==k].median(axis=0)
            centroids.loc[k]['cluster_size'] = np.sum(labels==k)
            centroids.loc[k]['centroid_disp'] = np.sum( cdist(centroids.loc[k][self.Antigen_list].values[np.newaxis],\
                                                        self.data[labels==k], \
                                                        parameters['distance_metric']) \
                                                        /centroids.loc[k]['cluster_size'])
        return centroids,labels,Antigen_to_cluster

    def _apply_gates(self,gates):
        gate_indicies=self.data['CD45']==-10e5  #intialize a blank boolean array
        for values in gates.values():
            gate_indicies+=self._gating_box(self.data,values[0],values[1],values[2])
        return self.data[~gate_indicies]


    def _gating_box(self,DF_array_data,x_ax,y_ax,c):
        '''
        corners should be [x1,y1,x2,y2]
        '''
        coords=[(c[0],c[1]),(c[2],c[1]),(c[2],c[3]),(c[0],c[3]),(c[0],c[1])]
        gate=Path(coords,closed=True)
        projection=np.array(DF_array_data[[x_ax,y_ax]])
        index=gate.contains_points(projection)
        return index

if __name__ == '__main__':
    #Filename = '/home/ngdavid/Desktop/MDS_Plates/09_27_11_Normal_BM_1/SynData_v3.hdf5'
    #Filename = '/home/ngdavid/Desktop/MDS_Plates/11-19244/SynData_v3.hdf5'
    #Filename = '/home/ngdavid/Desktop/MDS_Plates/12-02810/SynData_v3.hdf5'
    Filename = '/home/ngdavid/Desktop/MDS_Plates/12-02606/SynData_v3.hdf5'
    Filename = '/home/ngdavid/Desktop/MDS_Plates/11-17390/SynData_v3.hdf5'
    Gates={'low36' : ['CD36','CD45 KO',[0,0.3,0.3,0.4]],
           'high45': ['CD71 FITC','CD45',[0,0.7,1,1]],
           'low71' : ['CD71 FITC','CD45',[0,0.25,0.5,0.5]],
           '36_38' : ['CD36','CD38',[0,0.3,0.3,0.5]],
           '36_34' : ['CD36','CD34 PE-Cy7',[0,0.3,0.4,0.6]],
           'CD66B' : ['CD66B','CD45',[0.3,0.0,1,1]],
           'CD66'  : ['CD66','CD45',[0.3,0.0,1,1]],
           'CD15'  : ['CD15','CD45',[0.3,0.0,1,1]],
           'CD24'  : ['CD24','CD45',[0.3,0.0,1,1]],
           'Bcell' : ['CD19','CD45',[0.3,0.0,1,1]],
           'CD10'  : ['CD10','CD45',[0.35,0.0,1,1]],
           'Tcell' : ['CD2','CD45',[0.3,0.0,1,1]],
           'CD3'   : ['CD3','CD45',[0.35,0.0,1,1]],
           'CD64'  : ['CD64','CD45',[0.3,0.0,1,1]]}
           
    Parameters = {'display_pts' :60000,
                  'clusters_per_file':50,
                  'erythroid_clusters':25,
                  'centroid_filter':2,
                  'clustering': 'MiniBatchKmeans',
                  'distance_metric':'cityblock'}
                  #'ery_cluster_antigens' : ['CD45','CD36','CD98','CD38','CD49d']}
                  
    Axes=[('SSC-H','CD45'),('CD38','CD34'),('CD71 FITC','CD34'),('CD59','CD45'),('CD58','CD45'),('CD123','CD45')]
    
    Source={'CD59':0.5,'CD58':0.5,'CD45':0.1,'CD34':0.1,'CD38':0.1,'CD117':0.1,'CD36':0.1,'CD71 FITC':0.1,'CD4':0.1}
    
    Target={'CD38':0.6,'CD34':0.6,'CD45':0.6,'CD36':0.2}
    
    Removal_list=['CD71','CD34 PE-Cy7','CD38 A594','CD117 PE','APC-H','CD15 V450','CD19 PE-TR','CD123 PerCP-Cy55','CD235a']
    
    Isotypes=['mIgG1', 'mIgG2a', 'mIgG2b', 'mIgG3', 'mIgM', 'rIgG1', 'rIgG2a', 'rIgG2b', 'rIgM']
    
    test=SpadeGating(filename=Filename,
                     gates=Gates,
                     parameters=Parameters,
                     source=Source,
                     target=Target,
                     removal_list=Removal_list,
                     isotypes=Isotypes)
    #print test.centroids
    Antigen_list = ['CD45','CD34','CD38',
                    'CD47', 'CD108', 'SSEA-4', 'CD105', 'CD227', 'CD220', 'CD44',
                    'CD321', 'CD63', 'CD147', 'CD49d', 'CD317', 'CD81', 'CD59',
                    'CD58', 'CD95', 'CD98', 'CD164', 'CD36', 'CD29']
    test.dotplot(axes=Axes,display_pts=60000)
    test.waterfall_plot(title='RAEB-2',Antigen_list=Antigen_list)
"""
Defines enviromental and runtime variables for spade gating
"""
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

Source={'CD59':0.5,'CD58':0.5,'CD45':0.1,'CD34':0.1,'CD38':0.1,'CD117':0.1,'CD36':0.1,'CD71 FITC':0.1,'CD4':0.1}

Target={'CD38':0.6,'CD34':0.6,'CD45':0.6,'CD36':0.2}

Removal_list=['CD71','CD34 PE-Cy7','CD38 A594','CD117 PE','APC-H','CD15 V450','CD19 PE-TR','CD123 PerCP-Cy55','CD235a']

Isotypes=['mIgG1', 'mIgG2a', 'mIgG2b', 'mIgG3', 'mIgM', 'rIgG1', 'rIgG2a', 'rIgG2b', 'rIgM']

Disp_Antigen_list = ['CD45','CD34','CD38','CD71 FITC'
                'CD47', 'CD108', 'SSEA-4', 'CD105', 'CD227', 'CD220', 'CD44',
                'CD321', 'CD63', 'CD147', 'CD49d', 'CD317', 'CD81', 'CD59',
                'CD58', 'CD95', 'CD98', 'CD164', 'CD36', 'CD29']
                
common_antigens = ['SSC-H','FSC-A','CD45','CD34','CD36','CD38','CD71 FITC']

backbone_antigens = ['SSC-H','FSC-A','CD45','CD34','CD117','CD13','CD38','CD71 FITC']
NormNames = ['Normal1','Normal2','Normal3']
AbnNames = ['FLT3+AML','tMDS','RAEB-1','RCMD','RCUD']

table_of_contents = 'MDS_Plate_Samples.txt'
root_dir = '/home/ngdavid/FCS_Data/MDS_Plates/'

output_file = 'complete.csv'

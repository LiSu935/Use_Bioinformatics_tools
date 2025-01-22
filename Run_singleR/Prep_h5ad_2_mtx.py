# Assuming input data is in following directories:
# - major_dir 
# - - mouse_filtered_5fold_independent
# - - - 0
# - - - 1
# - - - 2
# - - - 3
# - - - 4
# - - - - "*_train*.h5ad", "*_test*.h5ad"




import os
import scanpy as sc

import glob
import numpy as np

import scanpy as sc
from scipy import sparse, io
import pandas as pd
    

os.listdir()



major_dir = '/cluster/pixstor/xudong-lab/suli/Alg_others/scBERT_jobs_Fei/data/'
data_directory_name = 'mouse_filtered_5fold_independent'
INPUT_obs_key = 'cell_subclass'
INPUT_var_key = 'gene_name'

# 20240724 after rename the independent adata and unifying the var name obs name:


tu_list = [(data_directory_name, INPUT_obs_key, INPUT_var_key)]

for tu in tu_list:
    
    ORI_INPUT_DIR = major_dir + tu[0]+'/'
    obs_key = tu[1]
    var_key = tu[2]
    
    # Get a list of subdirectory names
    subdirectories = [name for name in os.listdir(ORI_INPUT_DIR) if os.path.isdir(os.path.join(ORI_INPUT_DIR, name))]
    
    for sub in subdirectories:    
        Train_INPUT_FILE = glob.glob(ORI_INPUT_DIR+sub+"/"+"*_train*.h5ad")[0]
        Test_INPUT_FILE = glob.glob(ORI_INPUT_DIR+sub+"/"+"*_test*.h5ad")[0]
        
        file_paths = [Train_INPUT_FILE, Test_INPUT_FILE]
        
        file_names = list(map(lambda path: os.path.basename(path), file_paths))
        file_names_without_extension = list(map(lambda name: os.path.splitext(name)[0], file_names))
        
        # Assigning file names to specific variables
        train_name, test_name = file_names
        train_name_without_extension, test_name_without_extension = file_names_without_extension
    
        print(Train_INPUT_FILE)
        print(train_name_without_extension)
        
    
        train_ad = sc.read_h5ad(Train_INPUT_FILE)
        test_ad = sc.read_h5ad(Test_INPUT_FILE)
        
        print(train_ad)
        sc.pp.normalize_total(train_ad, target_sum=1e4)
        sc.pp.log1p(train_ad)
        
        print(type(train_ad.X))
    
            
        # Test obs_key, var_key
        #train_ad.obs["cell_type"] = train_ad.obs["celltype"]
        meta_df = train_ad.obs[[obs_key]]
        
        #meta_df['cell_type']
        
        meta_df.rename(columns={obs_key: "cell_type"}, inplace=True)    
    
            
        if isinstance(train_ad.X, np.ndarray):
            print("dense matrix input:")
            io.mmwrite(ORI_INPUT_DIR+sub+"/"+train_name_without_extension+".mtx", sparse.csr_matrix(train_ad.X))
        else:
            print('sparse matrix input')
            io.mmwrite(ORI_INPUT_DIR+sub+"/"+train_name_without_extension+".mtx", train_ad.X)
    
        meta_df.to_csv(ORI_INPUT_DIR+sub+"/"+train_name_without_extension+"_meta.csv", index=True)
        meta_df.to_csv(ORI_INPUT_DIR+sub+"/"+train_name_without_extension+"_label.csv", index=False, header=False)
        
        #train_ad.var['gene_name'] = train_ad.var_names.tolist()
        train_ad.var[[var_key]].to_csv(ORI_INPUT_DIR+sub+"/"+train_name_without_extension+"_genes.csv", index=True)
    
        # Test data
        sc.pp.normalize_total(test_ad, target_sum=1e4)
        sc.pp.log1p(test_ad)
        print(type(test_ad.X))
        meta_df = test_ad.obs[[obs_key]]
        meta_df.rename(columns={obs_key: "cell_type"}, inplace=True)    
      
        if isinstance(test_ad.X, np.ndarray):
            print("dense matrix input:")
            io.mmwrite(ORI_INPUT_DIR+sub+"/"+test_name_without_extension+".mtx", sparse.csr_matrix(test_ad.X))
        else:
            print('sparse matrix input')
            io.mmwrite(ORI_INPUT_DIR+sub+"/"+test_name_without_extension+".mtx", test_ad.X)
    
        meta_df.to_csv(ORI_INPUT_DIR+sub+"/"+test_name_without_extension+"_meta.csv", index=True)
        
        #train_ad.var['gene_name'] = train_ad.var_names.tolist()
        test_ad.var[[var_key]].to_csv(ORI_INPUT_DIR+sub+"/"+test_name_without_extension+"_genes.csv", index=True)
        
    
        
        ### Prep for the singleR
        
        sc.pp.highly_variable_genes(train_ad, flavor="seurat_v3", n_top_genes=2000, subset=True)
        
        sc.pp.highly_variable_genes(test_ad, flavor="seurat_v3", n_top_genes=2000, subset=True)
        
        test_ad.var_names.name = "gene_name"
        train_ad.var_names.name = "gene_name"
        
    
        if isinstance(train_ad.X, np.ndarray):
            train_df = pd.DataFrame(data=train_ad.X, index=train_ad.obs_names, columns=train_ad.var_names)
        else: 
            train_df = pd.DataFrame(data=train_ad.X.todense(), index=train_ad.obs_names, columns=train_ad.var_names) 
        
        if isinstance(test_ad.X, np.ndarray):
            test_df = pd.DataFrame(data=test_ad.X, index=test_ad.obs_names, columns=test_ad.var_names)
        else: 
            test_df = pd.DataFrame(data=test_ad.X.todense(), index=test_ad.obs_names, columns=test_ad.var_names)

        
        train_df_t = train_df.T
        train_df_t.to_csv(ORI_INPUT_DIR+sub+"/"+train_name_without_extension+"_df.csv")
        
        
        test_df_t = test_df.T
        test_df_t.to_csv(ORI_INPUT_DIR+sub+"/"+test_name_without_extension+"_df.csv")
        

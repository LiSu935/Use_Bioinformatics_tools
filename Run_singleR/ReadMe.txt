# Read me:

# Assuming input data is in following directories:
# - major_dir 
# - - mouse_filtered_5fold_independent
# - - - 0
# - - - 1
# - - - 2
# - - - 3
# - - - 4
# - - - - "*_train*.h5ad", "*_test*.h5ad"

Step 1: use Prep_h5ad_2_mtx.py to process from h5ad to all input files required by SingleR.
Step 2: use the run_singleR_batch_version.sh to run singleR. 
        Make sure change the path of `SCRIPT` where you save your 'run_singleR.R' file; change the `major_dir`, `data_directory_name`

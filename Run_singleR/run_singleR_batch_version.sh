#!/bin/bash -l
#SBATCH --partition general     # use partition Gpu
#SBATCH --account joshitr-lab
#SBATCH --cpus-per-task 1
#SBATCH --mem 96G          # memory (up to 256G)
#SBATCH --nodes 1           # number of nodes (1x node has 4x GPUs)
#SBATCH --ntasks-per-node 1
#SBATCH -t 2-00:00                      # max time limit is 7 days
#
## labels and outputs
#SBATCH -J geneformer # job name - shows up in sacct and squeue
#SBATCH -o run_singleR_BatchVersion-%j.out  # filename for the output from this job (%j = job#)#
echo "### Starting at: $(date) ###"
#SBATCH --gres gpu:A100:1       # request generic GPU resources

module load miniconda3/4.10.3_gcc_9.5.0


source activate r-singlecell-r4.1-py3.10
SCRIPT='run_singleR.R'

major_dir="/cluster/pixstor/xudong-lab/suli/Alg_others/scBERT_jobs_Fei/data/"
data_directory_name='mouse_filtered_5fold_independent'

declare -a StringArray=(${data_directory_name}) # 'mouse_unified_5fold'  #"mouse_filtered_5fold_independent") # "COVID_5fold" "NSCLC_5fold" # "MergedMonkey_5fold" "NSCLC_newsplit" "elegans_5fold" "MergedHuman_5fold" "MergedHuman_5fold" "mouse_5fold" ) # 'MergedHuman' 'mouse_115746' "ms_5fold"  "elegans_filtered_5fold"
for dir in "${StringArray[@]}"; do
    parent_dir="${major_dir}${dir}"
    # Use find command to list subdirectories
    # Execute a loop to process each subdirectory
    find "$parent_dir" -mindepth 1 -maxdepth 1 -type d | while read -r subdir; do
        # Process each subdirectory here
        echo "Processing directory: $subdir"
        # Add your processing commands here
        input_dir="$subdir/"
        output_dir="$parent_dir/"
        train_file=$(find "$input_dir" -type f -name "*train*df.csv")
        train_file_input=$(basename "$train_file")
        test_file=$(find "$input_dir" -type f -name "*test*df.csv")
        test_file_input=$(basename "$test_file")
        ref_label_file=$(find "$input_dir" -type f -name "*train*_label_new.csv")
        ref_label_file_input=$(basename "$ref_label_file")
        Rscript ${SCRIPT} -i ${input_dir} -o ${output_dir} -r ${train_file_input} -q ${test_file_input} -l ${ref_label_file_input};
    done
done



echo "### Ending at: $(date) ###"

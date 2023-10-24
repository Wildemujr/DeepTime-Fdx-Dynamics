#!/bin/zsh

# Define the necessary variables
pdb_id=$1
chain_id=$2
other_args="2fdn All"

# Dynamically define the directory containing your scripts based on the script's location
script_dir="$(dirname "$0")"

# Define the log file in the current execution directory
log_file="$PWD/execution_log.txt"

# Clear previous log file and write a new header
echo "Execution Log - $(date)" | tee $log_file

# Reset SECONDS variable to zero at the start
SECONDS=0
total_execution_time=0

# Order of Execution (Files Only)
echo -e "\nExecuting ProteoCorrAnalyser.py..." | tee -a $log_file
python3.11 $script_dir/ProteoCorrAnalyser.py $pdb_id 2fdn All $chain_id A >> $log_file 2>&1
echo "Execution time: $SECONDS seconds." | tee -a $log_file
total_execution_time=$(($total_execution_time + $SECONDS))
SECONDS=0  # reset SECONDS

echo -e "\nExecuting BayesianWeightOptimization.py..." | tee -a $log_file
python3.11 $script_dir/BayesianWeightOptimization.py All $pdb_id.chain$chain_id.nowa >> $log_file 2>&1
echo "Execution time: $SECONDS seconds." | tee -a $log_file
total_execution_time=$(($total_execution_time + $SECONDS))
SECONDS=0  # reset SECONDS

echo -e "\nExecuting ProximityBasedCysFeSelector.py..." | tee -a $log_file
python3.11 $script_dir/ProximityBasedCysFeSelector.py All $pdb_id $chain_id >> $log_file 2>&1
echo "Execution time: $SECONDS seconds." | tee -a $log_file
total_execution_time=$(($total_execution_time + $SECONDS))
SECONDS=0  # reset SECONDS

echo -e "\nExecuting ANMModeAveragesCompiler.py..." | tee -a $log_file
python3.11 $script_dir/ANMModeAveragesCompiler.py All >> $log_file 2>&1
echo "Execution time: $SECONDS seconds." | tee -a $log_file
total_execution_time=$(($total_execution_time + $SECONDS))
SECONDS=0  # reset SECONDS

echo -e "\nExecuting island_dynamic_programming_V9.py..." | tee -a $log_file
python3.11 $script_dir/island_dynamic_programming_V9.py >> $log_file 2>&1
echo "Execution time: $SECONDS seconds." | tee -a $log_file
total_execution_time=$(($total_execution_time + $SECONDS))
SECONDS=0  # reset SECONDS

# For R scripts, using Rscript to execute
echo -e "\nExecuting TileBasedSimilarityVisualization.R..." | tee -a $log_file
Rscript --vanilla $script_dir/TileBasedSimilarityVisualization.R $pdb_id "Chain $chain_id" All >> $log_file 2>&1
echo "Execution time: $SECONDS seconds." | tee -a $log_file
total_execution_time=$(($total_execution_time + $SECONDS))
SECONDS=0  # reset SECONDS

echo -e "\nExecuting ProteinTileSimViz.R..." | tee -a $log_file
Rscript --vanilla $script_dir/ProteinTileSimViz.R $pdb_id "Chain $chain_id" All >> $log_file 2>&1
echo "Execution time: $SECONDS seconds." | tee -a $log_file
total_execution_time=$(($total_execution_time + $SECONDS))
SECONDS=0  # reset SECONDS

echo -e "\nExecuting ProteinChainOptimalPathAnalyzer.R..." | tee -a $log_file
Rscript --vanilla $script_dir/ProteinChainOptimalPathAnalyzer.R --pdb-id=$pdb_id --chain-id=$chain_id >> $log_file 2>&1
echo "Execution time: $SECONDS seconds." | tee -a $log_file
total_execution_time=$(($total_execution_time + $SECONDS))
SECONDS=0  # reset SECONDS

# # At the end, echo the total execution time
# echo -e "\nTotal Execution Time: $total_execution_time seconds." | tee -a $log_file

# Calculate hours, minutes, and seconds
hours=$((total_execution_time / 3600))
minutes=$(( (total_execution_time % 3600) / 60 ))
seconds=$((total_execution_time % 60))

# At the end, echo the total execution time
echo -e "\nTotal Execution Time: $hours hours → $minutes minutes → $seconds seconds." | tee -a $log_file

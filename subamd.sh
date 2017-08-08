#!/bin/bash

# To submit to the AMD nodes (nodes 01-10), use OpenLava.
# Submit the job writing: bsub < subamd.sh

# Select AMD queue [nodes 01-10]
#SUB -q normal

# Select Intel queue [nodes 11-14] - faster, but only 4 nodes x 20 CPUs with only 64GB memory
##BSUB -q intel

# Job name
#BSUB -J ResA01_5

# Redirect screen output to output.txt
#BSUB -o server_script_bsub.err

# Error File
#BSUB -e server_script_bsub.err

# Mail notification
#BSUB -u rafnuss@gmail.com

# Number of processes. <- number of tasks=cpu_cores to reserve for the job

#BSUB -n 48

#BSUB -m node08

# Launch the batch job
nohup /soft/matlab/r2016a/bin/matlab -nodisplay -nosplash -r Server_script  <dummy>  server_script.log

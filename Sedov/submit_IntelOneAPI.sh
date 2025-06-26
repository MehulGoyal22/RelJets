#!/bin/sh
#SBATCH -N 2 ## Number of nodes
#SBATCH --ntasks-per-node=24 ## Number of cores for node
#SBATCH --time=1-00:00:00 ## Time required to execute the program
#SBATCH --job-name=flash ## Name of application
#SBATCH --error=job.%J.err_2_node_48 ## Name of the output file
#SBATCH --output=job.%J.out_2_node_48 ## Name of the error file
#SBATCH --partition=small
#SBATCH --mem-per-cpu=2gb
source /home/apps/spack/share/spack/setup-env.sh
spack load gcc-runtime@14.2.0 /7q6alwo
spack load intel-oneapi-mpi@2021.14.1 /w5s33c3
spack load hdf5@1.14.5 /uipptrn
spack load zlib@1.3.1 /wcqgvku
spack load python@3.11.9 /c4fojts
cd /home/IITB/CompnalAstrphyNRelvity/rahulkashyap/scratch/flash_runs/run_03
## export OMP_NUM_THREADS=8
mpiexec.hydra -np $SLURM_NTASKS ./flash4

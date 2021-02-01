#!/bin/bash

#SBATCH -N 48 
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=4
#SBATCH --mem=182000
#SBATCH --time=24:00:00
#SBATCH --account=Pra20_5206
#SBATCH --gres=gpu:4       
#SBATCH --partition=m100_usr_prod
#SBATCH --qos=m100_qos_bprod
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=s.silvestri@tudelft.nl 

###salloc -N 24 --ntasks-per-node=12 --cpus-per-task=4 --mem=182000 --time=01:00:00 --account=Pra20_5206 --gres=gpu:4 --partition=m100_usr_prod --qos=m100_qos_bprod


mpirun ./RadChan > out_s 

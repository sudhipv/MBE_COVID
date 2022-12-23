#!/bin/bash
#SBATCH --account=account_name
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:15:00
##SBATCH --mem-per-cpu=3700M
#SBATCH --job-name=ff_mpi
#SBATCH --output=%x-%j.out
#SBATCH --mail-user="email"
#SBATCH --mail-type=ALL

#### Only for Niagara ###
module load CCEnv
module load StdEnv/2020
module load openmpi/4.0.3

#### for Graham, Beluga and Cedar
# module load intelmpi/2019.7.217


which mpiexec
which mpirun


mpirun -np 800 /install_location/FreeFem-sources/src/mpi/FreeFem++-mpi covid_scalar2L_southON.edp -nw -log_view


exit


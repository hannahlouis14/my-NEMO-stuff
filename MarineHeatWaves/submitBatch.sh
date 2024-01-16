#!/bin/bash
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^##
##         run python in a job on graham with Slurm                   ##
##            (May07, 2021: weissgib@ualberta.ca)                     ##
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^##
#SBATCH -A rrg-pmyers-ad
#SBATCH -J plt_horizontal_advection
#SBATCH --ntasks=1
#SBATCH --mem=8000M
#SBATCH -t 0-03:30        ## 0 day, 1 hour, 0 minutes
#SBATCH -o slurm-mem-%j.out
#SBATCH -e slurm-mem-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hlouis@ualberta.ca

module restore MYMODS310
module load geos

#source /home/hlouis/otters311/bin/activate
module load python/3.10.2 hdf5/1.10.6 netcdf/4.9.0 scipy-stack/2023a
#source ~/env-serial-workflow/bin/activate

source /home/hlouis/otters310/bin/activate
python /project/rrg-pmyers-ad/hlouis/scripts/MarineHeatWaves/calc_mhw.py

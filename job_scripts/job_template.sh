#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=200:mem=50gb
#PBS -N Linear_0.2_0.5

cd /rds/general/user/hc1017/projects/hda-22-23/live/Summer_projects/hc1017

module load anaconda3/personal
source activate r413

numrep=1000
n=1000
pk=500
nu_xy=0.2
ev_xy=0.5
nchunks=200

Rscript Simulation_Linear.R $numrep $n $pk $nu_xy $ev_xy $nchunks


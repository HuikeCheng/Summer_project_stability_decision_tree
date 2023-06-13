#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=20:mem=50gb
#PBS -N Sig_0.1_0.7_1

cd /rds/general/user/hc1017/projects/hda-22-23/live/Summer_projects/hc1017

module load anaconda3/personal
source activate r413

numrep=100
association=Sigmoidal
n=1000
pk=500
nu_xy=0.1
ev_xy=0.7
proportion=1
nchunks=20

Rscript Simulation_nl.R $numrep $association $n $pk $nu_xy $ev_xy $proportion $nchunks


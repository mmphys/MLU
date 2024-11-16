#!/bin/bash

# Set max wallclock time to one hour
#PBS -l walltime=24:00:00

# Set name of job shown in showq
#PBS -N GPU3000

# Email
#PBS -m abe
#PBS -M Michael.Marshall@ed.ac.uk

# Redirect stdout and stderr to log directory
#PBS -o ../log/
#PBS -e ../log/

# Specify the number of nodes and placement
#PBS -l select=8:ncpus=24:mpiprocs=4
#PBS -l place=scatter

# Budget code
#PBS -A dp008

# RUn on GPUs
#PBS -q QGPU

set -e

NUM_NODES=8
TASKS_PER_NODE=4

# I assume files are submitted from sub subdirectory
cd $PBS_O_WORKDIR
cd ..

######################################################
JOBID=${PBS_JOBID%%.*}
INFO_FILE=log/$PBS_JOBNAME.info.$JOBID

cat $PBS_NODEFILE > $INFO_FILE
numactl -H >> $INFO_FILE
cat /proc/cpuinfo >> $INFO_FILE

NUMLOG=`grep "processor" $INFO_FILE | wc -l`
if [ $NUMLOG -eq 24 ]
then
    echo "HYPERTHREADING IS NOT ENABLED" >> $INFO_FILE
elif [ $NUMLOG -eq 48 ]
then
    echo "HYPERTHREADING IS ENABLED" >> $INFO_FILE
else
    echo "SOMETHING IS WRONG WITH HYPERTHREADING" >> $INFO_FILE
fi
######################################################

mpi_tasks=$[NUM_NODES*TASKS_PER_NODE]

module purge
module load packages-tesseract
module load gcc
module load cuda

# Choose which version of Grid I want to buld against
export GridPre=/tessfs1$HOME/.localgpu
export GridPkg=/home/dp008/dp008/dc-fila1/utils/opa-psm2-11.2.77/usr
       GridPkgLib=$GridPkg/lib64
export GridBuild=/tessfs1$HOME/src/Grid/buildGPU
# These are for Open MPI ... not sure they're really needed?
export LIBRARY_PATH=$GridPre/lib:$LIBRARY_PATH
export FPATH=$GridPre/include:$FPATH
export CPATH=$GridPre/include:$CPATH
export C_INCLUDE_PATH=$GridPre/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$GridPre/include:$CPLUS_INCLUDE_PATH
export PATH=$GridPre/bin:$PATH
export LD_LIBRARY_PATH=$GridPre/lib:$GridPkgLib:$LD_LIBRARY_PATH

export PSM2_CUDA=1
export PSM2_GPUDIRECT=1
export TMPDIR=/dev/shm/

export OMP_NUM_THREADS=6

JOBNAME=${PBS_JOBNAME}
NODEFILE=${PBS_NODEFILE}

#workdir="${PBS_O_WORKDIR}"
application="./GPURun"
args=sub/P3000rest.xml
args="${args} --grid 24.24.24.64"
args="${args} --mpi 2.2.2.4"
#args="${args} --log Error,Warning,Message,Debug"
CMD="mpirun -x PSM2_MULTIRAIL=0 -x OMP_NUM_THREADS -x PATH -x PSM2_GPUDIRECT=1 -x PSM2_CUDA=1 -x LD_LIBRARY_PATH -machinefile $PBS_NODEFILE -map-by node -np $mpi_tasks $application ${args} --gpu-thread &> log/${JOBNAME}.l${JOBID}"

###############################################################
### You should not have to change anything below this line ####
###############################################################

#cd ${workdir}
echo -e "Job running in `pwd`.\n"

echo -e "JobID: ${PBS_JOBID}\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

cat ${NODEFILE} | uniq > log/${JOBNAME}.machine.${JOBID}
echo -e "\nNodes allocated:\n=================="
echo `cat log/${JOBNAME}.machine.${JOBID} | sed -e 's/\..*$//g'`

# echo $SLURM_JOB_NODELIST

# if [ "$SLURM_JOB_NODELIST" ]; then
#         #! Create a machine file:
#         export NODEFILE=`generate_pbs_nodefile`
# 	cat $NODEFILE
#         cat $NODEFILE | uniq > machine/$JOBID
#         echo -e "\nNodes allocated:\n================"
#         echo `cat machine/$JOBID | sed -e 's/\..*$//g'`
# fi

# echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nCheck links:\n==================\n"
ldd ${application}

echo -e "\nExecuting command:\n==================\n${CMD}\n"

eval ${CMD}

echo -e "\nEND ==================\n"

export PROCS_ON_EACH_NODE=12

# ************* SGE qsub options ****************
#Export env variables and keep current working directory
#$ -V -cwd
#Select parallel environment and number of parallel queue slots (nodes)
#$ -pe mpi-verbose 1
#$ -P mccombes-sgm.prj 
#Combine STDOUT/STDERR
#$ -j y
#Specify output file
#$ -o out.$JOB_ID
#Request resource reservation (reserve slots on each scheduler run until enough have been gathered to run the job
#$ -R y

# Add runtime indication
#$ -ac runtime=""
# ************** END SGE qsub options ************

export NCORES=`expr $PROCS_ON_EACH_NODE \* $NSLOTS`
export OMPI_MCA_btl=openib,self
export OSTYPE=linux
cat /proc/meminfo > dump
cat /proc/cpuinfo >> dump
./main < case_files/heli.cas > dump

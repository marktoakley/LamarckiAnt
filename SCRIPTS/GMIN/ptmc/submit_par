#!/bin/bash
##############################################################################
# js850> this is a submit script for a parallel GMIN PTMC or BSPT simulation.
# The output of the job will be written to directory /scratch/$USER/$PBS_JOBID/
# on each node.  If the job ends without error the files will be coppied back
# to the submit directory.

#
# usage: you must pass the number of nodes.  e.g. for 40 processers on 10 nodes
# qsub -l nodes=10:ppn=4  submit_par
#


# It starts off with a long comment block that contains several PBS/Torque
# directives. PBS/Torque directives start like this: 

#PBS 

# This is so that they appear as comments to the shell. To
# comment out a directive you must add another # to the start: ##PBS.
#
# PBS/Torque directives set various options that you can also set on the qsub
# command line. Doing it this way saves typing.
#
# This script is suitable for serial jobs.  Please see the other scripts in
# /info/pbs for parallel examples.
#
##############################################################################
# start of PBS directives
##############################################################################
#PBS -N ptmcGMIN
#PBS -q h32
###PBS -l nodes=2:ppn=2 
#
# Use this to choose a node if needed. Generally a very bad idea but
# sometimes vital for testing.
##PBS -l nodes=comp017
#
# Use this to adjust required walltime , up to the maximum for the queue
# PLEASE make an effort to estimate this and set it accordingly, plus a
# generous safety margin, rather than simply going for the maximum 
# allowed. Jobs with shorter walltimes get higher priority and are more
# likely to be backfilled. Accurate walltime estimates help the scheduler
# do a good job, keeping everyone on the machine happy.
##PBS -l walltime=163:01:00
#
##############################################################################
# Start of shell script proper. Do not put PBS directives after this point.
# PBS uses your login shell to interpret the script- you cannot write it in
# csh if your login shell is bash.
##############################################################################
#
# If you want to submit your job from the directory you want it to run in
# then do this (note that this has to exist on the host the job is running
# on; beware when using non-shared /scratch)
#
#
# Write out some helpful info
#

#filenames
userscratch=/scratch/$USER
localscratchdir=/scratch/$USER/${PBS_JOBID}
outputmaster=${PBS_O_WORKDIR}/output.master.${PBS_JOBID}
infofile=${PBS_O_WORKDIR}/info${PBS_JOBID}
nodefile=${PBS_O_WORKDIR}/nodes.${PBS_JOBID}
nodescratchfile=${PBS_O_WORKDIR}/nodescratch.${PBS_JOBID}

#variables
node=`cat $PBS_NODEFILE`
headnode=`uname -n`
nodes=`cat $PBS_NODEFILE | sed 's/.cm.cluster//' | sort | uniq`
numproc=`cat $PBS_NODEFILE | wc -w`
if [ "$numproc" -eq 48 ]; then
  echo "numproc is 48, you probably forgot to specify -l nodes=#" >&2
  exit
fi

#node list files
echo "$nodes" > $nodefile
for n in $nodes; do
  echo /nodescratch/$n/${USER}/$PBS_JOBID >> $nodescratchfile
done

#info file
echo "PBS assigned me these nodes:" > $infofile
wc $PBS_NODEFILE >> $infofile
cat $PBS_NODEFILE >> $infofile
echo First task running on: >> $infofile
echo $headnode >> $infofile
echo "qsub working directory" >> $infofile
echo ${PBS_O_WORKDIR} >> $infofile
echo "node working directory" >> $infofile
echo "$localscratchdir"  >> $infofile
echo "Job started at" >> $infofile
echo `date` >> $infofile

#output
echo "Starting job $PBS_JOBID"
echo "PBS assigned me these nodes:"
wc $PBS_NODEFILE
cat $PBS_NODEFILE
echo `date`
echo
#
# Here is where you should set any environment variables your job requires,
# because PBS won't read your shell startup files. The values of PATH and
# LD_LIBRARY_PATH are inherited from the environment you had when you 
# submitted the job, so most people won't have to do this. If you do need to
# manipulate the environment bear in mind that this script gets interpreted
# by your login shell, so you must use the correct syntax for it.
#
# Example of setting a variable for Bourne shell users.
#
#export LD_LIBRARY_PATH=/usr/local/pgi/linux86/lib:/usr/local/intel/mkl60/lib/32:/usr/local/intel_fc_80/lib
#echo "My LD_LIBRARY_PATH is:"
#echo $LD_LIBRARY_PATH
#echo
#

##########################################################################
#copy the startup files to the local scratch
##########################################################################
#########################################################################
#make a working directory in each of the nodes
#copy startup files to each of the nodes
#
#must manually create directories 1 to $numproc because the program creates the
#directories only in the headnode.  This is not optimal
#########################################################################
cd ${PBS_O_WORKDIR}
echo "nodes " $nodes
for n in $nodes; do
  echo "making directories now: from pwd" $PWD
  ssh $n "if [ ! -e $userscratch ]; then mkdir ${userscratch}; fi"
  ssh $n "if [ ! -e $localscratchdir ]; then mkdir ${localscratchdir}; fi"
  scp -p data coords $n:${localscratchdir}
  if [ -e temperatures.init ]; then
    scp -p temperatures.init $n:${localscratchdir}
  fi
  if [ -e overlap.input ]; then
    scp -p overlap.input $n:${localscratchdir}
  fi
done
#make directories and copy them over or copy directories
echo "numproc " $numproc
orderednodes=`cat $PBS_NODEFILE | sed 's/.cm.cluster//'`
echo "orderednodes " $orderednodes
count=1
for n in $orderednodes; do
  if [ ! -e $count ]; then mkdir $count; fi
  scp -pr $count $n:${localscratchdir}
  echo "scp -pr $count $n:${localscratchdir}"
  count=$((count+1))
done

#
#start the program
# you may need to change this line to run a different program
# 
mpirun -wdir $localscratchdir /home/$USER/git/GMIN/source/parbuild/GMIN > $outputmaster


echo "main program finished at"  >> $infofile
echo `date` >> $infofile
echo "starting cleanup"  >> $infofile
echo "main program finished at"
echo `date`
echo "starting cleanup"

###############################################################
# do cleanup:
#  mv the files back and delete them from localscratch
###############################################################
#
#copy files back from each of the nodes
#
for n in $nodes; do
  scp -rp $n:${localscratchdir}/* ${PBS_O_WORKDIR}/
  #if copy was successfull then delete folder
  if [ "$?" -eq 0 ]; then 
    ssh $n "rm -r $localscratchdir"
    #echo "not executing: ssh $n \"rm -r $localscratchdir\""
  else
    echo "error copying from node $n $localscratchdir" >&2
  fi
done



#
# write end time to info file
#
echo "Finished at"  >> $infofile
echo `date` >> $infofile
# Write out the PBS details of the job. If it ran for any sensible length
# of time this will tell you the CPU, walltime, and memory usage.
echo
echo "Job finished." 
echo `date`
echo "PBS details are:"
echo
qstat -f ${PBS_JOBID}@$HOSTNAME.ch.cam.ac.uk
echo

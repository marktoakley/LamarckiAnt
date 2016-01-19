#!/bin/bash
# This script helps with managing jobs submitted with the PBS submit script
# submit_par.  The script does a specified action on each node in the parallel
# job defined by file nodes.<PBS_JOBID>.  The default command is 'ls'

ex="ls"
cflag=""
dflag=""
qflag="-v"
while getopts "x:cdq" opt; do
  case $opt in
    x)
      ex="$OPTARG"
      echo "will excecute script '$ex' "
      ;;
    c)
      cflag="-c"
      echo "cflag"
      ;;
    d)
      dflag="-d"
      echo "dflag"
      ;;
    q)
      qflag=""
      echo "qflag"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done
shift $((OPTIND-1))

fname=$1
jobid=`echo $fname | sed 's/nodes\.//'`
if [ $# -ne 1 -o ! "$fname" -o "$fname" = "$jobid" ]; then
    echo "This script helps with managing jobs submitted with the PBS submit script submit_par.  The script does a specified action on each node in the parallel job defined by file nodes.<PBS_JOBID>.  The default command is 'ls'"
    echo ""
    echo usage:
    echo "$0 [-x comm] [-c] nodes.PBS_JOBID"
    echo "  -c: copy all files from nodes to current directory"
    echo "  -d: delete all files from nodes directory /scratch/$USER/pbsjobid"
    echo "  -x comm: execute command comm on all the nodes"
    exit
fi

echo $fname
echo $jobid

localscratch=/scratch/$USER/$jobid
nodes=`cat $fname`

if [ -n "$cflag" ]; then
  #make a backup of bsptrestart so if something goes wrong we're not totally lost
  #get the time step of the bsptrestart file.  it will be the first field in the first line
  echo "Making a backup of the bsptrestart files"
  for n in `seq 300`; do
    if [ ! -d $n ]; then
      break;
    fi
    oldfile=$n/bsptrestart
    mcstep=`head -n 1 $oldfile | awk '{print $1}'`
    newfile=$oldfile.$mcstep
    echo "cp -p $oldfile $newfile"
    cp -p $oldfile $newfile
  done

  #now copy all data from nodes
  for n in $nodes; do
    echo "copying files from $n:$localscratch"
    #scp -rp $n:$localscratch/* ./
    rsync -az $qflag $n:$localscratch/ ./
  done
  exit
fi

if [ -n "$dflag" ]; then
  echo "really delete directory: "
  echo "$localscratch"
  echo "on nodes:"
  echo "$nodes"
  echo "? (y/n)"
  read yn
  if [ "$yn" != y ]; then
    exit
  fi
  for n in $nodes; do
    echo "deleting directory $n:$localscratch"
    ssh $n "rm -r $localscratch"
    #rsync -avz $n:$localscratch/ ./
  done
  exit
fi

for n in $nodes; do
  echo 
  echo $n
  echo 
  ssh $n "cd $localscratch; $ex"
done

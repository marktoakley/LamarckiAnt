#!/bin/bash
# CONTROL SCRIPT FOR ITERATIVE RESP CHARGE GENERATION
# Original script by Kyle Sutherland-Cash, modified by Chris Whittleston (csw34)

# Useage: ./control_script.sh <name> <charge> <number of iterations>
# name = the name of the .pdb file i.e. OTV for OTV.pdb
# charge = the total charge for the molecule 
# number of iterations = how many charge generation cycles to run (it is recommended that you use 4+)

# Modify this if you have the scripts in a non-standard location
scriptpath=~/svn/SCRIPTS/AMBER/iterative_resp

###########################################################
# You shouldn't need to modify anything below this point! #
###########################################################

# Sanity check for number of arguements
expected_args=3
if [ $# -lt $expected_args ]
then
   echo "Missing arguements! :("
   echo "./control_script.sh <name> <charge> <number of iterations>"
   exit
fi

# Assign arguements to variables
name=$1
charge=$2
niterations=$3

# Run through the cycles for specified number of iterations
for ((i=1; i <= $niterations ; i++))
do
    mkdir iteration_${i}
    mkdir iteration_${i}/gmin
    mkdir iteration_${i}/gamess
    ${scriptpath}/run_gmin.sh ${name} ${i}
    ${scriptpath}/run_gamess.sh ${name} ${i} ${charge}
    penultimate=${i}
    final=$((${i} + 1))
done

# For the final iteration, just run GMIN
mkdir iteration_${final}
mkdir iteration_${final}/gmin
${scriptpath}/run_gmin.sh ${name} ${final}

# Generate a prepin file from the charges and a matching PDB
sed 's/   [1-9]   /   1   /g' iteration_${final}/gmin/lowest1.1.pdb > ./tmp.pdb
cat iteration_${penultimate}/gamess/${name}.log | ${scriptpath}/gmstoresp.sh ./tmp.pdb > ./${name}_final.prep

# Clean up
rm ./tmp.pdb
mv ./NEWPDB.PDB ./${name}_final.pdb

#!/bin/bash
# COMPONENT SCRIPT FOR ITERATIVE RESP CHARGE GENERATION (SPAWNS GMIN RUNS)
# Original script by Kyle Sutherland-Cash, modified by Chris Whittleston (csw34)

# The following files should be present in the working directory:
# data = the GMIN data file
# atomgroups = file to specify GROUPROTATION groups for the molecule
# min.in = AMBER parameter file
# gmin_job = qsub script to run AMBGMIN

#############################################################################################
# You should not need to modify any of this script. It is called by control_script.sh only! #
#############################################################################################

# Sanity check to see if GMIN data file is present
if [ -e ./data ] 
then
    echo "Running GMIN"
else
    echo "ERROR: Missing GMIN files! Make sure required files are present (data, min.in, atomgroups gmin_job)"  
    exit
fi
 
# Assign arguement to variable
name=$1

# For the first iteration, copy the files from build_structure
if [ $2 -eq 1 ]
then
    iteration=./iteration_${2}

    cp ./build_structure/${name}.inpcrd ${iteration}/gmin/coords.inpcrd
    cp ./build_structure/${name}.prmtop ${iteration}/gmin/coords.prmtop
    cp ./min.in ${iteration}/gmin/
    cp ./atomgroups ${iteration}/gmin/
    cp ./data ${iteration}/gmin/
    cp ./gmin_job ${iteration}/gmin/gmin_${name}_${2}
fi

# For subsequent iterations, copy files from the main directory and propagate...
if [ $2 -ne 1 ]
then
    prev_iteration=./iteration_$((${2} - 1))
    iteration=./iteration_${2}
        
    cp ./min.in ${iteration}/gmin/
    cp ./atomgroups ${iteration}/gmin/
    cp ./data ${iteration}/gmin/
    cp ./gmin_job ${iteration}/gmin/gmin_${name}_${2}
    cp ${prev_iteration}/gmin/min1.1.rst ${iteration}/gmin/coords.inpcrd
    cp ${prev_iteration}/gamess/new_coords.prmtop ${iteration}/gmin/coords.prmtop
fi

# Run the GMIN job
cd ${iteration}/gmin/
qsub ./gmin_${name}_${2}
# Wait until it's finished
while [ ! -e ./lowest1.1.pdb ]
do
    sleep 120
done
# Change back to the other directory
cd ../../

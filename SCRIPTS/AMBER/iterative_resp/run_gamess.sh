#!/bin/bash
# COMPONENT SCRIPT FOR ITERATIVE RESP CHARGE GENERATION (RUNS GAMESS)
# Original script by Kyle Sutherland-Cash, modified by Chris Whittleston (csw34)

# Modify this if you have the scripts in a non-standard location
scriptpath=~/svn/SCRIPTS/AMBER/iterative_resp

#############################################################################################
# You shouldn't need to modify anything in this script as it is called by control_script.sh #
#############################################################################################

# Get the iteration
name=$1
iteration=./iteration_${2}
charge=$3

# Copy across the lowest energy structure
cp ${iteration}/gmin/lowest1.1.pdb ${iteration}/gamess/${name}.pdb
# Write GAMESS input
${scriptpath}/write_gamess_input.py ${iteration}/gamess/${name}.pdb \
                                    ${iteration}/gamess/${name}.inp \
                                    ${iteration}/gamess/${name}.pdb \
                                    ${iteration}/gamess/gamess_${name}_${2} \
                                    ${charge}
# Run the GAMESS job
cd ${iteration}/gamess/
chmod u+x ./gamess_${name}_${2}
qsub ./gamess_${name}_${2}
# Monitor the GAMESS job for completion
while [ ! -e ./gamess_${name}_${2}.e* ]
do
    sleep 120
done

# Change back to the other directory
cd ../../
# When it's done, write out a new coords file for the next GMIN run 
${scriptpath}/gamess_to_charges.sh ${iteration}/gamess/${name}.pdb \
                                   ${iteration}/gamess/${name}.log \
                                   ${charge} \
                                   ${iteration}/gmin/coords.prmtop \
                                   ${iteration}/gamess/new_coords.prmtop

# Remove scratch files to allow the next iteration to run with the same name
rm ~/scr/${name}.dat

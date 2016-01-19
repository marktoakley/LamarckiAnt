#!/bin/bash
# Script to calculate the interaction enthalpy between a ligand and protein during a GMIN run
# by Chris Whittleston (8/12/09) - csw34@cam.ac.uk
# This script is intended for use with the A9INTE keyword in AMBGMIN
# You need to place a copy of this script in the directory where you are running GMIN 
# There are only two changes you need to make below:

# 1. Make sure you have 'sander' in your path. If not, see the comment below
# 2. Change the residue specifying the ligand to match your system in the 'grep' line below

# echo the SANDER input file 
# For more info on how this works, look up the IDECOMP keyword in the AMBER manual
cat > intres.in <<EOF
Interaction energy per residue input for SANDER
 &cntrl
  imin   = 1,
  idecomp = 3,
  ncyc = 1,
  maxcyc = 1,
  igb = 1, saltcon=0.1,
  ntb    = 0,
  cut    = 999.9
  rgbmax = 8.22,
 /
First set, numbers here are irrelevent as we do a pairwise decomposition!
RES -1 318
END
Second set
RES 319
END
END 
EOF

# run SANDER (assumes you have it in your PATH, if not, see below)
sander -O -i intres.in -c coords.intres -p coords.prmtop -o intres.out
# *** YOU MIGHT NEED TO UNCOMMENT SOMETHING HERE ***
# If you do not already have AMBER9+ in your path (try typing 'sander' on the command line), you can
# uncomment the below 
#/home/csw34/amber/exe/sander -O -i intres.in -c coords.intres -p coords.prmtop -o intres.out

# *** YOU NEED TO SET SOMETHING HERE ***
# This is where we grep the interaction enthalpy out from the AMBER output file
# You need to specify which residue in your system is the ligand, for example - if the
# ligand is residue 220, the line below should look like this:
# grep 'TDC     220->' intres.out | sed '$d' | awk '{ SUM += ($5+$6+$7+$8)} END { print SUM }' > intE
# *** MODIFY THE LINE BELOW FOR YOUR SYSTEM ***
grep 'TDC     220->' intres.out | sed '$d' | awk '{ SUM += ($5+$6+$7+$8)} END { print SUM }' > intE

# The energy in intE is then read back into GMIN  

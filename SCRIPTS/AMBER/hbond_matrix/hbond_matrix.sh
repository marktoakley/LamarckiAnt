#!/bin/bash
# Script to produce a hydrogen-bond matrix for the structure provided using AMBER ptraj
# The matrix for all hydrogen-bonds is output to standard out. Matricies for all,
# backbone only and sidechain only hydrogen-bonds are created as files hbond_all.mat,
# hbond_bb.mat and hbond_sc.mat respectively.
# Chris Whittleston April 2014 (csw34@cam.ac.uk)

# Required arguements:
# $1 = the AMBER topology file (e.g. coords.prmtop)
# $2 = the AMBER coordinate file (e.g. coords.inpcrd, min1.rst, quench1.pdb etc. Use .rst if possible.)
# $3 = ptraj donor/acceptor list
# $4 = file containing residues of interest (1 per line)

# Optional arguements:
# $5 = the hydrogen-bond distance cutoff (default 3.00A) 
# $6 = the hydrogen-bond angle cutoff (default 120 degrees)
# Check for the correct number of arguements
expected_args=4

if [ $# -lt $expected_args ]
then
   echo "Missing arguements! :( <required> (optional) e.g."
   echo "./hbond_matrix.sh <topology file> <rst/pdb file> <ptraj donor/acceptors> <residues of interest> (distance cutoff) (angle cutoff)"
   exit
fi

# Set hydrogen-bond cutoffs if specified (ptraj defaults are 3.00 and 120.00)
distance_cutoff=3.00
angle_cutoff=120.00

if [ $# -gt $expected_args ]
then
   echo "hbond out hbond.out distance $5 angle $6" >> hbond.in.bottom 
else
   echo "hbond out hbond.out distance ${distance_cutoff} angle ${angle_cutoff}" >> hbond.in.bottom 
fi

# Construct full ptraj input file
echo "trajin ${2}" > hbond.in.top
cat hbond.in.top $3 hbond.in.bottom > hbond.in

# Run ptraj
ptraj $1 < hbond.in > ptraj.tmp

# Calculate matrix elements
# OUTER RESIDUE LOOP
for res1 in `cat $4`
do
# INNER RESIDUE LOOP
   for res2 in `cat $4`
   do
# If not on the diagonal (self)
      if [ "${res1}" != "${res2}" ]; then
# All hydrogen-bonds
         nhbonds_all=`grep ":${res1}@" hbond.out | grep ":${res2}@" | wc | awk '{print $1}'`
         echo $nhbonds_all >> ${res1}_all.tmp
# Backbone hydrogen-bonds (involving O or H atoms)
         nhbonds_bb=`grep ":${res1}@" hbond.out | grep ":${res2}@" | grep -E "@O |@H " | wc | awk '{print $1}'`
         echo $nhbonds_bb >> ${res1}_bb.tmp
# Sidechain hydrogen-bonds
         nhbonds_sc=$((nhbonds_all - nhbonds_bb))
         echo $nhbonds_sc >> ${res1}_sc.tmp
      else
# If on the diagonal, just echo 0
         echo "0" >> ${res1}_all.tmp
         echo "0" >> ${res1}_bb.tmp
         echo "0" >> ${res1}_sc.tmp
      fi
   done
done

# Construct the paste commands to assemble matricies
command_all="paste "
command_bb="paste "
command_sc="paste "
for res in `cat $4`
do
   command_all="${command_all} ${res}_all.tmp "
   command_bb="${command_bb} ${res}_bb.tmp "
   command_sc="${command_sc} ${res}_sc.tmp "
done

# Execute paste commands
$command_all > mat_all.tmp
$command_bb > mat_bb.tmp
$command_sc > mat_sc.tmp

# Nicely format the output and print it
column -t mat_all.tmp
# Output to files (all, just backbone, just sidechain) 
column -t mat_all.tmp > hbond_all.mat
column -t mat_bb.tmp > hbond_bb.mat
column -t mat_sc.tmp > hbond_sc.mat

# Remove temporary files 
rm hbond.out hbond.in hbond.in.top hbond.in.bottom *.tmp

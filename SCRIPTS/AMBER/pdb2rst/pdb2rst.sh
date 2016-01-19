#!/bin/bash
# Script to take a PDB file (from AMBGMIN or A9OPTIM) and using the appropriate topology file, generate an
# AMBER .rst file

# Check that ptraj is in the users $PATH
check=`which ptraj | awk '{print $2}'`
if [ "$check" = "no" ]
   then
   echo "You don't have ptraj in your PATH. Do you have AMBERHOME set in your .bashrc file?"
   exit
fi

# Check to make sure the user has specified the correct number of arguements 
EXPECTED_ARGS=2
if [ $# -ne $EXPECTED_ARGS ]
   then 
       echo 'Usage: ./pdb2rst.sh <prmtop file> <pdb file>'
       echo 'e.g. ./pdb2rst.sh coords.prmtop lowest1.pdb' 
       exit
fi

# Set up some variables to make the following clearer 
prmtop=$1
pdbfile=$2
rstfile=`echo $pdbfile | sed 's/pdb/rst/'`

echo "Creating AMBER restart file (inpcrd format) "$rstfile" from "$pdbfile" using "$prmtop

# Use a temporary pdb file name as if the name is too long, it breaks ptraj
# Also remove all information after the atom positions as bad formating can break ptraj
cut -b1-54 $pdbfile > temp.pdb

# Create the ptraj input file
cat > ptraj.in <<EOF
trajin temp.pdb
trajout temp.rst restart title title
EOF

# Run ptraj
ptraj $prmtop < ptraj.in > ptraj.out

# Remove all the 0.0000000 entries used for velocities
sed '/0\.0000000/d' temp.rst.1 | sed 's/title//' > $rstfile

# Clean up :)
rm $initialrst temp.rst temp.pdb ptraj.in

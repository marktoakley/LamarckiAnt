#!/bin/bash
# Script to add N and C to the appropriate residues in an AMBER topology file before it is symmetrised and to restore it after.
# You need a PDB file corresponding to the same sequence as the topology file you are altering and containing OXT atoms for terminal residues.
# The easiest way to create this file is to do 'savepdb' in your LeaP script where you make the topology file itself.

# Requires the following files to be in your path:
# perm-prmtop.ff03.py (for the ff03 and ff99SB forcefields)
# perm-prmtop.ff03ua.py (for the united atom ff03 forcefield)
# perm-prmtop.ff02.py (only needed if using ff02)
# NOTE: in each of these files, you need to set the 'path' at the top to point to a working AMBER(9/10) directory containing our modified AMBER files
# Chris Whittleston (csw34@cam.ac.uk) 2009 

# Deal with arguements 
TOPFILEIN=$1
TOPFILEOUT=$2
PDBFILE=$3
FORCEFIELD=$4

EXPECTED_ARGS=4

# Check correct usage (i.e. arguements provided to the script)
if [ $# -ne $EXPECTED_ARGS ]
   then NAME='Usage: A9prmtoptermini.sh <input prmtop> <output prmtop> <input PDB file (contains OXTs)> <forcefield of ff99SB/ff03/ff03ua>'
        echo 'e.g. A9prmtoptermini.sh coords.prmtop.old coords.prmtop input.pdb ff99SB'
       exit
fi

# extract the original residue block to be added back in later
sed -n '/FLAG RESIDUE_LABEL/,/FLAG RESIDUE_POINTER/p' $TOPFILEIN | sed -n '3,$p' | sed '$d' > prmtop_origresidues

# extract everything up to the residue block
sed -n '/VERSION/,/FLAG RESIDUE_LABEL/p' $TOPFILEIN > prmtop_oldtop
echo '%FORMAT(20a4)' >> prmtop_oldtop

# extract everything below the residue block
sed -n '/%FLAG RESIDUE_POINTER/,$p' $TOPFILEIN > prmtop_oldbottom

# awk line to print out a residue block in AMBER prmtop format from the PDB file given to this script. You need to have the termini marked by OXT atoms.
# Thanks go to David Burke for this piece of hackage :)
awk ' BEGIN { oxt=1}  /ATOM/ && ( $3=="N" || $3=="C1") && length(reslabel)>0{ n++; printf"%-4s", reslabel; reslabel=sprintf("%-3s",$4)} /ATOM/ && ( $3=="N" ) && oxt==1 { reslabel=sprintf("N%-3s",$4); oxt=0} /OXT/ { reslabel=sprintf("C%-3s",$4); oxt=1} n==20 { printf"\n"; n=0 }  END  { printf"%-4s\n", reslabel } ' $PDBFILE > prmtop_termresidues

# construct the modified prmtop file
cat prmtop_oldtop > working.prmtop
cat prmtop_termresidues >> working.prmtop 
cat prmtop_oldbottom >> working.prmtop

# do the symmetrisation
# this command should change with forcefield!
# i.e. for ff03 and ff99SB, you should run perm-prmtop.ff03.py 
# for ff03ua (united atom) you should run perm-prmtop.ff03ua.py

if [[ "$FORCEFIELD" == "ff99SB" ]]
   then perm-prmtop.ff03.py working.prmtop workingsym.prmtop
elif [[ "$FORCEFIELD" == "ff03" ]]
   then perm-prmtop.ff03.py working.prmtop workingsym.prmtop
elif [[ "$FORCEFIELD" == "ff03ua" ]]
   then perm-prmtop.ff03ua.py working.prmtop workingsym.prmtop
elif [[ "$FORCEFIELD" == "ff02" ]]
   then perm-prmtop.ff02.py working.prmtop workingsym.prmtop
else
   echo "Invalid forcefield specified :s Use ff99SB/ff03/ff03ua only"
fi

# Save a copy of the prmtop file with termini
cp workingsym.prmtop plustermini.prmtop

# now need to replace the original residue block
# extract everything up to the residue block from symmetrised prmtop
sed -n '/VERSION/,/FLAG RESIDUE_LABEL/p' workingsym.prmtop > prmtop_top
echo '%FORMAT(20a4)' >> prmtop_top

# extract everything below the residue block for symmetrised prmtop
sed -n '/%FLAG RESIDUE_POINTER/,$p' workingsym.prmtop > prmtop_bottom

# construct the final output prmtop

cat prmtop_top > $TOPFILEOUT
cat prmtop_origresidues >> $TOPFILEOUT
cat prmtop_bottom >> $TOPFILEOUT

# remove temporary files

rm prmtop* working*

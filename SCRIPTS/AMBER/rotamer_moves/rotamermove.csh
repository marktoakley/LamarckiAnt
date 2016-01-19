#!/bin/csh
##########################################
# ROTAMER MOVES FOR GMIN - DRIVER SCRIPT #
# Chris Whittleston October 2010         #
# Uses PdbRotamerSearch by David Burke   #
##########################################

######################################
# STUFF YOU CAN MODIFY IN THE SCRIPT #
######################################
# ROTAMER LIBRARY 
set sctype = -bbdlovell
set sclib = ./penultimate.lib

# RMSD THRESHOLDS
set minrmsd = 0.5
set rmsdtries = 10

# EXAMPLE USE outside GMIN:
# csh rotamermove.csh 1 1.0 8 104 0.004

########################################
# VALUES SET BY GMIN (ROTAMER keyword) #
########################################

set pdb = beforerotamer.pdb
# ROTMAXCHANGE
# Maximum number of residues to change sidechain rotamers for
set maxrot = $1
echo 'maxrotchange = '$maxrot

# ROTPSELECT
# Probability of selecting a sidechain to change
set rotp = $2
echo 'rotthreshold = '$rotp

# ROTCUTOFF
# Cutoff for rotamer moves (measured from residuecentre)
set cutoff = $3
echo 'distthreshold = '$cutoff

# ROTCENTRE
# Residue to use as centre for cutoff 
@ lig = $4 - 1
echo 'residuecentre = '$4

# ROTOCCUW
# The amount the occupation frequency affects the rotamer selection i.e.
# rotoccuw = 0.0 -> all rotamers possible
# rotoccuw = 0.004 -> only rotamers with a > 0.4% occupation selected
set occuw = $5
echo 'bbdminprob = '$occuw

#####################
# MAKE ROTAMER MOVE #
#####################

# CLEAN UP OLD FILES
rm afterrotamer.pdb afterrotamer.xyz reschanged.dat 

# RUN DAVID BURKE'S ROTAMER CODE
./PdbRotamerSearch  -maxrotchange $maxrot -rotthreshold $rotp -distthreshold $cutoff -residuecentre $lig -use_water -pdb $pdb -bbdminprob $occuw -bbdlib $sclib $sctype -outpdb -v -maxrmsdtries $rmsdtries -rmsdthreshold $minrmsd -tryharder -hydrogens -randgroup > ${pdb:r}.SRLOG

# OTHER FLAGS USED (and what they do!)
# -hydrogens : specifies that hydrogens should be added back to changed residues using David's code to avoid a LEaP bug
# -randgroup : specifies a RANDOM rotamer within -distthreshold should be changed rather than than the closest first (very important when changing 1)

# UNUSED OPTIONS FOR DAVID'S CODE (could be added to the call above)
# -tweakchi 20 : tweak rotamer chi angles randomly by up to 20 degrees 
# -chithreshold 0.4 : 40% chance to tweak rotamer chi angle
# -tweakxyz 0.2 : as above but with xyz displacements (in this case 0.2A)
# -xyzthreshold 0.4 : as above but with xyz displacements

# THE BELOW IS NOW UNUSED - IT REMAINS IN CASE IT IS EVER USEFUL
# RUN TLEAP TO ADD BACK H (protonate could also be used and is maybe more reliable!)
#echo COMPLEX = loadpdb ${pdb:r}_newrot1.pdb > run.tleap
#echo savepdb COMPLEX afterrotamer.pdb >> run.tleap
#echo quit >> run.tleap
#$AMBERHOME/exe/tleap -s -f  $AMBERHOME/dat/leap/cmd/leaprc.glycam06  -f  $AMBERHOME/dat/leap/cmd/leaprc.ff99SB  -f run.tleap
# END OF UNUSED CODE

# EXTRACT XYZ AND LIST OF CHANGED RESIDUES
cp ${pdb:r}_newrot1.pdb afterrotamer.pdb
awk '{print $7,$8,$9} '  afterrotamer.pdb | sed '/^ *$/d' >  afterrotamer.xyz
sed -n '/### LIST OF CHANGED SIDECHAINS/,$p' beforerotamer.SRLOG | awk '/REBUILDING/ {print $4}' | awk -F= '{print $2}' > reschanged.dat

# CLEAN UP INTERMEDIATE FILES
rm beforerotamer_newrot1.pdb beforerotamer.SRLOG   

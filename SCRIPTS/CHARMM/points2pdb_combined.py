#!/usr/bin/python

# Python script to make a PDB file containing a concatenation of all the coordinates frames from a points file.
#
# NOTE this script works for a "points" file, NOT an xyz file. You must remove the number-of-atoms line,
# the following line, and the atom symbols from an xyz file before running this script.
#
# The syntax is
# <directory path>/points2pdb_combined.py <name of reference .crd file> <name of points file>
#
# Output is <name of points file>.pdb
#
# Original script by Mey Khalili
# Edits by Joanne Carr (jmc49@cam.ac.uk) to include a chain ID in the output PDB file. The first character 
# of the 8th column of data is taken as the chain ID if this column is present; otherwise you'll get 'A's.

import sys
from string import ljust

def prettyprint(fst,snd):
	global out
	a=int(fst[0])
	print>> out,'ATOM',
	el = snd.split()
	print>>out, "%6d%5s%5s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"%(int(fst[0]), ljust(fst[3],3),ljust(fst[2],4), ljust(fst[4][0],1), int(fst[1]), float(el[0]), float(el[1]), float(el[2]), 1.0, 0.0)

seg=[]
chainid=[]
inp1 = open(sys.argv[1],'r')
out=open(sys.argv[2]+".pdb",'w')

for each in inp1:
	if each.startswith('*'):
                # ignore comment lines
		pass
	else:
		els=each.split()
		if len(els) > 7 :
                # there is a segid in the crd file
			seg.append(els[0:4]+els[7:8])
		elif len(els) < 8 and len(els) > 1 :
                # there isn't a segid in the crd file, so use a dummy character
			seg.append(els[0:4]+['A'])
		else :
                # it's the number of atoms line
			seg.append(els[0:4])
inp1.close()

# get the number of atoms
cntr=map(int,seg[0])

inp2 = open(sys.argv[2],'r')
count=0

for each in inp2:
	count+=1    	
	prettyprint(seg[count],each)
	if count >=cntr[0]:
		print >> out, "END"	
		count=0
out.close()
inp2.close()

#!/usr/bin/python

# Little script to convert output files from charmm/OPTIM to strict CRD format for use in VMD.
# Note that the .cor file extension is best for use with VMD, so that it doesn't get confused 
# between charmm and amber. Open as "vmd newoutput.cor"

# Usage: 
# <directory path>/standardize_crd.py <input filename>

import string,sys,math
from string import ljust

if len(sys.argv)<1:
        print 'Need to give the input filename as the argument'
        sys.exit()

inp=open(sys.argv[1],'r')
out=open('newoutput.cor','w')

for each in inp:
	els=each.split()
	if each.startswith('*') :
		pass
	elif len(els) == 1 :
		print>> out, "%5d"%(int(els[0]))
	else:
		print>> out, "%5d%5d%5s%5s%10.5f%10.5f%10.5f%5s%5s%10.5f"%(int(els[0]), int(els[1]), ljust(' '+els[2],5), ljust(' '+els[3],5), float(els[4]), float(els[5]), float(els[6]), ljust(' '+els[7],5), ljust(' '+els[1],5), 0.0)

# From the charmm documentation:
#         ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
#          I5    I5  1X A4 1X A4 F10.5 F10.5 F10.5 1X A4 1X A4 F10.5

inp.close()
out.close()

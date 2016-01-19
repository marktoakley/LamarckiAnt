#!/usr/bin/env python

import os
import os.path
import sys
import string

###############################################################
##                                                            #
## Edyta Malolepsza                                           #
## David Wales' group, University of Cambridge                #
## in case of problems please send email: em427@cam.ac.uk     #
##                                                            #
###############################################################
##                                                            #
## program finds in prmtop file from LEaP wrong defined order #
## of atoms in IMPROPER, permutes appropriate atoms and write #
## new prmtop file                                            #
##                                                            #
## how to use:                                                #
## ./perm-top.py NAME_OF_OLD_PRMTOP NAME_OF_NEW_PRMTOP        #
##                                                            #
## IMPORTANT:                                                 #
## 1. please change names of terminal amino acid residues     #
##    according to warnings below                             #
## 2. please change path to libraries                         #
## 3. program changes the atom order ONLY for amino acid and  #
##    nucleic residues                                        #
##                                                            #
###############################################################

# khs26> changed the path to use the $AMBERHOME environment variable
amberhome = os.environ["AMBERHOME"]
path = os.path.join(amberhome, "dat/leap/lib")

#########################
## some useful functions
#########################

def exchange_atoms(atom_type, a, aa, residue, dihedrals, currentAtomNumber):
  find_atom = a[aa.index(residue)].index(atom_type)
  atomNumber = find_atom+currentAtomNumber
  atomNumberIndex = atomNumber*3
  for j in range(len(dihedrals)):
    if (dihedrals[j][1]==str(atomNumberIndex)):
      d1 = dihedrals[j][0]
      d2 = dihedrals[j][1]
      dihedrals[j][0] = d2
      dihedrals[j][1] = d1

def exchange_atoms_nt(atom_type, a, aa, residue, dihedrals):
  find_atom = a[aa.index(residue)].index(atom_type)
  for j in range(len(dihedrals)):
    if (dihedrals[j][1]==str(atomIndex[find_atom])):
      d1 = dihedrals[j][0]
      d2 = dihedrals[j][1]
      dihedrals[j][0] = d2
      dihedrals[j][1] = d1

def exchange_atoms_arg(a, aa, residue, dihedrals, currentAtomNumber):
      ## IMPROPER responsible for trouble with NH2 group permutation: 
      find_atom1 = a[aa.index(residue)].index('NE')
      atomNumber1 = find_atom1+currentAtomNumber
      atomNumberIndex1 = atomNumber1*3
      find_atom2 = a[aa.index(residue)].index('NH1')
      atomNumber2 = find_atom2+currentAtomNumber
      atomNumberIndex2 = atomNumber2*3
      find_atom3 = a[aa.index(residue)].index('CZ')
      atomNumber3 = find_atom3+currentAtomNumber
      atomNumberIndex3 = atomNumber3*3
      find_atom4 = a[aa.index(residue)].index('NH2')
      atomNumber4 = find_atom4+currentAtomNumber
      atomNumberIndex4 = atomNumber4*3
      for j in range(len(dihedrals)):
        if ((dihedrals[j][0]==str(atomNumberIndex1)) and (dihedrals[j][1]==str(atomNumberIndex2))):
          d0 = dihedrals[j][0]
          d1 = dihedrals[j][1]
          dihedrals[j][0] = d1
          dihedrals[j][1] = d0

def exchange_atoms_ring1(a, aa, residue, dihedrals):
  find_atom1 = a[aa.index(residue)].index('CD1')
  atomNumber1 = find_atom1+currentAtomNumber
  atomNumberIndex1 = atomNumber1*3
  find_atom2 = a[aa.index(residue)].index('CD2')
  atomNumber2 = find_atom2+currentAtomNumber
  atomNumberIndex2 = atomNumber2*3
  for j in range(len(dihedrals)):
    if ((dihedrals[j][0]==str(atomNumberIndex1)) and (dihedrals[j][1]==str(atomNumberIndex2))):
      d0 = '-'+dihedrals[j][0]
      d1 = dihedrals[j][1]
      d3 = dihedrals[j][3][1:]
      dihedrals[j][0] = d1
      dihedrals[j][1] = d3
      dihedrals[j][3] = d0

def exchange_atoms_ring2(a, aa, residue, dihedrals):
  find_atom1 = a[aa.index(residue)].index('CG')
  atomNumber1 = find_atom1+currentAtomNumber
  atomNumberIndex1 = atomNumber1*3
  find_atom2 = a[aa.index(residue)].index('CE2')
  atomNumber2 = find_atom2+currentAtomNumber
  atomNumberIndex2 = atomNumber2*3
  find_atom3 = a[aa.index(residue)].index('CZ')
  atomNumber3 = find_atom3+currentAtomNumber
  atomNumberIndex3 = atomNumber3*3
  find_atom4 = a[aa.index(residue)].index('CD2')
  atomNumber4 = find_atom4+currentAtomNumber
  atomNumberIndex4 = atomNumber4*3
  find_atom5 = a[aa.index(residue)].index('CD1')
  atomNumber5 = find_atom5+currentAtomNumber
  atomNumberIndex5 = atomNumber5*3
  find_atom6 = a[aa.index(residue)].index('CE1')
  atomNumber6 = find_atom6+currentAtomNumber
  atomNumberIndex6 = atomNumber6*3
#  for j in range(len(dihedrals)): # this is ok
#    if ((dihedrals[j][0]==str(atomNumberIndex1)) and (dihedrals[j][1]==str(atomNumberIndex2))):
  for j in range(len(dihedrals)):
    if ((dihedrals[j][0]==str(atomNumberIndex3)) and (dihedrals[j][1]==str(atomNumberIndex4))):
      d1 = '-'+dihedrals[j][1]
      d2 = dihedrals[j][3][1:]
      dihedrals[j][1] = d2
      dihedrals[j][3] = d1
  for j in range(len(dihedrals)):
    if ((dihedrals[j][0]==str(atomNumberIndex5)) and (dihedrals[j][1]==str(atomNumberIndex3))):
      d1 = '-'+dihedrals[j][1]
      d2 = dihedrals[j][3][1:]
      dihedrals[j][1] = d2
      dihedrals[j][3] = d1
  for j in range(len(dihedrals)):
    if ((dihedrals[j][0]==str(atomNumberIndex1)) and (dihedrals[j][1]==str(atomNumberIndex6))):
      ## to compare IMPROPER before and after permutation
##test      a1 = (int(dihedrals[j][0])-currentAtomNumber)/3
##test      a2 = (int(dihedrals[j][1])-currentAtomNumber)/3
##test      a3 = (int(dihedrals[j][2][1:])-currentAtomNumber)/3
##test      a4 = (int(dihedrals[j][3][1:])-currentAtomNumber)/3
##test      print dihedrals[j], a[aa.index(residue)][a1], a[aa.index(residue)][a2], a[aa.index(residue)][a3], a[aa.index(residue)][a4]
      d1 = '-'+dihedrals[j][0]
      d2 = dihedrals[j][3][1:]
      dihedrals[j][0] = d2
      dihedrals[j][3] = d1
##test      a1 = (int(dihedrals[j][0])-currentAtomNumber)/3
##test      a2 = (int(dihedrals[j][1])-currentAtomNumber)/3
##test      a3 = (int(dihedrals[j][2][1:])-currentAtomNumber)/3
##test      a4 = (int(dihedrals[j][3][1:])-currentAtomNumber)/3
##test      print dihedrals[j], a[aa.index(residue)][a1], a[aa.index(residue)][a2], a[aa.index(residue)][a3], a[aa.index(residue)][a4]

def exchange_atoms_ring3(a, aa, residue, dihedrals):
  find_atom1 = a[aa.index(residue)].index('CE1')
  atomNumber1 = find_atom1+currentAtomNumber
  atomNumberIndex1 = atomNumber1*3
  find_atom2 = a[aa.index(residue)].index('CE2')
  atomNumber2 = find_atom2+currentAtomNumber
  atomNumberIndex2 = atomNumber2*3
  for j in range(len(dihedrals)):
    if ((dihedrals[j][0]==str(atomNumberIndex1)) and (dihedrals[j][1]==str(atomNumberIndex2))):
      d0 = '-'+dihedrals[j][0]
      d1 = dihedrals[j][1]
      d3 = dihedrals[j][3][1:]
      dihedrals[j][0] = d1
      dihedrals[j][1] = d3
      dihedrals[j][3] = d0

####################################
## reading all_amino02.lib library
####################################

print '\nDear user, please notice that only residues from the following libraries are taken into account:'
print '  ions94.lib'

print '  all_amino02.lib'
aalib = open("%s/all_amino02.lib" % path).read()
aa = string.split(aalib, "\n")
q1 = aa.index("!!index array str")
q2 = aa.index("!entry.ALA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg")

aaNames = [] # amino acid names
aTypes = []  # atom types
aNames = []  # atom names

for i in range(q2-q1-1):
  aaNames.append(aa[q1+1+i][2:5])

for i in range(len(aaNames)):
  q1 = aa.index("!entry.%s.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" % aaNames[i])
  q2 = aa.index("!entry.%s.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" % aaNames[i])
  aT = []
  aN = []
  for j in range(q2-q1-1):
    aT.append((string.split(aa[q1+1+j])[0]).replace('"',''))
    aN.append((string.split(aa[q1+1+j])[1]).replace('"',''))
  aTypes.append(aT)
  aNames.append(aN)

######################################
## reading all_aminont02.lib library
######################################

print '  all_aminont02.lib'
aalib = open("%s/all_aminont02.lib" % path).read()
aa = string.split(aalib, "\n")
q1 = aa.index("!!index array str")
q2 = aa.index("!entry.ACE.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg")

aantNames = [] # N terminus amino acid names
antTypes = []  # N terminus atom types
antNames = []  # N terminus atom names

for i in range(q2-q1-1):
  aantNames.append((aa[q1+1+i].replace('"','')).replace(' ',''))

for i in range(len(aantNames)):
  q1 = aa.index("!entry.%s.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" % aantNames[i])
  q2 = aa.index("!entry.%s.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" % aantNames[i])
  aT = []
  aN = []
  for j in range(q2-q1-1):
    aT.append((string.split(aa[q1+1+j])[0]).replace('"',''))
    aN.append((string.split(aa[q1+1+j])[1]).replace('"',''))
  antTypes.append(aT)
  antNames.append(aN)

######################################
## reading all_aminoct02.lib library
######################################

print '  all_aminoct02.lib'
aalib = open("%s/all_aminoct02.lib" % path).read()
aa = string.split(aalib, "\n")
q1 = aa.index("!!index array str")
q2 = aa.index("!entry.CALA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg")

aactNames = [] # C terminus amino acid names
actTypes = []  # C terminus atom types
actNames = []  # C terminus atom names

for i in range(q2-q1-1):
  aactNames.append((aa[q1+1+i].replace('"','')).replace(' ',''))

for i in range(len(aactNames)):
  q1 = aa.index("!entry.%s.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" % aactNames[i])
  q2 = aa.index("!entry.%s.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" % aactNames[i])
  aT = []
  aN = []
  for j in range(q2-q1-1):
    aT.append((string.split(aa[q1+1+j])[0]).replace('"',''))
    aN.append((string.split(aa[q1+1+j])[1]).replace('"',''))
  actTypes.append(aT)
  actNames.append(aN)

#####################################
## reading all_nucleic02.lib library
#####################################

print '  all_nucleic02.lib'
aalib = open("%s/all_nucleic02.lib" % path).read()
aa = string.split(aalib, "\n")
q1 = aa.index("!!index array str")
q2 = aa.index("!entry.DA.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg")

nucNames = [] # nucleic names
nucTypes = []  # nucleic atom types
nucaNames = []  # nucleic atom names

for i in range(q2-q1-1):
  nucNames.append((aa[q1+1+i].replace('"','')).replace(' ',''))

for i in range(len(nucNames)):
  q1 = aa.index("!entry.%s.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" % nucNames[i])
  q2 = aa.index("!entry.%s.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" % nucNames[i])
  aT = []
  aN = []
  for j in range(q2-q1-1):
    aT.append((string.split(aa[q1+1+j])[0]).replace('"',''))
    aN.append((string.split(aa[q1+1+j])[1]).replace('"',''))
  nucTypes.append(aT)
  nucaNames.append(aN)

#################################
## reading original prmtop file
#################################

prmtop = open(sys.argv[1]).read()
f = string.split(prmtop, "\n")

q0 = f.index("%FLAG POINTERS                                                                  ")
q1 = f.index("%FLAG ATOM_NAME                                                                 ")
q2 = f.index("%FLAG CHARGE                                                                    ")
q3 = f.index("%FLAG RESIDUE_LABEL                                                             ")
q4 = f.index("%FLAG RESIDUE_POINTER                                                           ")
q5 = f.index("%FLAG DIHEDRALS_INC_HYDROGEN                                                    ")
q6 = f.index("%FLAG DIHEDRALS_WITHOUT_HYDROGEN                                                ")
q7 = f.index("%FLAG EXCLUDED_ATOMS_LIST                                                       ")

## names of tables are related to names in prmtop file

atomNumber = int(string.split(f[q0+2])[0])

atomName = []
residueLabel = []
dihedralsIncHydrogen = []
dihedralsWithoutHydrogen = []
atomIndex = []

an = 0
line = 0
while (an<atomNumber):
  for j in range(20):
    if (an<atomNumber):
      an = an+1
      atomName.append(f[q1+2+line][j*4:(j+1)*4].strip())
    else:
      break
  line = line+1

for i in range(q4-q3-2):
  for j in range((len(f[q3+2+i])+1)/4):
    residueLabel.append(string.strip(f[q3+2+i][j*4:4*(j+1)]))

caa = 0
naa = 0
for i in range(len(residueLabel)):
  if (aactNames.count(residueLabel[i])>0): caa = caa+1
  if (aantNames.count(residueLabel[i])>0): naa = naa+1

if (caa==0):
  print "\n-----------------------------------------------------------------------------"
  print 'There is no C terminus amino acid in topology file!'
  print 'If system does not contain amino acids - continue.'
  print 'Otherwise please rename each C terminal residue by adding \'C\' to the name, e.g. ALA -> CALA'
  print 'Remember to follow format of topology file!'
  print "-----------------------------------------------------------------------------\n"
  ## the only exception for C terminus amino acids is NME

if (naa==0):
  print "-----------------------------------------------------------------------------"
  print 'There is no N terminus amino acid in topology file!'
  print 'If system does not contain amino acids - continue.'
  print 'Otherwise please rename each N terminal residue by adding \'N\' to the name, e.g. ALA -> NALA'
  print 'Remember to follow format of topology file!'
  print "-----------------------------------------------------------------------------"
  ## the only exception for N terminus amino acids is ACE

for i in range(q6-q5-2):
  for j in range(len(string.split(f[q5+2+i]))/5):
    dihedralsIncHydrogen.append(string.split(f[q5+2+i][j*40:40*(j+1)]))

for i in range(q7-q6-2):
  for j in range(len(string.split(f[q6+2+i]))/5):
    dihedralsWithoutHydrogen.append(string.split(f[q6+2+i][j*40:40*(j+1)]))

for i in range(len(atomName)):
  atomIndex.append(i*3)

############################################################
## groups of amino acids according to permutation behaviour
############################################################

## group0 and group0n: nothing to do
## group1: problem only with C terminus COO- group: CA-O-C-OXT -> O-CA-C-OXT

group0 = ['GLY','ALA','VAL','LEU','MET','ILE','SER','THR','CYS','LYS','HIP','HIE','HID']
group0n = ['NGLY','NALA','NVAL','NLEU','NMET','NILE','NSER','NTHR','NCYS','NLYS','NHIP','NHIE','NHID']
group1 = ['CGLY','CALA','CVAL','CLEU','CMET','CILE','CSER','CTHR','CCYS','CLYS','CPRO','CTRP','CHIP','CHIE','CHID']

#################################################################
## groups of nucleic residues according to permutation behaviour
#################################################################

group2 = ['DA', 'DA3', 'DA5', 'DAN', 'RA', 'RA3', 'RA5', 'RAN']
group3 = ['DC', 'DC3', 'DC5', 'DCN', 'RC', 'RC3', 'RC5', 'RCN']
group4 = ['DG', 'DG3', 'DG5', 'DGN', 'RG', 'RG3', 'RG5', 'RGN']
group5 = ['DT', 'DT3', 'DT5', 'DTN', 'RU', 'RU3', 'RU5', 'RUN'] ## nothing to do

#####################################################
## groups of species without any IMPROPER to correct
#####################################################

group7 = ['WAT', 'CIO', 'Cl-', 'Cs+', 'IB', 'K+', 'Li+', 'MG2', 'Na+', 'Rb+']

## groupUknown: residue not counted in library
groupUknown = []

for i in range(len(residueLabel)):
  if ((aaNames.count(residueLabel[i])==0) and (aactNames.count(residueLabel[i])==0) and (aantNames.count(residueLabel[i])==0) and (nucNames.count(residueLabel[i])==0) and (group7.count(residueLabel[i])==0)):
    groupUknown.append(residueLabel[i])

if (len(groupUknown)!=0): 
  print '\nThere are some residues missing in considered libraries:', groupUknown
  print 'Program just skips them\n'

currentAtomNumber = 0

######################################################################
## main part - permutation of atom positions in appropriate IMPROPERs
######################################################################

for i in range(len(residueLabel)):
#  print '----', i+1, '----', residueLabel[i]
  if (group7.count(residueLabel[i])>0): continue

  elif (groupUknown.count(residueLabel[i])>0): continue
    ##################################
    ## residue not counted in library
    ##################################

  elif (aantNames.count(residueLabel[i])>0):
    #########################
    ## N terminus amino acid
    #########################

    if (group0n.count(residueLabel[i])>0): 
      currentAtomNumber = currentAtomNumber+len(antTypes[aantNames.index(residueLabel[i])])
      continue
    elif (residueLabel[i]=='NASN'):
      exchange_atoms_nt('HD21', antTypes, aantNames, residueLabel[i], dihedralsIncHydrogen)
    elif (residueLabel[i]=='NGLN'):
      exchange_atoms_nt('HE21', antTypes, aantNames, residueLabel[i], dihedralsIncHydrogen)
    elif (residueLabel[i]=='NARG'):
      exchange_atoms_nt('HH11', antTypes, aantNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms_nt('HH21', antTypes, aantNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms_arg(antTypes, aantNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='NASP'):
      exchange_atoms_nt('OD1', antTypes, aantNames, residueLabel[i], dihedralsWithoutHydrogen)
    elif (residueLabel[i]=='NGLU'):
      exchange_atoms_nt('OE1', antTypes, aantNames, residueLabel[i], dihedralsWithoutHydrogen)
    elif (residueLabel[i]=='NPHE'):
      exchange_atoms_ring1(antTypes, aantNames, residueLabel[i], dihedralsWithoutHydrogen)
      exchange_atoms_ring2(antTypes, aantNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms_ring3(antTypes, aantNames, residueLabel[i], dihedralsIncHydrogen)
    elif (residueLabel[i]=='NTYR'):
      exchange_atoms_ring1(antTypes, aantNames, residueLabel[i], dihedralsWithoutHydrogen)
      exchange_atoms_ring2(antTypes, aantNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms_ring3(antTypes, aantNames, residueLabel[i], dihedralsWithoutHydrogen)
#        if (dihedralsWithoutHydrogen[j].count(str(atomNumberIndex1))>0):
#          print '"""', dihedralsWithoutHydrogen[j], atomNumberIndex1, dihedralsWithoutHydrogen[j].count(str(atomNumberIndex1))
#          print '!!!', dihedralsWithoutHydrogen[j], atomNumberIndex1, dihedralsWithoutHydrogen[j].index(str(atomNumberIndex1))

    currentAtomNumber = currentAtomNumber+len(antTypes[aantNames.index(residueLabel[i])])

  elif (aactNames.count(residueLabel[i])>0):
    #########################
    ## C terminus amino acid
    #########################

    if (group1.count(residueLabel[i])>0): ## res belongs to group1
      exchange_atoms('O', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='CASP'):
      exchange_atoms('OD1', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='CGLU'):
      exchange_atoms('OE1', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
      exchange_atoms('O', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='CGLN'):
      exchange_atoms('HE21', actTypes, aactNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
      exchange_atoms('O', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='CASN'):
      exchange_atoms('HD21', actTypes, aactNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
      exchange_atoms('O', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='CARG'):
      exchange_atoms('HH11', actTypes, aactNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
      exchange_atoms('HH21', actTypes, aactNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
      exchange_atoms('O', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
      exchange_atoms_arg(actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='CPHE'):
      exchange_atoms_ring1(actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen)
      exchange_atoms_ring2(actTypes, aactNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms_ring3(actTypes, aactNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms('O', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='CTYR'):
      exchange_atoms_ring1(actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen)
      exchange_atoms_ring2(actTypes, aactNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms('O', actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
      exchange_atoms_ring3(actTypes, aactNames, residueLabel[i], dihedralsWithoutHydrogen)

    currentAtomNumber = currentAtomNumber+len(actTypes[aactNames.index(residueLabel[i])])

  elif (aaNames.count(residueLabel[i])>0):
    ###########################
    ## not terminal amino acid
    ###########################

    if (group0.count(residueLabel[i])>0): 
      currentAtomNumber = currentAtomNumber+len(aTypes[aaNames.index(residueLabel[i])])
      continue
    elif (residueLabel[i]=='GLU'):
      exchange_atoms('OE1', aTypes, aaNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='ASP'):
      exchange_atoms('OD1', aTypes, aaNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='ASN'):
      exchange_atoms('HD21', aTypes, aaNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='GLN'):
      exchange_atoms('HE21', aTypes, aaNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='ARG'):
      exchange_atoms('HH11', aTypes, aaNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
      exchange_atoms('HH21', aTypes, aaNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
      exchange_atoms_arg(aTypes, aaNames, residueLabel[i], dihedralsWithoutHydrogen, currentAtomNumber)
    elif (residueLabel[i]=='PHE'):
      exchange_atoms_ring1(aTypes, aaNames, residueLabel[i], dihedralsWithoutHydrogen)
      exchange_atoms_ring2(aTypes, aaNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms_ring3(aTypes, aaNames, residueLabel[i], dihedralsIncHydrogen)
    elif (residueLabel[i]=='TYR'):
      exchange_atoms_ring1(aTypes, aaNames, residueLabel[i], dihedralsWithoutHydrogen)
      exchange_atoms_ring2(aTypes, aaNames, residueLabel[i], dihedralsIncHydrogen)
      exchange_atoms_ring3(aTypes, aaNames, residueLabel[i], dihedralsWithoutHydrogen)

    currentAtomNumber = currentAtomNumber+len(aTypes[aaNames.index(residueLabel[i])])

  elif (nucNames.count(residueLabel[i])>0):
    ###################
    ## nucleic residue
    ###################

    if (group5.count(residueLabel[i])>0):
      currentAtomNumber = currentAtomNumber+len(nucTypes[nucNames.index(residueLabel[i])])
      continue
    elif (group2.count(residueLabel[i])>0):
      exchange_atoms('H61', nucTypes, nucNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
    elif (group3.count(residueLabel[i])>0):
      exchange_atoms('H41', nucTypes, nucNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)
    elif (group4.count(residueLabel[i])>0):
      exchange_atoms('H21', nucTypes, nucNames, residueLabel[i], dihedralsIncHydrogen, currentAtomNumber)

    currentAtomNumber = currentAtomNumber+len(nucTypes[nucNames.index(residueLabel[i])])

  else:
    print 'Something strange happened... residue %d is neither in libraries and nor out of libraries' % residueLabel[i]
    sys.exit()

##################################
## preparation of new prmtop file
##################################

newprmtop = open(sys.argv[2],'w')

for i in range(q5+2):
  newprmtop.write("%s\n" % f[i])

an = 0
line = 0
while (an<len(dihedralsIncHydrogen)):
  for j in range(2):
    if (an<len(dihedralsIncHydrogen)):
      for k in range(5):
        newprmtop.write("%8s" % (dihedralsIncHydrogen[an][k]))
      an = an+1
    else:
      break
  newprmtop.write("\n")
  line = line+1

for i in range(2):
  newprmtop.write("%s\n" % f[q6+i])

an = 0
line = 0
while (an<len(dihedralsWithoutHydrogen)):
  for j in range(2):
    if (an<len(dihedralsWithoutHydrogen)):
      for k in range(5):
        newprmtop.write("%8s" % (dihedralsWithoutHydrogen[an][k]))
      an = an+1
    else:
      break
  newprmtop.write("\n")
  line = line+1

for i in range(len(f)-q7-1):
  newprmtop.write("%s\n" % f[q7+i])

newprmtop.close()

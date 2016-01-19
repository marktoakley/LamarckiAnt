#! /usr/bin/env python
import sys

# Generate group rotations for cyclic peptides
# Should work for all the naturally occurring amino acids (maybe not PRO)
# If optional 2nd argument supplied, this works out groups for multi-ring systems

if (len(sys.argv)<2 or len(sys.argv)>3):
    print"Usage: CyclicPeptideGroupRotate.py <filename> <rings>"
    print"    <filename> a .pdb file"
    print"    <rings> number of cyclic peptide rings (optional, default=1)"
    sys.exit("Wrong number of arguments supplied")

inp = open(sys.argv[1])

if len(sys.argv)==2: rings=1
else: rings=int(sys.argv[2])

#out=open('atomgroups','w')

class Atom:
  name=""
  index =0
  acidname=""
  def __str__(this):
    return this.name
  def __init__(this):
    this.name = ''
  def __str__(this):
    return this.name+this.acidname+" %d "%(this.index)
  def __init__(this):
    this.name = ''
    this.index = 0
    this.acidname=''
 
  def copy(this):
    n = Atom()
    n.name=this.name
    n.acidname=this.acidname
    n.atom=this.index
    return n

class RigidBody:
  name = ""
  start = 0
  finish = 0
  atoms = []
  scale = 1.0
  prob = 0.5

# This function reads and stores atom name, its corresponding residue name and its index
def readatom(s):
  a = Atom()
  a.acidname=s[17:20].strip()
  a.name = s[12:16].strip()
  a.index = s[6:11].strip()
  return a
  
def makeRigidBody(i,size):
  GRP=GRPlist[i]
  rb = RigidBody()
  rb.atoms=[]
  rb.start=GRPlist[i-1][0].index
  for j in range (i, i+size):
    if (i//ringsize==j//ringsize):
      k=j
    else:
      k=j-ringsize
    if (i//ringsize==(j+1)//ringsize):
      l=j+1
    else:
      l=j+1-ringsize      
    GRP=GRPlist[k]
    rb.atoms.extend(GRP)    
  rb.finish=GRPlist[l][0].index
  return rb

# Set up the list of groups
# The first item is an atom in the backbone
# All others are the non-backbone atoms appended to it

global GRPlist
global ringsize
global PheList
GRPlist=[]
GRP=[]
PheList=[]
for each in inp:
  els = each.split()
  if len(els)==1: continue
  if (els[0]=='TER'): continue
  if (els[0]=='REMARK'): continue
  AM=readatom(each)
  if AM.name in ['CA']:
    if AM.acidname in ['PHE']:
       PheList.append(AM.index)
  if AM.name in ['CA','C','N']:
    if len(GRP)>0:
      GRPlist.append(GRP)
      GRP=[]
  GRP.append(AM)
GRPlist.append(GRP)
ringsize=len(GRPlist)/rings
  
# We now have a list of groups
# Now we need to combine into rigid bodies
RBlist=[]

# Rotations of peptide bonds
for i in range(len(GRPlist)):
  GRP=GRPlist[i]
  if GRP[0].name=='C':
    rb = makeRigidBody(i,2)
    rb.name="PEP"
    rb.name += str(i)
    rb.scale=1.0
    rb.prob=0.5
    RBlist.append(rb)    

# Cis/trans isomers
for i in range(len(GRPlist)):
  GRP=GRPlist[i]
  if GRP[0].name=='C':
    rb = makeRigidBody(i,1)
    rb.name="CT"
    rb.name += str(i)
    rb.scale=1.0
    rb.prob=0.25
    RBlist.append(rb)

# Corner flips
for i in range(len(GRPlist)):
  GRP=GRPlist[i]
  if GRP[0].name=='C':
    rb = makeRigidBody(i,5)
    rb.name="CORNER"
    rb.name+=str(i)
    rb.scale=1.0
    rb.prob=0.2  
    RBlist.append(rb)

for RB in RBlist:
  print "GROUP ",RB.name,RB.start,RB.finish,len(RB.atoms),RB.scale,RB.prob
  for AM in RB.atoms:
    print AM.index
#out.close()

#Side chains
for i in range(len(PheList)):
  ii=int(PheList[i])
  name ="PHEA"+str(i)
  print "GROUP",name,ii,ii+2,14,1.0,0.2
  for j in range (ii+2,ii+16):
    print j
  name ="PHEB"+str(i)
  print "GROUP",name,ii+2,ii+5,11,1.0,0.2
  for j in range (ii+5,ii+16):
    print j

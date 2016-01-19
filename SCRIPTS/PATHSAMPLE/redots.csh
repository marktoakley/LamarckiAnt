#!/bin/csh

# rm OPTIM.* vector.dump.* path.info*
set ts=$1

mkdir $1
rm $1/eigmin
cp -r input.crd perm.allow $1
cd $1
ln -s ../points.ts points.ts
ln -s ../ts.data ts.data
ln -s ../min.A min.A
ln -s ../min.B min.B
ln -s ../min.data min.data
cp ../odata.ts odata

while ($ts <= $2 )

   echo CHARMM > pathdata
   echo NATOMS 215 >> pathdata
   echo COPYFILES  input.crd perm.allow >> pathdata
   echo COPYOPTIM >> pathdata
   echo EXTRACTTS $ts >> pathdata
   ../pathsample.2.1 >& pathsample.out
   cp extractedts start
   (time ~/bin/COPTIM.4.0) >& OPTIM.ts
   grep Smallest OPTIM.ts | tail -1 >> eigmin
   echo $ts `grep -c CONV OPTIM.ts`

   @ ts +=1

end

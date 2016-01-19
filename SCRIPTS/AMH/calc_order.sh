#!/bin/bash

numb=`wc min.data | awk '{print$1}'`
dumb=`cat pathdata.back | grep NATOMS | awk '{print$2}'`

echo $numb
echo $dumb

for  ((i=1; i<=$numb; i+=1))
  do
   echo $i out of $numb
   echo AMH            > pathdata
   echo NATOMS $dumb  >> pathdata
   echo SYSTEM AMH    >> pathdata
   echo EXTRACTMIN $i >> pathdata
   echo AMHQ $i       >> pathdata
   ./pathsample.2.1   >& ORDER_PARAM/pathsample.out.$i
  done


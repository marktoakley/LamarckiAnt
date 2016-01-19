#!/bin/sh

count=0
cores=200
start=1
finish=42000

while [ $count -lt $cores ] ; do

   count=$((count+1))
   mkdir order.$count
   cd order.$count
   echo core start and finish are $core $start $finish
   cp ../qsub.template qsub.calcorder
   cp ../pathdata.template pathdata
   echo CALCORDER $start $finish >> pathdata 
   qsub qsub.calcorder
   cd ..
   start=$((start+42000))
   finish=$((finish+42000))

done


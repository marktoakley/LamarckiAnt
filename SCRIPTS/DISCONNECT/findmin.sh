#!/bin/bash

#Usage: To get the identities of minima of interest from the disconnectivity graph.

#Required input: tree.ps, the identities of two minima corresponding to the 
#very left and very right minima. All the minima lying in between these two
#minima will be printed out

#./findmin.sh min_up min_down min_left min_right

if [ "$2" == "" ]
then
  echo ""
  echo 'Please give the identities of minima. Usage: ./findmin.sh min_up min_down min_left min_right'
  echo ""
  exit
fi
temp=`grep "mt ($1) show" tree.ps|awk '{print $2}'`
up=`echo "$temp*100/1"|bc`
temp=`grep "mt ($2) show" tree.ps|awk '{print $2}'`
down=`echo "$temp*100/1"|bc`
temp=`grep "mt ($3) show" tree.ps|awk '{print $1}'`
left=`echo "$temp*100/1"|bc`
temp=`grep "mt ($4) show" tree.ps|awk '{print $1}'`
right=`echo "$temp*100/1"|bc`
j=1
k=0
for i in `grep show tree.ps|grep mt|awk '{print $1 " " $2 " " $4}'`
do
  if [ "$j" -eq "1" ]
  then
    j=2
    a=`echo "$i*100/1"|bc`
    if [ "$a" -ge "$left" ] && [ "$a" -le "$right" ]
    then
      k=1
    fi
  elif [ "$j" -eq "2" ]
  then
    j=3
    if [ "$k" -eq 1 ]
    then
      a=`echo "$i*100/1"|bc`
      if [ "$a" -le "$up" ] && [ "$a" -ge "$down" ]
      then
        k=2
      fi
    fi
  else
    j=1
    if [ "$k" -eq 2 ]
    then
      echo $i|sed 's/(//g'|sed 's/)//g'
    fi
    k=0
  fi
done

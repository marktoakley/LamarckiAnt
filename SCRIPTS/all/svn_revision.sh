#!/bin/bash

# directory where this script resides
export shd="`dirname $(readlink -f $0)`"
export wg_dir="$shd/../../GMIN/source"
#export inc_dir="$shd/../../INCLUDE"

dir="."
[ ! -z $1 ] && dir="$1" 
#svnversion $dir | sed 's/.*://' | sed 's/M//' > version.tmp

svn info $wg_dir >& n  
s=`cat n | awk '/Revision:/ { print $2 }'`
rm n
[ ! -z "$s" ] && echo "$s" > SVNREV

if [ -f SVNREV ]; then
 	 cat SVNREV
else
 	exit 1
fi

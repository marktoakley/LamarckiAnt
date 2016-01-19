#!/bin/bash
# Chris Whittleston 2012 (csw34@cam.ac.uk)
# visualise_UNTRAP.sh
# Script to add IDMIN lines to the dinfo file from PATHSAMPLE output. This may be 
# useful when tuning the UNTRAP options to select pairs of minima to connect.
#
# Using the script:
# 1. Choose UNTRAP options (see documentation) and run a dummy PATHSAMPLE run (use DUMMYRUN)
# 2. Kill the job once it has produced a list of minima to reconnect (look for 'getupair> sorted list of')
# 3. Run this script
# 4. Run disconnectionDPS
# 5. Look at the tree using gv and make sure the minima highlighted are sensible!
# 6. Remove DUMMYRUN and start UNTRAPping!

# Arguements:
# $1 = PATHSAMPLE output file name
# $2 = number of minima to visualise

# Check correct number of arguements given
EXPECTED_ARGS=2
if [ $# -lt $EXPECTED_ARGS ]
then
   echo "Usage: ./visualise_UNTRAP.sh <PATHSAMPLE output file> <number of minima to visualise>"
   exit
fi

# Extract the first column of minima to be connected 
grep -A $2 'getupair> sorted list of' $1 | tail -n $2 | awk '{print $1}' > minima.tmp

# Concatonate on the second column 
grep -A $2 'getupair> sorted list of' $1 | tail -n $2 | awk '{print $2}' >> minima.tmp

# OPTIONAL - UNCOMMENT IF YOU WANT TO RETAIN EXISTING IDMIN LINES 
# Find the minima already identified in dinfo and add them also
#grep 'IDMIN' dinfo | awk '{print $2}' >> minima.tmp

# Sort minima numerically and remove duplicates (quite likely!)
sort -n minima.tmp | uniq > minima.dat
rm minima.tmp

# Remove all existing IDMIN lines in dinfo
cp dinfo dinfo.old
sed '/IDMIN/d' dinfo.old > dinfo

# Loop over minima, adding lines to the dinfo file as we go. 
for min in `cat minima.dat`
do 
   echo "IDMIN ${min}" >> dinfo
done

# Done!
echo 'Done! Now run disconnectionDPS :)'

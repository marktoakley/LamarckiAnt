#!/bin/bash
# for all running jobs print the jobid and the working directory

joblist=`qstat -u $USER | awk '(/^[0-9]/){print $1}' | sed 's/[^0-9].*$//'`


for job in $joblist; do 
  printf $job 
  qstat -f $job | awk '(s==1){line2 = $1; printf("  %s%s\n", line1, line2); exit(1)}(/init_work_dir/){s=1; line1 = $3}'


done

#!/bin/bash
#===============================================================================
#
#          FILE:  vf_fill_files.sh
# 
#         USAGE:  ./vf_fill_files.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 21/11/10 15:57:32 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

for f in ${base_files[@]}
	do
		files[$i]=$d/$f
		i=$(($i+1))
done


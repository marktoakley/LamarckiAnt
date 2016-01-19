#!/bin/bash
#===============================================================================
#
#          FILE:  bu_get_mav.sh
# 
#         USAGE:  ./bu_get_mav.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 24/11/10 16:16:46 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

set_mod_home

source $mod_home/init/bash
module av >& mav
cat mav | sed '/^$/d' | sed '/^---/d' | sed '/^_M/d' >& mav.n; mv mav.n mav
modules=( ` cat mav ` )
echo "" >& mav.n
for m in ${modules[@]}; do
		echo "$m" >> mav.n 
done
mv mav.n mav


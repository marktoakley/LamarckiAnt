#!/bin/bash
#===============================================================================
#
#          FILE:  bu_split_compiler_name.sh
# 
#         USAGE:  ./bu_split_compiler_name.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 23/11/10 12:04:13 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

# short compiler name
export compiler_name=` echo $compiler | awk -F "/" '{ print $1 } ' ` 
# bits? - 32 or 64
export compiler_bits=` echo $compiler | awk -F "/" '{ print $2 } ' ` 
# compiler version
#
export compiler_version=` echo $compiler | awk -F "/" '{ print $3 } ' ` 
#

#!/bin/bash
#===============================================================================
#
#          FILE:  bu_set_mod_home.sh
# 
#         USAGE:  ./bu_set_mod_home.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 25/11/10 17:27:45 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

case "$HOSTNAME" in
  clust) mod_home=/usr/share/modules ;;
  *) mod_home=$MODULESHOME ;;
esac

export mod_home

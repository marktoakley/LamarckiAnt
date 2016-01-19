#!/bin/bash
#===============================================================================
#
#          FILE:  bu_load_compiler.sh
# 
#         USAGE:  ./bu_load_compiler.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: Oleksandr Poplavskyy (OP), op226@gmail.com
#       COMPANY: Univ. Cambridge
#       CREATED: 24/11/10 13:13:05 GMT
#      REVISION:  ---
#===============================================================================

# Treat unset variables as an error

#use set_mod_home
set_mod_home
		# play with modules {{{
		source $mod_home/init/bash
                       #compiler=$default_compiler
		mload $compiler
		export FULL_COMPILER_NAME=$compiler
		export MAKE_OPTS=$make_opts
		# }}}


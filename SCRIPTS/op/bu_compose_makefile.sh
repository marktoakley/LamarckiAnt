#!/bin/bash

# syntax: bu_compose_makefile.sh Makefile_name bu_opts

# Commandline arguments to this script:
# 	first argument - the name for the Makefile
# 	other arguments => bu_opts , ie,  bu script options
# 		(don't mix up with the original Makefile options! )

m=$1
shift

cat $m.intro
cat $m.def
cat $m.obj

for bu_opt in $@
	do
	  case "$bu_opt" in
	  	"--with-charmm") cat "$m.charmm" ;;
	  	"--with-amber9") cat "$m.amber9" ;;
	  esac
done

cat $m.$compiler_name
cat $m.rules
cat $m.dep

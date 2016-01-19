#!/bin/bash

export shd="`dirname $(readlink -f $0)`"

w0="write(file,10)"

# (A bit changed by op226) Chris's code
## for retrieving SVN version number {{{

#svnversion . | sed 's/.*://' | sed 's/M//' > version.tmp

##svnversion returns 'exported' if you have it installed but run it in a non svn directory.
##If svnversion is not installed at all, the file will be empty. In either case, we want to
##use the version number recorded previously.

#if [ "`cat version.tmp`" != "exported" ] && [ "`cat version.tmp`" != "" ]; then

## If we are working with svn however, 
## we want to overwrite the version number with the current one.

#cp version.tmp VERSION
#fi
##Remove the temporary file
#rm version.tmp
##Return the version from the VERSION file
#cat VERSION
##}}}

svn_rev=` $shd/svn_revision.sh `

vars=( "prog" "fflags" "fc_full_name" "fc_exec" "make_opts" )

# get parameters, e.g, compiler etc., from command line 
# {{{
while [ ! -z "$*"  ]; do
  	var="$1"
	for known_var in ${vars[@]}
		do
			case "$var" in
			  	"$known_var") export $var="$2" ;;
			esac
        done
	shift
done

copyright="Copyright (C) 1999-2010 David J. Wales"

case "$prog" in
  	GMIN) prog_full="A program for finding global minima" ;;
	OPTIM) prog_full="A program for optimizing geometries and calculating reaction pathways" ;;
	PATHSAMPLE) prog_full="A driver for OPTIM to create stationary point databases and perform kinetic analysis"
	;;
esac
# }}}

# send to output: subroutine display_version(file)
# {{{
cat << EOF      
subroutine display_version(file)

integer file

$w0 "==========================================="
$w0
$w0 "$prog - $prog_full"
$w0
$w0 "$copyright"
$w0
$w0 "SVN revision: $svn_rev" 
$w0 
$w0 "Compilation time: ` date `"
$w0 "Compiled by $USER@$HOSTNAME on ` uname -o` ` uname -m`"
$w0
$w0 "Compiler name:  $fc_full_name"
$w0 "Compiler executable:  $fc_exec"
$w0 "Compiler flags: $fflags"
$w0 "Command-line options passed to makefile: $make_opts "
EOF

# what does the compiler say about itself? - to be completed 
#{{{
#cat <<EOF
#$w0 
#$w0 "Compiler version info from the compiler itself:" 
#$w0
#cmd=""

#vcmds=( "-v" "-V" "--version" )

#for version_cmd in "${vcmds[@]}"
	#do
		#cmd="$cmd || $fc_exec $version_cmd" 
#done

#` $cmd ` | sed 's/^/"/g; s/$/"/g' | sed "s/^/$w0  /g" 
#EOF
#}}}

cat << EOF
$w0
$w0 "==========================================="

10 format(a)

end
EOF
# }}}

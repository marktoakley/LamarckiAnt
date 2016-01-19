#!/bin/bash

pref_eoo="$this_script"

# latex/pdf stuff {{{
# variables and options {{{
export ps2pdf_opts="-dOptimize=true -dUseFlateCompression=true \
	       -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 \
		-dSubsetFonts=true -dEmbedAllFonts=true \
		-dAutoFilterColorImages=false \
		-dAutoFilterGrayImages=false \
		-dColorImageFilter=/FlateEncode \
		-dGrayImageFilter=/FlateEncode \
		-dModoImageFilter=/FlateEncode "	

export dvips_opts="-Pamz -Pcmz"

export pdfv="gv"
export psv="gv"
# }}}

# functions {{{
ps_merge(){ # {{{

gs -dNOPAUSE -sDEVICE=pswrite -sOutputFile=$1 $2 $3 -c quit

} # }}}

pdf_merge(){ # {{{

pdftk $2 $3 cat output $1.new
mv $1.new $1

} # }}}

pdf_add(){ # {{{

if [ -f $1 ]; then
pdf_merge $1 $1 $2
else
cp $2 $1
fi

} # }}}

ps_add(){ # {{{

  if [ -f $1 ]; then
	ps_merge $1 $1 $2
      else
	cp $2 $1
fi 

} # }}}

dvitopdf(){ # {{{

#dvips $dvips_opts -o $1.ps $1.dvi
#ps2pdf $ps2pdf_opts $1.ps $1.pdf
dvips  -o $1.ps $1.dvi
ps2pdf $1.ps $1.pdf

} # }}}
# }}}
# }}}
# date/time functions {{{

date_dmy_hms(){ # {{{

date +"%D   %H:%M:%S" 

} # }}}
date_dm_hm(){ # {{{

date +"%d%m-%H%M" 

} # }}}
date_hms(){ # {{{

date +\%H:\%M:\%S 

} # }}}
# time_dhms: decompose time in days, hours, minutes, and seconds, and display this
time_dhms(){ # {{{

te=$1
if [ ! -z "$2" ]; then 
	tb=$2
else
	tb=0
fi

let t=$(($te-$tb))
# number of days
let days=$(($t/86400))

# number of seconds in the last (incomplete) day

let t=$(($t%86400))
let hrs=$(($t/3600))

let t=$(($t%3600))
let mins=$(($t/60))

let t=$(($t%60))
let secs=$t
let gts=0

days_s=""
if [ $days -gt 0 ]; then
 	days_s="$days (days)"
	gts=1
fi	

hrs_s=""
if [ $hrs -gt 0 ]; then
 	hrs_s="$hrs (hrs)"
	gts=1
fi	

mins_s=""
if [ $mins -gt 0 ]; then
 	mins_s="$mins (mins)"
	gts=1
fi	

secs_s=""
if [ $secs -gt 0 ]; then
 	secs_s="$secs (secs)"
	gts=1
fi	

res="$days_s $hrs_s $mins_s $secs_s"

if [ $gts -eq 0 ]; then
  res=" < 1 sec"
fi

echo "$res"

} # }}}

# return date in seconds

date_in_secs(){ # {{{

secs=` date +"%S" `
mins=` date +"%M" `
hrs=` date +"%H" `
days=` date +"%d" `

echo " $( echo " $secs+60*$mins+3600*$hrs+86400*$days " | bc ) "

} # }}}
# }}}
# output {{{

printvar(){ # {{{

while [ ! -z "$1" ]; do 
	case "$1" in
	  	[a-z]*) 
		case "$1" in
			"force") format_string="%1.1e\n"  ;; 	
			"force_long") format_string="%1.4e\n"  ;; 	
			"nsteps") format_string="%1.0e\n"  ;; 	
			"float") format_string="%f\n"  ;; 	
		esac
		printf "$format_string" $2
		;;
	esac
        shift
done

} # }}}
eo(){ # {{{

  case "$1" in
   	"-n") echo "" >& $log_file ;;
        *) echo "$pref> $1" >> $log_file    ;;
  esac
} # }}}

# eoo {{{

eoos(){

pref_eoo="$this_script.code"
eoo "$*"
pref_eoo="$this_script"

}

eoof(){
eoo $* >& /dev/null
}

eooe(){
pref_eoo="$this_script.err"
eoo "$*"
pref_eoo="$this_script"
}

eoo(){
echo "$pref_eoo> $*"
echo "$pref_eoo> $*" >> $log_file
}

# }}}

esl(){ 
echo "**************************************************"
}

eo_s(){
echo "@@ $1" >> $outf
}
# }}}
# copying/synchronizing {{{
cpl(){

while [ ! -z "$1" ]; do
  	dr=` date_dm_hm `
	#scp -r $1 op226@leonov:~/gmGof/$dr/
	#scp -r $1 op226@leonov:~/gmGof/
	shift
done

} # }}}

# general {{{

dir_size(){

du -hc $1 | awk '/total/{ print $1 }'  

}

use(){

source $shd/"$this_script"_$1.sh

}

#}}}

exists(){
# Check whether a program exists 
type -P $1 &>/dev/null 
echo $?

}

tti(){

ec=1
[ $1 ] && ec=0
return $ec

}

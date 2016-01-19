proc plotGMINms {sf1} {
    global scalefac
    set scalefac $sf1  
    enabletrace
}

proc enabletrace {args} {
     global vmd_frame numatoms

      trace variable vmd_frame([molinfo top]) w do_drawellipse
}

proc disabletrace {args} {
     global vmd_frame numatoms

      trace vdelete vmd_frame([molinfo top]) w do_drawellipse
}

proc initGMINms {filename} {
   global adata bdata cdata a11d a12d a13d a21d a22d a23d a31d a32d a33d xdata ydata zdata
   global vmd_frame numParamsRead
   global molid numatoms
   # assume we are in the script directory and filename is in the data directory.

   disabletrace

   set wholeFileName [expr {"$filename"}]

   # open the data file and assign a fileId.
   set fileId [open $wholeFileName r]

   # extract num of atoms from file
   
   mol load xyz $filename
   set molid [molinfo top] 
   set frame [molinfo $molid get numframes]
   set cf 0 
      #   set vmd_frame($molId) 
   set numatoms [gets $fileId]
   # scan until the end of the file comes up
   while { [eof $fileId] == 0 } {
     # read next line in file
     set energy [gets $fileId]
  
     # parse line
     set numParamsRead [scan $energy " Energy of minimum %d =%f  first found at step %d" minIndex minEnergy minStep]

     # check the number of parameters read is correct before proceeding, else it is not a start of frame line.
     if {$numParamsRead == 3} {


#	mol rename $molId "minimum_$minIndex _$minEnergy"

        # output energy line by way of a heart beat [ie so user knows program is running].
        puts $energy
        puts $cf 
        puts $molid 

        # loop through $numAtoms lines obtaining information about each atom in molecule and then plot the data in VMD.
        for {set atom 0} {$atom < $numatoms} {set atom [expr $atom + 1]} {

           set atom_detail [gets $fileId]
           set numParamsRead [scan $atom_detail " %c %f %f %f  ellipse %f %f %f %f %f %f %f %f %f %f %f %f  atom_vector %f %f %f" \
                atomType xdata($cf,$atom) ydata($cf,$atom) zdata($cf,$atom) \
                adata($cf,$atom) bdata($cf,$atom) cdata($cf,$atom) \
                a11d($cf,$atom) a12d($cf,$atom) a13d($cf,$atom) \
                a21d($cf,$atom) a22d($cf,$atom) a23d($cf,$atom) \
                a31d($cf,$atom) a32d($cf,$atom) a33d($cf,$atom) n1 n2 n3]
        }

       set cf [expr $cf + 1]
     }
  } 
}





## ellipse stuff

proc f {a u v} {
  expr {$a*sin($u) * sin($v)}
}

proc g {b u v} {
  expr {$b*cos($u) * sin($v)}
}

proc h {c u v} {
  expr {$c*cos($v)}
}

# a = semi-axis in x direction
# b = semi-axis in y direction
# c = semi-axis in z direction
# gamma = first rotation about z axis in xy plane (body rotation of gamma - for ellipse of revolution has no effect) 
# (  cos(gamma) sin(gamma) 0)
# ( -sin(gamma) cos(gamma) 0)
# (     0        0         1)
#
# alpha = second rotation - normal vector starts parallel with z axis - tilt it alpha degrees in xz plane. 
# This makes the normal at alpha degrees to z-axis
# (  cos(alpha)  0  sin(alpha) )
# (     0        1      0      )
# ( -sin(alpha)  0  cos(alpha) )
#
# beta = third rotation about z axis in xy plane. - Rotates plane containing normal to an angle beta away from x-axis
# (  cos(beta) sin(beta) 0)
# ( -sin(beta) cos(beta) 0)
# (     0        0       1)
# 
#  End result of rotations: 
#  The ellipse normal is pointing at angle (alpha,beta) - corresponding to typical spherical polars), 
#  body has rotated an angle gamma about it's normal axis.
# x = translation in x direction.
# y = translation in y direction
# z = translation in z direction
# eColour is the colour of the ellipse in VMD colour table.

proc draw_ellipse {index atom eColour nSteps} {
  global adata bdata cdata hhdata iidata jjdata xdata ydata zdata
  global a11d a12d a13d a21d a22d a23d a31d a32d a33d
  
  set a [expr $adata($index,$atom)/2]
  set b [expr $bdata($index,$atom)/2]
  set c [expr $cdata($index,$atom)/2]
#  set alpha [expr $hhdata($index,$atom)]
#  set beta [expr $iidata($index,$atom)]
#  set gamma [expr $jjdata($index,$atom)]
  set x [expr $xdata($index,$atom)]
  set y [expr $ydata($index,$atom)]
  set z [expr $zdata($index,$atom)]
# rotation matrix is printed out in the file
 set a11 [expr $a11d($index,$atom)]
 set a12 [expr $a12d($index,$atom)]
 set a13 [expr $a13d($index,$atom)]
 set a21 [expr $a21d($index,$atom)]
 set a22 [expr $a22d($index,$atom)]
 set a23 [expr $a23d($index,$atom)]
 set a31 [expr $a31d($index,$atom)]
 set a32 [expr $a32d($index,$atom)]
 set a33 [expr $a33d($index,$atom)]

  

  set PI 3.14159265358979
#  set printFlag 0
#  # convert euler angles from degrees to radians
#  set alpha [expr $alpha*$PI/180]
#  set beta [expr $beta*$PI/180]
#  set gamma [expr $gamma*$PI/180]
#  puts "in here $index $atom $hhdata($index,$atom), $idata($index,$atom) $jdata($index,$atom) $cdata($index,$atom)"
# create individual parameters for matrix to reduce lookup operations

#  matrix[0] = cg*(cb*ca) - sg*sa;
#  matrix[1] = -(sg*(cb*ca)) - cg*sa;
#  matrix[2] = sb*ca;
#
#  matrix[3] = cg*(cb*sa) + sg*ca;
#  matrix[4] = cg*ca - sg*(cb*sa);
#  matrix[5] = sb*sa;
#
#  matrix[6] = -(cg*sb);
#  matrix[7] = sg*sb;
#  matrix[8] = cb;

#  set a11 [expr cos($gamma)*cos($beta)*cos($alpha) - sin($gamma)*sin($alpha)]
#  set a12 [expr -1*sin($gamma)*cos($beta)*cos($alpha) - cos($gamma)*sin($alpha)]
#  set a13 [expr    sin($beta)*cos($alpha)]
#  set a21 [expr cos($gamma)*cos($beta)*sin($alpha) + sin($gamma)*cos($alpha)]
#  set a22 [expr -1*sin($gamma)*cos($beta)*sin($alpha) + cos($gamma)*cos($alpha)]
#  set a23 [expr sin($alpha)*sin($beta)]
#  set a31 [expr -1*cos($gamma)*sin($beta)]
#  set a32 [expr sin($gamma)*sin($beta)]
#  set a33 [expr cos($beta)]
  
#  puts "calculated matrix elements: $a11 $a12 $a13"

  #define limits: u is analagous to azimuthal phi in conventional spherical polars 
  set minu 0
  set maxu [expr 2*$PI]
  set stepu [expr $maxu / $nSteps]

  #extend limits by one extra step to cope with edge effects.
  set maxu [expr $maxu+$stepu]
  set minu [expr $minu-$stepu]

  #define limits: V is analagous to theta in conventional spherical polars 
  set minv [expr 0]
  set maxv [expr $PI]
  set stepv [expr $maxv / $nSteps]

  #extend limits by one extra step to cope with edge effects.
  set maxv [expr $maxv+$stepv]
  set minv [expr $minv-$stepv]

  # loop through the u and v indices from -1 to nSteps+1 defining data, gradient and rotations and translations thereof.
  set u [expr $minu]
  for {set uIndex -1} {$uIndex <= $nSteps+1} {set uIndex [expr $uIndex + 1]} {
    set v [expr  $minv]
    for {set vIndex -1} {$vIndex <= $nSteps+1} {set vIndex [expr $vIndex + 1]} {
	#compute f,g and h (conventionally x,y and z) positions of unrotated ellipsoid at origin.
	set fdata($uIndex,$vIndex) [f $a $u $v]
        set gdata($uIndex,$vIndex) [g $b $u $v]
        set hdata($uIndex,$vIndex) [h $c $u $v]
        
	#compute normal vector at each position of unrotated ellipsoid at origin.
	  set fdata_g($uIndex,$vIndex) [expr {$fdata($uIndex,$vIndex)/($a*$a)}]
        set gdata_g($uIndex,$vIndex) [expr {$gdata($uIndex,$vIndex)/($b*$b)}]
        set hdata_g($uIndex,$vIndex) [expr {$hdata($uIndex,$vIndex)/($c*$c)}]
        set norm_constant [expr sqrt($fdata_g($uIndex,$vIndex)*$fdata_g($uIndex,$vIndex)+$gdata_g($uIndex,$vIndex)*$gdata_g($uIndex,$vIndex)+$hdata_g($uIndex,$vIndex)*$hdata_g($uIndex,$vIndex))]
	  set fdata_g($uIndex,$vIndex) [expr {$fdata_g($uIndex,$vIndex)/$norm_constant }]
        set gdata_g($uIndex,$vIndex) [expr {$gdata_g($uIndex,$vIndex)/$norm_constant }]
        set hdata_g($uIndex,$vIndex) [expr {$hdata_g($uIndex,$vIndex)/$norm_constant }]


        #perform rotation and translation of ellipsoid data
        set fdata_r($uIndex,$vIndex) [expr {($fdata($uIndex,$vIndex)*$a11+$gdata($uIndex,$vIndex)*$a12+$hdata($uIndex,$vIndex)*$a13) + $x}]
        set gdata_r($uIndex,$vIndex) [expr {($fdata($uIndex,$vIndex)*$a21+$gdata($uIndex,$vIndex)*$a22+$hdata($uIndex,$vIndex)*$a23) + $y}]
        set hdata_r($uIndex,$vIndex) [expr {($fdata($uIndex,$vIndex)*$a31+$gdata($uIndex,$vIndex)*$a32+$hdata($uIndex,$vIndex)*$a33) + $z}]
      
        #perform rotation and translation of ellipsoid gradient data
        set fdata_r_g($uIndex,$vIndex) [expr {($fdata_g($uIndex,$vIndex)*$a11+$gdata_g($uIndex,$vIndex)*$a12+$hdata_g($uIndex,$vIndex)*$a13) + $x}]
        set gdata_r_g($uIndex,$vIndex) [expr {($fdata_g($uIndex,$vIndex)*$a21+$gdata_g($uIndex,$vIndex)*$a22+$hdata_g($uIndex,$vIndex)*$a23) + $y}]
        set hdata_r_g($uIndex,$vIndex) [expr {($fdata_g($uIndex,$vIndex)*$a31+$gdata_g($uIndex,$vIndex)*$a32+$hdata_g($uIndex,$vIndex)*$a33) + $z}]

        #increment v
        set v [expr $v + $stepv]
    }
    #increment u
    set u [expr $u + $stepu]
  }
  # make another pass through to plot it; this time going from 0 to nsteps-1
  for {set uIndex 0} {$uIndex < $nSteps} {set uIndex [expr $uIndex + 1]} {
    for {set vIndex 0} {$vIndex < $nSteps} {set vIndex [expr $vIndex + 1]} {

      # get the next two corners
      set u2Index [expr $uIndex + 1]
      set v2Index [expr $vIndex + 1]
      draw color $eColour
      draw trinorm "$fdata_r($uIndex,$vIndex)  $gdata_r($uIndex,$vIndex)  $hdata_r($uIndex,$vIndex)" \
                    "$fdata_r($u2Index,$vIndex)  $gdata_r($u2Index,$vIndex)  $hdata_r($u2Index,$vIndex)" \
                    "$fdata_r($u2Index,$v2Index) $gdata_r($u2Index,$v2Index) $hdata_r($u2Index,$v2Index)" \
                    "$fdata_r_g($uIndex,$vIndex)  $gdata_r_g($uIndex,$vIndex)  $hdata_r_g($uIndex,$vIndex)" \
                    "$fdata_r_g($u2Index,$vIndex)  $gdata_r_g($u2Index,$vIndex)  $hdata_r_g($u2Index,$vIndex)" \
                    "$fdata_r_g($u2Index,$v2Index) $gdata_r_g($u2Index,$v2Index) $hdata_r_g($u2Index,$v2Index)"
                  
	      draw trinorm "$fdata_r($u2Index,$v2Index) $gdata_r($u2Index,$v2Index) $hdata_r($u2Index,$v2Index)" \
                    "$fdata_r($uIndex,$v2Index) $gdata_r($uIndex,$v2Index) $hdata_r($uIndex,$v2Index)" \
                    "$fdata_r($uIndex,$vIndex)  $gdata_r($uIndex,$vIndex)  $hdata_r($uIndex,$vIndex)" \
                    "$fdata_r_g($u2Index,$v2Index) $gdata_r_g($u2Index,$v2Index) $hdata_r_g($u2Index,$v2Index)" \
                    "$fdata_r_g($uIndex,$v2Index) $gdata_r_g($uIndex,$v2Index) $hdata_r_g($uIndex,$v2Index)" \
                    "$fdata_r_g($uIndex,$vIndex)  $gdata_r_g($uIndex,$vIndex)  $hdata_r_g($uIndex,$vIndex)"
    }
  }
}

proc draw_ljsite {index atom dist radius eColour} {
  global adata bdata cdata hhdata iidata jjdata xdata ydata zdata

  set alpha [expr $hhdata($index,$atom)]
  set beta [expr $iidata($index,$atom)]
  set gamma [expr $jjdata($index,$atom)]
  set x [expr $xdata($index,$atom)]
  set y [expr $ydata($index,$atom)]
  set z [expr $zdata($index,$atom)]

  set dist [expr $dist*$adata($index,$atom)/2]

  set PI 3.14159265358979
  set printFlag 0
  # convert euler angles from degrees to radians
  set alpha [expr $alpha*$PI/180]
  set beta [expr $beta*$PI/180]
  set gamma [expr $gamma*$PI/180]

  set a11 [expr cos($gamma)*cos($beta)*cos($alpha) - sin($gamma)*sin($alpha)]
  set a12 [expr -1*sin($gamma)*cos($beta)*cos($alpha) - cos($gamma)*sin($alpha)]
  set a13 [expr    sin($beta)*cos($alpha)]
  set a21 [expr cos($gamma)*cos($beta)*sin($alpha) + sin($gamma)*cos($alpha)]
  set a22 [expr -1*sin($gamma)*cos($beta)*sin($alpha) + cos($gamma)*cos($alpha)]
  set a23 [expr sin($alpha)*sin($beta)]
  set a31 [expr -1*cos($gamma)*sin($beta)]
  set a32 [expr sin($gamma)*sin($beta)]
  set a33 [expr cos($beta)]

        set xlj [expr {$dist*$a11 + $x}]
        set ylj [expr {$dist*$a21 + $y}]
        set zlj [expr {$dist*$a31 + $z}]
        draw color green

        draw sphere "$xlj $ylj $zlj" radius $radius
        draw cylinder "$x $y $z" "$xlj $ylj $zlj" radius 0.05
}


proc do_drawellipse {name element op} {
     global vmd_frame numatoms scalefac

      draw delete all
     for {set atom 0} {$atom < $numatoms} {set atom [expr $atom + 1]} {

        draw_ellipse [expr $vmd_frame([molinfo top])] $atom 110 10
#        draw_ljsite [expr $vmd_frame([molinfo top])] $atom [expr -1.0*$scalefac] 0.1 100

     }

#draw_ellipse [expr $a/2] [expr $b/2] [expr $c/2] [expr $h] [expr $i] $j $x $y $z 100 20
#draw_ljsite 1.0 0.2 [expr $h] [expr $i] $j $x $y $z 100

}

# for plotting chiropoles
# takes an .xyz file with molecule centers and mu_hats
# produces red/blue cylinders

# usage:
# open vmd in command line console (or Tk console), enter:
#  source ~/wherever/plotchiro.tcl
# then to load a molecule,
#  chiro chiro.xyz
# After the filename you can put an radius and length for the cylinders. Default is 0.05 0.5.

proc chiro {filename {rad 0.05} {len 0.5}} {
    global cyl_rad cyl_len

    set cyl_rad $rad
    set cyl_len $len
    
    initchiro $filename
    enabletrace

    animate goto start
}

proc enabletrace {args} {
    global vmd_frame natoms
    trace variable vmd_frame([molinfo top]) w drawcyls
}

proc disabletrace {args} {
    global vmd_frame natoms
    trace vdelete vmd_frame([molinfo top]) w drawcyls
}

proc initchiro {filename} {
    global vmd_frame molid natoms nread
    global xdata ydata zdata adata bdata cdata

    disabletrace

    # open file for reading
    set whole_filename [expr {"$filename"}]
    set file_id [open $whole_filename r]

    # read number of atoms
    mol load xyz $filename
    set molid [molinfo top]
    set frame [molinfo $molid get numframes]
    set cf 0
    set natoms [gets $file_id]

    # scan until EOF
    while { [eof $file_id] == 0 } {
        # read next line
        set energy [gets $file_id]
        set nread [scan $energy " energy of minimum %d =%f   first found at step %d" minind mine minstep]

        # double check that we read the three parameters. otherwise this is not the start of a frame
        if {$nread == 3} {
            # loop through the atoms
            for {set atom 0} {$atom < $natoms} {set atom [expr $atom + 1]} {
                set atom_detail [gets $file_id]
                set nread [scan $atom_detail " %c %f %f %f atom_vector %f %f %f" \
                    atom_type xdata($cf,$atom) ydata($cf,$atom) zdata($cf,$atom) \
                    adata($cf,$atom) bdata($cf,$atom) cdata($cf,$atom)]
            }

           set cf [expr $cf + 1]
        }
    }
}

proc drawcyls {args} {
    global cyl_rad cyl_len
    global vmd_frame natoms
    global xdata ydata zdata adata bdata cdata

    draw delete all

    mol delrep 0 top

    set ind [expr $vmd_frame([molinfo top])]

    for {set atom 0} {$atom < $natoms} {set atom [expr $atom + 1]} {
        set x [expr {$xdata($ind,$atom)}]
        set y [expr {$ydata($ind,$atom)}]
        set z [expr {$zdata($ind,$atom)}]
        set a [expr {$adata($ind,$atom)}]
        set b [expr {$bdata($ind,$atom)}]
        set c [expr {$cdata($ind,$atom)}]

        set xend [expr {$x+$cyl_len*$a}]
        set yend [expr {$y+$cyl_len*$b}]
        set zend [expr {$z+$cyl_len*$c}]
        draw color red
        draw cylinder "$x $y $z" "$xend $yend $zend" radius $cyl_rad

        set xend [expr {$x-$cyl_len*$a}]
        set yend [expr {$y-$cyl_len*$b}]
        set zend [expr {$z-$cyl_len*$c}]
        draw color blue
        draw cylinder "$x $y $z" "$xend $yend $zend" radius $cyl_rad
    }
}

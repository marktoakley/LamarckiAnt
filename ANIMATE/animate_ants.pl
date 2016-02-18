#! /usr/bin/perl -w
# Animate LamarckiAnt paths
#
# This script generates an animated gif showing the structures generated during
# a LamarckiAnt run. With the default setting, every strucure in a generation
# is represented by a thin grey line. The best structure from each generation
# and the best overall structure are represented by thicker red and black
# lines.

# Count number of generations
$max=(split/\s+/,`ls ANT.* | wc`)[1]-1;
# Set format of plot
open(OUTPUT, "> plot.gpi");
print OUTPUT "clear\n";
print OUTPUT "reset\n";
print OUTPUT "set terminal gif animate delay 30\n";
print OUTPUT "set output \"animate.gif\"\n";
print OUTPUT "set xlabel 'Residue'\n";
print OUTPUT "set ylabel 'Torsion'\n";
print OUTPUT "set isosample 40\n";
print OUTPUT "set hidden3d\n";
print OUTPUT "set style data lines\n";
print OUTPUT "set yrange [0:360]\n";
print OUTPUT "unset key\n";
# Print one frame per generation
for ($i=1;$i<=$max;$i++){
	print OUTPUT "set title 'Generation $i'\n";
       print OUTPUT "plot 'ANT.$i' u 1:(\$2*180/3.14159) lw 1 linecolor rgb 'grey','ANT.$i' u 1:(\$3*180/3.14159) w linespoints lw 2 linecolor rgb 'red','ANT.$i' u 1:(\$4*180/3.14159) w linespoints lw 2 linecolor rgb 'black'\n" ;
}
# Pause on final generation for a few frames
for ($i=1;$i<=5;$i++){
	print OUTPUT "set title 'Generation $max'\n";
       print OUTPUT "plot 'ANT.$max' u 1:(\$2*180/3.14159) lw 1 linecolor rgb 'grey','ANT.$max' u 1:(\$3*180/3.14159) w linespoints lw 2 linecolor rgb 'red','ANT.$max' u 1:(\$4*180/3.14159) w linespoints lw 2 linecolor rgb 'black'\n" ;
}
close (OUTPUT);

system ("gnuplot plot.gpi");

system ("firefox animate.gif");

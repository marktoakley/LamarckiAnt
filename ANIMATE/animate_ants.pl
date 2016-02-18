#! /usr/bin/perl -w

open(OUTPUT, "> plot.gpi");
print OUTPUT "clear\n";
print OUTPUT "reset\n";
print OUTPUT "set terminal gif animate delay 30\n";
print OUTPUT "set output \"animate.gif\"\n";

print OUTPUT "set isosample 40\n";
print OUTPUT "set hidden3d\n";
print OUTPUT "set style data lines\n";
print OUTPUT "set yrange [0:360]\n";
print OUTPUT "unset key\n";
for ($i=1;$i<=100;$i++){
	print OUTPUT "set title 'Generation $i'\n";
       print OUTPUT "plot 'ANT.$i' u 1:(\$2*180/3.14159) lw 1 linecolor rgb 'grey','ANT.$i' u 1:(\$3*180/3.14159) w linespoints lw 2 linecolor rgb 'red','ANT.$i' u 1:(\$4*180/3.14159) w linespoints lw 2 linecolor rgb 'black'\n" ;
}
close (OUTPUT);
system ("gnuplot plot.gpi");
system ("firefox animate.gif");

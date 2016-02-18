#! /usr/bin/perl -w

$max=128;

open(OUTPUT, "> plot.gpi");
print OUTPUT "clear\n";
print OUTPUT "reset\n";
print OUTPUT "set terminal gif animate delay 20\n";
print OUTPUT "set output \"animate.gif\"\n";
print OUTPUT "set isosample 40\n";
print OUTPUT "unset key\n";
#print OUTPUT "set pm3d map\n";
print OUTPUT "set yrange [0:360]\n";
print OUTPUT "set palette grey\n";
print OUTPUT "set contour base\n";
print OUTPUT "set nosurface\n";
print OUTPUT "set view 0,0\n";
print OUTPUT "set cntrparam levels incremental 0.001 ,0.01, 0.201\n";
for ($i=1;$i<=$max;$i++){
	print OUTPUT "set title '$i'\n";
	print OUTPUT "splot \"PHEROMONE\.$i\" u 2:(\$1*10):3 w dots\n";
}
print OUTPUT "splot \"PHEROMONE\.$max\" u 2:(\$1*10):3\n";
print OUTPUT "splot \"PHEROMONE\.$max\" u 2:(\$1*10):3\n";
print OUTPUT "splot \"PHEROMONE\.$max\" u 2:(\$1*10):3\n";
print OUTPUT "splot \"PHEROMONE\.$max\" u 2:(\$1*10):3\n";
print OUTPUT "splot \"PHEROMONE\.$max\" u 2:(\$1*10):3\n";
system ("gnuplot plot.gpi");
system ("firefox animate.gif");

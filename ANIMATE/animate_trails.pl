#! /usr/bin/perl -w
# Animate LamarckiAnt pheromone trails
#
# This script generates an animated gif showing the evolotion of the pheromone
# trail during a LamarckiAnt run.

# Count number of generations
$max=(split/\s+/,`ls ANT.* | wc`)[1]-1;

open(OUTPUT, "> plot.gpi");
print OUTPUT "clear\n";
print OUTPUT "reset\n";
print OUTPUT "set terminal gif animate delay 20\n";
print OUTPUT "set output \"animate.gif\"\n";
print OUTPUT "set xlabel 'Residue'\n";
print OUTPUT "set ylabel 'Torsion'\n";
print OUTPUT "set isosample 40\n";
print OUTPUT "unset key\n";
print OUTPUT "set yrange [0:360]\n";
print OUTPUT "unset surface\n";
print OUTPUT "set view map\n";

#print OUTPUT "set cntrparam levels incremental 0.001 ,0.01, 0.201\n";
#print OUTPUT "set palette grey\n";

#print OUTPUT "set contour base\n";
print OUTPUT "set style data pm3d\n";
print OUTPUT "set style function pm3d\n";
print OUTPUT "set pm3d implicit at b\n";

# Print one frame per generation
for ($i=1;$i<=$max;$i++){
	print OUTPUT "set title 'Generation$i'\n";
	print OUTPUT "splot \"PHEROMONE\.$i\" u 2:(\$1*10):3\n";
}
#Pause for a few frames at final generation
for ($i=1;$i<=5;$i++){
	print OUTPUT "splot \"PHEROMONE\.$max\" u 2:(\$1*10):3\n";
}

system ("gnuplot plot.gpi");
system ("firefox animate.gif");

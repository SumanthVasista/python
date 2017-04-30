#! /bin/bash

exists1=$(ls | grep -c -w "cluster.dat")
if [ $exists1 -gt 0 ]; then
	printf "#!/usr/bin/gnuplot\nunset key\nset datafile separator \" \"\nset object rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb 'white' fillstyle solid noborder\nset grid\nset style fill transparent solid 0.5\nset style circle radius 0.01\nset palette model RGB defined ( 0 'orange', 1 'blue' )\nif (!exists(\"filename\")) filename='default.dat'\nplot[0:5][0:6] \"cluster.dat\" u 2:3:( \$4 == 0 ? 0 : 1 ) with circles palette, $1 * x + $4, $2 * x + $4, $3 * x + $4\npause -1\n" > scatterplot.plg && gnuplot -e 'filename="cluster.dat"' scatterplot.plg
else
	echo "File not found"
fi

reset
set size 1,1

set multiplot

#first
set size 1,0.5
set origin 0,0.5
set arrow from graph 0,first 7.815 to graph 1,first 7.815 nohead lc rgb "red" lw 3 front

plot 'NIS_lidar' with lines

#second
set origin 0,0
set arrow from graph 0,first 7.815 to graph 1,first 7.815 nohead lc rgb "red" lw 3 front
plot 'NIS_radar' with lines

unset multiplot
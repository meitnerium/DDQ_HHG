set multiplot
set size 1,0.5
set origin 0.0,0.5
plot exp((-(x+0)**2)/1.5)*cos(x*3)
set origin 0.0,0.0
set xrange [-40:40]
plot exp((-(x+20)**2)/45)+exp((-(x-20)**2)/45)
set nomultiplot
pause -1

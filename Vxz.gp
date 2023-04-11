reset
set terminal qt font ',20'

set view 60,115

#set size square
set style data lines
#set title "{/Symbol D}t = 0.3 s,  {/Symbol D}x = 0.01 m"
set border lc rgb 'black' lw 1
set cbrange [-1:1]
set cblabel "malla"
set style fill   solid 2.00 border lt -1

set xlabel 'malla y' 
set ylabel 'malla z'
set xtics font ",12" rotate by 45
set ytics  font ",12"
set ztics   font ",10"
set mxtics 4 
set mytics 4
#set xrange [0:1]
#set contour base
#set clabel 
#set colorbox vert user origin .9,.05 size .04,.5
#set cntrparam levels incremental -1, 0.1, 1

set pm3d at s
splot "Vxz.dat" u 1:2:3  w l ls 7 palette notitle

set terminal pdf color
set output 'grafico_Vxz.pdf'
replot 
exit

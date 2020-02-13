reset
palette=system("echo ${HOME}/KNOSOS/PLT/palette.plt")
l palette
dz=`tail -1 B.map |cut -f5 -d" "`
N=int(2*pi/dz)
f=100
MOD(x,y)=x-int(x/y)*y
p "ph1_2d.re" u ((MOD($2,2*pi/N))/pi*N/2):(MOD($3,2*pi)/2./pi):($4*f) palette pt 5 ps 2 notitle
set terminal epslatex color standalone
set output 'varphi1oT.tex'
set size square
set xrange [0:1]
set yrange [0:1]
set cbrange [GPVAL_DATA_CB_MIN:GPVAL_DATA_CB_MAX]
set xtics 0,1
set ytics 0,1
set cbtics GPVAL_DATA_CB_MIN,GPVAL_DATA_CB_MAX-GPVAL_DATA_CB_MIN
set xlabel offset 0,1.2
set ylabel offset 3,0
set cblabel offset -3,0
set xlabel '$\zeta/(2\pi/N)$'
set ylabel '$\theta/(2\pi)$'
set cblabel '$e\varphi_1/T_i [\times 10^{-2}]$'
set format cb '$%1.1f$'
dz=`tail -1 B.map |cut -f5 -d" "`
N=int(2*pi/dz)
f=100
MOD(x,y)=x-int(x/y)*y
p "ph1_2d.re" u ((MOD($2,2*pi/N))/pi*N/2):(MOD($3,2*pi)/2./pi):($4*f) palette pt 5 ps 2 notitle
set output
!latex varphi1oT
#!dvipdf varphi1oT.dvi
#!cp -p varphi1oT.pdf /home/u6156/

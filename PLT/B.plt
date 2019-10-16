reset
N=`grep nzperiod ddkes2.data|cut -f2 -d=|cut -f1 -d,`
B00=`grep borbi\(0\,0\) ddkes2.data|cut -f2 -d=|cut -f1 -d,`
p "B.map" u ((MOD($2,2*pi/N))/pi*N/2):(MOD($3,2*pi)/2./pi):($4/B00) palette pt 5 ps 2 notitle
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
set cblabel '$B/B_{00}$'
set format cb '$%1.2f$'
MOD(x,y)=x-int(x/y)*y
N=`cat ddkes2.data |grep nzperiod|cut -f2 -d=|cut -f1 -d,`
B00=`grep borbi\(0\,0\) ddkes2.data|cut -f2 -d=|cut -f1 -d,`
p "B.map" u ((MOD($2,2*pi/N))/pi*N/2):(MOD($3,2*pi)/2./pi):($4/B00) palette pt 5 ps 2 notitle

reset
N=`grep nzperiod ddkes2.data|cut -f2 -d=|cut -f1 -d,`
B00=`grep borbi\(0\,0\) ddkes2.data|cut -f2 -d=|cut -f1 -d,`
p "B.map" u ((MOD($2,2*pi/N))/pi*N/2):(MOD($3,2*pi)/2./pi):($4/B00) palette pt 5 ps 2 notitle
set size square
set xrange [-0.2:1.2]
set yrange [-0.7:0.7]
set cbrange [GPVAL_DATA_CB_MIN:GPVAL_DATA_CB_MAX]
set xtics 0,1
set ytics -1,1
set cbtics GPVAL_DATA_CB_MIN,GPVAL_DATA_CB_MAX-GPVAL_DATA_CB_MIN
set xlabel '$\zeta/(2\pi/N)$'
set ylabel '$\theta/(2\pi)$'
set cblabel '$B/B_{00}$'
set cblabel offset -3,0
set format cb '$%1.2f$'
MOD(x,y)=x-int(x/y)*y
N=`grep nzperiod ddkes2.data|cut -f2 -d=|cut -f1 -d,`
B00=`grep borbi\(0\,0\) ddkes2.data|cut -f2 -d=|cut -f1 -d,`
siota=-1
p   "B.map" u ($2/2/pi*N-1):($3/2/pi+0*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+0):($3/2/pi+0*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+1):($3/2/pi+0*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+2):($3/2/pi+0*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+3):($3/2/pi+0*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+4):($3/2/pi+0*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+5):($3/2/pi+0*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+6):($3/2/pi+0*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N-1):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+0):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+1):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+2):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+3):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+4):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+5):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+5):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "B.map" u ($2/2/pi*N+6):($3/2/pi+1*siota):($4/B00) pt 5 palette notitle\
, "LOG/fort.3000" u ($3/2/pi*N):($4/2/pi) pt 5 ps 0.5 lt 2 notitle\
, "LOG/fort.3000" u ($1/2/pi*N):($2/2/pi) pt 7 ps 0.5 lt 3 notitle\
, "LOG/fort.2000" u ($2/2/pi*N+0):($3/2/pi+0*siota) pt 7 ps 1 lt 4 notitle
#, "LOG/fort.3000" u ($3/2/pi*N):($4/2/pi) pt 5 ps 1 lt 2 notitle\
#, "LOG/fort.3000" u ($1/2/pi*N):($2/2/pi) pt 7 ps 1 lt 3 notitle\
#, "LOG/fort.2000" u ($2/2/pi*N+0):($3/2/pi+0*siota) pt 7 ps 2 lt 4 notitle

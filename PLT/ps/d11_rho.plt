reset
set terminal epslatex color standalone
set output 'd11_rho.tex'
set log y
set key center center
set xlabel offset 0,0.0
set ylabel offset 1,0
set xlabel '$\psi/\psi_{LCFS}$'
set ylabel '${\hat D_{11}^*}$'
set xrange [0:1]
#set yrange [1E-3:3E+3]
set format y '$10^{%T}$'
set xtics 0.2,0.2
chip=`grep chip s_03??/ddkes2.data|cut -f2 -d=|cut -f1 -d,`
psip=`grep chip s_03??/ddkes2.data|cut -f3 -d=|cut -f1 -d,`
iota=abs(chip/psip)
B00=`grep borbi\(0\,0\) s_03??/ddkes2.data|cut -f2 -d=|cut -f1 -d,`
bzeta=`grep bzeta s_03??/ddkes2.data|cut -f3 -d=|cut -f1 -d,`
R=abs(bzeta)/B00
f=2*iota/(pi*R)
l "./local.plt"
p "<awk '{if(($1>0.05)&&($1<0.9)&&($2==3e-6)&&($3==0e-0)) print}' results.data" u 1:(f*0.5*($6+$7)):(f*0.5*($6-$7)) w e lc 3 pt 5 ps 2 title '$v_{E*}=0        $, {\ttfamily DKES}'\
, "<awk '{if(($1>0.05)&&($1<0.9)&&($2==3e-6)&&($3==3e-4)) print}' results.data" u 1:(f*0.5*($6+$7)):(f*0.5*($6-$7)) w e lc 5 pt 5 ps 2 title '$v_{E*}B_{0,0}=3\times 10^{-4}\,$T, {\ttfamily DKES}'\
, "<awk '{if(($1>0.05)&&($1<0.9)&&($2>2e-6)&&($2<4e-6)&&($3>-2e-4)&&($3<4e-9)) print}' results.knosos" u 1:(f*$7) w lp lt 1 lw 2 lc 3 pt 6 ps 1 title '$v_{E*}=0        $, {\ttfamily KNOSOS}'\
, "<awk '{if(($1>0.05)&&($1<0.9)&&($2>2e-6)&&($2<4e-6)&&($3 >2e-4)&&($3<4e-4)) print}' results.knosos" u 1:(f*$7) w lp lt 1 lw 2 lc 5 pt 6 ps 1 title '$v_{E*}B_{0,0}=3\times 10^{-4}\,$T, {\ttfamily KNOSOS}' 
set output
!latex d11_rho
!dvipdf d11_rho.dvi
!cp -p d11_rho.pdf /home/u6156/

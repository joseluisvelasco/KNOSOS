set terminal epslatex color standalone
set output 'd11.tex'
set log
set key top right
set xlabel offset 0,0.0
set ylabel offset 1,0
set xlabel '$\nu_*$'
set ylabel '${\hat D_{11}^*}$'
set format x '$10^{%T}$'
set format y '$10^{%T}$'
set ylabel offset -0.5,0
chip=`grep chip ddkes2.data|cut -f2 -d=|cut -f1 -d,`
psip=`grep chip ddkes2.data|cut -f3 -d=|cut -f1 -d,`
iota=abs(chip/psip)
B00=`grep borbi\(0\,0\) ddkes2.data|cut -f2 -d=|cut -f1 -d,`
bzeta=`grep bzeta ddkes2.data|cut -f3 -d=|cut -f1 -d,`
R=abs(bzeta)/B00
f=2*iota/(pi*R)
l "./local.plt"
p "<awk '{if(($2==0.00000E-00)&&($1<1e6)) print}' results.data" u ($1*R/iota):(f*($5+$6)/2):(f*($5-$6)/2) w e notitle lc 3 pt 5 ps 2\
, "<awk '{if(($2==1.00000E-05)&&($1<1e6)) print}' results.data" u ($1*R/iota):(f*($5+$6)/2):(f*($5-$6)/2) w e notitle lc 4 pt 5 ps 2\
, "<awk '{if(($2==3.00000E-05)&&($1<1e6)) print}' results.data" u ($1*R/iota):(f*($5+$6)/2):(f*($5-$6)/2) w e notitle lc 2 pt 5 ps 2\
, "<awk '{if(($2==1.00000E-04)&&($1<1e6)) print}' results.data" u ($1*R/iota):(f*($5+$6)/2):(f*($5-$6)/2) w e notitle lc 7 pt 5 ps 2\
, "<awk '{if(($2==3.00000E-04)&&($1<1e6)) print}' results.data" u ($1*R/iota):(f*($5+$6)/2):(f*($5-$6)/2) w e notitle lc 5 pt 5 ps 2\
, "<awk '{if(($2==1.00000E-03)&&($1<1e6)) print}' results.data" u ($1*R/iota):(f*($5+$6)/2):(f*($5-$6)/2) w e notitle lc 1 pt 5 ps 2\
, "<awk '{if(($2==3.00000E-03)&&($1<1e6)) print}' results.data" u ($1*R/iota):(f*($5+$6)/2):(f*($5-$6)/2) w e notitle lc 9 pt 5 ps 2\
, "<awk '{if(($2==0.00000E-00)&&($1<1e6)) print}' results.knosos" u ($1*R/iota):(f*$5) w lp title '$v_{E*}=0$' lt 1 pt 6 lw 2 lc 3\
, "<awk '{if(($2==1.00000E-05)&&($1<1e6)) print}' results.knosos" u ($1*R/iota):(f*$5) w lp title '$v_{E*}B_{0,0}=1\times 10^{-5}\,$T' lt 1 pt 6 lw 2 lc 4\
, "<awk '{if(($2==3.00000E-05)&&($1<1e6)) print}' results.knosos" u ($1*R/iota):(f*$5) w lp title '$v_{E*}B_{0,0}=3\times 10^{-5}\,$T' lt 1 pt 6 lw 2 lc 2\
, "<awk '{if(($2==1.00000E-04)&&($1<1e6)) print}' results.knosos" u ($1*R/iota):(f*$5) w lp title '$v_{E*}B_{0,0}=1\times 10^{-4}\,$T' lt 1 pt 6 lw 2 lc 7\
, "<awk '{if(($2==3.00000E-04)&&($1<1e6)) print}' results.knosos" u ($1*R/iota):(f*$5) w lp title '$v_{E*}B_{0,0}=3\times 10^{-4}\,$T' lt 1 pt 6 lw 2 lc 5\
, "<awk '{if(($2==1.00000E-03)&&($1<1e6)) print}' results.knosos" u ($1*R/iota):(f*$5) w lp title '$v_{E*}B_{0,0}=1\times 10^{-3}\,$T' lt 1 pt 6 lw 2 lc 1\
, "<awk '{if(($2==3.00000E-03)&&($1<1e6)) print}' results.knosos" u ($1*R/iota):(f*$5) w lp title '$v_{E*}B_{0,0}=3\times 10^{-3}\,$T' lt 1 pt 6 lw 2 lc 9
set output
!latex d11
!dvipdf d11.dvi
!cp -p d11.pdf /home/u6156/

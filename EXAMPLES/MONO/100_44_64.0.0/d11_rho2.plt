set log y
set key top center
set xlabel offset 0,0.0
set ylabel offset 1,0
set xlabel '$\psi/\psi_{LCFS}$'
set ylabel '${\hat D_{11}^*}$'
set xrange [0:0.9]
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
p "<awk '{if(($1>0.05)&&($1<0.9)&&($2==3e-5)&&($3==0e-0)) print}' results.data" u 1:(f*0.5*($6+$7)):(f*0.5*($6-$7)) w e lt 3 pt 5 ps 2 title '{\ttfamily DKES}, $v_{E*}=0        $'\
, "<awk '{if(($1>0.05)&&($1<0.9)&&($2==3e-5)&&($3==3e-4)) print}' results.data" u 1:(f*0.5*($6+$7)):(f*0.5*($6-$7)) w e lt 5 pt 5 ps 2 title '{\ttfamily{DKES}, $v_{E*}B_{0,0}=3\times 10^{-4}\,$T'\
, "<awk '{if(($1>0.05)&&($1<0.9)&&($2>2e-5)&&($2<4e-5)&&($3>-2e-4)&&($3<4e-9)) print}' results.knosos" u 1:(f*$7) w lp lt 3 pt 6 ps 1 title '{\ttfamily KNOSOS}, $v_{E*}=0        $'\
, "<awk '{if(($1>0.05)&&($1<0.9)&&($2>2e-5)&&($2<4e-5)&&($3 >2e-4)&&($3<4e-4)) print}' results.knosos" u 1:(f*$7) w lp lt 5 pt 6 ps 1 title '{\ttfamily KNOSOS}, $v_{E*}B_{0,0}=3\times 10^{-4}\,$T' 

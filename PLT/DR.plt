set term x11 title "`pwd`"


set xlabel '$r/a$'
set ylabel '$E_r/|\nabla r|\,$[kV/m]'
set key left
set xrange [0.47:]
#set yrange [0:3]
#set xtics 0.5,0.05
#set ytics 0.0.5
DR1(x)=a1+x*(b1+x*(c1+x*(d1+x*e1)))
fit DR1(x) "<awk '{if($5==1) print}' Er.map" u (sqrt($1)):($9/$8) via a1,b1,c1,d1,e1
DR2(x)=a2+x*(b2+x*(c2+x*(d2+x*e2)))
fit DR2(x) "<awk '{if($5==2) print}' Er.map" u (sqrt($1)):($9/$8) via a2,b2,c2,d2,e2
!rm fit.log
p "<awk '{if($1~/l/) print}' er_lr.dat" u 3:($5/DR1($3)):4:($6/DR1($3)) w xyerror pt 1 lt 1 lc 1 title '$E_r/|\nabla r|$, ~~left, DR'\
, "<awk '{if($1~/r/) print}' er_lr.dat" u 3:($5/DR2($3)):4:($6/DR2($3)) w xyerror pt 1 lt 1 lc 3 title '$E_r/|\nabla r|$, right, DR'\
, "Er.map" u (sqrt($1)):($10/$9*$8/1e3) w d lt 0 title '$-d_r(\varphi_0+\varphi_1)$, ~~~~~~~ {\ttfamily KNOSOS}'\
, "<awk '{if($5==1) print}' Er.map" u (sqrt($1)):($10/$9*$8/1e3) pt 7 ps 1 lt 1 title '$-d_r(\varphi_0+\varphi_1)$, ~~left, {\ttfamily KNOSOS}'\
, "<awk '{if($5==2) print}' Er.map" u (sqrt($1)):($10/$9*$8/1e3) pt 7 ps 1 lt 3 title '$-d_r(\varphi_0+\varphi_1)$, right, {\ttfamily KNOSOS}'


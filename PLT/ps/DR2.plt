set terminal epslatex color standalone
set output 'DR.tex'

set xlabel '$r/a$'
set ylabel '$E_r/|\nabla r|\,$[kV/m]'
set key left
set xrange [0.47:]
#set yrange [0:3]
#set xtics 0.5,0.05
#set ytics 0.0.5
DR1(x)=a1+x*(b1+x*(c1+x*(d1+x*e1)))
fit DR1(x) "<awk '{if($10==1) print}' Er.map" u (sqrt($1)):($8/$7) via a1,b1,c1,d1,e1
DR2(x)=a2+x*(b2+x*(c2+x*(d2+x*e2)))
fit DR2(x) "<awk '{if($10==1) print}' Er.map" u (sqrt($1)):($8/$7) via a2,b2,c2,d2,e2
!rm fit.log
p "<awk '{if($10>-1) print}' Er.map" u (sqrt($1)):($9/$8*$7/1e3) w d lt 0 title '$-d_r(\varphi_0+\varphi_1)$, ~~~~~~~ {\ttfamily KNOSOS}'\
, "<awk '{if($10==1) print}' Er.map" u (sqrt($1)):($9/$8*$7/1e3) pt 7 ps 1 lt 1 title '$-d_r(\varphi_0+\varphi_1)$, ~~left, {\ttfamily KNOSOS}'\
, "<awk '{if($10==2) print}' Er.map" u (sqrt($1)):($9/$8*$7/1e3) pt 7 ps 1 lt 3 title '$-d_r(\varphi_0+\varphi_1)$, right, {\ttfamily KNOSOS}'
#rep "<awk '{if($10==1) print}' Er.map" u (sqrt($1)):($7/1e3) pt 5 ps 1 lt 6 title '$-d_r\varphi_0$, {\ttfamily KNOSOS}'

set output
!latex DR
!dvipdf DR.dvi
!cp -p DR.pdf /home/u6156/



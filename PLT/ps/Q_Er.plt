reset
set terminal epslatex color standalone
set output 'Q_Er.tex'
set xlabel '$E_r\,$[kV/m]'
set ylabel '$Q_b\,$[MW/m$^2$]'
set xrange [-20:20]
set yrange [-0.02:0.4]
set xtics -20,10
set ytics 0,0.2
#set log y
#set format y '$10^{%T}$'
eVtoMJ=1.602e-6
p "KNOSOS+DKES/flux.amb.0" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w l lt 1 lc 9 lw 4 notitle\
, "KNOSOS+DKES/flux.amb.0" u ($2/1e3):(eVtoMJ*$13*$16*$18) w l lt 1 lc 9 lw 4 notitle\
, "KNOSOSlDKES/flux.amb.0" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w p pt 4 ps 1 lc 1 title '$\hat Q_e$, {\ttfamily KNOSOS}'\
, "KNOSOSlDKES/flux.amb.0" u ($2/1e3):(eVtoMJ*$13*$16*$18) w p pt 4 ps 1 lc 3 title '$\hat Q_i$, {\ttfamily KNOSOS}'\
,      "KNOSOS/flux.amb.0" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w p pt 5 ps 1 lc 1 title '$Q_e$, {\ttfamily KNOSOS}'\
,      "KNOSOS/flux.amb.0" u ($2/1e3):(eVtoMJ*$13*$16*$18) w p pt 5 ps 1 lc 3 title '$Q_i$, {\ttfamily KNOSOS}'\
, "KNOSOS+DKES/flux.amb.0" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w l lt 1 lc 9 lw 4 title ' {\ttfamily KNOSOS+DKES}'\
, "KNOSOS+DKES/flux.amb.0" u ($2/1e3):(eVtoMJ*$13*$16*$18) w l lt 1 lc 9 lw 4 notitle
#, "KNOSOS+DKES/flux.amb.0" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w l lt 7 lw 4 title '$Q_e$, {\ttfamily KNOSOS+DKES}'\
#, "KNOSOS+DKES/flux.amb.0" u ($2/1e3):(eVtoMJ*$13*$16*$18) w l lt 7 lw 4 title '$Q_i$, {\ttfamily KNOSOS+DKES}'
#p "KNOSOS+DKES/flux.amb.0" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w l lt 7 lw 4 title '$Q_e$, {\ttfamily KNOSOS+DKES}'\
#, "KNOSOSlDKES/flux.amb.0" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w p pt 4 ps 1 lt 1 title '$\hat Q_e$, {\ttfamily KNOSOS}'\
#,      "KNOSOS/flux.amb.0" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w p pt 5 ps 1 lt 1 title '$Q_e$, {\ttfamily KNOSOS}'\
#, "KNOSOS+DKES/flux.amb.0" u ($2/1e3):(eVtoMJ*$13*$16*$18) w l lt 7 lw 4 title '$Q_i$, {\ttfamily KNOSOS+DKES}'\
#, "KNOSOSlDKES/flux.amb.0" u ($2/1e3):(eVtoMJ*$13*$16*$18) w p pt 4 ps 1 lt 3 title '$\hat Q_i$, {\ttfamily KNOSOS}'\
#,      "KNOSOS/flux.amb.0" u ($2/1e3):(eVtoMJ*$13*$16*$18) w p pt 5 ps 1 lt 3 title '$Q_i$, {\ttfamily KNOSOS}'
set output
!latex Q_Er
!dvipdf Q_Er.dvi
!cp -p Q_Er.pdf /home/u6156/

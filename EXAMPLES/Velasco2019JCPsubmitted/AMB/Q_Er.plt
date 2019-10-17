reset
set xlabel '$E_r\,$[kV/m]'
set ylabel '$Q_b\,$[MW/m$^2$]'
set xrange [-20:20]
set yrange [-0.02:0.4]
set xtics -20,10
set ytics 0,0.2
#set log y
#set format y '$10^{%T}$'
eVtoMJ=1.602e-6
p      "KNOSOS/flux.amb" u ($2/1e3):(eVtoMJ* $4 *$7 *$9) w p pt 5 ps 1 lc 1 title '$Q_e$, {\ttfamily KNOSOS}'\
,      "KNOSOS/flux.amb" u ($2/1e3):(eVtoMJ*$13*$16*$18) w p pt 5 ps 1 lc 3 title '$Q_i$, {\ttfamily KNOSOS}'\

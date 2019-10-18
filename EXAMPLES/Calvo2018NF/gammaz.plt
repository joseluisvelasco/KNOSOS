!cat nophi/Er-*/imp.knosos       >       nophi/imp.knosos
!cat euterpe/Er-*/imp.knosos     >     euterpe/imp.knosos
!cat euterpe+cla/Er-*/imp.knosos > euterpe+cla/imp.knosos
!cat knosos/Er-*/imp.knosos      >      knosos/imp.knosos
set xrange [-15.5:-5.5]
set yrange [-20:1]
set xtics -14,2
set ytics -20,5
set xlabel "E_r [kV/m]"
set ylabel "Gamma_z/eta [m/s]"
set grid
set key bottom left
Ti=1.30374514099716 
p   "<awk '{if($2==24) print}' euterpe/imp.knosos" u ($9*Ti):5 w lp pt 6 ps 2 lt 3 title "w/o varphi_1, Z_z=24"
rep "<awk '{if($2==34) print}' euterpe/imp.knosos" u ($9*Ti):5 w lp pt 6 ps 2 lt 1 title "w/o varphi_1, Z_z=34"
rep "<awk '{if($2==44) print}' euterpe/imp.knosos" u ($9*Ti):5 w lp pt 6 ps 2 lt 2 title "w/o varphi_1, Z_z=44"
rep "<awk '{if($2==44) print}' euterpe+cla/imp.knosos" u ($9*Ti):5 w lp pt 6 ps 2 lt 0 notitle
rep "<awk '{if($2==24) print}' nophi/imp.knosos" u ($9*Ti):5 w lp pt 4 ps 2 lt 3 title "w.  varphi_1, Z_z=24"
rep "<awk '{if($2==34) print}' nophi/imp.knosos" u ($9*Ti):5 w lp pt 4 ps 2 lt 1 title "w.  varphi_1, Z_z=34"
rep "<awk '{if($2==44) print}' nophi/imp.knosos" u ($9*Ti):5 w lp pt 4 ps 2 lt 2 title "w.  varphi_1, Z_z=44"
pause -1
set xrange [-11.5:-8.5]
set yrange [-0.5:0.1]
set xtics -11,1
set ytics -0.4,0.2
set key top right
rep

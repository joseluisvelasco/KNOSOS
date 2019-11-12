set yrange [0:0.4]

set xlabel 'R [m]'
set ylabel 'Z [m]'

y0=0.872596422412738 
m=0.872355600749459
d=0.0275

set cblabel 'e\varphi_1/T_i'
p "Er.map" u (sqrt($13*$13+$12*$12)):((($13-y0-$12*m>-d)&&($13-y0-$12*m<+d))?$14:1/0):($4/$6) palette pt 5 ps 2 notitle

pause -1

l "/home/u6156/KNOSOS/PLT/palette.plt" 
set palette defined (-3 "blue", 0 "white", 1 "red")
set cblabel '-d_r\varphi_1 [kV/m]'
p "Er.map" u (sqrt($13*$13+$12*$12)):((($13-y0-$12*m>-d)&&($13-y0-$12*m<+d))?$14:1/0):(($10-$9)*$8/$9/1e3) palette pt 5 ps 2 notitle


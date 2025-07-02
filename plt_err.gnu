############################
#### File Selection #######

i = int(ARG1)
j = int(ARG2)

set title 'SCF Error per Iteration with Pulay = 0.3'
set xlabel 'Iteration'
set ylabel 'Error'
set key inside top right
#set yrange[0:1]

plot sprintf("err_V_n%d.dat",i)with lines title sprintf("Voltage = -%d V", i), sprintf("err_V_n%d.dat",j)with lines title sprintf("Voltage = -%d V", j)
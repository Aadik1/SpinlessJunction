
############################
#### Iteration Range #######

min = int(ARG1) 
max = int(ARG2)


##############################
#### Plots ###################

set title 'Error for SCF at Voltages'
set xlabel 'Iterations' 
set ylabel 'Error'
set yrange[0:10]
set xrange[0:150]

plot for [i=min:max] sprintf("err_V_n%d.dat",i) using 1:2 with lines title sprintf("Voltage: -%d eV",i), 'err_V_n0.dat' using 1:2 with lines title "Voltage: 0 eV"


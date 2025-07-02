
set title 'Current and dI/dV vs Voltage, Hub = 0.1eV'
set xlabel 'Voltge (eV)'
set ylabel 'Current (arb. units)'
set key inside bottom right
#set xrange[-1:1]

#plot 'Volt_Current_0.dat' u 1:2 w l title 'No Interactions', 'Volt_Current_1_sf.dat' u 1:2 w l title 'Hartree-Fock, Hub = 3 eV'

#plot 'Volt_Current_0.dat' u 1:2 w l title 'No Interactions', 'Volt_Current_1_sf.dat' u 1:2 w l title 'Hartree-Fock, Hub = 2 eV', 'Volt_Current_1_3sf.dat' u 1:2 w l title 'Hartree-Fock, Hub = 3 eV'

#, 'Volt_Current_2_sf.dat' u 1:2 w l title 'Second Born Approximation, Hub = 0.1 eV'

plot 'Volt_Current_2.dat' u 1:2 title 'Current' , 'dI_dV.dat' u 1:2 title 'dI/dV'


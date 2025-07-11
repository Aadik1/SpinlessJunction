
############################
#### Iteration Range #######

min = int(ARG1) 
max = int(ARG2)

#############################
##### Index of GF ###########

j = int(ARG3)	

##############################
#### Plots ###################

set title sprintf("Spectral Function vs Frequencies for G^{r}_{%d%d}",j,j)
set xlabel 'Frequency (w)' 
set ylabel 'A(w)'

plot for [i=min:max] sprintf("sf_%d.dat",i) using 1:(column(j+1)) with lines  title sprintf("Iteration %d",i)


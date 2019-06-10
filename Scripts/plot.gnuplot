set xrange [1:30]
plot 'light_100.phiL_0_rhoL.mass3pt.dat' using 1:2:($2-$3):($2+$4) with yerrorbars, \
     'light_75.phiL_0_rhoL.mass3pt.dat' using 1:2:($2-$3):($2+$4) with yerrorbars, \
     'light_50.phiL_0_rhoL.mass3pt.dat' using 1:2:($2-$3):($2+$4) with yerrorbars

set xrange [1:30]
plot 'phiL_0_rhoL.mass3pt.dat' using 1:2:($2-$3):($2+$4) with yerrorbars, \
    '' using ($1+0.5):($2+0.1):($2+0.1-$3):($2+0.1+$4) with yerrorbars, \
    '' using 1:(0.3 + 3 * exp(-0.8 * ($1+2))) smooth csplines

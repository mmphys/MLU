
set term pdf

set pointsize 0.3
set xlabel 'time'

set output "a_pseudoscalar_mass.pdf"
set title 'Pseudoscalar-Pseudoscalar ( {/Times:Italic π} = {/Times:Italic γ}_5 ) on 24^3 ensemble'
set xrange [1:30]
set ylabel 'Effective mass m(t)'
plot \
'light_100.phiL_0_rhoL.mass.dat' using 1:2:($2-$3):($2+$4) with yerrorbars title "100 ev/25 conf", \
'light_75.phiL_0_rhoL.mass.dat' using ($1+0.1):2:($2-$3):($2+$4) with yerrorbars title "75 ev/25 conf", \
'light_50.phiL_0_rhoL.mass.dat' using ($1+0.2):2:($2-$3):($2+$4) with yerrorbars title "50 ev/25 conf", \
'phiL_0_rhoL.mass.dat' using ($1+0.3):2:($2-$3):($2+$4) with yerrorbars title "Mike 100 ev/10 conf", \
'felix_100_10.phiL_0_rhoL.mass.dat' using ($1+0.4):2:($2-$3):($2+$4) with yerrorbars title "Felix 100 ev/10 conf", \
(0 * x + 0.1915) lt rgb "red" title "m_{/Symbol p}^{ref}=0.1915(8)"

set output "a_pseudoscalar_corr.pdf"
set xrange [0:63]
my_es=5 # error bar scale
set ylabel sprintf("Correlator C(t) x %d",my_es)
set logscale y
set key center top center #title 'Legend' box 3
plot \
'light_100.phiL_0_rhoL.corr.dat' using 1:2:($2-$3*my_es):($2+$4*my_es) with yerrorbars title "100 ev/25 conf", \
'light_75.phiL_0_rhoL.corr.dat' using ($1+0.1):2:($2-$3*my_es):($2+$4*my_es) with yerrorbars title "75 ev/25 conf", \
'light_50.phiL_0_rhoL.corr.dat' using ($1+0.2):2:($2-$3*my_es):($2+$4*my_es) with yerrorbars title "50 ev/25 conf", \
'phiL_0_rhoL.corr.dat' using ($1+0.3):2:($2-$3*my_es):($2+$4*my_es) with yerrorbars title "Mike 100 ev/10 conf", \
'felix_100_10.phiL_0_rhoL.corr.dat' using ($1+0.4):2:($2-$3*my_es):($2+$4*my_es) with yerrorbars title "Felix 100 ev/10 conf"

set output "a_axial_corr.pdf"
set title 'Axial-Axial ( {/Times:Italic γ}_0 {/Times:Italic γ}_5 ) on 24^3 ensemble'
plot \
'light_100.phiL_ax_0_rhoL_ax.corr.dat' using 1:2:($2-$3):($2+$4) with yerrorbars title "100 ev/25 conf", \
'light_75.phiL_ax_0_rhoL_ax.corr.dat' using ($1+0.1):2:($2-$3):($2+$4) with yerrorbars title "75 ev/25 conf", \
'light_50.phiL_ax_0_rhoL_ax.corr.dat' using ($1+0.2):2:($2-$3):($2+$4) with yerrorbars title "50 ev/25 conf"

set output "a_axial_mass.pdf"
set xrange [2:30]
set ylabel 'Effective mass m(t)'
unset logscale y
set key right top center
plot \
'light_100.phiL_ax_0_rhoL_ax.mass.dat' using 1:2:($2-$3):($2+$4) with yerrorbars title "100 ev/25 conf", \
'light_75.phiL_ax_0_rhoL_ax.mass.dat' using ($1+0.1):2:($2-$3):($2+$4) with yerrorbars title "75 ev/25 conf", \
'light_50.phiL_ax_0_rhoL_ax.mass.dat' using ($1+0.2):2:($2-$3):($2+$4) with yerrorbars title "50 ev/25 conf", \
(0 * x + 0.1915) lt rgb "red" title "m_{/Symbol p}^{ref}=0.1915(8)"

set output "a_vector_mass.pdf"
set title 'Vector-Vector ( {/Symbol:Italic r} = {/Times:Italic γ}_i ) on 24^3 ensemble'
plot \
'light_100.phiL_gx_0_rhoL_gx.mass.dat' using 1:2:($2-$3):($2+$4) with yerrorbars title "{/Times:Italic γ}_x 100 ev/25 conf", \
'light_100.phiL_gy_0_rhoL_gy.mass.dat' using ($1+0.1):2:($2-$3):($2+$4) with yerrorbars title "{/Times:Italic γ}_y 100 ev/25 conf", \
'light_100.phiL_gz_0_rhoL_gz.mass.dat' using ($1+0.2):2:($2-$3):($2+$4) with yerrorbars title "{/Times:Italic γ}_z 100 ev/25 conf", \
(0 * x + 0.1915) lt rgb "red" title "m_{/Symbol p}^{ref}=0.1915(8)"

set output "a_vector_corr.pdf"
set xrange [0:63]
my_es=1 # error bar scale
set ylabel sprintf("Correlator C(t) x %d",my_es)
set logscale y
set key center top center #title 'Legend' box 3
plot \
'light_100.phiL_gx_0_rhoL_gx.corr.dat' using 1:2:($2-$3):($2+$4) with yerrorbars title "{/Times:Italic γ}_x 100 ev/25 conf", \
'light_100.phiL_gy_0_rhoL_gy.corr.dat' using ($1+0.1):2:($2-$3):($2+$4) with yerrorbars title "{/Times:Italic γ}_y 100 ev/25 conf", \
'light_100.phiL_gz_0_rhoL_gz.corr.dat' using ($1+0.2):2:($2-$3):($2+$4) with yerrorbars title "{/Times:Italic γ}_z 100 ev/25 conf"

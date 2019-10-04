set term pdf

set pointsize 0.3
set xlabel 'time'

set output "a_" . prefix . ".meff.pdf"
suffix = sprintf(".cosh.%s.txt", seq)
file1 = prefix . "_ps_ps" . suffix
file2 = prefix . "_ps_ax" . suffix
file3 = prefix . "_ax_ps" . suffix
file4 = prefix . "_ax_ax" . suffix

set title "V" . inVec . ", ΔT=" . inDeltaT . ", n^2=" . inP2
set xrange [1:10]
set yrange [0.35:0.85]
set ylabel 'a m_{eff}'
plot \
  file1 using ($1-0.1):2:($2-$3):($2+$4) with yerrorbars title "ps - ps" lc rgb "red", \
  file2 using 1:2:($2-$3):($2+$4) with yerrorbars title "ps - ax" lc rgb "blue", \
  file3 using ($1+0.1):2:($2-$3):($2+$4) with yerrorbars title "ax - ps" lc rgb "green", \
  file4 using ($1+0.2):2:($2-$3):($2+$4) with yerrorbars title "ax - ax" lc rgb "magenta", #\
  (0 * x + 0.1915) lt rgb "red" title "m_{/Symbol p}^{ref}=0.1915(8)"

#----------------------------------------------

#set output "a_" . prefix . ".corr.nonlog.pdf"
suffix = sprintf(".corr.%s.txt", seq)
file1 = prefix . "_ps_ps" . suffix
file2 = prefix . "_ps_ax" . suffix
file3 = prefix . "_ax_ps" . suffix
file4 = prefix . "_ax_ax" . suffix

set title "V" . inVec . ", ΔT=" . inDeltaT . ", n^2=" . inP2
set xrange [0:14]
set autoscale y
set ylabel 'Correlator C_3(t)'

#plot \
file1 using ($1-0.1):2:($2-$3):($2+$4) with yerrorbars title "ps - ps" lc rgb "red", \
file2 using 1:2:($2-$3):($2+$4) with yerrorbars title "ps - ax" lc rgb "blue", \
file3 using ($1+0.1):2:($2-$3):($2+$4) with yerrorbars title "ax - ps" lc rgb "green", \
file4 using ($1+0.2):2:($2-$3):($2+$4) with yerrorbars title "ax - ax" lc rgb "magenta"

set output "a_" . prefix . ".corr.pdf"
set logscale y
set format y "10^{%L}"
plot \
  file1 using ($1-0.1):(abs($2)):($2-$3):($2+$4) with yerrorbars title "ps - ps" lc rgb "red", \
  file2 using 1:(abs($2)):($2-$3):($2+$4) with yerrorbars title "ps - ax" lc rgb "blue", \
  file3 using ($1+0.1):(abs($2)):($2-$3):($2+$4) with yerrorbars title "ax - ps" lc rgb "green", \
  file4 using ($1+0.2):(abs($2)):($2-$3):($2+$4) with yerrorbars title "ax - ax" lc rgb "magenta"

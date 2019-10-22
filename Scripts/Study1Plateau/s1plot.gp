
set term pdf

set pointsize 0.6
set xlabel 'initial fit time'

MyField="E0"

set output "a_fit_".MyField.".pdf"
set title MyField.' from 2-exponential fit ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'

set xrange [4.8:11.2]

set ylabel MyField.'(t)'
#set palette defined ( 15 "blue", 16 "red", 17 "green")
plot for [idx=0:*] FileParams index idx using ($1-0.1+0.05*idx):MyField:(column(MyField)-column(MyField."ErLow")):(column(MyField)+column(MyField."ErHigh")) with yerrorbars title columnhead(2), \
  (0.99656) with lines lt rgb "gray" title "E_0=0.99656(95) ArXiv:1812.08791"

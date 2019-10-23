#show colornames
set term pdf

set pointsize 0.6
set xlabel 'initial fit time'
set xrange [4.8:11.2]
#set palette defined ( 15 "blue", 16 "red", 17 "green")

do for [MyFieldNum = 0:1] {
  MyField="E".MyFieldNum
  set output "a_fit_".MyField.".pdf"
  if (MyFieldNum==0) {
    set arrow from 4.8,0.99565 to 11.2,0.99565 nohead front lc rgb "gray40" lw 0.25  dashtype "-"
    set label "E_0=0.99656(95), ArXiv:1812.08791" at 4.9,0.99656 font "Arial,8" front textcolor "grey40" offset character 0,0.2
    TitleFieldName=MyField
  } else {
    unset arrow
    unset label
    TitleFieldName="Δ".MyField
  }
  set title TitleFieldName.' from 2-exponential fit ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
  set ylabel TitleFieldName.'(t)'
  plot for [idx=0:*] FileParams index idx using ($1-0.1+0.05*idx):MyField:(column(MyField)-column(MyField."ErLow")):(column(MyField)+column(MyField."ErHigh")) with yerrorbars title columnhead(2)
}

#set pointsize 0.6
#set xlabel 'initial fit time'

set output "a_fit_chisq.pdf"
set title '{/Times:Italic χ}^2 per d.o.f. dependence on initial/final fit times ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
set xrange [4.8:12.2]
#set yrange [0.95:1.1]
set ylabel 'Chi squared per degree of freedom'
plot for [idx=0:*] FileParams index idx using ($1-0.10):"ChiSqPerDof":(column(-2)+1):(column(-2)+1) with linespoints pt variable lc variable title columnhead(2)

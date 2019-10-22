set term pdf

set pointsize 0.6
set xlabel 'initial fit time'

MyTitle(n)="tf=".n

set output "a_fit_chisq.pdf"
set title '{/Times:Italic χ}^2 per d.o.f. dependence on initial/final fit times ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
set xrange [5.8:12.2]
#set yrange [0.95:1.1]
set ylabel 'Chi squared per degree of freedom'
#set palette defined ( 15 "blue", 16 "red", 17 "green")
plot for [idx=0:*] FileParams index idx using ($1-0.10):"ChiSqPerDof":(column(-2)+1):(column(-2)+1) with linespoints pt variable lc variable title columnhead(2)

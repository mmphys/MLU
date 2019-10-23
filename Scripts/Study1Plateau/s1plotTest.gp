set term pdf

#Function to return operator name
OpName(n) = (n==0) ? "g5" : "gT5"
OpText(n) = (n==0) ? "Γ_5" : "Γ_4 Γ_5"
YOffset(snk,src) = (snk==0 && src==1) ? -.002 : 0

Seq="4147798751"
ti=6
tf=16
prefix="h1_l_p_0_0_0"
suffix=".txt"
InBasename=prefix.".corr_".ti."_".tf."."

set pointsize 0.6
set xlabel 'initial fit time'
set xrange [4.8:11.2]

set output "a_reconstruct_".ti."_".tf.".pdf"
set title 'Masses reconstructed from fit on ['.ti.','.tf.'] ( {/Times:Italic Γ}_5, {/Times:Italic Γ}_4 {/Times:Italic Γ}_5 on 24^3 ensemble)'
set xrange [4.8:12.2]
set ylabel 'Synthesised correlator mass from fit'

if (YOffset(0,1)!=0) {
  set label OpText(0)."-".OpText(1)." offset to avoid overlap with ".OpText(1)."-".OpText(0) \
    at graph 0.99, 0.025 right font "Arial,8" front textcolor "grey40"
}

set style fill transparent solid 0.2 noborder

#plot for [src=0:1] for [snk=0:1] InBasename.OpName(snk)."_".OpName(src).".mass.".Seq.suffix \
  using 1:2:($2-$3):($2+$4) with yerrorbars title OpText(snk)." - ".OpText(src)

plot for [snk=0:1] for [src=0:1] FileParams.InBasename.OpName(snk)."_".OpName(src).".mass.".Seq.suffix \
  using 1:($2-$3+YOffset(snk,src)):($2+$4+YOffset(snk,src)) with filledcurves title "Error ".OpText(snk)." - ".OpText(src), \
  for [snk=0:1] for [src=0:1] FileParams.InBasename.OpName(snk)."_".OpName(src).".mass.".Seq.suffix \
  using 1:($2+YOffset(snk,src)) with lines title "Midpoint ".OpText(snk)." - ".OpText(src)

#using ($1-0.10):"ChiSqPerDof":(column(-2)+1):(column(-2)+1)

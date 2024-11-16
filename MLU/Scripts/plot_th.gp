#!/usr/bin/env gnuplot

###################################################
# Make a plot of simulated 3pt function
###################################################

NumHeavy=4
array HeavyName[NumHeavy]
do for [i=1:NumHeavy] { HeavyName[i]=sprintf("h%d",i-1) }
array Energy0[NumHeavy] = [ 0.96698, 1.05284, 1.11479, 1.16423 ]
array ZPoint0[NumHeavy] = [ 34.62, 34.83, 34.70, 34.28 ]
array ZWall0[NumHeavy]  = [ 12424., 12147., 11900., 11646. ]
array Energy1[NumHeavy] = [ 1.4241, 1.5035, 1.5570, 1.5948 ]
array ZPoint1[NumHeavy] = [ 47.3, 51.7, 55.0, 57.5 ]
array ZWall1[NumHeavy]  = [ -7703., -6961., -5915., -4641. ]

if( 0 ) {
do for [i=1:NumHeavy] {
  print i, ": ".HeavyName[i].", E0=", Energy0[i], ", ZPoint0=", ZPoint0[i], ", ZWall0=", ZWall0[i], \
          Energy1[i], ", ZPoint0=", ZPoint1[i], ", ZWall0=", ZWall1[i]
}}

set term pdf
set pointsize 0.6

NumDeltaT=7
array DeltaTArray[NumDeltaT] = [ 12, 14, 16, 20, 24, 28, 32 ]
do for [DeltaTIndex=1:NumDeltaT] {
DeltaT=DeltaTArray[DeltaTIndex]

# 0 = Point, 1 = Wall
do for [SnkT=0:1] {
do for [SrcT=0:1] {
# Which heavy are we interested in at source and sink?
snk = 0
src = 3

sTypSnk=(SnkT == 0 ? "P" : "W")
sTypSrc=(SrcT == 0 ? "P" : "W")
EndPoint(epNum,epType) = sprintf("h%d", epNum).(epType == 0 ? "P" : "W")
SnkSrc = "Effective (exp) mass ".EndPoint(snk,SnkT)." <- ".EndPoint(src,SrcT)

set output EndPoint(snk,SnkT)."_".EndPoint(src,SrcT).sprintf("_dt_%d",DeltaT).".pdf"

# Set Variables
iSnk = snk + 1
iSrc = src + 1
Ef0 = Energy0[iSnk]
Ef1 = Energy1[iSnk]
Ei0 = Energy0[iSrc]
Ei1 = Energy1[iSrc]
Af0 = SnkT == 0 ? ZPoint0[iSnk] : ZWall0[iSnk]
Af1 = SnkT == 0 ? ZPoint1[iSnk] : ZWall1[iSnk]
Ai0 = SrcT == 0 ? ZPoint0[iSrc] : ZWall0[iSrc]
Ai1 = SrcT == 0 ? ZPoint1[iSrc] : ZWall1[iSrc]
Exp_f0dt = exp(-Ef0 * DeltaT)
Exp_f1dt = exp(-Ef1 * DeltaT)
Prefactor_f0i0 = Af0 * Ai0 / ( 4 * Ef0 * Ei0 )
Prefactor_f1i0 = Af1 * Ai0 / ( 4 * Ef1 * Ei0 )
Prefactor_f0i1 = Af0 * Ai1 / ( 4 * Ef0 * Ei1 )
Prefactor_f1i1 = Af1 * Ai1 / ( 4 * Ef1 * Ei1 )

array CGnd[DeltaT + 1]
array CGndEx[DeltaT + 1]
array CEx[DeltaT + 1]
do for [t=0:DeltaT] {
  Exp_f0i0 = exp(t * (Ef0 - Ei0))
  Exp_f1i0 = exp(t * (Ef1 - Ei0))
  Exp_f0i1 = exp(t * (Ef0 - Ei1))
  Exp_f1i1 = exp(t * (Ef1 - Ei1))
  CGnd[t + 1] = Prefactor_f0i0 * Exp_f0i0 * Exp_f0dt
  CGndEx[t + 1] = CGnd[t + 1] + Prefactor_f1i0 * Exp_f1i0 * Exp_f1dt + Prefactor_f0i1 * Exp_f0i1 * Exp_f0dt
  CEx[t + 1] = CGndEx[t + 1] + Prefactor_f1i1 * Exp_f1i1 * Exp_f1dt
}

if( 0 ) {
set logscale y
  plot CGnd
}

array aEffGnd[DeltaT]
array aEffGndEx[DeltaT]
array aEffEx[DeltaT]
do for [t=1:DeltaT] {
  aEffGnd[t] = t - 0.6 + {0,1} * log(CGnd[t] / CGnd[t+1])
  aEffGndEx[t] = t - 0.5 + {0,1} * log(CGndEx[t] / CGndEx[t+1])
  aEffEx[t] = t - 0.4 + {0,1} * log(CGndEx[t] / CGndEx[t+1])
}

set title SnkSrc
plot aEffGnd using 2:3 title "Gnd", aEffGndEx using 2:3 title "+Gnd-Ex", aEffEx using 2:3 title "+Ex-Ex"
}
}
}

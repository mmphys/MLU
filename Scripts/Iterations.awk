#!/opt/local/bin/gawk -f

# Tell me about the iterations for each CG

# Environment variables for maximum inner/outer iteration counts/quark
# ${q}in  Maximum inner iteration count. Default light/other=15000/5000
# ${q}out Maximum outer iteration count. Default 100

# Remember which quark we are processing
/Computing quark propagator/ {
  PropName=$11
  PropName=substr(PropName,2,length(PropName)-2)
  #print $0
  split(PropName, Prop, "_")
  switch( Prop[1] ) {
    case "SeqProp":
      q=Prop[4]
      break
    default:
      q=Prop[3]
      break
  }
}

# Validate iteration counts for the current quark
/Inner CG iterations/ {
  Inner=$12
  Outer=$14
  Final=$18
  #print q, Inner, Outer, Final

  if (q in InMin)
  {
    if ( Inner < InMin[q] )
      InMin[q] = Inner
    if ( Inner > InMax[q] )
      InMax[q] = Inner
    if ( Outer < OutMin[q] )
      OutMin[q] = Outer
    if ( Outer > OutMax[q] )
      OutMax[q] = Outer
    if ( Final < FinalMin[q] )
        FinalMin[q] = Final
    if ( Final > FinalMax[q] )
        FinalMax[q] = Final
  }
  else
  {
    # This is the first iteration count for this quark
    InMin[q]=Inner
    InMax[q]=Inner
    OutMin[q]=Outer
    OutMax[q]=Outer
    FinalMin[q]=Final
    FinalMax[q]=Final
    # See whether max iteration counts have been passed in
    MaxInner[q]=q"in"  in ENVIRON ? ENVIRON[q"in" ] + 0 :
                q == "l" ? 15000 : 5000
    MaxOuter[q]=q"out" in ENVIRON ? ENVIRON[q"out"] + 0 : 100
    print "MaxInner[" q "]=" MaxInner[q]
    print "MaxOuter[" q "]=" MaxOuter[q]
  }
  # If I'm outside the allowed iteration counts, retrun error
  if ((MaxInner[q] > 0 && (Inner >= MaxInner[q] || Final >= MaxInner[q])) ||
      (MaxOuter[q] > 0 && Outer >= MaxOuter[q]))
  {
    MyErr = 1
    print PropName
  }
}

END {
  print "Filename Quark InnerMin InnerMax OuterMin OuterMax FinalMin FinalMax"
  for (q in InMin) {
    print FILENAME, q, InMin[q], InMax[q], OutMin[q], OutMax[q], FinalMin[q], FinalMax[q]
  }
  exit MyErr
}

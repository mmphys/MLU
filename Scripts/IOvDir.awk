#!/opt/local/bin/gawk -f

# Parse the output of Benchmark_IO_vs_dir

BEGIN { print "Stripe Size SizeSuffix IOType IODir Rate RateUnits Num NumUnits Time TimeUnits Dir" }

/-- Directory/ {
Dir=$10
Num=split(Dir,DirArray,"/")
split(DirArray[Num],StripeSize,"_")
Stripe=0 + substr(StripeSize[1],7)
Size  =0 + substr(StripeSize[2],5,length(StripeSize[2]) - 5)
SizeSuffix=substr(StripeSize[2],length(StripeSize[2]))
}

/ bytes in / {
  if ($8=="IOobject:")
  {
    Field=9
    IOType="C-Lime"
  }
  else
  {
    Field=11
    IOType="std-IO"
  }
  ParseIO(Field,IOType)
}

function ParseIO(Field,IOType)
{
  # Read/Write sometimes initial capped, and sometimes followed by ':'
  if( substr($Field,2,3) == "ead" )
    IODir="Read"
  else
    IODir="Write"
  # Time units sometimes has a trailing ','
  TimeUnits=$(Field+5)
  if( substr(TimeUnits,length(TimeUnits)) == "," )
  {
    TimeUnits=substr(TimeUnits,1,length(TimeUnits)-1)
  }
  print Stripe, Size, SizeSuffix, IOType, IODir, $(Field+6), $(Field+7), $(Field+1), $(Field+2),
        $(Field+4), TimeUnits, Dir
}

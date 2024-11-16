# This is an example FitRanges.sh
# See Utility/Scripts/MakeFit.sh

# These are the fit ranges for F1M as at Tue 19 July 2022,
# however see the F1M directory on Tursa for the most up-to-date

function FitRanges()
{
  case $1 in
    # Can also set options, e.g.: Options="--Hotelling 0.005";;
    hh) Range=6:18:10:10,R0:-2:0:5:1; Range2="1,";;
    hs) Range=6:16:10:10,R0:-2:0:5:1; Range2="1,";;  #5:10:6:7;;
    hl) Range=6:14:8:10,R0:-2:0:5:1 ; Range2="1,";;
    ss) Range=4:10:6:15,R0:-2:0:5:1 ; Range2="1,";;
    sl) Range=4:10:7:14,R0:-2:0:5:1 ; Range2="1,";;
    ll) Range=2:6:6:18,R0:-2:0:5:1  ; Range2="1,";;
  esac
}

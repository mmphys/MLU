HL_Contract_Num="$1"
HL_Contract_QL="$2"
HL_Contract_QR="$3"
HL_Contract_P="$4"

if [[ "$HL_Contract_Num" == "" || "$HL_Contract_QL" == "" \
   || "$HL_Contract_QR" == "" || "$HL_Contract_P" == "" ]]
then
    echo "Usage: ${0##*/} prefix quarkL quarkR momentum"
    exit 1
fi

case "$HL_Contract_Num" in
  0) HL_Contract="c30p";;
  1) HL_Contract="c3p0";;
  2) HL_Contract="c2p0";;
  *) echo "Error: Contraction type \"$HL_Contract_Num\" should be 0, 1 or 2"; exit 1;;
esac

HL_Contract_Base=${HL_Contract}_${HL_Contract_QL}_${HL_Contract_QR}_p2_${HL_Contract_P}
HL_Contract_XML=${HL_Contract_Base}.xml
if [[ ! -r "$HL_Contract_XML" ]]
then
    echo "Error: \"$HL_Contract_XML\" does not exist"
    exit 1
fi

# Create the output directory
HL_Contract_Output="${HL_Contract_Base}"
mkdir -p "${HL_Contract_Output}"
if [[ ! -d "$HL_Contract_Output" ]]
then
    echo "Error: Can't create ouput directory \"$HL_Contract_Output\""
    exit 1
fi

# Clean up the tmp directory. Optional - contractor creates this
#HL_Contract_Tmp="tmp_${HL_Contract_Base}"
#rm -rf ${HL_Contract_Tmp}
#mkdir -p ${HL_Contract_Tmp}

#create the PBS script
HL_Contract_PBS="${HL_Contract_Base}.pbs"
export HL_Contract_XML
export HL_Contract_JobName="${HL_Contract_QL}${HL_Contract_QR}c${HL_Contract_Num}p${HL_Contract_P}"
cat hlcontract.pbs | envsubst '$HL_Contract_XML $HL_Contract_JobName' > "$HL_Contract_PBS"

#now submit the command
cmd="qsub \"$HL_Contract_PBS\""
echo "$cmd"
eval "$cmd"

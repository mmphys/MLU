#!/bin/bash

step=40
#configs=`seq 3040 $step 3440`
configs=3040
echo configs=$configs
#tplt_xml=pion.template
tplt_xml=pion_load_peramb.template
#tplt_xml=pion_LapEvec.template

case $tplt_xml in
  pion*.template) tplt_pbs=pion.template;;
  *)              tplt_pbs=$tplt_xml;;
esac

mkdir -p sub
mkdir -p log
for i in $configs
do
    j=$((i + step))
    this_pbs=sub/${tplt_pbs/.template/}.$i.pbs
    this_xml=sub/${tplt_xml/.template/}.$i.xml
    sed -e "s|{cfg}|$i|g" -e "s|{end}|$j|g" -e "s|{step}|$step|g" -e "s|{xml}|$this_xml|g" $tplt_pbs.pbs > $this_pbs
    sed -e "s|{cfg}|$i|g" -e "s|{end}|$j|g" -e "s|{step}|$step|g" $tplt_xml.xml > $this_xml
    cd log
    qsub ../$this_pbs
    cd ..
done

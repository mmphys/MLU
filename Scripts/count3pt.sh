#!/usr/bin/env bash

for f in ${1:-[ZG][F2][PW]*.xml}
do
echo "$f"
echo "   Modules: $(grep -c '<module>' $f)"
echo "       2pt: $(grep -c '<name>meson_2' $f)"
echo "       3pt: $(grep -c '<name>meson_3' $f)"
echo "    Prop  : $(grep -c '<name>Prop' $f)"
echo "    Prop l: $(egrep -c '<name>Prop_.*_l_p_' $f)"
echo "    Prop s: $(egrep -c '<name>Prop_.*_s_p_' $f)"
echo "    Prop h: $(egrep -c '<name>Prop_.*_h[0-9]_p_' $f)"
echo "Seq prop  : $(grep -c '<name>SeqProp_*_' $f)"
echo "Seq prop l: $(egrep -c '<name>SeqProp_.*_l_g' $f)"
echo "Seq prop s: $(egrep -c '<name>SeqProp_.*_s_g' $f)"
echo "Seq prop h: $(egrep -c '<name>SeqProp_.*_h[0-9]_g' $f)"
done

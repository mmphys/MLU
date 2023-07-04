#!/usr/bin/env bash

export Ensemble=${Ensemble:-C2}
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit 2
fi

#set -x
set -e

# Initialise variables from user

MELFit=${MELFit:-MELFit}

# Derived variables

OutBase=$Ensemble/$MELFit/Fit_sp2_

############################################################

# Make the fit selection files I use

############################################################

for f in disp dispind; do
[ -e ${OutBase}$f.txt ] || cat > ${OutBase}$f.txt <<- EndFitChoices
	# Ensemble: $Ensemble
	# This is my favoured set of fit choices
	# D_s
	h6413_s_p2_0  2ptp2/h6413_s/h6413_s_p2_0.corr_8_24_13_24.g5P_g5W.model.694228835.h5
	# K - using dispersion relation
	s_l_p2_0  2ptp2/s_l/s_l.corr_5_20_6_22_5_19_6_20_5_17.g5P.model.694228835.h5
	s_l_p2_1  2ptp2/s_l/s_l.corr_5_20_6_22_5_19_6_20_5_17.g5P.model.694228835.h5
	s_l_p2_2  2ptp2/s_l/s_l.corr_5_20_6_22_5_19_6_20_5_17.g5P.model.694228835.h5
	s_l_p2_3  2ptp2/s_l/s_l.corr_5_20_6_22_5_19_6_20_5_17.g5P.model.694228835.h5
	s_l_p2_4  2ptp2/s_l/s_l.corr_5_20_6_22_5_19_6_20_5_17.g5P.model.694228835.h5
EndFitChoices
done

############################################################

# Make the form factors

############################################################

series='disp dispind' PlotFormFactor.sh
series='disp dispind' UnCorr= PlotFormFactor.sh

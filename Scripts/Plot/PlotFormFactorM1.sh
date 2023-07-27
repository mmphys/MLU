#!/usr/bin/env bash

export Ensemble=${Ensemble:-M1}
if ! [ -d $Ensemble ]; then
  echo "Ensemble $Ensemble doesn't exist. Change directory?"
  exit 2
fi

#set -x
set -e

# Initialise variables from user

MELFit=${MELFit:-MELFit}

# Derived variables

OutDir=$Ensemble/$MELFit
OutBase=Fit_sp2_

############################################################

# Make the fit selection files I use

############################################################

(
mkdir -p $OutDir
cd $OutDir

[ -e ${OutBase}disp.txt ] || cat > ${OutBase}disp.txt <<- EndFitChoices
	# Ensemble: $Ensemble
	# This is my favoured set of fit choices
	# D_s
	h447_s_p2_0  2ptp2/h447_s/h447_s_p2_0.corr_10_28_15_28.g5P_g5W.model.1835672416.h5
	# K - using dispersion relation
	s_l_p2_0  2ptp2/s_l/s_l.corr_6_18_6_18_6_18_6_16_7_19.g5P.model.1835672416.h5
	s_l_p2_1  2ptp2/s_l/s_l.corr_6_18_6_18_6_18_6_16_7_19.g5P.model.1835672416.h5
	s_l_p2_2  2ptp2/s_l/s_l.corr_6_18_6_18_6_18_6_16_7_19.g5P.model.1835672416.h5
	s_l_p2_3  2ptp2/s_l/s_l.corr_6_18_6_18_6_18_6_16_7_19.g5P.model.1835672416.h5
	s_l_p2_4  2ptp2/s_l/s_l.corr_6_18_6_18_6_18_6_16_7_19.g5P.model.1835672416.h5
EndFitChoices

[ -e ${OutBase}renorm.txt ] || ln -s ${OutBase}disp.txt ${OutBase}renorm.txt
)

############################################################

# Make the form factors

############################################################

series=disp PlotFormFactor.sh
series=disp UnCorr= PlotFormFactor.sh
series=renorm ZV= PlotFormFactor.sh

#!/usr/bin/env bash

export Ensemble=${Ensemble:-M2}
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
	h447_s_p2_0  2ptp2/h447_s/h447_s_p2_0.corr_10_26_19_26.g5P_g5W.model.2475361470.h5
	# K - using dispersion relation
	s_l_p2_0  2ptp2/s_l/s_l.corr_6_20_6_19_7_20_7_18_6_15.g5P.model.2475361470.h5
	s_l_p2_1  2ptp2/s_l/s_l.corr_6_20_6_19_7_20_7_18_6_15.g5P.model.2475361470.h5
	s_l_p2_2  2ptp2/s_l/s_l.corr_6_20_6_19_7_20_7_18_6_15.g5P.model.2475361470.h5
	s_l_p2_3  2ptp2/s_l/s_l.corr_6_20_6_19_7_20_7_18_6_15.g5P.model.2475361470.h5
	s_l_p2_4  2ptp2/s_l/s_l.corr_6_20_6_19_7_20_7_18_6_15.g5P.model.2475361470.h5
EndFitChoices

for series in renorm AltZV Jan24
do
  [ -e ${OutBase}$series.txt ] || ln -s ${OutBase}disp.txt ${OutBase}$series.txt
done
)

############################################################

# Make the form factors

############################################################

series=disp PlotFormFactor.sh
series=disp UnCorr= PlotFormFactor.sh
# renorm=(mostly NPR) renormalised. No need for ZV
series='renorm renormC AltZV Jan24' Suffix=_mostly ZV= PlotFormFactor.sh
# renorm=(mostly NPR) renormalised. Apply Fully NPR correction
series='renorm renormC AltZV Jan24' FullyNP= PlotFormFactor.sh

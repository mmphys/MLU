#!/usr/bin/env bash

export Ensemble=${Ensemble:-F1M}
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

[ -e ${OutBase}old.txt ] || cat > ${OutBase}old.txt <<- EndFitChoices
	# Ensemble: $Ensemble
	# This is my favoured set of fit choices
	# D_s
	h385_s_p2_0  2ptp2/h385_s/h385_s_p2_0.corr_6_29_5_29.g5P_g5W.model.3285139232.h5
	# K
	s_l_p2_0  2ptp2/s_l/s_l_p2_0.corr_10_26_7_26.g5P_g5W.model.3285139232.h5
	s_l_p2_1  2ptp2/s_l/s_l_p2_1.corr_8_28_9_28.g5P_g5W.model.3285139232.h5
	s_l_p2_2  2ptp2/s_l/s_l_p2_2.corr_8_26_8_28.g5P_g5W.model.3285139232.h5
	s_l_p2_3  2ptp2/s_l/s_l_p2_3.corr_8_28_8_28.g5P_g5W.model.3285139232.h5
	s_l_p2_4  2ptp2/s_l/s_l_p2_4.corr_10_28_8_28.g5P_g5W.model.3285139232.h5
	s_l_p2_5  2ptp2/s_l/s_l_p2_5.corr_10_28_12_28.g5P_g5W.model.3285139232.h5
	s_l_p2_6  2ptp2/s_l/s_l_p2_6.corr_10_28_12_28.g5P_g5W.model.3285139232.h5
EndFitChoices

[ -e ${OutBase}better.txt ] || cat > ${OutBase}better.txt <<- EndFitChoices
  # Ensemble: $Ensemble
  # Work in progress
  # These first couple of fits reselected
  # D_s
  h385_s_p2_0  2ptp2/h385_s/h385_s_p2_0.corr_13_28_13_29.g5P_g5W.model.3285139232.h5
  # K
  s_l_p2_0  2ptp2/s_l/s_l_p2_0.corr_7_31_8_20.g5P_g5W.model.3285139232.h5
  s_l_p2_1  2ptp2/s_l/s_l_p2_1.corr_10_26_9_21.g5P_g5W.model.3285139232.h5
  s_l_p2_2  2ptp2/s_l/s_l_p2_2.corr_9_21_8_21.g5P_g5W.model.3285139232.h5
  s_l_p2_3  2ptp2/s_l/s_l_p2_3.corr_9_21_8_21.g5P_g5W.model.3285139232.h5
  s_l_p2_4  2ptp2/s_l/s_l_p2_4.corr_9_19_10_17.g5P_g5W.model.3285139232.h5
  s_l_p2_5  2ptp2/s_l/s_l_p2_5.corr_7_28_12_28.g5P_g5W.model.3285139232.h5
  s_l_p2_6  2ptp2/s_l/s_l_p2_6.corr_7_23_12_20.g5P_g5W.model.3285139232.h5
EndFitChoices

[ -e ${OutBase}disp.txt ] || cat > ${OutBase}disp.txt <<- EndFitChoices
  # Ensemble: $Ensemble
  # Work in progress
  # These first couple of fits reselected
  # D_s
  h385_s_p2_0  2ptp2/h385_s/h385_s_p2_0.corr_13_28_13_29.g5P_g5W.model.3285139232.h5
  # K
  s_l_p2_0  2ptp2/s_l/s_l.corr_7_30_10_26_9_19_9_19_9_19_7_26_7_23.g5P.model.3285139232.h5
  s_l_p2_1  2ptp2/s_l/s_l.corr_7_30_10_26_9_19_9_19_9_19_7_26_7_23.g5P.model.3285139232.h5
  s_l_p2_2  2ptp2/s_l/s_l.corr_7_30_10_26_9_19_9_19_9_19_7_26_7_23.g5P.model.3285139232.h5
  s_l_p2_3  2ptp2/s_l/s_l.corr_7_30_10_26_9_19_9_19_9_19_7_26_7_23.g5P.model.3285139232.h5
  s_l_p2_4  2ptp2/s_l/s_l.corr_7_30_10_26_9_19_9_19_9_19_7_26_7_23.g5P.model.3285139232.h5
  s_l_p2_5  2ptp2/s_l/s_l.corr_7_30_10_26_9_19_9_19_9_19_7_26_7_23.g5P.model.3285139232.h5
  s_l_p2_6  2ptp2/s_l/s_l.corr_7_30_10_26_9_19_9_19_9_19_7_26_7_23.g5P.model.3285139232.h5
EndFitChoices

for series in renorm AltZV
do
  [ -e ${OutBase}$series.txt ] || ln -s ${OutBase}disp.txt ${OutBase}$series.txt
done
)

############################################################

# Make the form factors

############################################################

series='old' PlotFormFactor.sh
series='old' FitSeries='std' UnCorr= PlotFormFactor.sh
series='better' PlotFormFactor.sh
#series='disp' PlotFormFactor.sh
# renorm=(mostly NPR) renormalised. No need for ZV
series='renorm AltZV' Suffix=_mostly ZV= PlotFormFactor.sh
# renorm=(mostly NPR) renormalised. Apply Fully NPR correction
series='renorm AltZV' FullyNP= PlotFormFactor.sh

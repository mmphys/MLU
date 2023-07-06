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

OutBase=$Ensemble/$MELFit/Fit_sp2_

############################################################

# Make the fit selection files I use

############################################################

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
  h385_s_p2_0  2ptp2/h385_s/h385_s_p2_0.corr_16_40_16_29.g5P_g5W.model.3285139232.h5
  # K
  s_l_p2_0  2ptp2/s_l/s_l_p2_0.corr_8_31_8_20.g5P_g5W.model.3285139232.h5
  # These still to be chosen (for now same as old)
  s_l_p2_1  2ptp2/s_l/s_l_p2_1.corr_8_28_9_28.g5P_g5W.model.3285139232.h5
  s_l_p2_2  2ptp2/s_l/s_l_p2_2.corr_8_26_8_28.g5P_g5W.model.3285139232.h5
  s_l_p2_3  2ptp2/s_l/s_l_p2_3.corr_8_28_8_28.g5P_g5W.model.3285139232.h5
  s_l_p2_4  2ptp2/s_l/s_l_p2_4.corr_10_28_8_28.g5P_g5W.model.3285139232.h5
  s_l_p2_5  2ptp2/s_l/s_l_p2_5.corr_10_28_12_28.g5P_g5W.model.3285139232.h5
  s_l_p2_6  2ptp2/s_l/s_l_p2_6.corr_10_28_12_28.g5P_g5W.model.3285139232.h5
EndFitChoices

############################################################

# Make the form factors

############################################################

series='old' PlotFormFactor.sh
series='old' FitSeries='std' UnCorr= PlotFormFactor.sh
series='better' PlotFormFactor.sh

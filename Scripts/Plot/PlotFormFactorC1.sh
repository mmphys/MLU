#!/usr/bin/env bash

export Ensemble=${Ensemble:-C1}
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
	h6413_s_p2_0  2ptp2/h6413_s/h6413_s_p2_0.corr_6_27_14_27.g5P_g5W.model.2263212701.h5
	# K - using dispersion relation
	s_l_p2_0  2ptp2/s_l/s_l.corr_6_23_6_23_5_20_5_20_5_18.g5P.model.2263212701.h5
	s_l_p2_1  2ptp2/s_l/s_l.corr_6_23_6_23_5_20_5_20_5_18.g5P.model.2263212701.h5
	s_l_p2_2  2ptp2/s_l/s_l.corr_6_23_6_23_5_20_5_20_5_18.g5P.model.2263212701.h5
	s_l_p2_3  2ptp2/s_l/s_l.corr_6_23_6_23_5_20_5_20_5_18.g5P.model.2263212701.h5
	s_l_p2_4  2ptp2/s_l/s_l.corr_6_23_6_23_5_20_5_20_5_18.g5P.model.2263212701.h5
EndFitChoices

[ -e ${OutBase}old.txt ] || cat > ${OutBase}old.txt <<- EndFitChoices
	# Ensemble: $Ensemble
	# This is the first set of fit choices I made on this ensemble
	# D_s
	h6413_s_p2_0  2ptp2/h6413_s/h6413_s_p2_0.corr_6_27_14_27.g5P_g5W.model.2263212701.h5
	# K
	s_l_p2_0  2ptp2/s_l/s_l_p2_0.corr_6_23_7_23.g5P_g5W.model.2263212701.h5
	s_l_p2_1  2ptp2/s_l/s_l_p2_1.corr_6_23_7_23.g5P_g5W.model.2263212701.h5
	s_l_p2_2  2ptp2/s_l/s_l_p2_2.corr_5_20_5_20.g5P_g5W.model.2263212701.h5
	s_l_p2_3  2ptp2/s_l/s_l_p2_3.corr_5_20_5_20.g5P_g5W.model.2263212701.h5
	s_l_p2_4  2ptp2/s_l/s_l_p2_4.corr_5_18_5_18.g5P_g5W.model.2263212701.h5
EndFitChoices

[ -e ${OutBase}priorP.txt ] || cat > ${OutBase}priorP.txt <<- EndFitChoices
	# Ensemble: $Ensemble
	# Fit n^2=0, then get overlap coefficients at each momentum enforcing dispersion
	# Not such a good idea
	# D_s
	h6413_s_p2_0  2ptp2/h6413_s/h6413_s_p2_0.corr_6_27_14_27.g5P_g5W.model.2263212701.h5
	# K
	# Each momentum fitted independently, but using dispersion relation for non-zero momenta
	s_l_p2_0  2ptp2/s_l/s_l_p2_0.corr_6_23_7_23.g5P_g5W.model.2263212701.h5
	s_l_p2_1  2ptp2/s_l/s_l_p2_1.priorP_6_23_7_23.corr_6_23.g5P_g5W.model.2263212701.h5
	s_l_p2_2  2ptp2/s_l/s_l_p2_2.priorP_6_23_7_23.corr_5_20.g5P_g5W.model.2263212701.h5
	s_l_p2_3  2ptp2/s_l/s_l_p2_3.priorP_6_23_7_23.corr_5_20.g5P_g5W.model.2263212701.h5
	s_l_p2_4  2ptp2/s_l/s_l_p2_4.priorP_6_23_7_23.corr_5_18.g5P_g5W.model.2263212701.h5
EndFitChoices

[ -e ${OutBase}priorPW.txt ] || cat > ${OutBase}priorPW.txt <<- EndFitChoices
	# Ensemble: $Ensemble
	# Fit n^2=0, then get overlap coefficients at each momentum enforcing dispersion
	# Not such a good idea
	# D_s
	h6413_s_p2_0  2ptp2/h6413_s/h6413_s_p2_0.corr_6_27_14_27.g5P_g5W.model.2263212701.h5
	# K
	# Each momentum fitted independently, but using dispersion relation for non-zero momenta
	s_l_p2_0  2ptp2/s_l/s_l_p2_0.corr_6_23_7_23.g5P_g5W.model.2263212701.h5
	s_l_p2_1  2ptp2/s_l/s_l_p2_1.priorPW_6_23_7_23.corr_6_23_7_23.g5P_g5W.model.2263212701.h5
	s_l_p2_2  2ptp2/s_l/s_l_p2_2.priorPW_6_23_7_23.corr_5_20_5_20.g5P_g5W.model.2263212701.h5
	s_l_p2_3  2ptp2/s_l/s_l_p2_3.priorPW_6_23_7_23.corr_5_20_5_20.g5P_g5W.model.2263212701.h5
	s_l_p2_4  2ptp2/s_l/s_l_p2_4.priorPW_6_23_7_23.corr_5_18_5_18.g5P_g5W.model.2263212701.h5
EndFitChoices

[ -e ${OutBase}betterP.txt ] || cat > ${OutBase}betterP.txt <<- EndFitChoices
	# Ensemble: $Ensemble
	# Better fit choices than priorP, but still not such a good idea
	# D_s
	h6413_s_p2_0  2ptp2/h6413_s/h6413_s_p2_0.corr_6_27_14_27.g5P_g5W.model.2263212701.h5
	# K
	# Each momentum fitted independently, but using dispersion relation for non-zero momenta
	s_l_p2_0  2ptp2/s_l/s_l_p2_0.corr_6_23_7_23.g5P_g5W.model.2263212701.h5
	s_l_p2_1  2ptp2/s_l/s_l_p2_1.betterP_6_23_7_23.corr_6_23.g5P_g5W.model.2263212701.h5
	s_l_p2_2  2ptp2/s_l/s_l_p2_2.betterP_6_23_7_23.corr_6_20.g5P_g5W.model.2263212701.h5
	s_l_p2_3  2ptp2/s_l/s_l_p2_3.betterP_6_23_7_23.corr_6_20.g5P_g5W.model.2263212701.h5
	s_l_p2_4  2ptp2/s_l/s_l_p2_4.betterP_6_23_7_23.corr_6_18.g5P_g5W.model.2263212701.h5
EndFitChoices

[ -e ${OutBase}betterPW.txt ] || cat > ${OutBase}betterPW.txt <<- EndFitChoices
	# Ensemble: $Ensemble
	# D_s
	# Better fit choices than priorPW, but still not such a good idea
	h6413_s_p2_0  2ptp2/h6413_s/h6413_s_p2_0.corr_6_27_14_27.g5P_g5W.model.2263212701.h5
	# K
	# Each momentum fitted independently, but using dispersion relation for non-zero momenta
	s_l_p2_0  2ptp2/s_l/s_l_p2_0.corr_6_23_7_23.g5P_g5W.model.2263212701.h5
	s_l_p2_1  2ptp2/s_l/s_l_p2_1.betterPW_6_23_7_23.corr_6_23_7_23.g5P_g5W.model.2263212701.h5
	s_l_p2_2  2ptp2/s_l/s_l_p2_2.betterPW_6_23_7_23.corr_6_20_5_20.g5P_g5W.model.2263212701.h5
	s_l_p2_3  2ptp2/s_l/s_l_p2_3.betterPW_6_23_7_23.corr_6_20_6_20.g5P_g5W.model.2263212701.h5
	s_l_p2_4  2ptp2/s_l/s_l_p2_4.betterPW_6_23_7_23.corr_6_18_6_18.g5P_g5W.model.2263212701.h5
EndFitChoices

for series in renorm renormold AltZV Jan24
do
  [ -e ${OutBase}$series.txt ] || ln -s ${OutBase}disp.txt ${OutBase}$series.txt
done
)

############################################################

# Make the form factors

############################################################

series='disp old priorP priorPW betterP betterPW' PlotFormFactor.sh
series=disp FitSeries=std UnCorr= PlotFormFactor.sh
# renorm=(mostly NPR) renormalised. No need for ZV
series='renormold renorm AltZV Jan24' Suffix=_mostly ZV= PlotFormFactor.sh
# renorm=(mostly NPR) renormalised. Apply Fully NPR correct.
series='renormold renorm AltZV Jan24' FullyNP= PlotFormFactor.sh

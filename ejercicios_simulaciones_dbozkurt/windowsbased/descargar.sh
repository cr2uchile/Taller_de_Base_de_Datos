#!/bin/bash
###Taller de Bases de Datos, CR2, Lunes 3 de Septiembre, 2018###
###Procesamiento y visualización del dato del modelo climático regional con NCL, CDO y NCO. Deniz Bozkurt###

#PARTE-1: Descargar datos del servidor#

export CR2_DATASITE=ftp://qhawayra2.dgf.uchile.cl/mma_simulations/Local/RegCM4
typeset -i a
a=1980
while [ "$a" -lt "2051" ]
do
if [ $a -lt 2006 ]; then
curl -o tas_CL-09_MPI-M-MPI-ESM-MR_historical_r1i1p1_CR2-RegCM4-6_v4_mon_${a}0101-${a}1231.nc ${CR2_DATASITE}/historical/mon/tas/CR2-RegCM4-6/MPI-M-MPI-ESM-MR/r1i1p1/tas_CL-09_MPI-M-MPI-ESM-MR_historical_r1i1p1_CR2-RegCM4-6_v4_mon_${a}0101-${a}1231.nc
fi
if [ $a -gt 2005 ]; then
curl -o tas_CL-09_MPI-M-MPI-ESM-MR_rcp85_r1i1p1_CR2-RegCM4-6_v4_mon_${a}0101-${a}1231.nc ${CR2_DATASITE}/rcp85/mon/tas/CR2-RegCM4-6/MPI-M-MPI-ESM-MR/r1i1p1/tas_CL-09_MPI-M-MPI-ESM-MR_rcp85_r1i1p1_CR2-RegCM4-6_v4_mon_${a}0101-${a}1231.nc
fi
a=a+1
done

mkdir data_cr2_curso
mv tas_*nc data_cr2_curso/


#!/bin/bash
###Taller de Bases de Datos, CR2, Lunes 3 de Septiembre, 2018###
###Procesamiento y visualización del dato del modelo climático regional con NCL, CDO y NCO. Deniz Bozkurt###

#PARTE-2: Procesamiento del dato#

cd data_cr2_curso/
echo ">> Copiar nco para ejecutarlo sobre data_cr2_curso"
cp -r /cygdrive/c/nco/* /cygdrive/c/ejercicios_simulaciones_dbozkurt/windowsbased/data_cr2_curso/
echo ">> Concatenar los datos"
cdo copy tas_CL-09*nc tas_1980_2050_CL-09.nc #Concatenar los datos para tener un solo archivo que consta de la temperatura mensual de 1980 a 2050
echo ">> Seleccionar un box y calcular temperatura anual promedio"
cdo yearavg -sellonlatbox,-80,-62,-57,-15 tas_1980_2050_CL-09.nc tas_1980_2050_CL-09_box_yearly.nc #Seleccionar un box y calcular temperatura anual promedio de cada año en este box
echo ">> Extraer series de tiempo para Santiago"
cdo -outputtab,year,lon,lat,value -remapnn,lon=-70.67_lat=-33.45 tas_1980_2050_CL-09_box_yearly.nc > santiago_yearly_temp.txt #Extraer series de tiempo de temperatura anualpara Santiago

echo ">> Regresión lineal"
cdo regres tas_1980_2050_CL-09_box_yearly.nc tas_1980_2050_CL-09_box_yearly_regres.nc #Regresión lineal
echo ">> Cambiar el nombre del variable"
./ncrename.exe -v tas,tastrend tas_1980_2050_CL-09_box_yearly_regres.nc #Cambiar el nombdre del variable

#PARTE-3: Visualización del dato#
cd ../
ncl tas_trend_cdo.ncl

ncl time_series_cdo.ncl

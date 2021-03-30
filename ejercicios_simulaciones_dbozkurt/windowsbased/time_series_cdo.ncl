begin


nrows = 71
ncols = 4
;Location of the files
data  = asciiread("data_cr2_curso/santiago_yearly_temp.txt",(/nrows,ncols/),"float")

;For time axis
year = data(:,0)
ndim = dimsizes(year)
syear = year(0)
eyear = year(ndim-1)

;Extract the data array
data_temp = data(:,3)-273.15
data_temp_runave = runave(data_temp,5,0)


;PLOT and DISPLAY ARRANGEMENTS
wks = gsn_open_wks ("png","temp_santiago")
res=True
res@gsnMaximize = True
res@tiYAxisFontHeightF = 0.03
res@tmYLLabelFontHeightF = 0.03
res@tmXBLabelAngleF = 45
res@tiXAxisFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.03
res@trXMinF = syear-.5
res@trXMaxF = eyear+1.5
res@vpWidthF              = 1.4
res@vpHeightF             = 0.6
res@tmXBMinorOn   = False 
res@tmYLMinorOn   = False
res@tmXTOn        = False 
res@tmYROn        = False 
res@tmXTBorderOn  = False 
res@tmYRBorderOn  = False 
res@tiYAxisString = "Temperature (~S~o~N~C)"

ldash   = (/0,1/)
lthick   = (/4,4/)
lcolors = (/"black","red"/)

res@xyDashPatterns = ldash
res@xyLineColors = lcolors
res@xyLineThicknesses = lthick
;print (year)
;print (eyear)
;print (ndim)
years = fspan(syear,eyear,ndim)
print (years)
plot = gsn_csm_xy(wks,years,(/data_temp,data_temp_runave/),res)
end

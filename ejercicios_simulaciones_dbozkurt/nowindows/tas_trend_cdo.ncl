begin

latsouth = -57
latnorth = -15
lonwest = -80
loneast = -62

;Location of the files
model = addfile("data_cr2_curso/tas_1980_2050_CL-09_box_yearly_regres.nc","r")

;Selection of variables
model_pre = model->tastrend(time|0,{lat|latnorth:latsouth},{lon|lonwest:loneast})
printVarSummary(model_pre)
model_pre = model_pre*10 ;(per decade)
title = "TASTREND_1980_2050"


;PLOT and DISPLAY ARRANGEMENTS
wks = gsn_open_wks ("ps",title)

gsn_define_colormap (wks,"MPL_YlOrRd")
;gsn_reverse_colormap (wks)

res                        = True               ; plot mods desired
res@cnFillOn               = True               ; color fill
res@cnLineLabelsOn      = False       ; Turn on contour line labels
res@cnLinesOn             = False

res@mpPerimOn = False
res@lbLabelBarOn         = False
res@mpFillOn = False
res@mpDataBaseVersion = "MediumRes"
res@gsnCenterString = ""
res@mpGridLineDashPattern  = 2                  ; lat/lon lines as dashed
res@pmTickMarkDisplayMode  = "Always"           ; turn on tickmarks
res@gsnSpreadColors = True
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF        =  0           ; set min contour level
res@cnMaxLevelValF        =  0.5            ; set max contour level
res@cnLevelSpacingF       =  0.05           ; set contour interval

res@mpLimitMode           = "LatLon"     ; required
res@mpMinLatF             = latsouth         ; min lat
res@mpMaxLatF             = latnorth         ; max lat
res@mpMinLonF             = lonwest        ; min lon
res@mpMaxLonF             = loneast        ; max lon

res@mpOutlineBoundarySets       = "National"
res@mpGeophysicalLineThicknessF = 0.5
res@mpNationalLineThicknessF    = 0.5

res@cnMissingValPerimOn     = True            ; turn on the missing value perimeter
res@cnMissingValPerimColor  = "grey80"           ; outline the perimeter red
res@cnMissingValFillPattern = 0               ; choose fill pattern 0 (solid fill) 
res@cnMissingValFillColor   = "grey80"          ; color fill missing areas blue
res@tmXBLabelFontHeightF  = 0.021
res@tmYLLabelFontHeightF  = 0.021
res@gsnRightString = ""
res@gsnLeftString = ""
res@gsnRightStringFontHeightF = 0.012
res@gsnCenterStringFontHeightF = 0.02
res@gsnLeftStringFontHeightF = 0.012
res@gsnDraw              = False              ; don't draw yet
res@gsnFrame             = False             ; don't advance frame tee
res@gsnAddCyclic           = False              ; regional data
;res@tfDoNDCOverlay = True
res@tiMainOn = True
res@tiXAxisOn = False
res@tiMainString = ""
;res@gsnMaximize = True


plots  = new(1,graphic)

plots(0) = gsn_csm_contour_map(wks,model_pre,res)

;---Create resources for paneling
resP                            = True
resP@gsnMaximize                = True
resP@gsnPanelLabelBar           = True
resP@txFontHeightF              = .20
resP@lbLabelFontHeightF         = 0.015
resP@pmLabelBarWidthF = 0.04
resP@pmLabelBarHeightF = 0.4
resP@lbOrientation = "Vertical"
resP@lbTitleOn        =  True                ; turn on title
resP@lbTitleString    = "~S~o~N~C/10yrs"                ; title string (Pa*100/100)
resP@lbTitlePosition  = "Bottom"              ; title position
resP@lbTitleFontHeightF= .012                ; make title smaller
resP@lbTitleDirection = "Across"             ; title direction
resP@lbLabelStride        = 2
resP@mpPerimOn = False

resP@pmLabelBarParallelPosF = 0.25
resP@pmLabelBarOrthogonalPosF = 0.0
;---Panel the plot.
resP@gsnPanelRowSpec = True
;resP@gsnPanelCenter  = False
gsn_panel(wks,plots,(/1,1/),resP)

end

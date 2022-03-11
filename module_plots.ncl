;======================================================================
; Author : Bappaditya
; This module contains plotting functions
;======================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./shapefile_utils.ncl"





;======================================================================
; Procedure to check the ROC score (xy plot)
;======================================================================
procedure plot_xy_var(x, y, legend, plottype, plotname)

;local 
begin

	wks	= gsn_open_wks (plottype,plotname)

	res							= True
	res@gsnMaximize 	        = True
	res@gsnFrame				= False
	res@gsnDraw					= False

	res@tiYAxisString			= "Hit Rate"
	res@tiXAxisString			= "False Alarm"

	res@xyLineThicknessF		= 15
	res@tiXAxisFontHeightF		= 0.035
	res@tiYAxisFontHeightF		= 0.035
	res@tmXBLabelFontHeightF	= 0.030
	res@tmYLLabelFontHeightF	= 0.030

	res@xyDashPatterns			= (/15,0,16/)
	res@xyLineColors			= (/"red","black","royalblue3"/)

	res@pmLegendDisplayMode		= "Always"
	res@lgPerimOn				= False		; Legend perimeter.
	res@pmLegendParallelPosF	= 0.6
	res@pmLegendOrthogonalPosF	= -0.5
	res@pmLegendWidthF			= 0.14		; Change width of legend
	res@pmLegendHeightF			= 0.14

	res@xyExplicitLegendLabels = (/"Bel. Norm.   : " + sprintf("%5.2f", legend(0)), \
                                   "Norm.          : " + sprintf("%5.2f", legend(1)), \
                                   "Abv. Norm. : " + sprintf("%5.2f", legend(2))/)

	res@lgLabelFont				= 25
	res@lgLabelFontHeightF		= 0.03	; set the legend label font thickness
	res@lgPerimFill            	= "HollowFill"		; Fill legend box w/white
	res@lgPerimFillColor       	= "white"		; so it masks XY curves

	plot	= gsn_csm_xy (wks,x,y,res)

	lnres						= True
	lnres@gsLineDashPattern		= 2
	lnres@gsLineThicknessF		= 5
	gsn_polyline(wks,plot,fspan(0,1,2),fspan(0,1,2),lnres)

	draw(plot)
	frame(wks)

	if (plottype .eq. "png")
		system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
	end if

;======================================================================
end
;======================================================================





;======================================================================
; Procedure to check the spatial pattern
; Check for the Indian region by plotting
;======================================================================

procedure plot_2D_var(var,plottype,region)

local min_lat, max_lat, min_lon, max_lon
begin

;======================================================================
; Specifications
;======================================================================

	if (region .eq. "Tropics")
		min_lat 		= -30
		max_lat 		=  30
		min_lon 		=  20
		max_lon 		= 300
	end if
	if (region .eq. "SouthAsia")
		min_lat 		= -30
		max_lat 		=  40
		min_lon 		=  30
		max_lon 		= 120
	end if
	if (region .eq. "India")
		min_lat 		=   5
		max_lat 		=  40
		min_lon 		=  65
		max_lon 		= 100
	end if

;======================================================================
; Resources for plotting
;======================================================================

	plotname	= "fig_"+var@long_name
	wks		= gsn_open_wks(plottype, plotname)

	res	= True
	res@gsnMaximize 	        = True

	res@mpFillOn				= False			; turn off gray fill over land
	if (region .eq. "World")
		res@tmXBTickSpacingF			= 60.
		res@gsnAddCyclic				= True
		res@gsnStringFontHeightF		= 0.025
		res@tmXBLabelFontHeightF		= 0.025
		res@tmYLLabelFontHeightF		= 0.025
		res@lbLabelFontHeightF			= 0.018
		res@lbTitleFontHeightF			= 0.020
		res@pmLabelBarOrthogonalPosF	= 0.15
	end if
	if (region .ne. "World")
		res@mpMinLatF					= min_lat
		res@mpMaxLatF					= max_lat
		res@mpMinLonF					= min_lon
		res@mpMaxLonF					= max_lon
		res@gsnAddCyclic				= False
	end if
	if (region .eq. "SouthAsia")
		res@gsnStringFontHeightF		= 0.025
		res@tmXBLabelFontHeightF		= 0.025
		res@tmYLLabelFontHeightF		= 0.025
		res@lbLabelFontHeightF			= 0.018
		res@tmYLTickSpacingF			= 20.
		res@pmLabelBarOrthogonalPosF	= 0.15
	end if
	if (region .eq. "Tropics")
		res@gsnStringFontHeightF		= 0.025
		res@tmXBLabelFontHeightF		= 0.018
		res@tmYLLabelFontHeightF		= 0.018
		res@lbLabelFontHeightF			= 0.018
		res@lbTitleFontHeightF			= 0.018
		res@pmLabelBarOrthogonalPosF	= 0.30
	end if
	if (region .eq. "India")
		res@mpShapeMode					= "FreeAspect"
		res@vpWidthF					= 0.6
		res@vpHeightF					= 0.4
		res@tmXBTickSpacingF			= 30.
		res@tmYLTickSpacingF			= 15.
		res@gsnStringFontHeightF		= 0.023
		res@tmXBLabelFontHeightF		= 0.023
		res@tmYLLabelFontHeightF		= 0.023
		res@lbLabelFontHeightF			= 0.018
		res@lbTitleFontHeightF			= 0.018
		res@pmLabelBarOrthogonalPosF	= 0.15
	end if

	res@gsnSpreadColors			= True		; use full range of colors

	res@cnLevelSelectionMode	= "ManualLevels"
;	res@cnMinLevelValF			=  -0.1
;	res@cnMaxLevelValF			=   0.8
;	res@cnLevelSpacingF			=   0.1
	res@cnMinLevelValF			=  0.1
	res@cnMaxLevelValF			=  0.9
	res@cnLevelSpacingF			=  0.1
;	res@cnLevelSelectionMode	= "ExplicitLevels"
;	res@cnLevels				= (/ -1,-0.5,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.5,1 /)   ; set levels

	res@tiMainString			= ""
	res@gsnLeftString			= ""
	res@gsnRightString			= ""

	if (region .eq. "World")
		res@tmXBTickSpacingF		= 60.
	end if
	if (region .eq. "India")
		res@tmXBTickSpacingF		= 5.
		res@tmYLTickSpacingF		= 5.
	end if
	if (region .eq. "SouthAsia")
		res@tmYLTickSpacingF		= 20.
	end if

;	gsn_define_colormap(wks,"BlueWhiteOrangeRed")
	gsn_define_colormap(wks,"BlueYellowRed")
;	gsn_define_colormap(wks,"sunshine_9lev")
;	gsn_define_colormap(wks,"BlAqGrWh2YeOrReVi22")
;	gsn_define_colormap(wks,"hotres")
;;	gsn_define_colormap(wks,"NCV_blu_red")
; Flux and Temp Related
;	gsn_define_colormap(wks,"MPL_hot")
;	gsn_define_colormap(wks,"BlueWhiteOrangeRed")
;	gsn_define_colormap(wks,"BlueDarkOrange18")
;	gsn_define_colormap(wks,"temp_19lev")
; Precip Related
;;	gsn_define_colormap(wks,"precip4_diff_19lev")
;	gsn_define_colormap(wks,"precip_diff_12lev")
;	gsn_define_colormap(wks,"precip_11lev")
;	gsn_define_colormap(wks,"CBR_drywet")
;	gsn_define_colormap(wks,"MPL_BrBG")
;	gsn_define_colormap(wks,"MPL_gist_yarg")

;	gsn_reverse_colormap(wks)

; Plotting options for Contours
	res@cnFillOn								= True		;turn on color
	res@cnLinesOn								= False
	res@cnLineLabelsOn							= False
	res@gsnContourLineThicknessesScale			= 0.1
	res@gsnContourZeroLineThicknessF			= 10
; To switch on Raster Mode, set cnLinesOn = False, else set True
	res@cnFillMode								= "RasterFill"

	res@mpOutlineOn								= True		; turn on map outlines
	res@mpGeophysicalLineColor					= "black"
	if (plottype .eq. "eps")
		res@mpGeophysicalLineThicknessF			= 2
	end if
	if (plottype .eq. "png")
		res@mpGeophysicalLineThicknessF			= 4
	end if

	if(isdefined(var@units)) then
		lbTitleString							= var@units
		res@lbTitleOn							= True
		res@lbTitleFont							= 25
		res@lbTitlePosition						= "Right"
		res@lbTitleDirection					= "Across"
		res@lbTitleString						= lbTitleString
		res@lbBoxMinorExtentF					= 0.2
	end if

;	symMinMaxPlt (var,20,False,res)
	contour		= gsn_csm_contour_map_ce(wks,var,res)

	if (plottype .eq. "png")
		system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
	end if

;======================================================================
end
;======================================================================





;======================================================================
; Procedure to plot the spatial pattern of a variable from a file
;======================================================================

procedure plot_2D_file(filerd,plottype,region)

local filein, vNames, nNames
begin

	filein			= addfile(filerd,"r")

	vNames	= getfilevarnames (filein)
	nNames	= dimsizes (vNames)
	print (vNames)
	var		= filein->$vNames(2)$

	plot_2D_var(var,plottype,region)

;======================================================================
end
;======================================================================


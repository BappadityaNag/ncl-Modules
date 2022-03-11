;======================================================================
; Author : Bappaditya
; This module contains functions to compute statistics (RPSS/RPSSd)
; Files should have data formatted as  (year|29, cat|3)
;======================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./shapefile_utils.ncl"





;======================================================================
; Function to calculate the RPSS
;======================================================================

procedure score_rpss_2D( modelfile, obvsfile, stat, no_of_ens, fileout )

; local 
begin

;======================================================================
; Read Data
;======================================================================

	filein		= addfile(modelfile,"r")
	vNames		= getfilevarnames(filein)
;	print(vNames)
	var			= filein->$vNames(4)$

	filein		:= addfile(obvsfile,"r")
	varo		= filein->$vNames(4)$

	ncat		= dimsizes(var&cat)

;======================================================================
; Calculation of RPS Forecast
;======================================================================

; Dimension of terciles should be as follows
;	terc	= (cat,time,lat,lon)

	var		:= var(cat|:,time|:,lat|:,lon|:)
	varo	:= varo(cat|:,time|:,lat|:,lon|:)

	pfcst			= var
	obvs			= varo
	pclim			= pfcst			; Casting pclim in the mould of pfcst
	pclim(0,:,:,:)	= 0.33333
	pclim(1,:,:,:)	= 0.66667
	pclim(2,:,:,:)	= 1.00000
	
	do icat = 1,ncat-1
		pfcst(icat,:,:,:)	= pfcst(icat-1,:,:,:) + var(icat,:,:,:)
		obvs(icat,:,:,:)	= obvs(icat-1,:,:,:) + varo(icat,:,:,:)
	end do

;	printVarSummary(pfcst)
;	printMinMax(pfcst,1)
;	printVarSummary(obvs)
;	printMinMax(obvs,1)

	sq_diff_fcst	= ( pfcst - obvs )^2
	sq_diff_clim	= ( pclim - obvs )^2
	copy_VarMeta(pfcst,sq_diff_fcst)
	copy_VarMeta(pclim,sq_diff_clim)

	rps_fcst	= dim_sum_n_Wrap( sq_diff_fcst,0 )
	rps_clim	= dim_sum_n_Wrap( sq_diff_clim,0 )
	copy_VarMeta(sq_diff_fcst(0,:,:,:),rps_fcst)
	copy_VarMeta(sq_diff_clim(0,:,:,:),rps_clim)

;======================================================================
; RPSS_2D
;======================================================================

	D	= 1.0/toint(no_of_ens) * (ncat^2 - 1) / 6*ncat
	if (stat .eq. "RPSS")
		score	= 1.0 - dim_avg_n(rps_fcst,0) / ( dim_avg_n(rps_clim,0) + D )
		copy_VarMeta(rps_fcst(0,:,:),score)
	end if

;	score	= 1.0 - dim_avg_n(rps_fcst,0) / dim_avg_n(rps_clim,0)
;	copy_VarMeta(rps_fcst(0,:,:),score)

;======================================================================
; RPSSd_2D
;======================================================================

	if (stat .eq. "RPSSd")

		iterN			= 100000
		dimvar			:= dimsizes(var)
		ntime			= dimvar(1)
		rps_clim_exp	= new( (/iterN, dimvar(2), dimvar(3)/), typeof(rps_clim), rps_clim@_FillValue )

		do iter = 0, iterN-1
			print( "Bootstrapping : " + iter )
			iw 			= generate_sample_indices( ntime, 1 )
			rps_clim_exp(iter,:,:)	= dim_avg_n_Wrap(rps_clim(iw,:,:),0)
		end do

		score	= 1.0 - dim_avg_n(rps_fcst,0) / dim_avg_n(rps_clim_exp,0)
		copy_VarMeta(rps_fcst(0,:,:),score)

	end if

;======================================================================
; Writing to file
;======================================================================
	score@long_name	= stat

	if (fileexists(fileout)) then
		print("Deleting "+fileout)
		system("rm "+fileout)
	end if
	write_2D( fileout, score )

;======================================================================
end
;======================================================================





;======================================================================
; Function to calculate the ROC
;======================================================================
procedure score_roc_2D( modelfile, obvsfile, stat, fileout )

; local 
begin

	nbin			= 10
;======================================================================
; Read Data
;======================================================================

	filein		= addfile(modelfile,"r")
	vNames		= getfilevarnames(filein)
	var			= filein->$vNames(4)$

	filein		:= addfile(obvsfile,"r")
	varo		= filein->$vNames(4)$

	ncat		= dimsizes(var&cat)
	nyrs		= dimsizes(var&time)
	nlat		= dimsizes(var&lat)
	nlon		= dimsizes(var&lon)

;======================================================================
; Calculate the ROC
;======================================================================

; var ( cat, lat, lon )

; Array Size is taken as nbin to include (0,0) point
	hr				= new((/ncat,nbin+1,nlat,nlon/),float,var@_FillValue)
	fa 				= new((/ncat,nbin+1,nlat,nlon/),float,var@_FillValue)
	tot_hits		= new((/ncat,nlat,nlon/),float,var@_FillValue)
	tot_falalrm		= new((/ncat,nlat,nlon/),float,var@_FillValue)
	hr(:,:,:,:)			= 0
	fa(:,:,:,:)			= 0
	tot_hits(:,:,:)		= 0
	tot_falalrm(:,:,:)	= 0

	binsize = 1.0/nbin
	do bin = nbin-1, 0, 1
		prob		= (tofloat(bin))/(nbin)
;		print( "                     " )
		print( "Probability >= " + prob + " and Probability < " + prob+"+"+binsize )
		do cat	= 0, ncat-1
;			print( "Cat = " + cat )
			do t = 0, nyrs-1
				print( "Cat = " + cat + " ... Year : " + var&time(t) )
				do ilat = 0, nlat-1
;					print( "Latitude : " + var&lat(ilat) )
					do ilon = 0, nlon-1
						if ( .not.ismissing(var(cat,t,ilat,ilon)) .and. .not.ismissing(varo(cat,t,ilat,ilon)) )
							if ( var(cat,t,ilat,ilon) .ge. prob .and. var(cat,t,ilat,ilon) .lt. prob+binsize )
;									print( "Condition Satisfied : Year = " + year(t) )
								if ( varo(cat,t,ilat,ilon) .eq. 1 )
;									print( "Hit Occurred : Year = " + year(t) )
									hr(cat,bin,ilat,ilon)		= hr(cat,bin,ilat,ilon) + 1
									tot_hits(cat,ilat,ilon)		= tot_hits(cat,ilat,ilon) + 1
								else
									fa(cat,bin,ilat,ilon)		= fa(cat,bin,ilat,ilon) + 1
									tot_falalrm(cat,ilat,ilon)	= tot_falalrm(cat,ilat,ilon) + 1
								end if
							end if
						end if
					end do
				end do
			end do
;			print( "Hits = " + hr(cat,bin,ilat,ilon) )
;			print( "False Alarm = " + fa(cat,bin,ilat,ilon) )
		end do
	end do

	do cat = 0, ncat-1
		do bin = nbin-2, 0, 1
			hr(cat,bin,:,:)	= hr(cat,bin+1,:,:) + hr(cat,bin,:,:)
			fa(cat,bin,:,:)	= fa(cat,bin+1,:,:) + fa(cat,bin,:,:)
		end do
	end do

;	print( "                     " )
;	print( "Total Hits = " + tot_hits(cat) )
;	print( "Total False Alarm = " + tot_falalrm(cat) )

	do cat = 0, ncat-1
		do ilat = 0, nlat-1
			do ilon = 0, nlon-1
				if( tot_hits(cat,ilat,ilon) .eq. 0 .and. tot_falalrm(cat,ilat,ilon) .eq. 0 ) then
					hr(cat,:,ilat,ilon)	= var@_FillValue
					fa(cat,:,ilat,ilon)	= var@_FillValue
				else
					if( tot_hits(cat,ilat,ilon) .eq. 0 ) then
						hr(cat,:,ilat,ilon)	= 0
					else if ( tot_falalrm(cat,ilat,ilon) .eq. 0) then
						fa(cat,:,ilat,ilon)	= 0
					else
						hr(cat,:,ilat,ilon)	= hr(cat,:,ilat,ilon)/tot_hits(cat,ilat,ilon)
						fa(cat,:,ilat,ilon)	= fa(cat,:,ilat,ilon)/tot_falalrm(cat,ilat,ilon)
					end if
					end if
				end if
			end do
		end do
	end do

;	print( "                     " )
;	print( "Hit Rate = " + hr(cat,::-1) + "            False Alarm = " + fa(cat,::-1))

;======================================================================
; Area under ROC
;======================================================================

	do cat = 0, ncat-1
		hr(cat,:,:,:)		= hr(cat,::-1,:,:)
		fa(cat,:,:,:)		= fa(cat,::-1,:,:)
	end do

;	print( "                     " )
	score		= new((/ncat,nlat,nlon/),float,var@_FillValue)
	score(:,:,:)		= 0.0
	copy_VarMeta(var(:,0,:,:),score)
	score@long_name	= stat

; Upper limit is exceeded by 1 to account for (1,1) point which is pushed left for inclusion of (0,0)
	do cat = 0, ncat-1
		do bin = 0, nbin-1
			score(cat,:,:) = score(cat,:,:) + 0.5 * ( hr(cat,bin,:,:)+hr(cat,bin+1,:,:) ) * ( fa(cat,bin+1,:,:)-fa(cat,bin,:,:) )
		end do
		score(cat,:,:) = where( ismissing(var(0,0,:,:)), var@_FillValue, score(cat,:,:) )
	end do

;======================================================================
; Writing to file
;======================================================================

	if (fileexists(fileout)) then
		print("Deleting "+fileout)
		system("rm "+fileout)
	end if
	write_3D( fileout, score )

;======================================================================
end
;======================================================================





;======================================================================
; Function to calculate the RPSS
;======================================================================
function score_rpss( modelfile, obvsfile, stat, no_of_ens )

; local 
begin

;======================================================================
; Read Data for Observations
;======================================================================

	ncat		= dimsizes(str_split(readAsciiHead(obvsfile, 1), " "))-1
	varo			= readAsciiTable(obvsfile, ncat+1, "float", 1)
	year			= toint(varo(:,0))
	varo			:= varo(:,1:)
	varo@_FillValue	= -99.900
	varo!0			= "year"
	varo&year		= year
	varo!1			= "cat"
	varo&cat		= ispan(1,ncat,1)
;	printVarSummary(varo)

;======================================================================
; Read Data for Model
;======================================================================

	var				= readAsciiTable(modelfile, ncat+1, "float", 1)
	year			= toint(var(:,0))
	var				:= var(:,1:)
	var@_FillValue	= -99.900
	var!0			= "year"
	var!1			= "cat"
	var&year		= year
	var&cat			= ispan(1,ncat,1)
;	printVarSummary(var)

;======================================================================
; Calculation of RPS Forecast
;======================================================================

; Dimension of terciles should be as follows
;	terc	= (cat,year)

	var		:= var(cat|:,year|:)
	varo	:= varo(cat|:,year|:)

	pfcst	= var
	obvs	= varo
	pclim	= pfcst			; Casting pclim in the mould of pfcst
	pclim(0,:)	= 0.33333
	pclim(1,:)	= 0.66667
	pclim(2,:)	= 1.00000
	
	do icat = 1,ncat-1
		pfcst(icat,:)	= pfcst(icat-1,:) + var(icat,:)
		obvs(icat,:)	= obvs(icat-1,:) + varo(icat,:)
	end do

;	print ("The Cumulative Predicted Model Values are : ")
;	print(year + " : " + pfcst(0,:) + " " + pfcst(1,:) + " " + pfcst(2,:))
;	print ("The Observed Values are : ")
;	print(year + " : " + obvs(0,:) + " " + obvs(1,:) + " " + obvs(2,:))

	sq_diff_fcst	= ( pfcst - obvs )^2
	sq_diff_clim	= ( pclim - obvs )^2
	copy_VarMeta(pfcst,sq_diff_fcst)
	copy_VarMeta(pclim,sq_diff_clim)

	rps_fcst	= dim_sum_n_Wrap( sq_diff_fcst,0 )
	rps_clim	= dim_sum_n_Wrap( sq_diff_clim,0 )
	copy_VarMeta(sq_diff_fcst(0,:),rps_fcst)
	copy_VarMeta(sq_diff_clim(0,:),rps_clim)

;======================================================================
; RPSS
;======================================================================

	D	= 1.0/toint(no_of_ens) * (ncat^2 - 1) / 6*ncat
	if (stat .eq. "rpss")
		return( 1.0 - dim_avg(rps_fcst) / ( dim_avg(rps_clim) + D ) )
	end if
;	return( 1.0 - dim_avg(rps_fcst) / dim_avg(rps_clim) )

;======================================================================
; RPSSd
;======================================================================

	if (stat .eq. "rpssd")

		iterN			= 100000
		dimvar := dimsizes(var)
		ntime	= dimvar(1)
		rps_clim_exp	= new( iterN, typeof(rps_clim), rps_clim@_FillValue )

		do iter = 0, iterN-1
;			print( "Iterations : " + iter )
			iw 			= generate_sample_indices( ntime, 1 )
			rps_clim_exp(iter)	= dim_avg(rps_clim(iw))
		end do

		return( 1.0 - dim_avg(rps_fcst) / dim_avg(rps_clim_exp) )

	end if

;======================================================================
end
;======================================================================





;======================================================================
; Function to calculate the ROC
;======================================================================
function score_roc( modelfile, obvsfile )

; local 
begin

	nbin			= 10
	ncat		= dimsizes(str_split(readAsciiHead(obvsfile, 1), " "))-1
	var				= readAsciiTable( modelfile, ncat+1, "float", 1 )
	varo			= readAsciiTable( obvsfile, ncat+1, "float", 1)
	year			= toint(var(:,0))
	var				:= var(:,1:3)
	varo			:= varo(:,1:3)
	nyrs			= dimsizes(var(:,0))

;======================================================================
; Calculate the ROC
;======================================================================

; var ( cat1, cat2, cat3 )

; Array Size is taken as nbin to include (0,0) point
	hr				= new((/ncat,nbin+1/),float,var@_FillValue)
	fa 				= new((/ncat,nbin+1/),float,var@_FillValue)
	tot_hits		= new(ncat,float,var@_FillValue)
	tot_falalrm		= new(ncat,float,var@_FillValue)
	hr(:,:)			= 0
	fa(:,:)			= 0
	tot_hits(:)		= 0
	tot_falalrm(:)	= 0

	binsize = 1.0/nbin
	do bin = nbin-1, 0, 1
		prob		= (tofloat(bin))/(nbin)
;		print( "                     " )
;		print( "Probability >= " + prob + " and Probability < " + prob+"+"+binsize )
		do cat	= 0, ncat-1
;			print( "Cat = " + cat )
			do t = 0, nyrs-1
				if ( var(t,cat) .ge. prob .and. var(t,cat) .lt. prob+binsize )
;					print( "Condition Satisfied : Year = " + year(t) )
					if ( varo(t,cat) .eq. 1 )
;						print( "Hit Occurred : Year = " + year(t) )
						hr(cat,bin)			= hr(cat,bin) + 1
						tot_hits(cat)		= tot_hits(cat) + 1
					else
						fa(cat,bin)			= fa(cat,bin) + 1
						tot_falalrm(cat)	= tot_falalrm(cat) + 1
					end if
				end if
			end do
;			print( "Hits = " + hr(cat,bin) )
;			print( "False Alarm = " + fa(cat,bin) )
		end do
	end do

	do bin = nbin-2, 0, 1
		hr(:,bin)	= hr(:,bin+1) + hr(:,bin)
		fa(:,bin)	= fa(:,bin+1) + fa(:,bin)
	end do

;	print( "                     " )
;	print( "Total Hits = " + tot_hits(cat) )
;	print( "Total False Alarm = " + tot_falalrm(cat) )

	do cat = 0, ncat-1
		hr(cat,:)	= hr(cat,:)/tot_hits(cat)
		fa(cat,:)	= fa(cat,:)/tot_falalrm(cat)
	end do

;	print( "                     " )
;	print( "Hit Rate = " + hr(cat,::-1) + "            False Alarm = " + fa(cat,::-1))

;======================================================================
; Area under ROC
;======================================================================

	do cat = 0, ncat-1
		hr(cat,:)		= hr(cat,::-1)
		fa(cat,:)		= fa(cat,::-1)
	end do

;	print( "                     " )
	area		= new(ncat,float,var@_FillValue)
	area(:)		= 0.0
; Upper limit is exceeded by 1 to account for (1,1) point which is pushed left for inclusion of (0,0)
	do cat = 0, ncat-1
		do bin = 0, nbin-1
			area(cat) = area(cat) + 0.5 * ( hr(cat,bin)+hr(cat,bin+1) ) * ( fa(cat,bin+1)-fa(cat,bin) )
		end do
;		print( "Cat = " + cat + " : Area under ROC = " + area(cat) )
	end do

;======================================================================
; Plot
;======================================================================

	plot_xy_var( fa, hr, area, "eps", "fig_ROC" )
	plot_xy_var( fa, hr, area, "png", "fig_ROC" )

;======================================================================
	return (area)
end
;======================================================================





;======================================================================
; Function to calculate the acc
;======================================================================
function acc( modelfile, obvsfile )

; local 
begin

;======================================================================
; Read Data for Observations
;======================================================================

	varo			= readAsciiTable(obvsfile, 2, "float", 1)
	year			= toint(varo(:,0))
	varo			:= varo(:,1)
	varo@_FillValue	= -999.90
	varo!0			= "year"
	varo&year		= year
;	printVarSummary(varo)

;======================================================================
; Read Data for Model
;======================================================================

	nens	= dimsizes(str_split(readAsciiHead(modelfile, 1), " "))-1
	var				= readAsciiTable (modelfile, nens+1, "float", 1)
	year			:= toint(var(:,0))
	var				:= var(:,1:)
	var@_FillValue	= -999.90
	var!0			= "year"
	var!1			= "ens"
	var&year		= year
	var&ens			= ispan(0,nens-1,1)
;	printVarSummary(var)

; Detrend
;======================================================================
	var_dtrend	= dtrend_msg_n(var&year,var,True,False,0)
	copy_VarMeta(var,var_dtrend)
	var		:= var_dtrend
	delete(var_dtrend) 
	varo_dtrend	= dtrend_msg_n(varo&year,varo,True,False,0)
	copy_VarMeta(varo,varo_dtrend)
	varo		:= varo_dtrend
	delete(varo_dtrend) 

; Enemble mean of Model Data
;======================================================================
	var		:= dim_avg_n_Wrap(var,1)

; Calculation of Anomaly Correlation Coefficient
;======================================================================
	return( escorc(var,varo) )

end
;======================================================================


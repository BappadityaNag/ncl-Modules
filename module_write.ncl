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
; Procedure to write a var with 2 dimensions (dim0, dim1)
;======================================================================

procedure write_2D(fileout,var)
;======================================================================

local fout, dNames, dsize, dim0, dim1, dimNames, dimSizes, dimUnlim
begin

;======================================================================
; Variable var have associated meta data: 
; (a) coordinate variables time, lat, lon and
; (b) attributes
;======================================================================

	fout		= addfile (fileout, "c")			; open output file
	dNames		= getvardimnames(var)
	dsize		= dimsizes(var)
	dim0		= var&$var!0$
	dim1		= var&$var!1$

;======================================================================
; explicitly declare file definition mode. Improve efficiency.
;======================================================================

    setfileoption(fout,"DefineMode",True)

;======================================================================
; create global attributes of the file
;======================================================================

	fAtt				= True				; assign file attributes
	fAtt@title			= ""
	fAtt@source_file	= ""
	fAtt@Conventions	= ""
	fAtt@creation_date	= systemfunc("date")
	fileattdef( fout, fAtt )				; copy file attributes

;======================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;======================================================================

	dimNames	= (/ dNames(0),  dNames(1) /)
	dimSizes	= (/ dsize(0),   dsize(1) /) 
	dimUnlim	= (/ False,  False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

;======================================================================
; predefine the dimensionality of the variables to be written out
;======================================================================
; Here we are using NCL functions to facilitate defining
; each variable's dimension name(s) and type.
; The following could be replaced with explicit, user defined dimension
; names different from those associated with the variable in memory.
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/))
;======================================================================

	filevardef(fout, dNames(0),     typeof(dim0), getvardims(dim0))
	filevardef(fout, dNames(1),     typeof(dim1), getvardims(dim1))
	filevardef(fout, var@long_name, typeof(var),  getvardims(var))
; variable name on the file 

; different from name on script
;======================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;======================================================================

	filevarattdef(fout, var@long_name, var)				; copy var attributes
	filevarattdef(fout, dNames(0), dim0)				; copy dim0 attributes
	filevarattdef(fout, dNames(1), dim1)				; copy dim1 attributes

;======================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;======================================================================

	setfileoption(fout,"DefineMode",False)

;======================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;======================================================================

	fout->$var!0$			= (/dim0/)
	fout->$var!1$			= (/dim1/)
	fout->$var@long_name$	= (/var/)

;======================================================================
end
;======================================================================






;======================================================================
; Procedure to write a var with 3 dimensions (dim0, dim1, dim2)
;======================================================================

procedure write_3D(fileout,var)
;======================================================================

local fout, dNames, dsize, dim0, dim1, dim2, dimNames, dimSizes, dimUnlim
begin

;======================================================================
; Variable var have associated meta data: 
; (a) coordinate variables time, lat, lon and
; (b) attributes
;======================================================================

	fout		= addfile (fileout, "c")			; open output file
	dNames		= getvardimnames(var)
	dsize		= dimsizes(var)
	dim0		= var&$var!0$
	dim1		= var&$var!1$
	dim2		= var&$var!2$

;======================================================================
; explicitly declare file definition mode. Improve efficiency.
;======================================================================

    setfileoption(fout,"DefineMode",True)

;======================================================================
; create global attributes of the file
;======================================================================

	fAtt				= True				; assign file attributes
	fAtt@title			= ""
	fAtt@source_file	= ""
	fAtt@Conventions	= ""
	fAtt@creation_date	= systemfunc("date")
	fileattdef( fout, fAtt )				; copy file attributes

;======================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;======================================================================

	dimNames	= (/ dNames(0),  dNames(1), dNames(2) /)
	dimSizes	= (/ dsize(0),   dsize(1),  dsize(2)  /) 
	dimUnlim	= (/ False,  False, False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

;======================================================================
; predefine the dimensionality of the variables to be written out
;======================================================================
; Here we are using NCL functions to facilitate defining
; each variable's dimension name(s) and type.
; The following could be replaced with explicit, user defined dimension
; names different from those associated with the variable in memory.
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/))
;======================================================================

	filevardef(fout, dNames(0),     typeof(dim0), getvardims(dim0))
	filevardef(fout, dNames(1),     typeof(dim1), getvardims(dim1))
	filevardef(fout, dNames(2),     typeof(dim2), getvardims(dim2))
	filevardef(fout, var@long_name, typeof(var),  getvardims(var))
; variable name on the file 

; different from name on script
;======================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;======================================================================

	filevarattdef(fout, var@long_name, var)				; copy var attributes
	filevarattdef(fout, dNames(0), dim0)				; copy dim0 attributes
	filevarattdef(fout, dNames(1), dim1)				; copy dim1 attributes
	filevarattdef(fout, dNames(2), dim2)				; copy dim2 attributes

;======================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;======================================================================

	setfileoption(fout,"DefineMode",False)

;======================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;======================================================================

	fout->$var!0$			= (/dim0/)
	fout->$var!1$			= (/dim1/)
	fout->$var!2$			= (/dim2/)
	fout->$var@long_name$	= (/var/)

;======================================================================
end
;======================================================================






;======================================================================
; Procedure to write a var with 4 dimensions (dim0, dim1, dim2, dim3)
;======================================================================

procedure write_4D(fileout,var)
;======================================================================

local fout, dNames, dsize, dim0, dim1, dim2, dim3, dimNames, dimSizes, \
      dimUnlim
begin

;======================================================================
; Variable var have associated meta data: 
; (a) coordinate variables time, lat, lon and
; (b) attributes
;======================================================================

	fout		= addfile (fileout, "c")			; open output file
	dNames		= getvardimnames(var)
	dsize		= dimsizes(var)
	dim0		= var&$var!0$
	dim1		= var&$var!1$
	dim2		= var&$var!2$
	dim3		= var&$var!3$

;======================================================================
; explicitly declare file definition mode. Improve efficiency.
;======================================================================

    setfileoption(fout,"DefineMode",True)

;======================================================================
; create global attributes of the file
;======================================================================

	fAtt				= True				; assign file attributes
	fAtt@title			= ""
	fAtt@source_file	= ""
	fAtt@Conventions	= ""
	fAtt@creation_date	= systemfunc("date")
	fileattdef( fout, fAtt )				; copy file attributes

;======================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;======================================================================

	dimNames	= (/ dNames(0),  dNames(1), dNames(2), dNames(3) /)
	dimSizes	= (/ dsize(0),   dsize(1),  dsize(2),  dsize(3)  /) 
	dimUnlim	= (/ False, False, False, False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

;======================================================================
; predefine the dimensionality of the variables to be written out
;======================================================================
; Here we are using NCL functions to facilitate defining
; each variable's dimension name(s) and type.
; The following could be replaced with explicit, user defined dimension
; names different from those associated with the variable in memory.
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/))
;======================================================================

	filevardef(fout, dNames(0),     typeof(dim0), getvardims(dim0))
	filevardef(fout, dNames(1),     typeof(dim1), getvardims(dim1))
	filevardef(fout, dNames(2),     typeof(dim2), getvardims(dim2))
	filevardef(fout, dNames(3),     typeof(dim3), getvardims(dim3))
	filevardef(fout, var@long_name, typeof(var),  getvardims(var))
; variable name on the file 

; different from name on script
;======================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;======================================================================

	filevarattdef(fout, var@long_name, var)				; copy var attributes
	filevarattdef(fout, dNames(0), dim0)				; copy dim0 attributes
	filevarattdef(fout, dNames(1), dim1)				; copy dim1 attributes
	filevarattdef(fout, dNames(2), dim2)				; copy dim2 attributes
	filevarattdef(fout, dNames(3), dim3)				; copy dim3 attributes

;======================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;======================================================================

	setfileoption(fout,"DefineMode",False)

;======================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;======================================================================

	fout->$var!0$			= (/dim0/)
	fout->$var!1$			= (/dim1/)
	fout->$var!2$			= (/dim2/)
	fout->$var!3$			= (/dim3/)
	fout->$var@long_name$	= (/var/)

;======================================================================
end
;======================================================================






;======================================================================
; Procedure to write a var with 5 dimensions (dim0, dim1, dim2, dim3, dim4)
;======================================================================

procedure write_5D( fileout, var )
;======================================================================

local fout, dNames, dsize, dim0, dim1, dim2, dim3, dim4, dimNames, dimSizes, \
      dimUnlim
begin

;======================================================================
; Variable var have associated meta data: 
; (a) coordinate variables time, lat, lon and
; (b) attributes
;======================================================================

	fout		= addfile (fileout, "c")			; open output file
	dNames		= getvardimnames(var)
	dsize		= dimsizes(var)
	dim0		= var&$var!0$
	dim1		= var&$var!1$
	dim2		= var&$var!2$
	dim3		= var&$var!3$
	dim4		= var&$var!4$

;======================================================================
; explicitly declare file definition mode. Improve efficiency.
;======================================================================

    setfileoption(fout,"DefineMode",True)

;======================================================================
; create global attributes of the file
;======================================================================

	fAtt				= True				; assign file attributes
	fAtt@title			= ""
	fAtt@source_file	= ""
	fAtt@Conventions	= ""
	fAtt@creation_date	= systemfunc("date")
	fileattdef( fout, fAtt )				; copy file attributes

;======================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;======================================================================

	dimNames	= (/ dNames(0),  dNames(1), dNames(2), dNames(3), dNames(4) /)
	dimSizes	= (/ dsize(0),   dsize(1),  dsize(2),  dsize(3),  dsize(4)  /) 
	dimUnlim	= (/ False, False, False, False, False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

;======================================================================
; predefine the dimensionality of the variables to be written out
;======================================================================
; Here we are using NCL functions to facilitate defining
; each variable's dimension name(s) and type.
; The following could be replaced with explicit, user defined dimension
; names different from those associated with the variable in memory.
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/))
;======================================================================

	filevardef(fout, dNames(0),     typeof(dim0), getvardims(dim0))
	filevardef(fout, dNames(1),     typeof(dim1), getvardims(dim1))
	filevardef(fout, dNames(2),     typeof(dim2), getvardims(dim2))
	filevardef(fout, dNames(3),     typeof(dim3), getvardims(dim3))
	filevardef(fout, dNames(4),     typeof(dim4), getvardims(dim4))
	filevardef(fout, var@long_name, typeof(var),  getvardims(var))
; variable name on the file 

; different from name on script
;======================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;======================================================================

	filevarattdef(fout, var@long_name, var)				; copy var attributes
	filevarattdef(fout, dNames(0), dim0)				; copy dim0 attributes
	filevarattdef(fout, dNames(1), dim1)				; copy dim1 attributes
	filevarattdef(fout, dNames(2), dim2)				; copy dim2 attributes
	filevarattdef(fout, dNames(3), dim3)				; copy dim3 attributes
	filevarattdef(fout, dNames(4), dim4)				; copy dim4 attributes

;======================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;======================================================================

	setfileoption(fout,"DefineMode",False)

;======================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;======================================================================

	fout->$var!0$			= (/dim0/)
	fout->$var!1$			= (/dim1/)
	fout->$var!2$			= (/dim2/)
	fout->$var!3$			= (/dim3/)
	fout->$var!4$			= (/dim4/)
	fout->$var@long_name$	= (/var/)

;======================================================================
end
;======================================================================






;======================================================================
; Procedure to write a var with 6 dimensions (dim0, dim1, dim2, dim3, dim4, dim5)
;======================================================================

procedure write_6D( fileout, var )
;======================================================================

local fout, dNames, dsize, dim0, dim1, dim2, dim3, dim4, dim5, dimNames, \
      dimSizes, dimUnlim
begin

;======================================================================
; Variable var have associated meta data: 
; (a) coordinate variables time, lat, lon and
; (b) attributes
;======================================================================

	fout		= addfile (fileout, "c")			; open output file
	dNames		= getvardimnames(var)
	dsize		= dimsizes(var)
	dim0		= var&$var!0$
	dim1		= var&$var!1$
	dim2		= var&$var!2$
	dim3		= var&$var!3$
	dim4		= var&$var!4$
	dim5		= var&$var!5$

;======================================================================
; explicitly declare file definition mode. Improve efficiency.
;======================================================================

    setfileoption(fout,"DefineMode",True)

;======================================================================
; create global attributes of the file
;======================================================================

	fAtt				= True				; assign file attributes
	fAtt@title			= ""
	fAtt@source_file	= ""
	fAtt@Conventions	= ""
	fAtt@creation_date	= systemfunc("date")
	fileattdef( fout, fAtt )				; copy file attributes

;======================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;======================================================================

	dimNames	= (/ dNames(0),  dNames(1), dNames(2), dNames(3), dNames(4), dNames(5) /)
	dimSizes	= (/ dsize(0),   dsize(1),  dsize(2),  dsize(3),  dsize(4),  dsize(5)  /) 
	dimUnlim	= (/ False, False, False, False, False, False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

;======================================================================
; predefine the dimensionality of the variables to be written out
;======================================================================
; Here we are using NCL functions to facilitate defining
; each variable's dimension name(s) and type.
; The following could be replaced with explicit, user defined dimension
; names different from those associated with the variable in memory.
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/))
;======================================================================

	filevardef(fout, dNames(0),     typeof(dim0), getvardims(dim0))
	filevardef(fout, dNames(1),     typeof(dim1), getvardims(dim1))
	filevardef(fout, dNames(2),     typeof(dim2), getvardims(dim2))
	filevardef(fout, dNames(3),     typeof(dim3), getvardims(dim3))
	filevardef(fout, dNames(4),     typeof(dim4), getvardims(dim4))
	filevardef(fout, dNames(5),     typeof(dim5), getvardims(dim5))
	filevardef(fout, var@long_name, typeof(var),  getvardims(var))
; variable name on the file 

; different from name on script
;======================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;======================================================================

	filevarattdef(fout, var@long_name, var)				; copy var attributes
	filevarattdef(fout, dNames(0), dim0)				; copy dim0 attributes
	filevarattdef(fout, dNames(1), dim1)				; copy dim1 attributes
	filevarattdef(fout, dNames(2), dim2)				; copy dim2 attributes
	filevarattdef(fout, dNames(3), dim3)				; copy dim3 attributes
	filevarattdef(fout, dNames(4), dim4)				; copy dim4 attributes
	filevarattdef(fout, dNames(5), dim5)				; copy dim5 attributes

;======================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;======================================================================

	setfileoption(fout,"DefineMode",False)

;======================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;======================================================================

	fout->$var!0$			= (/dim0/)
	fout->$var!1$			= (/dim1/)
	fout->$var!2$			= (/dim2/)
	fout->$var!3$			= (/dim3/)
	fout->$var!4$			= (/dim4/)
	fout->$var!5$			= (/dim5/)
	fout->$var@long_name$	= (/var/)

;======================================================================
end
;======================================================================






;======================================================================
; Procedure to write a var with 3 dimensions (dim0, dim1, dim2)
; First Dimension (for eq. time) is unlimited
;======================================================================

procedure write_3D_unlimited(fileout,var)
;======================================================================

local fout, dNames, dsize, dim0, dim1, dim2, dimNames, dimSizes, dimUnlim
begin

;======================================================================
; Variable var have associated meta data: 
; (a) coordinate variables time, lat, lon and
; (b) attributes
;======================================================================

	fout		= addfile (fileout, "c")			; open output file
	dNames		= getvardimnames(var)
	dsize		= dimsizes(var)
	dim0		= var&$var!0$
	dim1		= var&$var!1$
	dim2		= var&$var!2$

;======================================================================
; explicitly declare file definition mode. Improve efficiency.
;======================================================================

    setfileoption(fout,"DefineMode",True)

;======================================================================
; create global attributes of the file
;======================================================================

	fAtt				= True				; assign file attributes
	fAtt@title			= ""
	fAtt@source_file	= ""
	fAtt@Conventions	= ""
	fAtt@creation_date	= systemfunc("date")
	fileattdef( fout, fAtt )				; copy file attributes

;======================================================================
; predefine the coordinate variables and their dimensionality
; Note: to get an UNLIMITED record dimension, we set the dimensionality
; to -1 (or the actual size) and set the dimension name to True.
;======================================================================

	dimNames	= (/ dNames(0),  dNames(1), dNames(2) /)
	dimSizes	= (/        -1,   dsize(1),  dsize(2) /) 
	dimUnlim	= (/ True , False, False /)
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

;======================================================================
; predefine the dimensionality of the variables to be written out
;======================================================================
; Here we are using NCL functions to facilitate defining
; each variable's dimension name(s) and type.
; The following could be replaced with explicit, user defined dimension
; names different from those associated with the variable in memory.
; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/))
;======================================================================

	filevardef(fout, dNames(0),     typeof(dim0), getvardims(dim0))
	filevardef(fout, dNames(1),     typeof(dim1), getvardims(dim1))
	filevardef(fout, dNames(2),     typeof(dim2), getvardims(dim2))
	filevardef(fout, var@long_name, typeof(var),  getvardims(var))
; variable name on the file 

; different from name on script
;======================================================================
; Copy attributes associated with each variable to the file
; All attributes associated with each variable will be copied.
;======================================================================

	filevarattdef(fout, var@long_name, var)				; copy var attributes
	filevarattdef(fout, dNames(0), dim0)				; copy dim0 attributes
	filevarattdef(fout, dNames(1), dim1)				; copy dim1 attributes
	filevarattdef(fout, dNames(2), dim2)				; copy dim2 attributes

;======================================================================
; explicitly exit file definition mode. **NOT REQUIRED**
;======================================================================

	setfileoption(fout,"DefineMode",False)

;======================================================================
; output only the data values since the dimensionality and such have
; been predefined. The "(/", "/)" syntax tells NCL to only output the
; data values to the predefined locations on the file.
;======================================================================

	fout->$var!0$			= (/dim0/)
	fout->$var!1$			= (/dim1/)
	fout->$var!2$			= (/dim2/)
	fout->$var@long_name$	= (/var/)

;======================================================================
end
;======================================================================







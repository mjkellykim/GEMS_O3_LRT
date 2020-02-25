PRO ds_CAL_AMV_STATS, Sonde_Windspeed, AMV_Windspeed $ ;; Input
						  , Sonde_UW=Sonde_UW	$ ;; Input
						  , Sonde_VW=Sonde_VW	$ ;; Input
						  , AMV_UW=AMV_UW,AMV_VW=AMV_VW	$ ;; Input
						  , bias=bias, rmse=rmse $ ;; Output
              , nbias=nbias $
						  , avg_Sonde_Windspeed=avg_sonde_ws $ ;; Output
						  , avg_AMV_Windspeed=avg_amv_ws $ ;; Output
              , stddev_Sonde=stddev_Sonde $ ;; Output
              , stddev_AMV=stddev_AMV $ ;; Output
						  , Correlation=r	$ ;; Output
              , MVD=MVD,RMSVD=RMSVD,NRMSVD=NRMSVD $ ;; Output
              , NMVD=NMVD $
              , NaN=NaN,Count=N_data ;; Option

;;-----------------------------------------------------------------------------
IF((SIZE(Sonde_UW))[0]+(SIZE(Sonde_VW))[0] $
  +(SIZE(AMV_UW))[0]+(SIZE(AMV_VW))[0] EQ 4)THEN BEGIN
	vector_check='yes'
ENDIF ELSE BEGIN
	vector_check='no'
ENDELSE
;;=============================================================================
;; 1. NaN Value Filter	  		  
;;=============================================================================
;; All Data
n_data_1=(SIZE(Sonde_Windspeed))[-1]
n_data_2=(SIZE(AMV_Windspeed))[-1]

;;-----------------------------------------------------------------------------
;; Without NaN Value
IF(KEYWORD_SET(NaN))THEN BEGIN
  sonde_finiteidx = where(finite(sonde_windspeed) eq 1)
  amv_finiteidx = where(finite(amv_windspeed) eq 1)
  finiteidx = where(finite(sonde_windspeed) eq 1 and $
    finite(amv_windspeed) eq 1, n_data_1)
  n_data_2 = n_data_1
	;Sonde_WindSpeed=Sonde_WindSpeed[WHERE(FINITE(Sonde_WindSpeed) EQ 1, n_data_1)]
	;AMV_WindSpeed  =AMV_WindSpeed[WHERE(FINITE(AMV_WindSpeed) EQ 1, n_data_2)]
  _Sonde_WindSpeed=Sonde_WindSpeed[finiteidx]
  _AMV_WindSpeed  =AMV_WindSpeed[finiteidx]

  _amv_uw = amv_uw[finiteidx]
  _amv_vw = amv_vw[finiteidx]

  _sonde_uw = sonde_uw[finiteidx]
  _sonde_vw = sonde_vw[finiteidx]
ENDIF

;;=============================================================================

;;=============================================================================
;; 2. Scalar Calculation
;;=============================================================================
IF(n_data_1 EQ n_data_2)THEN BEGIN
	n_data=DOUBLE(n_data_1)

    DIFF=_AMV_Windspeed - _Sonde_Windspeed

	SUM_1=DIFF.TOTAL()
	SUM_2=(DIFF^2).TOTAL()

	Sonde_SUM=_Sonde_Windspeed.TOTAL()
	AMV_SUM=_AMV_Windspeed.TOTAL()
;;-----------------------------------------------------------------------------
;; Results 1
	Bias=Sum_1/n_data
	RMSE=SQRT(SUM_2/n_data)
    avg_Sonde_WS=Sonde_SUM/n_data
	avg_AMV_WS	=AMV_SUM/n_data

    NBias=bias/avg_sonde_ws
;;-----------------------------------------------------------------------------
	IF(ARG_PRESENT(stddev_Sonde) $
	  +ARG_PRESENT(stddev_AMV)	  $
	  +ARG_PRESENT(r) GT 0)THEN BEGIN

		Sonde_sd=TOTAL(_Sonde_Windspeed - avg_Sonde_WS)
		AMV_sd  =TOTAL(_AMV_Windspeed - avg_AMV_WS)

		Sonde_std=TOTAL((_Sonde_Windspeed - avg_Sonde_WS)^2)
		AMV_std	=TOTAL((_AMV_Windspeed - avg_AMV_WS)^2)
		r1		=TOTAL((_Sonde_Windspeed - avg_Sonde_WS)*(_AMV_Windspeed - avg_AMV_WS))
		r2		=SQRT(TOTAL((_Sonde_Windspeed - avg_Sonde_WS)^2) $
                     *TOTAL((_AMV_Windspeed - avg_AMV_WS)^2))
;;-----------------------------------------------------------------------------
;; Results 2
		stddev_Sonde=SQRT(Sonde_std/(n_data-1))
  		stddev_AMV	=SQRT(AMV_std/(n_data-1))
  	 	r	        =r1/r2
;;-----------------------------------------------------------------------------	
	ENDIF
ENDIF ELSE BEGIN
  PRINT,'Error!! Check Your Data Size and Dimension.'
ENDELSE
;;=============================================================================


;;=============================================================================
;; 3. Vector Calculation
;;=============================================================================
IF((vector_check EQ 'yes'))THEN BEGIN
	VD=SQRT((_AMV_UW - _Sonde_UW)^2 + (_AMV_VW - _Sonde_VW)^2)
	MVD=MEAN(VD)
    NMVD=MVD/avg_sonde_ws

    RMSVD=SQRT(MEAN((_AMV_UW-_Sonde_UW)^2 + (_AMV_VW-_Sonde_VW)^2))
    NRMSVD=RMSVD/avg_sonde_ws
ENDIF
;;=============================================================================
END

;PRO mj_gchem_amv_v0
;; TOMV Calculation for 200808

  ;; 1. Library Setting
  ;; ===============================================================================;;
  COMPILE_OPT IDL2
  MOD_GCHEM_AMV
  CD, CURRENT=runpath
  !PATH=!PATH+':'+EXPAND_PATH('+'+runpath)  ;; For compiling the local library

  ;; 2. Read Location Information
  ;; ===============================================================================;;
  readme='/DB/Data/geoschem/readme'
  eline=FILE_LINES(readme)
  OPENR,lun,readme,/GET_LUN
  tmp=''
  line=0
  idxs=[]  
  hgt=[]
  pres=[]
  etas=[]
  WHILE line LT eline DO BEGIN
    READF,lun, tmp
    IF STRMID(tmp,0,7) EQ '2. Vert' THEN BEGIN
      FOR iline=1, eline-line-1 DO BEGIN 
        READF,lun, tmp
        IF iline GT 5 AND STRTRIM(STRMID(tmp,0,2)) NE '' THEN BEGIN
          out=STRSPLIT(tmp,/EXT)
          idx=FLOAT(out[0])
          idxs=[idxs, idx]
          pres=[pres, out[3]]
          hgt=[hgt, out[2]]
          etas=[etas, out[1]]
          IF idx EQ 1 THEN BREAK
        ENDIF
      ENDFOR
      line=line+iline
    ENDIF
    line=line+1
  ENDWHILE ; line 
  FREE_LUN, lun
  hgt=REVERSE(hgt)
  pres=REVERSE(pres)


  ;; 3. Plot Setting
  ;; ===============================================================================;;
  ; set flags for process
  amv_cal=1
  plot1=0  ; plot concentration and wind field
    wnd_on=1  ;
    po_idx=0 
    iano = 0 ; 1: anomaly open
      ;po_lonidx=108 & po_latidx=96
      po_lonidx=108 & po_latidx=93
  plot2=0   ; plot calculated flux field
  plot3=1   ; plot cross correlation image
    wnd_type=1   ; 1:wind speed, 2:wind direction
  plot4=1   ; time series
    dx=2   ; averaging area
    kst_on=1
    inc_clim=1
    clim_xdr='o3_clim2.xdr'
  plot5=1   ; clim   iano should be 0 for creating o3_clim variable
    sellev=0
    stidx=0 & etidx=22  ; 0 ~ 22(01utc ~ 23utc) , only for used in plotting
  plot6=0   ; column (du)  
    unit_conv=1  ;  ppbv to du
  
  plot7=1 ; wind speed correlation map between tomv and geoschem 
    plot7_wnd_on = 1
    plot7_wnd_type = 1 ; 1:wind speed, 2:wind direction
    plot7_o3 = 0

  plot8=1;
    plot8_wnd_on = 1
    plot8_wnd_type = 1 ; 1:wind speed, 2:wind direction
    plot8_corr = 1

  plot9=0 ; wind speed and direction correlation map between tomv and geoschem 
    plot9_wnd_on = 1
    plot9_o3 = 0

  ;; 4. Read Data
  ;; ===============================================================================;;
  ws_amvs=[]
  wd_amvs=[]
  ws_models=[]
  wd_models=[]

  t_sta=[]
  tjul_sta=[]
  o3_sta=[]
  o3_sta_clim=[]

  IF plot4 EQ 1 AND inc_clim EQ 1 THEN BEGIN
    RESTORE, clim_xdr
    time_clim=[]
  ENDIF

  ;files=FILE_SEARCH('/DB/Data/geoschem/201605/*.nc'
  files=FILE_SEARCH('/DB/Data/geoschem/200808/*.nc')
  files=files[14:15]
  nfiles=N_ELEMENTS(files)
  FOR ifile=0,nfiles-1 DO BEGIN
    file=files[ifile]
    print, ' :: ifile = ', ifile, ',   : ',file

    time0 = ncread (file, 'time')
    lon0 = ncread (file, 'lon')
    lat0 = ncread (file, 'lat')
    lev0 = ncread (file, 'lev')
    o30 = ncread (file, 'IJ_AVG_S__O3')
    uwnd0 = ncread (file, 'DAO_3D_S__UWND')
    vwnd0 = ncread (file, 'DAO_3D_S__VWND')
    temp0 = ncread (file, 'DAO_3D_S__TMPU')

    date=STRMID(file, 10,8,/rev)
    ;caldat,(julday(1,1,1985,0,0,0)*86400.0d0+time0*3600)/86400.0d0,mon,day,year,hour,min,sec
    ;time=STRING(year,'(I04)')+STRING(mon,'(I02)')+STRING(day,'(I02)')+STRING(hour,'(I02)')
    time=date+STRING(FINDGEN(23)+1,F='(I02)')

    nlon=N_ELEMENTS(lon0)
    nlat=N_ELEMENTS(lat0)
    ntime=N_ELEMENTS(time0)
    nlev=N_ELEMENTS(lev0)

    lon=FLTARR(nlon,nlat)
    lat=FLTARR(nlon,nlat)

    IF plot5 EQ 1 THEN BEGIN
      IF ifile EQ 0 THEN o3_clim=FLTARR(nlon,nlat,nlev,ntime,nfiles)
      o3_clim[*,*,*,*,ifile]=o30
      time_clim=[time_clim,time]
    ENDIF
    IF iano EQ 1 THEN RESTORE, 'o3_clim.xdr'
    FOR i=0,nlon-1 DO lon[i,*]=lon0[i] 
    FOR i=0,nlat-1 DO lat[*,i]=lat0[i]

    ;FOR tidx=0,N_ELEMENTS(time)-1 DO BEGIN
    ;FOR tidx=5,5 DO BEGIN
    tidxs=FINDGEN(ntime) 
    ;tidxs=[0, 11., 12.] 
    levidxs=[10] ; level index for 
    FOR ilevidx=0,N_ELEMENTS(levidxs)-1 DO BEGIN
    FOR itidx=0,N_ELEMENTS(tidxs)-1 -1 DO BEGIN
      tidx=tidxs[itidx]
      levidx=levidxs[ilevidx]  ; 10: 850, 21: 550, 28: 267

      IF plot4 EQ 1 THEN BEGIN
        t_sta=[t_sta,time[tidx]]
        tjul_tmp=julday(STRMID(time[tidx],4,2),STRMID(time[tidx],6,2),STRMID(time[tidx],0,4),STRMID(time[tidx],8,2))
        tjul_sta=[tjul_sta,tjul_tmp]
        o3_sta=[o3_sta,MEAN(o30[po_lonidx-dx:po_lonidx+dx,po_latidx-dx:po_latidx+dx,levidx,tidx])]
        IF inc_clim EQ 1 THEN BEGIN
          o3_sta_clim=[o3_sta_clim,MEAN(o3_clim2[po_lonidx-dx:po_lonidx+dx,po_latidx-dx:po_latidx+dx,levidx,tidx])]
        ENDIF
      ENDIF
      IF amv_cal EQ 1 THEN BEGIN
        o3f_1=REFORM(o30[*,*,levidx,tidx-1])
        o3f_2=REFORM(o30[*,*,levidx,tidx])
        o3f_3=REFORM(o30[*,*,levidx,tidx+1])

        sat_time=time[[tidx-1, tidx, tidx+1]]


  ;; 5. TOMV Calculate
  ;; ===============================================================================;;
      ;;========================================================================
      ;; Set Environment
      ;;========================================================================
        mpi_set='N'
          mpi_num=[0]
        Nested='Y'

        ;Target_Size = 6 
        ;Search_Size = 10   ; have to larger than target_size+1, target+3 in case, target+5~6 in korus
        Target_Size = 6 
        Search_Size = 14   ; have to larger than target_size+1, target+3 in case, target+5~6 in korus
        TS = string(target_size, format='(I02)')
        SS = string(search_size, format='(I02)')
        ;Grid_Size   =8 
        Grid_Size   = target_size 

        nImages=3   ;; 2: t=t0-dt, t0       , No Averaged
                    ;; 3: t=t0-dt, t0, t0+dt, Averaged

        Interval=60 ;; mins

        First_Guess='Y'

        nx=nlon
        ny=nlat
      ;;-----------------------------------------------------------------------
        Configure={MPI_SET:(mpi_set.SubString(0,0)).ToUpper()         $
                  ,MPI_num:mpi_num                                    $
                  ,Nested:(Nested.SubString(0,0)).ToUpper()           $
                  ,First_Guess:(First_Guess.SubString(0,0)).ToUpper() $
                  ,Target_Size:Target_Size,Search_Size:Search_Size    $
                  ,Grid_Size:Grid_Size                                $
                  ,nImages:nImages                                    $
                  ,nx:nx,ny:ny}

        PRINT,'Pre-Processing [Configure] .................................... '
      ;;========================================================================


      ;;========================================================================
      ;; Run Programs
      ;;========================================================================
        PRINT,'Run AMV Main Program .......................................... '
        PRINT,''

        imgs=FLTARR(nlon,nlat,3)
        imgs[*,*,1]=o3f_1
        imgs[*,*,0]=o3f_2  ;target
        imgs[*,*,2]=o3f_3

        finite_image=FINITE(TOTAL(imgs,3))

        half_search=Configure.Search_Size/2
        half_Grid  =Configure.Grid_Size/2
        half_Target=Configure.Target_Size/2
        CC_Buf_Size=Configure.Search_Size-Configure.Target_size+1


        STDDEV_Val=FLTARR(Configure.Target_Size+1,Configure.Target_Size+1)
        ;STDDEV_Val=FLTARR(Configure.Target_Size,Configure.Target_Size) ; geun modi (Configure.Target_Size+1 -> Configure.Target_Size)
        Cross_Corr=FLTARR(CC_Buf_Size,CC_Buf_Size)

        TIC

        WINDS = !NULL
        corrs = !null

        FOR iscene= 1,2 DO BEGIN  ; iscene:1 t0-dt -> t0,  iscene:2 t0 -> t0+dt
          GRID_X=!null
          GRID_Y=!null
          MAX_Corr=!null
          DEST_X=!null
          DEST_Y=!null
          DEST_LON=!null ; geun added
          DEST_LAT=!null ; geun added
          bflag=!null

          FOR ix=0+half_search+half_Target,Configure.nx-half_search-half_Target,Configure.Grid_Size DO BEGIN
          FOR iy=0+half_search+half_Target,Configure.ny-half_search-half_Target,Configure.Grid_Size DO BEGIN

            finite_Buf=finite_image[ix-half_search-half_Target:ix+half_search+half_Target  $
                                   ,iy-half_search-half_Target:iy+half_search+half_Target]

            IF(TOTAL(finite_Buf) $
               EQ ((Configure.Search_Size+Configure.Target_size+1)^2))THEN BEGIN ;; Filtering NaN Value

              STDDEV_Val[*,*]=!Values.f_nan  ;; Calc Standard Deviation
              FOR ix_in= -half_Target, half_Target DO BEGIN   ; geun mod  (half_Target -> half_Target-1)
              FOR iy_in= -half_Target, half_Target DO BEGIN   ; geun mod  (half_Target -> half_Target-1)
                std_ix=ix_in+half_Target
                std_iy=iy_in+half_Target

                STDDEV_Val[std_ix,std_iy]                                          $
                          =STDDEV(imgs[ix+ix_in-half_Target:ix+ix_in+half_Target    $
                                     ,iy+iy_in-half_Target:iy+iy_in+half_Target,0])
              ENDFOR ;; ix_Buf
              ENDFOR ;; iy_Buf

              max_stddev=MAX(STDDEV_Val,max_pos)
              max_pos=ARRAY_INDICES(STDDEV_VAL,max_pos)

        ;;    Get Grid Center
              GRID_X=[GRID_X,ix+max_pos[0]-half_Target]
              GRID_Y=[GRID_y,iy+max_pos[1]-half_Target]

              Target=imgs[GRID_X[-1]-half_Target:GRID_X[-1]+half_Target $
                        ,GRID_Y[-1]-half_Target:GRID_Y[-1]+half_Target,0]

              Search=imgs[GRID_X[-1]-half_Search:GRID_X[-1]+half_Search $
                        ,GRID_Y[-1]-half_Search:GRID_Y[-1]+half_Search,iscene] ;t+1
                      ; if iscene eq 1 then t-1, if iscene eq 1 then t+1

              Cross_Corr[*,*]=!values.f_nan ;; variable for Calc. Cross Correlation

              ; filling cross correlation for each search grid
              FOR ix_in=0,CC_Buf_size-1 DO BEGIN
              FOR iy_in=0,CC_Buf_size-1 DO BEGIN
                Search_Buf=Search[ix_in:ix_in+Configure.Target_Size $
                                 ,iy_in:iy_in+Configure.Target_Size]

                Cross_Corr[ix_in,iy_in]=Correlate(Target,Search_Buf)
              ENDFOR
              ENDFOR

              tmp_max=MAX(Cross_Corr,max_pos)
              MAX_Corr=[MAX_Corr,tmp_max]
              max_pos=ARRAY_INDICES(Cross_Corr,max_pos)

              ibuf=[max_pos[0]-1,max_pos[0],max_pos[0]+1]
              jbuf=[max_pos[1]-1,max_pos[1],max_pos[1]+1]

              IF (MIN(ibuf) GE 0) AND (MAX(ibuf) LE CC_Buf_size-1) AND $
                 (MIN(jbuf) GE 0) AND (MAX(jbuf) LE CC_Buf_size-1) THEN BEGIN
                bflag=[bflag,1] 
                Buf_CC=Cross_Corr[ibuf[0]:ibuf[2] $  ;; 3x3
                                 ,jbuf[0]:jbuf[2]]

                bot_x=2*(Buf_CC[0,1]+Buf_CC[2,1]-2.*Buf_CC[1,1])
                bot_y=2*(Buf_CC[1,0]+Buf_CC[1,2]-2.*Buf_CC[1,1])

                IF(ABS(bot_x) GT 6.10E-5)THEN BEGIN
                  del_x=(Buf_CC[0,1]-Buf_CC[2,1])/bot_x
                ENDIF ELSE BEGIN
                  del_x=0.
                ENDELSE

                IF(ABS(bot_y) GT 6.10E-5)THEN BEGIN
                  del_y=(Buf_CC[1,0]-Buf_CC[1,2])/bot_y
                ENDIF ELSE BEGIN
                  del_y=0.
                ENDELSE

                id=del_x
                jd=del_y
              ENDIF ELSE BEGIN
                bflag=[bflag,0]
                id=0  &  jd=0
              ENDELSE 

  ;;
              ;Buf_Fit_CC=SFIT(Buf_CC,2,MAX_Degree=2,KX=Coeff)
                ;Detrm=4*coeff[5]*coeff[2] - coeff[4]^2

              ;IF(Detrm GT 1.5D-16) THEN BEGIN
                ;id=(coeff[1]*coeff[4]-2*coeff[3]*coeff[2])/Detrm -1
                ;jd=(coeff[3]*coeff[4]-2*coeff[1]*coeff[5])/Detrm -1        
                ;IF id GT SQRT(2) OR id LT 0 THEN id=0   ; geun added
                ;IF jd GT SQRT(2) OR jd LT 0 THEN jd=0   ; geun added
              ;ENDIF ELSE BEGIN
                ;id=0.
                ;jd=0.
              ;ENDELSE

              DEST_X=[DEST_X,GRID_X[-1]-(CC_BUF_size/2)+ibuf[1]+1+id]
              DEST_Y=[DEST_Y,GRID_Y[-1]-(CC_BUF_size/2)+jbuf[1]+1+jd]
              ;DEST_X=[DEST_X,GRID_X[-1]+(CC_BUF_size/2)-ibuf[1]-id]   ; iscene 1 test
              ;DEST_Y=[DEST_Y,GRID_Y[-1]+(CC_BUF_size/2)-jbuf[1]-jd]   ; iscene 1 test

              dest_x1=FIX(dest_x[-1])  &   dest_x2=FIX(dest_x[-1]+1)
              dest_y1=FIX(dest_y[-1])  &   dest_y2=FIX(dest_y[-1]+1)
              dlon=lon[dest_x2,dest_y2]-lon[dest_x1,dest_y1]
              dlat=lat[dest_x2,dest_y2]-lat[dest_x1,dest_y1]
              xrat=dest_x[-1]-dest_x1
              yrat=dest_y[-1]-dest_y1
              dest_lon=[dest_lon,lon[dest_x1,dest_y1]+dlon*xrat]
              dest_lat=[dest_lat,lat[dest_x1,dest_y1]+dlat*yrat]       

              ;stop

            ENDIF ;; Total(Buf)
          ENDFOR ;; ix
          ENDFOR ;; iy

          inu=N_ELEMENTS(grid_x)
          grid_lon=lon[grid_x,grid_y]
          grid_lat=lat[grid_x,grid_y]
          dis_us_cmp = AMV_X_Distance(grid_lat, grid_lon,dest_lat,dest_lon )
          dis_vs_cmp = AMV_Y_Distance(grid_lat, dest_lat )

          ; Designate Direction AFTER LAT/LON Calculation
          u_sign = REPLICATE(1., inu)
          v_sign = u_sign
       
          uback = WHERE(grid_x GT dest_x, nuback)
          vback = WHERE(grid_y GT dest_y, nvback)

          ;if(clust_method EQ 'YES') THEN stop
          IF (nuback GT 0.) THEN u_sign[uback] = -1.
          IF (nvback GT 0.) THEN v_sign[vback] = -1.

          IF iscene EQ 1 THEN idts = -3600. ELSE idts = 3600
          us = REPLICATE(!Values.f_NAN, inu)
          vs = us
          ;ws = us
          ;wd = us
          good = WHERE(finite(dis_us_cmp) GT 0. AND finite(dis_vs_cmp) GT 0., ngood)

          corr = REPLICATE(!Values.f_NAN, inu)

          IF (ngood GT 0) THEN BEGIN
            us[good]  = (dis_us_cmp[good] / idts ) * u_sign[good]
            vs[good]  = (dis_vs_cmp[good] / idts ) * v_sign[good]
            dum       = AMV_Convert_UV_Wvector( us[good],vs[good],ws1,wd1 )
            ;ws[good]  = ws1
            ;wd[good]  = wd1
            corr[good] = MAX_Corr[good]
          ENDIF
          bad = WHERE(ABS(us) GE 1000. OR ABS(vs) GE 1000., nbad)
          IF (nbad GT 0.) THEN BEGIN
            us[bad] = !Values.f_NAN
            vs[bad] = !Values.f_NAN
            ;ws[bad] = !Values.f_NAN
            ;wd[bad] = !Values.f_NAN
          ENDIF
          ;WIND = {WS:ws,WD:wd,US:us,VS:vs}
          nan=WHERE(bflag EQ 0,nnan)
          us[nan]=!values.f_nan  &  vs[nan]=!values.f_nan
          WIND = {US:us,VS:vs}
          winds=[winds,wind]     
          corrs=[corrs, corr]
        ENDFOR ; iscene
        TOC
        PRINT,'----------------------------------------------------------------'
        us_avg=MEAN(winds.us,dim=2)
        vs_avg=MEAN(winds.vs,dim=2)

        ncorrs = N_ELEMENTS(corrs)
        corrs = reform(corrs, ncorrs/2, 2)
        corrs_mean = mean(corrs, dimension=2)

        dum       = AMV_Convert_UV_Wvector(us_avg,vs_avg,ws,wd )

;=====================================================================
; test for correlation between t-1 and t+1
        ;stop
        wus = winds.us
        wvs = winds.vs
        dum = AMV_Convert_UV_Wvector(wus[*, 0], wvs[*, 0], ws0, wd0)
        dum = AMV_Convert_UV_Wvector(wus[*, 1], wvs[*, 1], ws1, wd1)

        ;print, correlate(ws0, ws1)
        ;print

        ;stop
;=====================================================================

        ;dum       = AMV_Convert_UV_Wvector(us,vs,ws,wd )
      ;;========================================================================
      ENDIF ; amv_cal

      ;tidx=tidx+1

      the_time=time[tidx]
      the_hgt=hgt[levidx]
      the_pres=pres[levidx]

      ;limit=[-11,60,55,150]
      limit=[10,85,50,145]   ; sele
      ;limit=[29,115,43,140]

      ;limit=[15,70,55,140]

      o3_tmp=REFORM(o30[*,*,levidx,tidx])


  ;; 6. Start Plotting
  ;; ===============================================================================;;
      IF plot1 EQ 1 THEN BEGIN
        tail='d'+date+'_t'+ STRMID(the_time,1,2,/REV) + '_l'+STRING(the_pres, F='(I04)')
        outname='gchem_comp_plot1_TS' + TS + '_SS' + SS
        out_ps=outname+tail+'.ps'
        out_png=outname+tail+'.png'
        ;CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=14,/NOMATCH,/INCH
        CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=7,/NOMATCH,/INCH
        pos=[0.12,0.12,0.9,0.91]
        !P.POSITION = pos

        ct=72
        nclev=51
        clev=(FINDGEN(nclev))*2
        barcol=REVERSE(FINDGEN(nclev-1)*253/(nclev-1)+1)

        clev[0]=-9999  &  clev[nclev-1]=9999

        ;znan=WHERE(z LT MEAN(z)-3*STDDEV(z) AND z GT mean(z)+3*STDDEV(z), zval)
        ;z[znan]=!VALUES.F_NAN

        LOADCT, ct 
        TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
        TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
        !P.COLOR=255
        !P.BACKGROUND=254
        !P.CHARSIZE=2
        !P.CHARTHICK=2

        ;title='!6Ozone (t=' + STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
        title='!6Ozone (d='+ date +', t='+ STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
        form='(I2)' 

        z=o3_tmp
        IF iano EQ 1 THEN BEGIN
          ;clev=(FINDGEN(nclev))*0.9-30.0
          clev=(FINDGEN(nclev))-25
          clev[0]=-9999  &  clev[nclev-1]=9999
          z=o3_tmp-o3_clim[*,*,levidx,tidx]
          form='(I3)' 
        ENDIF
        
        MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,POS=pos,title='!6'
        CONTOUR, z, lon, lat, /CELL_FILL, C_COLORS=barcol, LEVELS=clev,/OVER, $
                 XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos

        MAP_CONTINENTS,/CONTINENT,/COUNT,MLINETHICK=1.8,/COUNTRY,/HIRES 
        xtlen=20
        ytlen=10 

        IF wnd_on EQ 1 THEN BEGIN

      ; nwp wind
          uwnd=REFORM(uwnd0[*,*,levidx,tidx])
          vwnd=REFORM(vwnd0[*,*,levidx,tidx])
          ;PARTVELVEC, uwnd, vwnd, lon, lat, COLOR=CGCOLOR('grey'), FRACTION=0.05,/OVER, noclip=0, length=0.025
          tmp=WHERE(finite(uwnd) EQ 1 AND finite(vwnd) EQ 1, ntmp)
          tmp=FIX(ntmp*RANDOMU(seed,1000))
          ;stop
          uwnd_nwp=uwnd[tmp]
          vwnd_nwp=vwnd[tmp]
          dum = AMV_Convert_UV_Wvector(uwnd_nwp,vwnd_nwp,ws_nwp,wd_nwp )
          WINDBARB,lon[tmp],lat[tmp],ws_nwp*1.94384,wd_nwp,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('CHARCOAL'), clip=pos

      ; AMV wind
          WINDBARB,grid_lon,grid_lat,ws*1.94384,wd,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('MAGENTA'), clip=pos
          ;WINDBARB,grid_lon,grid_lat,ws*1.94384,wd,ASPECT=1,LENGTH=0.05,THICK=5,COL=CGCOLOR('MAGENTA'), clip=pos
          tmp=WHERE(ws*1.94384 LT 2.5,COMP=ttmp,ntmp)
          oplot, grid_lon[ttmp],grid_lat[ttmp], psym=plotsym_fn(/cir,/fill,scale=0.5)
          IF ntmp NE 0 THEN oplot, grid_lon[tmp],grid_lat[tmp], psym=plotsym_fn(/cir,scale=0.5)

      ; Plot target position
          ;cgoplot, dest_lon,dest_lat, psym=1
          ;cgoplot, grid_lon,grid_lat, psym=1
        ENDIF
        IF po_idx EQ 1 THEN BEGIN
          xcoord=[lon[po_lonidx-dx,po_latidx-dx],lon[po_lonidx+dx,po_latidx-dx],lon[po_lonidx+dx,po_latidx+dx], $
                  lon[po_lonidx-dx,po_latidx+dx],lon[po_lonidx-dx,po_latidx-dx]]
          ycoord=[lat[po_lonidx-dx,po_latidx-dx],lat[po_lonidx+dx,po_latidx-dx],lat[po_lonidx+dx,po_latidx+dx], $
                  lat[po_lonidx-dx,po_latidx+dx],lat[po_lonidx-dx,po_latidx-dx]]
          POLYFILL,xcoord, ycoord, COLOR=CGCOLOR('MAGENTA')
          PLOTS,xcoord, ycoord
          ;PLOTS,lon[po_lonidx,po_latidx], lat[po_lonidx,po_latidx], PSYM=plotsym_fn(/box,scal=2,/fill), COLOR=CGCOLOR('MAGENTA')
          ;PLOTS,lon[po_lonidx,po_latidx], lat[po_lonidx,po_latidx], PSYM=plotsym_fn(/box,scal=2)
       
        ENDIF ; po_idx

        CONTOUR, z, z, z, /NODAT,/NOERA,CHARSIZE=2,CHARTHICK=2, $
                 XSTYLE=1, YSTYLE=1,YTICKFORMAT='(I3)',XTICKFORMAT='(i4)',XTICKINTERVAL=xtlen,YTICKINTERVAL=ytlen, $
                 XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]], $
                 TITLE=title, XTITLE = 'Longitude',YTITLE='Latitude', POS=pos

        loadct_dg,ct
        COLORBAR_2, clev, barcol, /COL, YS =0.77, XS = 0.015, CHARSIZE=2, CHARTHICK=2,levelind=findgen(100)*5,  $
                                 FORMAT=form,UNIT='ppbv',/IVE,/NOFIRST;,/NOFRAME

        cgPs_Close,density=800, /png
        ;file_trans_ds,out_png,/DEL
      ENDIF ; plot1

      

      IF plot2 EQ 1 THEN BEGIN
        tail='t'+ STRMID(the_time,1,2,/REV) + '_l'+STRING(the_pres, F='(I04)')
        outname='gchem_flux_plot2_TS' + TS + '_SS' + SS
        out_ps=outname+tail+'.ps'
        out_png=outname+tail+'.png'
        CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=7,/NOMATCH,/INCH
        pos=[0.12,0.12,0.9,0.91]
        !P.POSITION = pos

        ct=22
        nclev=51
        ;clev=(FINDGEN(nclev))*2
        ;clev=(FINDGEN(nclev))*0.005   ; total column analysis
        clev=(findgen(nclev))*0.01   ; total column analysis
        barcol=FIX(FINDGEN(nclev-1)*253/nclev)
        clev[0]=-999  &  clev[nclev-1]=999

        ;znan=WHERE(z LT MEAN(z)-3*STDDEV(z) AND z GT mean(z)+3*STDDEV(z), zval)
        ;z[znan]=!VALUES.F_NAN

        LOADCT, ct 
        TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
        TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
        !P.COLOR=255
        !P.BACKGROUND=254
        !P.CHARSIZE=2
        !P.CHARTHICK=2

        title='!6Ozone (t=' + STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'

        RCON = 8.314
        hpa2pa = (1E2)
        mo3 = 48
        g2Mg = (1E-6)  ; gram to Megagram
        m2km = (1E-3)
        m2cm = (1E2)
        Navo=6.0221409e+23

        z=o3_tmp*1E-9 ; ppbv to vmr
        t_tmp=temp0[*,*,levidx,tidx]
        o3p_mcon = (pres[levidx]*hpa2pa*mo3)/(RCON*t_tmp)*z/(m2km)^3*g2Mg   ; mass

        IF wnd_on EQ 1 THEN BEGIN

        ; nwp wind
          uwnd=REFORM(uwnd0[*,*,levidx,tidx])
          vwnd=REFORM(vwnd0[*,*,levidx,tidx])

        ENDIF


        MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,POS=pos,title='!6'


        ; model flux
        ;o3p_flux_model=SQRT(uwnd^2+vwnd^2)*o3p_mcon  ; velocity[km/h]*mass_con
        ;CONTOUR, o3p_flux_model, lon, lat, /CELL_FILL, C_COLORS=barcol, LEVELS=clev,/OVER, $
                 ;XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos
        ;PARTVELVEC, uwnd, vwnd, lon, lat, COLOR=CGCOLOR('grey'), FRACTION=0.25,/OVER, noclip=0, length=0.025


        ; amv flux
        val=where(grid_lon GT limit[1] AND grid_lon LT limit[3] AND $
          grid_lat GT limit[0] AND grid_lat LT limit[2],nval)
        grid_lon=grid_lon[val]
        grid_lat=grid_lat[val]
        grid_x=grid_x[val]
        grid_y=grid_y[val]
        us_avg=us_avg[val]
        vs_avg=vs_avg[val]
        ntmp=N_ELEMENTS(grid_x)
        zz=FLTARR(nx,ny)*!VALUES.F_NAN
        xx=FLTARR(nx,ny)*!VALUES.F_NAN
        yy=FLTARR(nx,ny)*!VALUES.F_NAN
        uwnd_amv=FLTARR(nx,ny)*!VALUES.F_NAN
        vwnd_amv=FLTARR(nx,ny)*!VALUES.F_NAN
        FOR itmp=0,ntmp-1 DO BEGIN
          ztmp=SQRT(us_avg[itmp]^2+vs_avg[itmp]^2)*o3p_mcon[grid_x[itmp],grid_y[itmp]]     
          zz[grid_x[itmp],grid_y[itmp]]=ztmp 
          uwnd_amv[grid_x[itmp],grid_y[itmp]]=us_avg[itmp]
          vwnd_amv[grid_x[itmp],grid_y[itmp]]=vs_avg[itmp]
          xx[grid_x[itmp],grid_y[itmp]]=grid_lon[itmp]
          yy[grid_x[itmp],grid_y[itmp]]=grid_lat[itmp]
          col=-999
          FOR ilev=0,nclev-2 DO BEGIN
            IF (ztmp GE clev[ilev] AND ztmp LT clev[ilev+1]) THEN col=barcol[ilev]
          ENDFOR
          IF col NE -999 THEN PLOTS, grid_lon[itmp], grid_lat[itmp], PSYM=plotsym_fn(/box,scale=3,/fill),COLOR=col
        ENDFOR
        PARTVELVEC, us_avg, vs_avg, grid_lon, grid_lat, COLOR=CGCOLOR('grey'), FRACTION=1,/OVER, noclip=0, length=0.025


        xtlen=20
        ytlen=10 

        MAP_CONTINENTS,/CONTINENT,/COUNT,MLINETHICK=1.8,/COUNTRY,/HIRES 

        CONTOUR, z, z, z, /NODAT,/NOERA,CHARSIZE=2,CHARTHICK=2, $
                 XSTYLE=1, YSTYLE=1,YTICKFORMAT='(I3)',XTICKFORMAT='(i4)',XTICKINTERVAL=xtlen,YTICKINTERVAL=ytlen, $
                 XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]], $
                 TITLE=title, XTITLE = 'Longitude',YTITLE='Latitude', POS=pos

        LOADCT_DG,ct
        COLORBAR_2, clev, barcol, /COL, YS =0.77, XS = 0.015, CHARSIZE=1.9, CHARTHICK=2, LEVELIND = FINDGEN(100)*5, $
                                 FORMAT='(F4.1)',UNIT='',/IVE,/NOFIRST;,/NOFRAME

        XYOUTS,pos[2]-0.02,pos[3]+0.04,'Avg Flux', /nor
        XYOUTS,pos[2]-0.02,pos[3]+0.01,'(Mg/km2/h)', /nor,charsize=1

        cgPs_Close,density=800, /png
        ;FILE_TRANS_ds,out_png,/DEL
      ENDIF  ; plot2

      IF plot6 EQ 1 THEN BEGIN
        tail='d'+date+'_t'+ STRMID(the_time,1,2,/REV) + '_l'+STRING(the_pres, F='(I04)')
        outname='gchem_comp_plot6_TS' + TS + '_SS' + SS
        out_ps=outname+tail+'.ps'
        out_png=outname+tail+'.png'
        CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=7,/NOMATCH,/INCH
        pos=[0.12,0.12,0.9,0.91]
        !P.POSITION = pos

        ct=72
        nclev=51
        clev=(FINDGEN(nclev))*0.5
        ;clev=(FINDGEN(nclev))*3
        ;clev=(FINDGEN(nclev))*
        barcol=REVERSE(FINDGEN(nclev-1)*253/(nclev-1)+1)

        clev[0]=-9999  &  clev[nclev-1]=9999

        ;znan=WHERE(z LT MEAN(z)-3*STDDEV(z) AND z GT mean(z)+3*STDDEV(z), zval)
        ;z[znan]=!VALUES.F_NAN

        LOADCT, ct 
        TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
        TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
        !P.COLOR=255
        !P.BACKGROUND=254
        !P.CHARSIZE=2
        !P.CHARTHICK=2


; ppbv to du
        print, ' :: Convert ppbv to du'
        pprof=FLOAT(pres)
        zprof=FLOAT(hgt)
        nl=N_ELEMENTS(pprof)
        o3_du=FLTARR(nlon,nlat,nlev)
        FOR ix=0,nlon-1 DO BEGIN
        FOR iy=0,nlat-1 DO BEGIN
          o3p=REFORM(o30[ix,iy,*,tidx])*1E-9 ; ppbv to vmr   
          tp=REFORM(temp0[ix,iy,*,tidx])  ;                  
          o3p_du=vmr2du_geun(REVERSE(o3p),REVERSE(pprof),REVERSE(tp),alt=REVERSE(FLOAT(hgt)),/INTER)
          o3p_du=REVERSE(o3p_du) ; convet to bot2top
          o3_du[ix,iy,*]=o3p_du
          ;IF iy EQ 0 THEN print, STRING(ix,F='(I4)')+' / '+STRING(nlon,F='(I4)')
          ;IF iy EQ 0 THEN print, STRING(FLOAT(ix)/(nlon-1)*100.,F='(F5.1)')+' %'
        ENDFOR ; iy
        ENDFOR ; ix
;----------------------------

        blev=28
        o3_low=TOTAL(o3_du[*,*,0:blev],3)
        z=o3_low
        
        ;title='!6Ozone (t=' + STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
        title='!6Ozone (d='+ date +', t='+ STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(pres[blev], F='(I4)')+'hPa)'
        form='(I3)' 

        MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,POS=pos,title='!6'
        CONTOUR, z, lon, lat, /CELL_FILL, C_COLORS=barcol, LEVELS=clev,/OVER, $
                 XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos

        MAP_CONTINENTS,/CONTINENT,/COUNT,MLINETHICK=1.8,/COUNTRY,/HIRES 
        xtlen=20
        ytlen=10 

        IF wnd_on EQ 1 THEN BEGIN

      ; nwp wind
          uwnd=REFORM(uwnd0[*,*,levidx,tidx])
          vwnd=REFORM(vwnd0[*,*,levidx,tidx])
          ;PARTVELVEC, uwnd, vwnd, lon, lat, COLOR=CGCOLOR('grey'), FRACTION=0.05,/OVER, noclip=0, length=0.025
          tmp=WHERE(finite(uwnd) EQ 1 AND finite(vwnd) EQ 1, ntmp)
          tmp=FIX(ntmp*RANDOMU(seed,1000))
          uwnd_nwp=uwnd[tmp]
          vwnd_nwp=vwnd[tmp]
          dum = AMV_Convert_UV_Wvector(uwnd_nwp,vwnd_nwp,ws_nwp,wd_nwp )
          WINDBARB,lon[tmp],lat[tmp],ws_nwp*1.94384,wd_nwp,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('CHARCOAL'), clip=pos
        ENDIF
        CONTOUR, z, z, z, /NODAT,/NOERA,CHARSIZE=2,CHARTHICK=2, $
                 XSTYLE=1, YSTYLE=1,YTICKFORMAT='(I3)',XTICKFORMAT='(i4)',XTICKINTERVAL=xtlen,YTICKINTERVAL=ytlen, $
                 XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]], $
                 TITLE=title, XTITLE = 'Longitude',YTITLE='Latitude', POS=pos

        loadct_dg,ct
        COLORBAR_2, clev, barcol, /COL, YS =0.77, XS = 0.015, CHARSIZE=2, CHARTHICK=2,levelind=findgen(100)*5,  $
                                 FORMAT=form,UNIT='du',/IVE,/NOFIRST;,/NOFRAME

        cgPs_Close,density=800, /png
        ;file_trans_ds,out_png,/DEL
      ENDIF  ;plot6


      IF amv_cal EQ 1 THEN BEGIN

        ; for scatter density plot
        uwnd=REFORM(uwnd0[*,*,levidx,tidx])
        vwnd=REFORM(vwnd0[*,*,levidx,tidx])
        uwnd_model=us_avg
        vwnd_model=vs_avg
        dum = AMV_Convert_UV_Wvector(uwnd_model,vwnd_model,ws_model,wd_model)

        ; amv flux
        ntmp=N_ELEMENTS(grid_x)
        uwnd_amv=FLTARR(ntmp)
        vwnd_amv=FLTARR(ntmp)
        FOR itmp=0,ntmp-1 DO BEGIN
          uwnd_amv[itmp]=uwnd[grid_x[itmp],grid_y[itmp]]
          vwnd_amv[itmp]=vwnd[grid_x[itmp],grid_y[itmp]]
        ENDFOR
        dum = AMV_Convert_UV_Wvector(uwnd_amv,vwnd_amv,ws_amv,wd_amv)
        ws_amvs=[ws_amvs,ws_amv]
        wd_amvs=[wd_amvs,wd_amv]
        ws_models=[ws_models,ws_model]
        wd_models=[wd_models,wd_model]

        print, tidx
      ENDIF ; amv_cal



      IF plot7 EQ 1 THEN BEGIN
        tail='d'+date+'_t'+ STRMID(the_time,1,2,/REV) + '_l'+STRING(the_pres, F='(I04)')
        if plot7_wnd_type eq 1 then BEGIN
          typestr = 'ws_'
        endif else BEGIN
          typestr = 'wd_'
        endelse

        outname= typestr + 'corr_gchem_tomv_plot7_TS' + TS + '_SS' + SS
        out_ps=outname+tail+'.ps'
        out_png=outname+tail+'.png'
        ;CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=14,/NOMATCH,/INCH
        CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=7,/NOMATCH,/INCH
        pos=[0.12,0.12,0.9,0.91]
        !P.POSITION = pos

        title='!6Ozone (d='+ date +', t='+ STRMID(the_time,1,2,/REV) + $
          'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
        if plot7_o3 then begin
          ct=73
          nclev=51
          clev=(FINDGEN(nclev))/(nclev -1)
          print, clev
          barcol=REVERSE(FINDGEN(nclev-1)*253/(nclev-1)+1)

          clev[0]=-9999  &  clev[nclev-1]=9999

          ;znan=WHERE(z LT MEAN(z)-3*STDDEV(z) AND z GT mean(z)+3*STDDEV(z), zval)
          ;z[znan]=!VALUES.F_NAN

          LOADCT, ct 
          TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
          TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
          !P.COLOR=255
          !P.BACKGROUND=254
          !P.CHARSIZE=2
          !P.CHARTHICK=2

          ;title='!6Ozone (t=' + STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
          form='(I2)' 

          z=o3_tmp
          IF iano EQ 1 THEN BEGIN
            ;clev=(FINDGEN(nclev))*0.9-30.0
            clev=(FINDGEN(nclev))-25
            clev[0]=-9999  &  clev[nclev-1]=9999
            z=o3_tmp-o3_clim[*,*,levidx,tidx]
            form='(I3)' 
          ENDIF

          MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,$
            POS=pos,title=title
          CONTOUR, z, lon, lat, /CELL_FILL, C_COLORS=barcol, LEVELS=clev,/OVER, $
                   XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos

          COLORBAR_2, clev, barcol, /COL, YS =0.77, XS = 0.015, $
            CHARSIZE=2, CHARTHICK=2,$
            ;levelind=findgen(100)*5,  $
            FORMAT=form, $
            ;UNIT='ppbv', $
            /IVE,/NOFIRST;,/NOFRAME
          ;COLORBAR_2, clev, barcol, YS =0.015, XS = 0.85, $
            ;CHARSIZE=2, CHARTHICK=2,$
            ;;levelind=findgen(100)*5,  $
            ;FORMAT=form, $
            ;;UNIT='ppbv', $
            ;/IVE,/NOFIRST;,/NOFRAME
        endif ; plot7_o3

    ; nwp wind
        uwnd=REFORM(uwnd0[*,*,levidx,tidx])
        vwnd=REFORM(vwnd0[*,*,levidx,tidx])
        tmp=WHERE(finite(uwnd) EQ 1 AND finite(vwnd) EQ 1, ntmp)
        tmp=FIX(ntmp*RANDOMU(seed,1000))
        uwnd_nwp=uwnd[grid_x, grid_y]
        vwnd_nwp=vwnd[grid_x, grid_y]
        
        dum = AMV_Convert_UV_Wvector(uwnd_nwp,vwnd_nwp,ws_nwp,wd_nwp )

    ; AMV wind

        ;=============================================================
        ; calculate local correlation over the map 
        ;=============================================================
        corr_map_xsize = 10 
        corr_map_ysize = 10
        xstep = 5
        ystep = 5
        ;limit = 10, 85, 50, 145
        corr_nx = (limit[3]-limit[1]-corr_map_xsize/2)/xstep
        corr_ny = (limit[2]-limit[0]-corr_map_ysize/2)/ystep

        corr_map = fltarr(corr_nx, corr_ny)
        corr_map[*] = !values.f_nan

        for ilon=0, corr_nx-1 do BEGIN
          for ilat=0, corr_ny-1 do BEGIN
            lowlon = limit[1] + ilon*xstep
            highlon = limit[1] + ilon*xstep + corr_map_xsize
            lowlat = limit[0] + ilat*ystep
            highlat = limit[0] + ilat*ystep + corr_map_ysize

            ingrid = where(lon[grid_x, grid_y] ge lowlon and $
              lon[grid_x, grid_y] le highlon and $
              lat[grid_x, grid_y] ge lowlat and $
              lat[grid_x, grid_y] le highlat, $
              ingridnum)
            if plot7_wnd_type eq 1 then begin
              _corr = correlate(ws_nwp[ingrid], ws[ingrid])
            endif ELSE BEGIN
              _corr = correlate(cos(wd_nwp[ingrid]*!dtor), cos(wd[ingrid]*!dtor))
            ENDELSE
            corr_map[ilon, ilat] = _corr
          ENDFOR; ilat
        ENDFOR; ilon
        print, corr_map
        makeXY, limit[1]+corr_map_xsize/2, limit[3]-corr_map_xsize/2, xstep, $
          limit[0]+corr_map_ysize/2, limit[2]-corr_map_ysize/2, ystep, $
          xlon, ylat
        ;; end calculate local correlation
        ;=============================================================

        clev = findgen(10)/10.
        nclev = n_elements(clev)
        clev[0]=-9999  &  clev[nclev-1]=9999
        nclev = n_elements(clev)
        barcol=FINDGEN(nclev-1)*253/(nclev-1)+1

        loadct_dg, 73
        MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,$
          POS=pos;,title=title

        CONTOUR, corr_map, xlon, ylat, $
          /cell_fill, $
          C_COLORS=barcol, $
          LEVELS=clev, /OVER, $
          XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos

        COLORBAR_2, clev, barcol, YS =0.015, XS = 0.015, /col, $
          CHARSIZE=2, CHARTHICK=2,$
          ;levelind=findgen(100)*5,  $
          FORMAT=form, $
          ;UNIT='ppbv', $
          /IVE,/NOFIRST;,/NOFRAME
        MAP_CONTINENTS,/CONTINENT,/COUNT,MLINETHICK=1.8,/COUNTRY,/HIRES 
        xtlen=20
        ytlen=10 
        print, clev

        IF plot7_wnd_on EQ 1 THEN BEGIN

      ; nwp wind
          uwnd=REFORM(uwnd0[*,*,levidx,tidx])
          vwnd=REFORM(vwnd0[*,*,levidx,tidx])
          tmp=WHERE(finite(uwnd) EQ 1 AND finite(vwnd) EQ 1, ntmp)
          tmp=FIX(ntmp*RANDOMU(seed,1000))
          uwnd_nwp=uwnd[grid_x, grid_Y]
          vwnd_nwp=vwnd[grid_x, grid_Y]
          dum = AMV_Convert_UV_Wvector(uwnd_nwp,vwnd_nwp,ws_nwp,wd_nwp )
          WINDBARB,lon[grid_x, grid_Y],lat[grid_x, grid_Y],$
            ws_nwp*1.94384,wd_nwp,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('CHARCOAL'), clip=pos

      ; AMV wind
          WINDBARB,grid_lon,grid_lat,$
            ws*1.94384,wd,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('MAGENTA'), clip=pos
          tmp=WHERE(ws*1.94384 LT 2.5,COMP=ttmp,ntmp)
          oplot, grid_lon[ttmp],grid_lat[ttmp], psym=plotsym_fn(/cir,/fill,scale=0.5)
          IF ntmp NE 0 THEN oplot, grid_lon[tmp],grid_lat[tmp], psym=plotsym_fn(/cir,scale=0.5)

        ENDIF ;plot7_wnd_on


        IF po_idx EQ 1 THEN BEGIN
          xcoord=[lon[po_lonidx-dx,po_latidx-dx],$
            lon[po_lonidx+dx,po_latidx-dx],$
            lon[po_lonidx+dx,po_latidx+dx],$
            lon[po_lonidx-dx,po_latidx+dx],$
            lon[po_lonidx-dx,po_latidx-dx]]
          ycoord=[lat[po_lonidx-dx,po_latidx-dx],$
            lat[po_lonidx+dx,po_latidx-dx],$
            lat[po_lonidx+dx,po_latidx+dx],$
            lat[po_lonidx-dx,po_latidx+dx],$
            lat[po_lonidx-dx,po_latidx-dx]]
          POLYFILL,xcoord, ycoord, COLOR=CGCOLOR('MAGENTA')
          PLOTS,xcoord, ycoord
        ENDIF ; po_idx
        z=o3_tmp
        ;ct = 72

        CONTOUR, z, z, z, /NODAT,/NOERA,CHARSIZE=2,CHARTHICK=2, $
                 XSTYLE=1, YSTYLE=1,YTICKFORMAT='(I3)',XTICKFORMAT='(i4)',XTICKINTERVAL=xtlen,YTICKINTERVAL=ytlen, $
                 XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]], $
                 TITLE=title, XTITLE = 'Longitude',YTITLE='Latitude', POS=pos


        cgPs_Close,density=800, /png
        ;file_trans_ds,out_png,/DEL
      ENDIF ; plot7

      IF plot8 EQ 1 THEN BEGIN
        tail='d'+date+'_t'+ STRMID(the_time,1,2,/REV) + '_l'+STRING(the_pres, F='(I04)')
        if plot8_wnd_type eq 1 then BEGIN
          typestr = 'ws_'
        endif else BEGIN
          typestr = 'wd_'
        endelse

        outname = typestr + 'corr_gchem_tomv_plot8_TS' + TS + '_SS' + SS
        out_ps = outname + tail+'.ps'
        out_png = outname + tail+'.png'
        ;CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=14,/NOMATCH,/INCH
        CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=7,/NOMATCH,/INCH
        pos=[0.12,0.12,0.9,0.91]
        !P.POSITION = pos

        if plot8_corr then begin
          ct=73
          nclev=51
          clev=(FINDGEN(nclev))/(nclev -1)/2+ 0.5
          print, clev
          barcol=REVERSE(FINDGEN(nclev-1)*253/(nclev-1)+1)

          clev[0]=-9999  &  clev[nclev-1]=9999

          ;znan=WHERE(z LT MEAN(z)-3*STDDEV(z) AND z GT mean(z)+3*STDDEV(z), zval)
          ;z[znan]=!VALUES.F_NAN

          LOADCT, ct 
          TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
          TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
          !P.COLOR=255
          !P.BACKGROUND=254
          !P.CHARSIZE=2
          !P.CHARTHICK=2

          ;title='!6Ozone (t=' + STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
          title='!6Correlation (d='+ date +', t='+ STRMID(the_time,1,2,/REV) + $
            'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
          form='(I2)' 

          z=o3_tmp
          IF iano EQ 1 THEN BEGIN
            ;clev=(FINDGEN(nclev))*0.9-30.0
            clev=(FINDGEN(nclev))-25
            clev[0]=-9999  &  clev[nclev-1]=9999
            z=o3_tmp-o3_clim[*,*,levidx,tidx]
            form='(I3)' 
          ENDIF

          MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0, $
            /CYLIN,/NOERA,/NOBOR,POS=pos
          CONTOUR, corrs_mean, grid_lon, grid_lat, $
            /CELL_FILL, $
            /IRREGULAR, $
            C_COLORS=barcol, LEVELS=clev,/OVER, $
            XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]], $
            POS=pos,title=title
          
          MAP_CONTINENTS,/CONTINENT,/COUNT,MLINETHICK=1.8,/COUNTRY,/HIRES 

        endif ; plot8_o3
        ;stop
        ;===========================================
        ; correlation plot
        ;===========================================
        ;lond = 1 
        ;latd = 1
        ;;corrs_new = sph_scat(grid_lon, grid_lat, corrs_mean, GS=[lond, latd])
        ;ct = colortable(73)
        ;m = map('Orthographic', $
          ;/buffer, $
          ;CENTER_LONGITUDE=115, CENTER_LATITUDE=30)
        ;levels = findgen(11)/10
        ;c = CONTOUR(corrs_mean, grid_lon, grid_lat, $
          ;/buffer, $
          ;;/fill, $
          ;;n_levels=20, $
          ;c_value=levels, $
          ;grid_units='degrees', $
          ;rgb_table=ct, $
          ;title='TOMV Correlation', $
          ;overplot=m) 
        ;mc = MAPcontinents()
        ;cb = colorbar(title = 'Mean Correlation', target=c, $
          ;orientation=1, $
          ;position=[0.89, 0.05, 0.93, 0.9], $
          ;tickinterval=0.025, $
          ;textpos=1)
        ;c.save, date+'_correlation.png'
        ;c.close



    ; nwp wind
        uwnd=REFORM(uwnd0[*,*,levidx,tidx])
        vwnd=REFORM(vwnd0[*,*,levidx,tidx])
        tmp=WHERE(finite(uwnd) EQ 1 AND finite(vwnd) EQ 1, ntmp)
        tmp=FIX(ntmp*RANDOMU(seed,1000))
        uwnd_nwp=uwnd[grid_x, grid_y]
        vwnd_nwp=vwnd[grid_x, grid_y]
        
        dum = AMV_Convert_UV_Wvector(uwnd_nwp,vwnd_nwp,ws_nwp,wd_nwp )

    ; AMV wind

        ;=============================================================
        ; calculate local correlation over the map 
        ;=============================================================
        ;corr_map_dlon = 10 
        ;corr_map_dlat = 10
        ;xstep = 2
        ;ystep = 2
        ;;limit = 10, 85, 50, 145
        ;corr_nx = (limit[3]-limit[1]-corr_map_dlon)/xstep
        ;corr_ny = (limit[2]-limit[0]-corr_map_dlat)/ystep

        ;corr_map = fltarr(corr_nx, corr_ny)
        ;corr_map[*] = !values.f_nan

        ;for ilon=0, corr_nx-1 do BEGIN
          ;for ilat=0, corr_ny-1 do BEGIN
            ;lowlon = limit[1] + ilon
            ;highlon = limit[1] + ilon + corr_map_dlon
            ;lowlat = limit[0] + ilat
            ;highlat = limit[0] + ilat + corr_map_dlat

            ;ingrid = where(lon[grid_x, grid_y] ge lowlon and $
              ;lon[grid_x, grid_y] le highlon and $
              ;lat[grid_x, grid_y] ge lowlat and $
              ;lat[grid_x, grid_y] le highlat, $
              ;ingridnum)
            ;if plot8_wnd_type eq 1 then begin
              ;_corr = correlate(ws_nwp[ingrid], ws[ingrid])
            ;endif ELSE BEGIN
              ;_corr = correlate(cos(wd_nwp[ingrid]*!dtor), cos(wd[ingrid]*!dtor))
            ;ENDELSE

            ;corr_map[ilon, ilat] = _corr
          ;ENDFOR; ilat
        ;ENDFOR; ilon
        ;;makexy, limit[1]+3, limit[3]-5+3, 5, limit[0]+3, limit[2]-5+3, 5, $
          ;;xlon, ylat
        ;makeXY, limit[1]+corr_map_dlon/2, limit[3]-corr_map_dlon/2-1, xstep, $
          ;limit[0]+corr_map_dlat/2, limit[2]-corr_map_dlat/2-1, ystep, $
          ;xlon, ylat
        ;;; end calculate local correlation
        ;clev_ = findgen(10)/10.
        ;nclev_ = n_elements(clev_)
        ;clev_[0]=-9999  &  clev_[nclev_-1]=9999
        ;nclev_ = n_elements(clev_)
        ;barcol=REVERSE(FINDGEN(nclev_-1)*253/(nclev_-1)+1)
        ;MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,POS=pos,title='!6'
        ;loadct_dg, 63

        ;;cellfill = 0
        ;;if plot8_o3 eq 0 then cellfill = 1

        ;CONTOUR, corr_map, xlon, ylat, /cell_fill, C_COLORS=barcol, LEVELS=clev_, /OVER, $
                 ;XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos

        ;MAP_CONTINENTS,/CONTINENT,/COUNT,MLINETHICK=1.8,/COUNTRY,/HIRES 
        ;xtlen=20
        ;ytlen=10 


        IF plot8_wnd_on EQ 1 THEN BEGIN

          ; nwp wind
          uwnd=REFORM(uwnd0[*,*,levidx,tidx])
          vwnd=REFORM(vwnd0[*,*,levidx,tidx])
          tmp=WHERE(finite(uwnd) EQ 1 AND finite(vwnd) EQ 1, ntmp)
          tmp=FIX(ntmp*RANDOMU(seed,1000))
          uwnd_nwp=uwnd[grid_x, grid_Y]
          vwnd_nwp=vwnd[grid_x, grid_Y]
          dum = AMV_Convert_UV_Wvector(uwnd_nwp,vwnd_nwp,ws_nwp,wd_nwp )
          WINDBARB,lon[grid_x, grid_Y],lat[grid_x, grid_Y],$
            ws_nwp*1.94384,wd_nwp,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('CHARCOAL'), clip=pos

          ; AMV wind
          WINDBARB,grid_lon,grid_lat,$
            ws*1.94384,wd,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('MAGENTA'), clip=pos
          tmp=WHERE(ws*1.94384 LT 2.5,COMP=ttmp,ntmp)
          oplot, grid_lon[ttmp],grid_lat[ttmp], psym=plotsym_fn(/cir,/fill,scale=0.5)
          IF ntmp NE 0 THEN oplot, grid_lon[tmp],grid_lat[tmp], psym=plotsym_fn(/cir,scale=0.5)

          ; Plot target position
              ;cgoplot, dest_lon,dest_lat, psym=1
              ;cgoplot, grid_lon,grid_lat, psym=1
        ENDIF ;plot8_wnd_on


        IF po_idx EQ 1 THEN BEGIN
          xcoord=[lon[po_lonidx-dx,po_latidx-dx],$
            lon[po_lonidx+dx,po_latidx-dx],$
            lon[po_lonidx+dx,po_latidx+dx],$
            lon[po_lonidx-dx,po_latidx+dx],$
            lon[po_lonidx-dx,po_latidx-dx]]
          ycoord=[lat[po_lonidx-dx,po_latidx-dx],$
            lat[po_lonidx+dx,po_latidx-dx],$
            lat[po_lonidx+dx,po_latidx+dx],$
            lat[po_lonidx-dx,po_latidx+dx],$
            lat[po_lonidx-dx,po_latidx-dx]]
          POLYFILL,xcoord, ycoord, COLOR=CGCOLOR('MAGENTA')
          PLOTS,xcoord, ycoord
        ENDIF ; po_idx
        z=o3_tmp
        ct = 73

        CONTOUR, z, z, z, /NODAT,/NOERA,CHARSIZE=2,CHARTHICK=2, $
                 XSTYLE=1, YSTYLE=1,YTICKFORMAT='(I3)',XTICKFORMAT='(i4)',$
                 XTICKINTERVAL=xtlen,YTICKINTERVAL=ytlen, $
                 XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]], $
                 TITLE=title, XTITLE = 'Longitude',YTITLE='Latitude', POS=pos

        loadct_dg,ct
        form = '(f4.2)'
        COLORBAR_2, clev, barcol, /COL, YS =0.77, XS = 0.015, CHARSIZE=2, $
          CHARTHICK=2,levelind=findgen(100)*5,  $
          FORMAT=form,UNIT='ppbv',/IVE,/NOFIRST;,/NOFRAME

        cgPs_Close,density=800, /png
        ;file_trans_ds,out_png,/DEL
      ENDIF ; plot8

      IF plot9 EQ 1 THEN BEGIN
        tail='d'+date+'_t'+ STRMID(the_time,1,2,/REV) + '_l'+STRING(the_pres, F='(I04)')

        typestr = 'winds_'

        outname= typestr + 'corr_gchem_tomv_plot9_TS' + TS + '_SS' + SS
        out_ps=outname+tail+'.ps'
        out_png=outname+tail+'.png'
        ;CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=14,/NOMATCH,/INCH
        CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=7,/NOMATCH,/INCH
        pos=[0.12,0.12,0.9,0.91]
        !P.POSITION = pos

        title='!6Ozone (d='+ date +', t='+ STRMID(the_time,1,2,/REV) + $
          'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
        if plot9_o3 then begin
          ct=73
          nclev=51
          clev=(FINDGEN(nclev))/(nclev -1)
          print, clev
          barcol=REVERSE(FINDGEN(nclev-1)*253/(nclev-1)+1)

          clev[0]=-9999  &  clev[nclev-1]=9999

          ;znan=WHERE(z LT MEAN(z)-3*STDDEV(z) AND z GT mean(z)+3*STDDEV(z), zval)
          ;z[znan]=!VALUES.F_NAN

          LOADCT, ct 
          TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
          TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
          !P.COLOR=255
          !P.BACKGROUND=254
          !P.CHARSIZE=2
          !P.CHARTHICK=2

          ;title='!6Ozone (t=' + STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
          form='(I2)' 

          z=o3_tmp
          IF iano EQ 1 THEN BEGIN
            ;clev=(FINDGEN(nclev))*0.9-30.0
            clev=(FINDGEN(nclev))-25
            clev[0]=-9999  &  clev[nclev-1]=9999
            z=o3_tmp-o3_clim[*,*,levidx,tidx]
            form='(I3)' 
          ENDIF

          MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,$
            POS=pos,title=title
          CONTOUR, z, lon, lat, /CELL_FILL, C_COLORS=barcol, LEVELS=clev,/OVER, $
                   XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos

          COLORBAR_2, clev, barcol, /COL, YS =0.77, XS = 0.015, $
            CHARSIZE=2, CHARTHICK=2,$
            ;levelind=findgen(100)*5,  $
            FORMAT=form, $
            ;UNIT='ppbv', $
            /IVE,/NOFIRST;,/NOFRAME
          ;COLORBAR_2, clev, barcol, YS =0.015, XS = 0.85, $
            ;CHARSIZE=2, CHARTHICK=2,$
            ;;levelind=findgen(100)*5,  $
            ;FORMAT=form, $
            ;;UNIT='ppbv', $
            ;/IVE,/NOFIRST;,/NOFRAME
        endif ; plot9_o3

    ; nwp wind
        uwnd=REFORM(uwnd0[*,*,levidx,tidx])
        vwnd=REFORM(vwnd0[*,*,levidx,tidx])
        tmp=WHERE(finite(uwnd) EQ 1 AND finite(vwnd) EQ 1, ntmp)
        tmp=FIX(ntmp*RANDOMU(seed,1000))
        uwnd_nwp=uwnd[grid_x, grid_y]
        vwnd_nwp=vwnd[grid_x, grid_y]
        
        dum = AMV_Convert_UV_Wvector(uwnd_nwp,vwnd_nwp,ws_nwp,wd_nwp )

    ; AMV wind

        ;=============================================================
        ; calculate local correlation over the map 
        ;=============================================================
        corr_map_xsize = 10 
        corr_map_ysize = 10
        xstep = 5
        ystep = 5
        ;limit = 10, 85, 50, 145
        corr_nx = (limit[3]-limit[1]-corr_map_xsize/2)/xstep
        corr_ny = (limit[2]-limit[0]-corr_map_ysize/2)/ystep

        corr_map = fltarr(corr_nx, corr_ny)
        corr_map[*] = !values.f_nan

        for ilon=0, corr_nx-1 do BEGIN
          for ilat=0, corr_ny-1 do BEGIN
            lowlon = limit[1] + ilon*xstep
            highlon = limit[1] + ilon*xstep + corr_map_xsize
            lowlat = limit[0] + ilat*ystep
            highlat = limit[0] + ilat*ystep + corr_map_ysize

            ingrid = where(lon[grid_x, grid_y] ge lowlon and $
              lon[grid_x, grid_y] le highlon and $
              lat[grid_x, grid_y] ge lowlat and $
              lat[grid_x, grid_y] le highlat, $
              ingridnum)
            _corr = (correlate(ws_nwp[ingrid], ws[ingrid]) * $
              correlate(cos(wd_nwp[ingrid]*!dtor), cos(wd[ingrid]*!dtor)))
            corr_map[ilon, ilat] = _corr
          ENDFOR; ilat
        ENDFOR; ilon
        print, corr_map
        makeXY, limit[1]+corr_map_xsize/2, limit[3]-corr_map_xsize/2, xstep, $
          limit[0]+corr_map_ysize/2, limit[2]-corr_map_ysize/2, ystep, $
          xlon, ylat
        ;; end calculate local correlation
        ;=============================================================

        clev = findgen(10)/10.
        nclev = n_elements(clev)
        clev[0]=-9999  &  clev[nclev-1]=9999
        nclev = n_elements(clev)
        barcol=FINDGEN(nclev-1)*253/(nclev-1)+1

        loadct_dg, 73
        MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,$
          POS=pos;,title=title

        CONTOUR, corr_map, xlon, ylat, $
          /cell_fill, $
          C_COLORS=barcol, $
          LEVELS=clev, /OVER, $
          XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos

        COLORBAR_2, clev, barcol, YS =0.015, XS = 0.015, /col, $
          CHARSIZE=2, CHARTHICK=2,$
          ;levelind=findgen(100)*5,  $
          FORMAT=form, $
          ;UNIT='ppbv', $
          /IVE,/NOFIRST;,/NOFRAME
        MAP_CONTINENTS,/CONTINENT,/COUNT,MLINETHICK=1.8,/COUNTRY,/HIRES 
        xtlen=20
        ytlen=10 
        print, clev

        IF plot9_wnd_on EQ 1 THEN BEGIN

      ; nwp wind
          uwnd=REFORM(uwnd0[*,*,levidx,tidx])
          vwnd=REFORM(vwnd0[*,*,levidx,tidx])
          tmp=WHERE(finite(uwnd) EQ 1 AND finite(vwnd) EQ 1, ntmp)
          tmp=FIX(ntmp*RANDOMU(seed,1000))
          uwnd_nwp=uwnd[grid_x, grid_Y]
          vwnd_nwp=vwnd[grid_x, grid_Y]
          dum = AMV_Convert_UV_Wvector(uwnd_nwp,vwnd_nwp,ws_nwp,wd_nwp )
          WINDBARB,lon[grid_x, grid_Y],lat[grid_x, grid_Y],$
            ws_nwp*1.94384,wd_nwp,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('CHARCOAL'), clip=pos

      ; AMV wind
          WINDBARB,grid_lon,grid_lat,$
            ws*1.94384,wd,ASPECT=1,LENGTH=0.02,THICK=3,COL=CGCOLOR('MAGENTA'), clip=pos
          tmp=WHERE(ws*1.94384 LT 2.5,COMP=ttmp,ntmp)
          oplot, grid_lon[ttmp],grid_lat[ttmp], psym=plotsym_fn(/cir,/fill,scale=0.5)
          IF ntmp NE 0 THEN oplot, grid_lon[tmp],grid_lat[tmp], psym=plotsym_fn(/cir,scale=0.5)

        ENDIF ;plot9_wnd_on


        IF po_idx EQ 1 THEN BEGIN
          xcoord=[lon[po_lonidx-dx,po_latidx-dx],$
            lon[po_lonidx+dx,po_latidx-dx],$
            lon[po_lonidx+dx,po_latidx+dx],$
            lon[po_lonidx-dx,po_latidx+dx],$
            lon[po_lonidx-dx,po_latidx-dx]]
          ycoord=[lat[po_lonidx-dx,po_latidx-dx],$
            lat[po_lonidx+dx,po_latidx-dx],$
            lat[po_lonidx+dx,po_latidx+dx],$
            lat[po_lonidx-dx,po_latidx+dx],$
            lat[po_lonidx-dx,po_latidx-dx]]
          POLYFILL,xcoord, ycoord, COLOR=CGCOLOR('MAGENTA')
          PLOTS,xcoord, ycoord
        ENDIF ; po_idx
        z=o3_tmp
        ;ct = 72

        CONTOUR, z, z, z, /NODAT,/NOERA,CHARSIZE=2,CHARTHICK=2, $
                 XSTYLE=1, YSTYLE=1,YTICKFORMAT='(I3)',XTICKFORMAT='(i4)',XTICKINTERVAL=xtlen,YTICKINTERVAL=ytlen, $
                 XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]], $
                 TITLE=title, XTITLE = 'Longitude',YTITLE='Latitude', POS=pos


        cgPs_Close,density=800, /png
        ;file_trans_ds,out_png,/DEL
      ENDIF ; plot9


    ENDFOR ; tidx
    ENDFOR ; levidx
  ENDFOR ; ifile


  IF plot3 EQ 1 THEN BEGIN
    tail='t'+ STRMID(the_time,1,2,/REV) + '_l'+STRING(the_pres, F='(I04)')
    if wnd_type eq 2 then typestr = 'wd_'
    if wnd_type eq 1 then typestr = 'ws_'

    outname='gchem_cross_correlation_image_' + typestr + '_TS' + TS + '_SS' + SS
    out_ps=outname+tail+'.ps'
    out_png=outname+tail+'.png'
    pos=[0.12,0.12,0.9,0.91]
    !P.POSITION = pos

    ct=20
    nclev=51
    clev=(findgen(nclev))*0.01   ; total column analysis
    barcol=FIX(FINDGEN(nclev-1)*253/nclev)
    clev[0]=-999  &  clev[nclev-1]=999

    ;znan=WHERE(z LT MEAN(z)-3*STDDEV(z) AND z GT mean(z)+3*STDDEV(z), zval)
    ;z[znan]=!VALUES.F_NAN

    LOADCT, ct 
    TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
    TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
    !P.COLOR=255
    !P.BACKGROUND=254
    !P.CHARSIZE=2
    !P.CHARTHICK=2

    title='!6Ozone (t=' + STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
    ;YONFORM_SDP_V2, ws_model, ws_amv, 'Model Wind','AMV Wind', minv, maxv, 'wind_cc.ps', hist_max=1000

    hist_max = 10.0
    hist_min = 0.0

    IF wnd_type EQ 1 THEN BEGIN
      ; wind speed scatterplot
      minv=0  &  maxv=25
      x_para = ws_models         & X_name = 'Model Wind'
      y_para = ws_amvs           & Y_name = 'AMV Wind'
      tail = tail + '_ws'
    ENDIF ELSE BEGIN
      ; wind direction scatterplot
      minv=0  &  maxv=360
      x_para = wd_models         & X_name = 'Model Wind'
      y_para = wd_amvs           & Y_name = 'AMV Wind'
      tail = tail + '_wd'
    ENDELSE ; wnd_type

    delta = (maxv-minv)/80.              ;; set histogram interval of x, y axis

    cgPS_Open, out_ps
    cgDisplay, aspect=1
    LoadCT, 20, Ncolors=254, Bottom=1

    position = [0.16, 0.15, 0.83, 0.9]
    thick = (!D.Name EQ 'PS') ? 6 :3

    cgPlot, x_para, y_para, /NoData, $
          XTitle =  X_name, YTitle = Y_name, $
          XRange = [minv, maxv], YRange = [minv, maxv], $
          Position=position, Charsize=2.0,XThick=3.0, YThick=3.0, Thick=6.0, /NoErase

    linfit_xdata=findgen(200)-100
    back_eq1 = linfit_xdata*1.
    oplot, linfit_xdata, back_eq1, color=0, linestyle=1
    ;============================
    ; plot histogram
    ;============================

    ;for j = minv, maxv, delta  do begin
      ;for i = minv, maxv, delta  do begin
    for jtmp = minv, (maxv-minv)/delta  do begin
      for itmp = minv, (maxv-minv)/delta  do begin
        i=minv+(itmp-minv)*delta
        j=minv+(jtmp-minv)*delta
        idx_histo= where( y_para ge j and y_para lt j+delta and x_para ge i and x_para lt i+delta, idx_histo_num)

        xbox = [i, i+delta, i+delta, i]
        ybox = [j, j, j+delta, j+delta]
        if idx_histo_num gt 0.0 then begin
          ;print, i, j, idx_histo_num
          polyfill, xbox, ybox, color = bytscl(idx_histo_num, min = hist_min , max = hist_max, top = 254)    ;; set color bar min / max
        endif
      endfor
    endfor

    LoadCT, 20, Ncolors=254, Bottom=1
    TVLCT, CGColor('BLACK', /Triple), 0 ; Drawing color.
    cgColorbar,format='(i4)',range=[hist_min,hist_max], position=[0.84, 0.15, 0.88, 0.9],color=0, /vertical,/right, $     
                          divisions=5,bottom=1,ncolor=254,minor=1, Charsize=2


    ; static inforamtion
    nt=N_ELEMENTS(wd_amvs)
    xy_finite_idx = where(finite(x_para) and finite(y_para), finite_num) 
    x = x_para[xy_finite_idx]
    y = y_para[xy_finite_idx]
    x1 = fltarr(1, finite_num)
    x1[0,*] = x

    if wnd_type eq 1 then begin
      slope = regress(x1, y, const=const, sigma=_sigma, /relative_weight)
      cgtext, 5, 23, 'R='+strcompress(string(correlate(x,y), format='(f7.4)'), $
        /remove_all), /data
    endif else begin
      slope = regress(cos(x1*!dtor), cos(y*!dtor), const=const, sigma=_sigma, /relative_weight)
      cgtext, 0.5, 0.1, 'R='+strcompress(string(correlate(x,y), format='(f7.4)'), $
        /remove_all), /data
    endelse

    NofS  = 'N= ' + strcompress(string(nt,format='(i6)'), /remove_all)
    R     = 'R= ' + strcompress(string(correlate(x,y), format='(f7.4)'), /remove_all)
    rmse  = 'RMSE= '+ strcompress(string(sqrt(mean(( slope[0]*x + const[0] - y)^2.0)), format='(f6.2)'), /remove_all)
    reg   = 'Y= ' + strcompress(string(slope,format='(f6.2)')+'X + ' +string(const,format='(f6.2)'), /remove_all)
    ;bias  = 'Abs.GEMS-SON = '   + strcompress(string(avg1,format='(f6.2)')+ '!9+!X' + string(std1,format='(f6.2, "(DU)")'),/remove_all)
    ;rbias = '          = ' + strcompress(string(avg2,format='(f6.2)')+ '!9+!X' + string(std2,format='(f6.2, "(%)")'), /remove_all)
    print, nofs
    print, r
    print, reg
    print, rmse

    cgPS_Close, density=1000,/png
    ;file_trans_ds,out_png,/DEL


  ENDIF  ; plot3


  IF plot4 EQ 1 THEN BEGIN
    ;save,file='tmp.xdr',o3_sta, t_sta,/xdr
    ;RESTORE, 'tmp.xdr'
    
    tail='t'+ STRMID(the_time,1,2,/REV) + '_l'+STRING(the_pres, F='(I04)')
    outname='gchem_sdp_plot4_' + typestr + '_TS' + TS + '_SS' + SS
    out_ps=outname+tail+'.ps'
    out_png=outname+tail+'.png'
    CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=7,/NOMATCH,/INCH
    pos=[0.12,0.12,0.9,0.91]
    !P.POSITION = pos

    LOADCT, 33 
    TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
    TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
    !P.COLOR=255
    !P.BACKGROUND=254
    !P.CHARSIZE=2
    !P.CHARTHICK=2

    ;title='!6Ozone (t=' + STRMID(the_time,1,2,/REV) + 'utc, lev='+STRING(the_pres, F='(I4)')+'hPa)'
    title='Ozone time series'

    ntime=N_ELEMENTS(o3_sta)
    x=findgen(ntime)

    IF kst_on EQ 1 THEN BEGIN
      tjul_sta=tjul_sta+0.04166667*9
      jtime_arr=TIMEGEN(ntime,START=tjul_sta,UNITS='hours',STEP_SIZE=1)
      XTITLE='Time [KST]' 
    ENDIF ELSE BEGIN
      jtime_arr=TIMEGEN(ntime,START=tjul_sta[0],UNITS='hours',STEP_SIZE=1)
      XTITLE='Time [UTC]' 
    ENDELSE
    CALDAT,jtime_arr,mon,day,year,hour
    ;tc_mdh=STRING(mon,F='(I02)') $
          ;+STRING(day,F='(I02)')+'_'+STRING(hour,F='(I02)')+'00'
    ;tc_mdh=STRING(day,F='(I02)')+STRING(hour,F='(I02)')
    ;tc_mdh=STRING(mon,F='(I02)')+STRING(day,F='(I02)')
    tc_mdh=STRING(hour,F='(I02)')

    ;tname_int=72
    tname_int=3
    tname_start=0   ; tname_start+1 utc
    PLOT,[0,0],[0.0],TITLE=title $
        ,XRANGE=minmax(tjul_sta),YRANGE=[0,100] $
        ,THICK=5,XSTYLE=1,YSTYLE=1,/NODATA    $
        ,YTITLE='O3 [ppbv]' $
        ,XTITLE=xtitle $
        ,XTICKV=jtime_arr[tname_start:*:tname_int],XTICKNAME=tc_mdh[tname_start:*:tname_int] $
        ,XMINOR=4,XTICKS=N_ELEMENTS(jtime_arr[tname_start:*:tname_int])-1
    OPLOT,tjul_sta,o3_sta,THICK=5
    OPLOT,tjul_sta,o3_sta,PSYM=plotsym_fn(/box,scale=1,/fill)

    IF inc_clim EQ 1 THEN BEGIN
      OPLOT,tjul_sta,o3_sta_clim,THICK=5,COLOR=CGCOLOR('RED')
      OPLOT,tjul_sta,o3_sta_clim,PSYM=plotsym_fn(/cir,scale=0.5),COLOR=CGCOLOR('RED')
    ENDIF


    ;PLOT, x, o3_sta, xrange=minmax(x), ytitle='O3 [ppbv]', xstyle=1

    cgPS_Close, density=1000,/png
    ;file_trans_ds,out_png,/DEL

  ENDIF ; plot4


  IF plot5 EQ 1 THEN BEGIN
    
    o3_clim2=MEAN(o3_clim,dim=5)

    ;save, file='o3_clim2.xdr', o3_clim2,/xdr
    ;print, '  :: Save "o3_clim2.xdr"'

    FOR tloop=stidx,etidx DO BEGIN
   
      tail='o3_climatology_t'+ STRMID(time[tloop],1,2,/REV)  +'_l'+STRING(pres[sellev], F='(I04)')
                           
      outname='gchem_comp_TS' + TS + '_SS' + SS
      out_ps=outname+tail+'.ps'
      out_png=outname+tail+'.png'
      CGPS_OPEN, out_ps ,XSIZE=10,YSIZE=7,/NOMATCH,/INCH
      pos=[0.12,0.12,0.9,0.91]
      !P.POSITION = pos


      ct=72
      nclev=51
      clev=(FINDGEN(nclev))*2
      barcol=REVERSE(FINDGEN(nclev-1)*253/(nclev-1)+1)

      clev[0]=-9999  &  clev[nclev-1]=9999

      LOADCT, ct 
      TVLCT, CGColor('WHITE', /Triple), 254 ; Background color.
      TVLCT, CGColor('BLACK', /Triple), 255 ; Drawing color.
      !P.COLOR=255
      !P.BACKGROUND=254
      !P.CHARSIZE=2
      !P.CHARTHICK=2

      title='!6Ozone (t=' + STRMID(time[tloop],1,2,/REV) + 'utc, lev='+STRING(pres[sellev], F='(I4)')+'hPa)'
      ;title='!6Ozone (lev='+STRING(pres[sellev], F='(I4)')+'hPa)'

      ;z=MEAN((REFORM(o3_clim[*,*,sellev,tloop,*])),dim=3)
      z=REFORM(o3_clim2[*,*,sellev,tloop])
      MAP_SET,  LIMIT=limit, XMARGIN=0, YMARGIN=0,/CYLIN,/NOERA,/NOBOR,POS=pos,title='!6'
      CONTOUR, z, lon, lat, /CELL_FILL, C_COLORS=barcol, LEVELS=clev,/OVER, $
               XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]],POS=pos

      MAP_CONTINENTS,/CONTINENT,/COUNT,MLINETHICK=1.8,/COUNTRY,/HIRES 
      xtlen=20
      ytlen=10 


      CONTOUR, z, z, z, /NODAT,/NOERA,CHARSIZE=2,CHARTHICK=2, $
               XSTYLE=1, YSTYLE=1,YTICKFORMAT='(I3)',XTICKFORMAT='(i4)',XTICKINTERVAL=xtlen,YTICKINTERVAL=ytlen, $
               XRANGE=[limit[1], limit[3]], YRANGE=[limit[0], limit[2]], $
               TITLE=title, XTITLE = 'Longitude',YTITLE='Latitude', POS=pos

      loadct_dg,ct
      COLORBAR_2, clev, barcol, /COL, YS =0.77, XS = 0.015, CHARSIZE=2, CHARTHICK=2,levelind=findgen(100)*5,  $
                               FORMAT='(I2)',UNIT='ppbv',/IVE,/NOFIRST;,/NOFRAME

      cgPs_Close,density=800, /png
      ;file_trans_ds,out_png,/DEL

    ENDFOR ; tloop
  ENDIF  ;plot5



END

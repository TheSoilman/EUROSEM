C-----------------------------------------------------------------------
C  EUROSEM -- European Soil Erosion Model, version 3.2L (metric, Lahey)
C  Event-based model of runoff and soil erosion for small catchments,
C  derived from KINEROS.
C
C  CATCHMENT REPRESENTATION
C    The catchment is a cascade of up to 60 ELEMENTS (NELE in use):
C    overland-flow planes, channels and ponds.  Each element is solved
C    separately, in the order given by NLOG(1..NELE).  Outflow from an
C    upstream element is stored (arrays QS and SCON, offset NB) and
C    up by the element(s) below it.  The LAST element in the computation
C    order is taken as the basin outlet.
C
C  PROGRAM FLOW (this MAIN program)
C    1. Banner, run title, and screen-graph question  (interactive)
C    2. INTSEL   choose / confirm file names, open all files
C    3. READER   read parameter + rainfall files, defaults, checks
C    4. Loop over elements in computation order:
C          OUTPUT (echo parameters), then POND, PLANE or CHANNL
C    5. Basin hydrograph totals and the global water balance
C    6. STOP 'NORMAL COMPLETION'
C
C  FILES (Fortran unit numbers set in MAIN)
C    15  IREAD   parameter input file         (read by READER/PAREAD)
C     7  JREAD   rainfall input file          (read by READER)
C    16  IWRITE  dynamic output (hydrographs, per-element results)
C    99          static output  (summaries, global balance)
C     8  IDIAGN  diagnostic / auxiliary output (messages, errors)
C     4  IPOND   pond input file, if any      (see PONDRD)
C    12          observed data for graph comparison, optional (INTSEL)
C     9          EUROM.FIL, remembers file names between runs (INTSEL)
C
C  UNITS
C    Metric working units: m, kg, s.  Rainfall is in mm/h and times are
C    in minutes in the input files; CONVER converts to internal units.
C
C  ROUTINES IN THIS FILE
C    MAIN, HYDwri, INTSEL, AUXTXT, READER, PAREAD, READTM, BLOCK DATA,
C    DEFAUL, PONDRD, INSPEC, ERRPO, CONVER, OUTPUT
C-----------------------------------------------------------------------
C*********************************************************************
C
      PROGRAM MAIN
C              EUROSEM file set:  1 of 7
C
C************************************************************************
C    Metric working units: M,KG,SEC: RAIN IN MM./H
C------------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C  COMMON BLOCKS
C    The same declarations are repeated in every routine; they must
C    exactly everywhere, or storage is silently misinterpreted.
C      /IO/      Fortran unit numbers (see list above)
C      /RILL1/   rill dimensions along the slope, NK nodes
C      /CNTRL/   run control: NELE (no. of elements), NLOG (computation
C                order), DELT (time step), NTIME (output time unit),
C                NEROS (erosion option), NPN (pond flag per element)
C      /GEOM/    element geometry and links: XL length, W width, S
C                slope,
C                NU/NR/NL/NC1/NC2 upstream-element links, NP print flag,
C                SUMA contributing area, rill and interrill slope data
C      /EVENT/   rainfall event: TFIN (run length), rainfall time-depth
C                data (TI, QI for the current element; TIDD, QIDD per
C                gage), RNDPTH net rainfall depth per element
C      /PLANE1/  flow solution: T time, Q outflow, QL lateral inflow,
C                QS/NB storage of element outflow for downstream use
C      /CHAN/    channel flow arrays (upstream-boundary flow, area)
C      /INFIL/   infiltration parameters (soil moisture, capillary
C                drive)
C      /SED/     sediment parameters and concentrations (D50, RHOS,
C                CONC, SCON, SPLTEX, COHE, ...)
C      /BALANC/  per-element water-balance totals (STORA storage,
C                EINF infiltration volume)
C      /ESEM/    EUROSEM vegetation / interception parameters
C      /labels/, /pldat/   graph labels and observed data (graph option)
C-----------------------------------------------------------------------
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      COMMON /CNTRL/ NRES, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), rfr(60), DEPNO(60)
     4       ,rs(60), dero(60), sir(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /CHAN/ A1(20),A2(20), QUB(1000), AUB(1000), CO1, CO2,B,NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /INFIL/ AL(60), THI(60), THMX(60), ROC(60), RECS(60),
     1        THR, THFC, tempw, stone(60), RVAL
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60), ITEX(60), COHE(60)
      COMMON /BALANC/ STORA(60),EINF(60),CINF2(20), STODEP(60)
      COMMON /ESEM/ COVER(60), CUINT(60), ISHAPE(60), PLHGT(60),
     1        VPANG(60),PBASE(60)
      common /labels/ title, tlabl, qlabl, slabl, ifsed
      Common /pldat/ ifad,ndat,tdat(500),qdat(500),sdat(500),smromm,
     1                sumskg
C
      character*6 tlabl(2), tuns(3)
      Character*12 qlabl, slabl
      character*80 title
C
      CHARACTER*1 ANS,DASH
C      CHARACTER*15 NAME
C
      LOGICAL IFSED, IFAD, ifgrf
      LOGICAL LHEAD
      CHARACTER*1 HANS
      CHARACTER*80 HTITLE,HFIL(5)
      COMMON /HLESSL/ LHEAD
      COMMON /HLESSC/ HANS,HTITLE,HFIL
C
      DIMENSION QD(1000),RGR(1000),TGR(1000),SGR(1000),TGQ(1000)
C
C  NTELE is the maximum number of elements (the dimension 60 in /GEOM/).
      DATA NTELE /60/
      DATA DASH /'-'/,tuns/' SEC. ',' MIN. ',' HRS. '/
c
C                       ASK USER
C  SECTION 1: banner, run title, graph question  (interactive: waits for
C  keyboard input at each READ(*,...) below)
C  Headless run: title comes from eurosem.ctl, no questions asked.
      CALL HEADLS
      IF(LHEAD) THEN
        TITLE = HTITLE
        IFGRF = .FALSE.
      ELSE
      WRITE(*,136)
136   FORMAT(//,6X,
     - '**************************************************',/,6X,
     - '    EEEE U  U RRRR OOOO SSSS EEEE M     M   ',/,6X,
     - '    E    U  U R  R O  O S    E    MM   MM   ',/,6X,
     - '    EE   U  U RRR  O  O SSSS EEEE M M M M   ',/,6X,
     - '    E    U  U R R  O  O    S E    M  M  M   ',/,6X,
     - '    EEEE UUUU R  R OOOO SSSS EEEE M     M   ',/,6X,
     - '--------------------------------------------------',/,6X,
     - '            RUNNING WITH',/,6X,
     - '--------------------------------------------------',//,6X,
     - ' PROGRAM KINEROS/  Metric Lahey VERSION OF 11/97  ',/,6X,
     - '--------------------------------------------------',//,6X,
     - ' VERSION 3.2L, 11-97, for Lahey LISK graph library',/,6X,
     - '                 USE WITH CARE!!!',//,6X,
     - '--------------------------------------------------',//,6X,
     - '  Press Carriage Return to Continue:   ')
      read(*,*)
      write(*,137)
137   FORMAT(//,6X,
     - '--------------------------------------------------',//,6X,
     - '   PLEASE REPORT ALL BUGS/PROBLEMS TO:',//,6X,
     - '          J.N.QUINTON',//,6X,
     - '          SILSOE COLLEGE',/,6X,
     - '          SILSOE',/,6X,
     - '          BEDFORD MK45 4DT',/,6X,
     - '          UNITED KINGDOM',//,6X,
     - '          TEL 0525-860428',/,6X,
     - '          FAX 0525-861527',/,6X,
     - '--------------------------------------------------',//,6X,
     - ' Enter a 1 to 80 char. title for the output file: ',/,6X,
     - '->')
C  Run title, up to 80 characters; stored in TITLE (COMMON /labels/).
	READ(*,277) TITLE
277     FORMAT(A80)
C
C  Screen-graph question.  In this version the graphics call (GRAFXY,
C  further down) is commented out, so IFGRF mainly controls whether a
C  file of observed data is requested in INTSEL and whether plot arrays
C  are prepared.
      write(*,*) 'Do you want screen graph of output hydrograph?'
      read(*,610) ans
  610 format(A1)
      Ifgrf = .false.
      if(ans.eq. 'Y' .or. ans .eq. 'y') Ifgrf = .true.
      END IF
C
C  Fortran unit numbers (see list at top of program).
      ans = 'n'
      IVOUT=1
      IPOND  = 4
      IREAD  = 15
      IWRITE = 16    !  dynamic output file.  unit 99 = static output
      JREAD  = 7
      IDIAGN = 8
C  Choose or confirm file names, and open all input/output files.
      CALL INTSEL(Ifgrf)                                                ! INTSEL
C**      AREA=0.0
C  Read the parameter file and the rainfall file; apply defaults,
C  convert units and check the input.
      CALL READER                                                       ! READER
C**      TODPTH=0.0
C     WRITE OUT PROCESSING HEADER
C  Column headings for the per-element volume-balance lines that follow.
      WRITE (IWRITE,95)
   95 FORMAT(//,
     -       ' ',9X,'     ',5X,'    ',7X,'VOL. BAL. ',7X,'SED. TOTAL',/,
     -       ' ',9X,'ELE #',5X,'TYPE',7X,' ERROR  % ',7X,'  (KGS.)  ',/,
     -       ' ',9X,'-----',5X,'----',7X,'----------',7X,'----------')
C
C  Zero the per-element storage and infiltration totals.
      DO 100 I=1,NTELE
        STORA(I) = 0.0
        EINF(I)  = 0.0
  100 CONTINUE
      AREA = 0.0
      TODPTH =0.0
C*********************************** MAIN LOOP PROCESSING ELEMENTS
C-----------------------------------------------------------------------
C  ELEMENT LOOP.  I = position in computation order, J = element number.
C    AREA   accumulates the total basin area (m2): W*XL of each element.
C    OUTPUT echoes the element's parameters.  Then, by element type:
C       NPN(J) >= 1                 pond     -> POND
C       W(J) > 0 and NPN(J) < 1     plane    -> PLANE
C       W(J) = 0 and NPN(J) < 1     channel  -> CHANNL
C    (a channel has no plane width).
C    TODPTH accumulates rainfall volume (net depth * element area).
C-----------------------------------------------------------------------
      DO 50 I=1,NELE
       J=NLOG(I)
       AREA = AREA  + W(J)*XL(J)
C       
       CALL OUTPUT(J)
       IF(NPN(J).GE.1) THEN
         CALL POND(J)
       ENDIF
       IF(W(J).NE.0 .AND. NPN(J).LT.1) THEN
C
         CALL PLANE(J)
       ENDIF
       IF(W(J).EQ.0 .AND. NPN(J).LT.1) THEN
C
         CALL CHANNL(J)
       ENDIF
C
       DPTHAV = RNDPTH(J)*W(J)*XL(J)
       TODPTH = DPTHAV+TODPTH
C  HYDROGRAPH PRINT OUTPUT LOGIC VIA SUBROUTINE: HYDwri
C  (HYETOGRAPH PRINT OUTPUT LOGIC IN SUBROUTINE PLANE)
C  Print the hydrograph for this element if NP(J) > 0 (HYDwri).
C  NRP(J) is 3 for planes and channels (set in PAREAD) and comes from
C  the pond file for ponds, so this block also runs for the last element
C  regardless.
       IF(NRP(J) .GT. 1 .OR. I .GE. NELE) THEN
         If(NP(j) .gt. 0) CALL HYDwri(J)
         WRITE(IWRITE,98) (DASH,K=1,75)
         IF(I .LT. NELE) THEN
C  THEN REWRITE THE  PROCESSING HEADER
           WRITE (IWRITE,95)
         ENDIF
       ENDIF
       WRITE(IDIAGN,98) (DASH,K=1,75)
   98  FORMAT(75A1/)
C
   50 CONTINUE
C
C*********************************** END MAIN LOOP TO PROCESS ELEMENTS
      WRITE(IDIAGN,235)
235   FORMAT(//)
C
C     SUM THE OUTFLOW
C-----------------------------------------------------------------------
C  BASIN OUTFLOW.  At this point Q(1..NI) still holds the outflow
C  hydrograph of the LAST element processed, i.e. NLOG(NELE), so that
C  element must be the basin outlet.  Q is in m3/s.
C    SUMRO  total outflow volume (m3)
C    QD     outflow as depth rate over the whole basin (mm/h);
C           3.6E6 = 1000 mm/m * 3600 s/h
C-----------------------------------------------------------------------
      SUMRO = 0.
      DO 75 L=1,NI
        SUMRO = SUMRO + Q(L)*delt
        QD(L)=Q(L)*3.6e6/AREA
   75 CONTINUE
C  CONVER(2): restore the rainfall arrays TIDD and QIDD to their input
C  units (time in the NTIME unit, intensity in mm/h) for the
C  hyetograph.
      CALL CONVER(2)                                                    *CONVERT
C
C      OUTPUT TO ASCII FILE FOR PLOTS
C  SUMRMM  total rainfall depth over the basin (mm, area-weighted)
C  SMROMM  total runoff depth over the basin (mm)
      SUMRmm =TODPTH/AREA*1000.
      SuMROm = SUMRO/AREA
      smromm = sumrom*1000. 
      WRITE(IWRITE,127) sumrmm
C
C     * IF PLOTFILE FOR HYDROGRAPH PLOTS OPENED
C-----------------------------------------------------------------------
C  GRAPH PREPARATION (only if IFGRF).  Builds arrays for the hydrograph,
C  the sediment discharge (SGR, kg/min) and the hyetograph (TGR, RGR) in
C  the requested output time unit.  The call that would draw them
C  (GRAFXY) is commented out in this version.
C  NOTE: the hyetograph part is only valid for a single rain gage.
C-----------------------------------------------------------------------
      IF( IFGRF ) THEN
C          THEN WRITE OUT VARIABLES TO THE HYDROGRAPH PLOT FILE
C          OBTAIN PROPER UNITS
        tlabl(2) = tuns(ntime)
        if(neros .ge. 1) ifsed = .true.
        DIV = 1.0
C assume metric graph output
        sdfac = 60.*1000.
        IF(NTIME .EQ. 2) DIV = 60.
        IF(NTIME .EQ. 3) DIV = 3600.
        ymax = 0.
        tmax = tfin/div
C          * LOOP WRITING TO PLOT FILE
        DO 905 L=1,NI
          TGQ(L) = T(L)/DIV
          IF(QD(L) .GT. YMAX) YMAX = QD(L)
C               IF(IHYDFL .EQ. 1) WRITE(40,900) TGQ(L),QD(L),Q(L)
          IF(NEROS .GE. 1) THEN
C  Sediment discharge in kg/min
            SGR(L) = Q(L)*CONC(L)*RHOS(J)*SDFAC
C  convert actual conc dat to sed discharge, kg/min:
            If(L .le.NDAT) sdat(l) = sdat(l)*area*RHOS(j)*qdat(l)/60.
          END IF
  900 FORMAT(' ',5X,3(F12.6,5X))
  905   CONTINUE
C          * LOOP SETTING UP HYETOGRAPH PLOT FILE
C            NOTE: THIS IS ONLY GOOD FOR A SINGLE RAIN GAGE -
C                  MULTIPLE GAGES MUST HAVE SOME TYPE OF AVERAGING
C                  SCHEME
C            THEN WRITE OUT VARIABLES TO THE HYDROGRAPH PLOT FILE
C            OBTAIN PROPER UNITS: NOTE TIDD ARRAY STORED IN MIN.
         DIV = 1.0
         IF(NTIME .EQ. 1) DIV = 0.016667
         IF(NTIME .EQ. 3) DIV = 60.
         IGAGE = 1
  910  FORMAT(' ',5X,2(F10.6,5X))
         TID = TIDD(IGAGE,1)/DIV
         QID = QIDD(IGAGE,1)
         L = 1
         TGR(L) = TID
         RGR(L) = QID/10.
C          IF(IHYDFL .EQ.1) WRITE(41,910) TID,QID
C  J is reused here as a loop counter over rainfall data pairs; this is
C  safe only because the element index J is not needed after this point.
         DO 500 J=2,ND
           L = L+1
           TID = TIDD(IGAGE,J)/DIV
           TGR(L) = TID
           QID = QIDD(IGAGE,J)
           QIDM1 = QIDD(IGAGE,J-1)
           RGR(L) = QIDM1/10.
           IF(RGR(L) .GT. YMAX) YMAX = RGR(L)
C   IF(IHYDFL .EQ.1) WRITE(41,910) TID,QIDM1
           L = L+1
           TGR(L) = TID
           RGR(L) = QID/10.
           IF(RGR(L) .GT. YMAX) YMAX = RGR(L)
C            IF(IHYDFL.EQ.1) WRITE(41,910) TID,QID
  500   CONTINUE
C	   TMAX = TID
        NYL = INT((YMAX+0.1)/0.1)
        YMX = 0.1*NYL
        TMIN = 0.
c        CALL GRAFXY(1,0,TMIN,TMAX,0.,YMX,TGQ,QD,SGR,TGR,RGR,1000,NI,L)
C
	ENDIF
C     OUTPUT INTERCEPTION FILES
C     OUTPUT FOR GLOBAL WATER BALANCE
C-----------------------------------------------------------------------
C  GLOBAL WATER BALANCE
C    rainfall = basin runoff + storage remaining + infiltration
C    Volumes (m3) accumulated per element type are converted to mm over
C    the whole basin area.  The percent error is written at the end.
C-----------------------------------------------------------------------
      WRITE(IWRITE,106)
  106 FORMAT(1x ,//24X,'**** EVENT SUMMARY ****',/,
     -       ' ',28X,'----- -------',//,
     -   ' ',24X,'GLOBAL VOLUME BALANCE',/,
     -   ' ',12X,'VALUES ARE IN UNITS OF LENGTH (VOL./BASIN AREA)',/)
      WRITE(IWRITE,108) AREA
  108 FORMAT(' ',17X,'BASIN AREA = ',G16.8,'  (M**2)',//)
C
  127 FORMAT(' ',17X,'TOTAL RAINFALL DEPTH = ',F7.3,' (MM)'//)
C
C  Totals for planes (P), channels (C) and ponds (PD).
      TPINF  = 0.0
      TCINF  = 0.0
      TPSTOR = 0.0
      TCSTOR = 0.0
      TPDST  = 0.0
C     * LOOP SUMMING STORAGES AND INFILTRATION VOL. (CUBIC M)
      DO 310 I = 1,NELE
        IP = NLOG(I)
C          * IF A PLANE ELEMENT
        IF(W(IP) .GT. 0.001) THEN
C             THEN ADD IN THE PLANE STORAGE
C
C    PLANE STORAGE AND TOTAL INFILTRATION HAS UNITS OF CUBIC M
          TPSTOR = TPSTOR + STORA(IP)
          TPINF  = TPINF  + EINF(IP)
        ELSE IF(W(IP) .LE. 0.001 .AND. NPN(IP) .EQ. 0) THEN
C    THEN ADD IN CHANNEL (INCLUDING CONDUIT) STOR. AND INF
          TCSTOR = TCSTOR + STORA(IP)
          TCINF  = TCINF  + EINF(IP)
        ELSE IF(NPN(IP) .EQ. 1) THEN
C    THEN ADD IN THE STORAGE FOR THE POND
          TPDST  = TPDST  + STORA(IP)
        ENDIF
  310 CONTINUE
C
C     CONVERT TO mm OVER THE BASIN
C
      TPSTOR = TPSTOR/AREA* 1000.0
      TCSTOR = TCSTOR/AREA *1000.0
      TPDST  = TPDST/AREA * 1000.0
C
      TPINF  = TPINF/AREA * 1000.0
      TCINF  = TCINF/AREA * 1000.0
C
C     TOTAL OF RUNOFF AND REMAINING STORAGE (mm)
      TRSI = SMROMM + TPSTOR + TCSTOR + TPDST + TPINF + TCINF
C     COMPUTE PERCENT ERROR
      ER = (sumrmm - TRSI)/sumrmm * 100.0
C  Write the balance to both the dynamic output file and the static
C  summary file (unit 99).
      WRITE(99,977)
      WRITE(99,978) sumrmm
      WRITE(IWRITE,115) TPSTOR,TCSTOR,TPDST,TPINF,TCINF,SMROmm,
     &                sumro,TRSI,ER
      WRITE(99,115) TPSTOR,TCSTOR,TPDST,TPINF,TCINF,SMROmm,
     &                     sumro,TRSI,ER
C
  977 FORMAT(//,' ',2X,'GLOBAL VOLUME BALANCE',/,
     -       ' ',2X,'=====================',/)
C
  978 FORMAT(
     - ' ',2X,'TOTAL RAINFALL DEPTH                    = ',F9.3,
     -                                                      ' (MM)')
C
  115 FORMAT(/,
     - ' ',2X,'STORAGE REMAINING ON ALL PLANES         = ',F9.5,
     -                                                      ' (MM)',/,
     - ' ',2X,'STORAGE REMAINING IN CHANNELS+CONDUITS  = ',F9.5,
     -                                                      ' (MM)',/,
     - ' ',2X,'STORAGE REMAINING IN PONDS              = ',F9.5,
     -                                                      ' (MM)',/,
     - ' ',2X,'TOTAL INFILTRATION FROM ALL PLANES      = ',F9.5,
     -                                                      ' (MM)',/,
     - ' ',2X,'TOTAL INFILTRATION FROM ALL CHANNELS    = ',F9.5,
     -                                                      ' (MM)',/,
     - ' ',2X,'TOTAL BASIN RUNOFF                      = ',F9.5,
     -                                 ' (MM)',F15.4,' CU.M.',/,
     - ' ',4X,'                                        ---------',/,
     - ' ',2X,'TOTAL OF STOR., INFIL. AND RUNOFF TERMS = ',F9.5,
     -                                                    ' (MM)',//,
     -/,' ',15X,'*** GLOBAL VOL. ERROR  = ',F10.4,' PERCENT ***',//)
C     
C  NOTE for a GUI or wrapper program: a normal run ends with this STOP.
      CALL ENDOK
      END
C***********************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE HYDwri(J)
C    Writes the hydrograph table and hydrology summary for element J.
C    Called from MAIN when NP(J) > 0.
C      - time series of rain, discharge, concentration and sediment
C        discharge to the dynamic output file (unit IWRITE)
C      - summary (net rain, peak rain rate, time to runoff, duration,
C        peak flow and time, peak sediment discharge) to unit 99
C      - for planes with rills, the rill dimensions down the slope
C-----------------------------------------------------------------------
      SUBROUTINE HYDwri(J)
C
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /CNTRL/ NRES, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60),ASR(60), rfr(60), DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1,CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1     T(1000), Q(1000),HUB(1000), DX,DT, INDEX, THETA, XNU, GRAV,
     2     NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
C
C      common/splice/ thinp(1000)
      DIMENSION TD(1000),Qmmh(1000),QS4(1000)
C
      CHARACTER*3 TYPE
      CHARACTER*10 IU1,IU2
      DATA IU1,IU2/' SQ. METER',' HECTARES '/
  610 FORMAT(A1)
C
C  Guard: more than 100 rainfall data pairs overflows the QI/TI arrays.
C  This STOPs the whole program (a GUI must report this separately).
      If(ND .ge. 101) Then
        Write(*,*)' dimension of output array exceeded. '
        CALL FATAL('Dimension of output array exceeded (ND>100)')
      End If
C  Contributing area of the element (m2), also in hectares.
      AREAP = SUMA(J)
      HECTA = areap/1.e4
C
C  DIV converts time in seconds to the output unit chosen by NTIME:
C  1 = seconds, 2 = minutes.  (3 = hours is written here, but INSPEC
C  rejects NTIME > 2, so only 1 and 2 are usable.)
      DIV = 1.
      IF(NTIME .EQ. 2) DIV = 60.
      IF(NTIME .EQ. 3) DIV = 3600.
C  Running values: peak discharge QDMAX at time TDMAX, peak sediment
C  discharge QSMAX at time TSMAX, start time and duration of runoff.
C    INIT FOR TIME TO PEAK AND PEAK RATE
      QDMAX = 0.0
      TDMAX = 0.0
      qsmax =0.0
      tsmax=0.0
C
C    INIT FOR TIMETO RUNOFF AND DURATION
      TDSTAR=0.0
      TDDUR=0.0
      rrmx = 0.
C-----------------------------------------------------------------------
C  Loop over time steps:
C    TD    time in output units
C    QMMH  discharge as depth rate over the contributing area (mm/h)
C    QS4   sediment discharge (kg/min) = Q*60*1000*CONC*RHOS
C  Also tracks the start of runoff (first step with QMMH > 0), its
C  duration (until QMMH returns to zero), and the two peaks.
C-----------------------------------------------------------------------
      DO 75 L=1,NI
        TD(L) = T(L)/DIV
C
C   FROM m**3/sec TO MM/HR
       Qmmh(L) = Q(L)*3.6e6/suma(j)
C   GIVES SOIL LOSS IN KG/MIN
       QS4(L) = (Q(L)*60.*1000.)*CONC(L)*RHOS(J)
C    FIND TIME TO PEAK AND PEAK RATE
       IF(QMMH(L).GT.0.0.AND.TDSTAR.LE.0.0)  TDSTAR=TD(L)
       IF (TDSTAR.GT.0.0)THEN
         IF (TDDUR.LE.0.0) THEN
           IF (QMMH(L).EQ.0.0) THEN
             IF (TD(L).GT.TDSTAR) TDDUR=TD(L)-TDSTAR      
           ENDIF
         ENDIF
       ENDIF
       IF( Qmmh(L) .GT. QDMAX ) THEN
         QDMAX = Qmmh(L)
         TDMAX = TD(L)
       ENDIF
       IF( Qs4(L) .GT. QsMAX ) THEN
         QsMAX = Qs4(L)
         TsMAX = TD(L)
       ENDIF
   75 CONTINUE
C  Net rainfall depth for this element (mm), and the peak rainfall
C  intensity (mm/h) from the rainfall data pairs (3.6E6 converts m/s to
C  mm/h).
      sumrn = rndpth(j)*1000.
      Do 80 n=1,nd
        rmmh = qi(n)*3.6e6
        If(rmmh .gt. rrmx) rrmx = rmmh
   80 Continue
C
C  If the discharge was still positive at the last step, the runoff
C  duration runs to the end of the simulation.
      IF (TDDUR.LE.0.0) THEN
        IF (QMMH(NI).GT.0.0)  TDDUR=TD(NI)-TDSTAR
      END IF
      If(NTIME .le. 1) Then
        TYPE = 'SEC'
      Else If(NTIME .eq. 2) Then
        TYPE='MIN'
      Else If(NTIME .ge. 3) Then
        TYPE = 'HRS'
      End If
C
C  Heading, then the time series table: time, rain intensity, Q (m3/min
C  and mm/h), sediment concentration and sediment discharge (kg/min).
C  The rain intensity shown is taken from rain gage 1.
	WRITE(IWRITE,81) J,AREAP,IU1,HECTA,IU2
   81 FORMAT(//,25X,'HYDROGRAPH FOR ELEMENT ',I3,/2X,'  CONTRIBUTING ARE
     -A =' ,F12.2,A10,' OR',F10.4,A10)
C
       WRITE(IWRITE,120) TYPE
      m = 1
      Do L=1,NI
        If(TIdd(1,m+1)/60. .lt. td(L)) m = m+1
        rmmh = qidd(1,m)*3.6e6
        WRITE(IWRITE,121) TD(L), rmmh, Q(L)*60., Qmmh(L),
     1       CONC(L), QS4(L)
C     1       thinp(L), QS4(L)
      End Do
C
  120 FORMAT(/3X,'Time(',A3,') ',2x'Rain(mm/h)',3X,'Q(M3/Min)',4X,'Q(MM/
     1H)',4X,'CONC.',4X,'QS(KG/MIN)'//)
  121 Format(3X,F9.2,4x,F6.2,4X,F9.5,2X,F9.3,3X,F10.8,3X,G10.4)
C  Peak flow line on the dynamic file; the full hydrology summary
C  goes to the static file (unit 99).
       WRITE(IWRITE,125) TDMAX,TYPE,QDMAX
       WRITE(99,925)J,sumrn,rrmx,TDSTAR,TYPE,TDDUR,TYPE,TDMAX,TYPE,
     -    QDMAX,tsmax,qsmax
  125 FORMAT(/,9X,'TIME TO PEAK FLOW RATE =',G12.5,' (',A3,')',
     2       /,9X,'PEAK FLOW RATE =',9X,G12.5,'(MM/H)',/)       
  925 FORMAT(/,3X,'HYDROLOGY SUMMARY, ELEMENT',I3,/,
     -         3X,'==============================',
     -  //,3x,'NET RAINFALL                            =',G12.5,
     -                                                   ' (MM)'
     -   /,3X,'PEAK RAINFALL RATE                      =',g12.5,
     -                                                  ' (MM/H)',   
     -  //,3X,'TIME TO RUNOFF                          =',G12.5,
     -                                                   ' (',A3,')',
     2   /,3X,'DURATION OF RUNOFF                      =',G12.5,
     -                                                   ' (',A3,')',
     -   /,3X,'TIME TO PEAK FLOW RATE                  =',G12.5, 
     -                                                   ' (',A3,')',
     2   /,3X,'PEAK FLOW RATE                          =',G12.5,
     -                                                    '(MM/H)',
     -  //,3x,'TIME TO PEAK SEDIMENT DISCHARGE         =',G12.5,
     -                                                    '(MIN)',
     -   /,3x,'PEAK SEDIMENT DISCHARGE                 =',G12.5,
     -                                                 '(kg/MIN)',/)
C  Rill dimension summary down the slope (planes with rills only).
C     TEST TO CHECK THAT ELEMENT IS A PLANE
      IF (W(J).GT.0.0 .and. RILLW(J) .GT. 0.001) Then
        WRITE(99,979)J         
        DO I=1,NK
          WRITE(99,980)FLOAT(I-1)*DX,1000.*RILLD1(I),
     &             1000.*RILLW1(I),1000.*(rilld1(I)-arrild(i)),
     &             1000.*(RILLW1(I)-ARRILW(I))
        END DO
      END IF
C                                                    rill output
  979 FORMAT(//,
     - ' ',4X,'RILL DIMENSION SPATIAL SUMMARY, ELEMENT',I3,/,
     - ' ',4X,'-------------------------------------------',//,
     - ' ',4X,'DISTANCE ',2X,'RILL DEPTH',2X,'RILL WIDTH',4X,
     - ' DEPTH ',4X,' WIDTH ',/,
     - ' ',4X,'DOWNSLOPE',28X,'INCREASE ',4X,'INCREASE ',/,
     - ' ',4X,'   M     ',2X,'   mm     ',2X,'   mm     ',2X,
     - '  mm   ',4X,'  mm   ',//)
  980 FORMAT(' ',2X,F9.2,3X,F9.2,3X,F9.2,3X,F9.2,3X,F9.2)
C
      WRITE(IWRITE,99)
  99  FORMAT(' ')
      RETURN
      END
C**********************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE INTSEL(IFGRF)
C    Interactive file selection.  Asks for (or confirms) the file names,
C    then opens all files.  The names are remembered between runs in
C    EUROM.FIL, in this order:
C        1 dynamic output   2 auxiliary (diagnostic) output
C        3 pond input
C        4 parameter input  5 rainfall input   6 static output
C        7 observed (recorded) data for graph comparison
C    Pressing RETURN at a prompt keeps the remembered name.
C    NOTE: all file names are CHARACTER*12, so a name, including any
C    directory path, is limited to 12 characters.
C-----------------------------------------------------------------------
      SUBROUTINE INTSEL(IFGRF)
C
C  INTERACTIVE FILE SELECTION ROUTINE
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
C
      Common /pldat/ ifad,ndat,tdat(500),qdat(500),sdat(500),smromm,
     1                sumskg
      CHARACTER*80 inpt,FIL1,FIL2,FIL3,FIL4,FIL5,fil6,FIL7,dtitl,BLANK
      LOGICAL LHEAD
      CHARACTER*1 HANS
      CHARACTER*80 HTITLE,HFIL(5)
      COMMON /HLESSL/ LHEAD
      COMMON /HLESSC/ HANS,HTITLE,HFIL
      CHARACTER*1 ANS
      logical new, ifgrf,ifad
      BLANK = '            '
      fil7 = blank
C
      FIL3 = '            '
C  Headless run: names come from eurosem.ctl (see HEADLS); no prompts
C  and no EUROM.FIL.
      IF(LHEAD) THEN
        FIL4 = HFIL(1)
        FIL5 = HFIL(2)
        FIL1 = HFIL(3)
        FIL2 = HFIL(4)
        FIL6 = HFIL(5)
        GO TO 3
      END IF
C  Read the remembered names from the last run, if EUROM.FIL has any.
      OPEN(9,FILE='EUROM.FIL',STATUS='UNKNOWN')
      new = .true.
      READ(9,1000,END=2)FIL1
      READ(9,1000,END=2)FIL2
      READ(9,1000,END=2)FIL3
      READ(9,1000,END=2)FIL4
      READ(9,1000,END=2)FIL5
      read(9,1000,END=2)fil6
      read(9,1000,END=2)FIL7
       new = .false.
      Write(*,1010)
      write(*,1006) FIL4
      WRITE(*,1007) FIL5
      If(fil3 .ne. blank) WRITE(*,1005) FIL3
      WRITE(*,1013) FIL6
      WRITE(*,1011) FIL1
      WRITE(*,1012) FIL2
      if(FIL7 .NE. BLANK .and. IFGRF) WRITE(*,1017) FIL7
C
C  Offer to reuse all remembered names; otherwise ask for each in turn.
      WRITE(*,1015)
1015  FORMAT(/,' ','  Use These I/O FILES? (Y or N) ')
      READ(*,1000) ANS
      IF(ANS.EQ.'y' .OR. ANS.EQ.'Y') go to 3
C                                              PARAMETER INPUT
   2   continue
       WRITE(*,1001)
       If(.not.new) write(*,1020)
       READ(*,1000) inpt
       If(inpt .ne. blank) FIL4 = inpt
C                                              RAINFALL INPUT
       WRITE(*,1002)
       If(.not.new) write(*,1020)
       READ(*,1000) inpt
       If(inpt .ne. blank) FIL5 = inpt
C                                              STAtic OUTPUT
       WRITE(*,1014)
       If(.not.new) write(*,1020)
       READ(*,1000) inpt
       If(inpt .ne. blank) FIL6 = inpt
C                                              dynamic OUTPUT
       WRITE(*,1003)
       If(.not. new) write(*,1020)
       READ(*,1000) inpt
       If(inpt .ne. blank) FIL1 = inpt
C                                              DIAGNOSTIC OUTPUT
       WRITE(*,1004)
       If(.not. new) write(*,1020)
       READ(*,1000) inpt
       If(inpt .ne. blank) FIL2 = inpt
C                                             POND INPUT
       WRITE(*,1008)
       READ(*,1000) ANS
       IF(ANS.EQ.'Y' .OR. ANS.EQ.'y') THEN
         WRITE(*,1009)
         READ(*,1000) FIL3
       Else
          fil3 = blank
       END IF
C                                             Recorded Data
      If(IFGRF) Then
        FIL7 = blank
        Write(*,1016)
        Read(*,1000) Fil7
      End If
C                                         
C  Open the files.  Output files are created if absent; the parameter
C  and rainfall files must already exist (STATUS='OLD').
C    16 dynamic output   8 diagnostic output   99 static output
C    15 parameter input   7 rainfall input
    3 Continue
      OPEN(16,FILE=FIL1,STATUS='UNKNOWN')
      OPEN(8,FILE=FIL2,STATUS='UNKNOWN')
      OPEN(99,FILE=FIL6,STATUS='UNKNOWN')
      OPEN(15,FILE=FIL4,STATUS='OLD')
      OPEN(7,FILE=FIL5,STATUS='OLD')
C  Optional observed-data file (only with the graph option): two header
C  lines, then up to 500 rows of  time(min)  discharge(mm/h)  sediment
C  concentration(by volume).
      If(fil7 .ne. blank .and. IFGRF ) Then
        OPEN(12,FILE=FIL7,STATUS='OLD')
C                                  Read Actual Data for Plot Comparison:
C        time in minutes, qdat in mm/h, and sdat as concentration by vol
        Read(12,902) Dtitl
        Read(12,902) dum
        IFAD = .true.
        Do 40 i= 1,500
          Read(12,*,end=43) tdat(i),qdat(i),sdat(i)
 40     Continue
 43     Continue
        NDAT = i-1
      Else
        IFAD = .false.
      End If
  902 Format(a)
C  
C  Write the header and key of the diagnostic file (see AUXTXT).
      Call AUXTXT(IDIAGN)
C
      If(fil3 .ne. blank) Then
        WRITE(IWRITE,1005)FIL3
        WRITE(IDIAGN,1005)FIL3
      End If
      WRITE(IWRITE,1006)FIL4
      WRITE(IDIAGN,1006)FIL4
      WRITE(IWRITE,1007)FIL5
      WRITE(IDIAGN,1007)FIL5
C  Save the names in EUROM.FIL for the next run.
c                               FOR NEXT TIME USE:
      IF(LHEAD) RETURN
      REWIND 9
      WRITE(9,1000)FIL1
      WRITE(9,1000)FIL2
      WRITE(9,1000)FIL3
      WRITE(9,1000)FIL4
      WRITE(9,1000)FIL5
      write(9,1000)fil6
      write(9,1000)fil7
      close(9)
      RETURN
1000  FORMAT(A)
1001  FORMAT(/,' ','NAME OF INPUT',
     +' PARAMETER FILE (UP TO 12 CHAR.):')
1002  FORMAT(/,' ','NAME OF INPUT',
     +' RAINFALL FILE (UP TO 12 CHAR.) :')
1003  FORMAT(/,' ','NAME OF DYNAMIC',
     +' OUTPUT FILE (UP TO 12 CHAR.) :')
 1004 FORMAT(/,' ','NAME OF AUXILIARY',
     +' OUTPUT FILE (UP TO 12 CHAR):')
 1005 FORMAT('   INPUT POND FILE:       ',3X,A)
 1006 FORMAT('   INPUT PARAMETER FILE:  ',3X,A)
 1007 FORMAT('   INPUT RAINFALL FILE:   ',3X,A/)
 1008 FORMAT(/'   ARE ANY ELEMENTS PONDS? (Y OR N) :')
 1009 FORMAT(/,' ','NAME OF INPUT',
     +' POND FILE (UP TO 12 CHAR.)     :')
 1010 Format(' Filename Assignments in Memory: '/)
 1011 Format('   Dynamic Output File:   ',3x,A)
 1012 Format('   Output Auxiliary File: ',3x,A)
 1013 Format('   Static Output File:    ',3x,a)
 1014 FORMAT(/,' ','NAME OF STATIC',
     +' OUTPUT FILE (UP TO 12 CHAR):')
 1016 Format(/1x,' Enter Name of Recorded Data File: (If any): ')
 1017 Format(/'   Recorded Data File:    ',3x,a)
 1020 Format(' (RETURN to accept existing file)  ')
      END
C-----------------------------------------------------------------------     
C  Writes the titles and header of the EUROSEM auxiliary (diagnostic)
C  output file; error messages are appended here during the run.
	subroutine auxtxt(idiagn)
C     provides titles and additional header information for the      
C     EUROSEM auxilary output
C    
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), rfr(60), DEPNO(60)
     4       ,rs(60), dero(60), sir(60)
C
      write(idiagn,100)
      write(idiagn,101)
 100  format(15x,'EUROSEM AUXILARY OUTPUT',/,15X,
     -'=======================',/)
 101  FORMAT (2X,'The following is a key to descriptors and error',/,
     - 'messages which may appear in this file.',/,
     - ' Errors',/)
C     - ' ======',/,
C     - 't,cx =  : ') 
      return
      end
C*********************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE READER
C    Reads the run set-up from the parameter file (unit IREAD) and the
C    rainfall file (unit JREAD), sets defaults, checks the input and
C    converts to internal units (CONVER(1)).
C
C    Both files are free-format records interleaved with template lines
C    that are skipped by line COUNT (READTM), so the number of header
C    lines in each section must be exactly as below.
C
C    PARAMETER FILE
C       skip 6   NELE NPART CLEN TFIN DELT THETA TEMPW      (1 record)
C       skip 5   NTIME NEROS                                (1 record)
C       skip 10  NELE records:  computation-order no., element no.
C       skip 7   then PAREAD reads, for each element: 6 records, each
C                preceded by 2 skipped lines (see PAREAD), and 1 skipped
C                terminator line after the last record.
C
C    RAINFALL FILE
C       skip 9   NGAGES MAXND                               (1 record)
C       skip 5   NELE records:  element no., gage no., weight
C       skip 8   then for each gage:
C           skip 4   gage no., ND                           (1 record)
C           skip 6   ND records:  time (min), cumulative depth (mm)
C           skip 1
C-----------------------------------------------------------------------
      SUBROUTINE READER
C
      COMMON /CNTRL/ nres, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60),ASR(60), rfr(60), DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000),Q(1000),HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /ESEM/ COVER(60), CUINT(60), ISHAPE(60), PLHGT(60),
     1        VPANG(60),PBASE(60)
      COMMON /PONDS/ AP(10,10,3), TWP(10,10,3), EL(10,3), QSO(20,3),
     1   ELQ(20,3),ELZQ(3),ELVZ(3),NOPQ(3),MEL(3),NSEC(3),XSEC(10,3),
     2   ELZX(10), JOP(3), ELST(3), ALV, DIFUS(3), SIGMAS(60), NPART
C
      common /labels/ title, tlabl, qlabl, slabl, ifsed
      Character*12 qlabl, tlabl, slabl
      CHARACTER TITLE*80
C
      DATA NTELE /60/
C
C************************************************************************
C*     Code inserted for run-time overrides of SI and NP parameters,
C*     and run-time input of recession factor RECS
C************************************************************************
C  Echo the run title to the dynamic and diagnostic output files.
      WRITE(IWRITE,278) TITLE
      WRITE(IDIAGN,278) TITLE
278     FORMAT(' ',20X,' === DESCRIPTIVE RUN TITLE === ',/,' ',A80,/)
C
************************************************************************
C  Array limits used by the solvers: LENQS (size of the outflow storage
C  QS/SCON), LENQ (time steps per element), LENH1 (nodes per element).
      LENQS = 4000
      LENQ  = 1000
      LENH1 = 20
      GRAV  = 9.806
C
C  RUN CONTROL RECORD:
C    NELE   number of elements    NPART  (set to 1 by DEFAUL if <= 0)
C    CLEN   characteristic length (m), used to set the node spacing
C    TFIN   run duration (min)         DELT   time step (min)
C    THETA  time-weighting factor of the implicit scheme
C    TEMPW  water temperature (deg C)
      CALL READTM(IREAD,6)                                                READTM
C            modified for EUROSEM:
      READ(IREAD,*,END=9001)NELE,NPART,CLEN,TFIN,DELT,THETA,TEMPW
C
C  Fix manning resistance law:
C  NRES = 1 selects the Manning resistance law.
      nres = 1
C  Check storage limit:
C  Storage limit: the run may have at most 999 time steps.
      If(TFIN/DELT .gt. 999.) Then
        CALL FATAL('TFIN/DELT must be less than 1000')
      End If
C  INSPEC(k,..) checks what was just read; the END= branches above land
C  here when the file ends early.  Code 1 = run control record.
 9001 CALL INSPEC(1,0,0)                                                 INSPEC
      CALL READTM(IREAD,5)                                               READTM
C  OPTIONS RECORD:  NTIME output time unit (1 = sec, 2 = min; INSPEC
C                   rejects values above 2)
C                   NEROS erosion option (0 = no sediment calculation)
      READ(IREAD,*,END=9002)NTIME,NEROS
 9002 CALL INSPEC(2,0,0)                                                 INSPEC
      CALL READTM(IREAD,10)                                              READTM
C  Defaults: THETA = 0.7, TEMPW = 10 deg C, NTIME = 2 (minutes),
C  NELE = 1.
C                                SET DEFAULTS
      IF(THETA.LE.0.0) THETA = 0.7
C             Allow for frozen soil, but limit viscosity calc.:
      IF(TEMPW.LE.0.0)   TEMPW  = 10.
      IF(NTIME.LE.0)   NTIME = 2
      IF(NELE.LE.0)    NELE  = 1
C  COMPUTATION ORDER: each record gives a position IORD and the element
C  number J to compute at that position, stored as NLOG(IORD) = J.
C  Upstream elements must come before the elements they feed.
C           SEQUENTIAL ORDER OF NLOG
      DO 108 L=1,NELE
        READ(IREAD,*,ERR=9003) IORD,NLOG(IORD)
C        IF(L .NE. IORDER) THEN
C          WRITE(*,900)
C          WRITE(IDIAGN,900)
C 900  FORMAT(//,' ',10X,' *** FATAL ERROR *** ',/,
C     -      ' ',' THE COMP. ORDER (NLOG) FROM THE INPUT PARAMETER ',/,
C     -      ' ',' MUST BE IN SEQUENTIAL ASCENDING ORDER ',//)
C          STOP
C        ENDIF
 108  CONTINUE
 9003 CALL INSPEC(3,0,0)                                                 INSPEC
C
C  Read the element-wise data (geometry, soil, cover, ...).
      CALL PAREAD                                                        PAREAD
C
C  RAINFALL FILE.  NGAGES = number of rain gages, MAXND = largest number
C  of time-depth pairs for any gage.
      CALL READTM(JREAD,9)                                               READTM
      READ(JREAD,*,END=10008) NGAGES,MAXND
      CALL READTM(JREAD,5)                                               READTM
C  Link each element to a rain gage, with a weight (<= 0 becomes 1.0).
      DO 105 I=1,NELE
        READ(JREAD,*) IEL,IGAGES(IEL),WTRAIN(IEL)
        IF(WTRAIN(IEL) .LE. 0.0) WTRAIN(IEL) = 1.0
 105  CONTINUE
C +++++++++++++++++++++++++++++++++++++++++++++++  START LOOP OVER GAGES      ++
C                                                                             ++
C                                                                             ++
      CALL READTM(JREAD,8)                                               READTM
C  LOOP OVER RAIN GAGES: read ND time-depth pairs (time in min,
C  cumulative depth in mm) for each gage.
10008 DO 230 IGAGE=1,NGAGES
       CALL READTM(JREAD,4)                                               READTM
       READ(JREAD,*) IGNUM,ND
       CALL READTM(JREAD,6)                                               READTM
C  Only the first ND pairs are read, ND being the count given in the
C  file; any further pairs listed are ignored (skipped by the next
C  READTM).
       DO 104 I=1,ND
        READ(JREAD,*) TI(I),QI(I)
104    CONTINUE
       CALL READTM(JREAD,1)                                               READTM
C-----------------------------------------------------------------------
C  If the last rainfall time is before TFIN, a closing point is
C  appended: time TFIN + 0.05 min, with the same cumulative depth (no
C  more rain).  The small offset avoids a divide-by-zero when intensity
C  is computed.  MAXND is raised if the extra point exceeds it.
C-----------------------------------------------------------------------
10009  IF(TI(ND) .GE. TFIN) GO TO 10010
       ND = ND + 1
C       ADD SMALL TIME TO TI(ND) TO AVOID DIVIDE BY ZERO IN
C                 INTENSITY CALCUALTION
       TI(ND) = TFIN + 0.05
C    SET FLAG FOR QRAT INTERPOLATION
       QI(ND) = QI(ND-1)
       IF(ND .gt. MAXND) Then
         MXAND=MAXND
         MAXND=MAXND+1
         IBACK=NGAGES-1
         DO 205 IGAGEX=1,IBACK
C     ADD SOME SMALL AMOUNT OF TIME TO TIDD(IGAGEX,MXAND)
C               SO DON;T GET DIVIDE BY ZERO WHEN COMP. INTENSITY
          TIDD(IGAGEX,MAXND)=TIDD(IGAGEX,MXAND) + 0.05
          QIDD(IGAGEX,MAXND)=QIDD(IGAGEX,MXAND)
          WRITE(IDIAGN,207)
  205    Continue
207         FORMAT(1X,'MAXND HAS BEEN INCREASED BY 1')
       End If
       WRITE(IDIAGN,203)
203   FORMAT( ' A POINT QI(ND+1),TI(ND+1) HAS BEEN ADDED TO RAIN DATA SO
     1 THAT TI(ND+1) = TFIN')
10010 CALL INSPEC(5,0,0)                                                 INSPEC
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                        NESTED LOOP         ++
C  Copy each gage's data into TIDD/QIDD (time, depth), padding short
C  gages up to MAXND points with the last depth at times 0.05 min
C  apart.
       DO 230 LND=1,MAXND
        IF(LND .le. ND .or. LND .gt. MAXND) Then
          TIDD(IGNUM,LND)=TI(LND)
C  cum depth in mm
          QIDD(IGNUM,LND)=QI(LND)
          cumd(ignum,LND) = QI(LND)
        Else 
C   FORM INCREASING TIDD ARRAY AND  MAINTAIN LAST DEPTH IN QIDD ARRAY
          TIDD(IGNUM,LND)=TI(ND) + 0.05*(LND-ND)
          QIDD(IGNUM,LND)=QI(ND)
          cumd(ignum,LND) = QI(ND)
        End If
230   CONTINUE
C                                                         END LOOPS           ++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  From here on ND = MAXND, the length of all gage records.
      ND=MAXND
C
C  Kinematic viscosity of water (m2/s) from the temperature TEMPW.
C     COMPUTE KINEMATIC VISCOSITY IN SQ. M/SEC.
C
      XNU   =(.0000017756/(1.+0.03368*TEMPW + 0.000221*TEMPW*TEMPW))
C
C  Clear the offsets that locate each element's outflow in QS/SCON.
      DO 40 I=1,NTELE
        NB(I) = 0
   40 CONTINUE

C  Convert times to seconds, rainfall to m/s, discharges to m3/s, etc.
CC     ---  CONVERT TIME INTO SECONDS, RAIN to M/SEC 
C             AND DISCHARGE INTO CU.M PER SEC.
      CALL CONVER(1)                                                    CONVERT
      RETURN
      END
C***************************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE PAREAD
C    Reads the element-wise records of the parameter file, for each of
C    the NELE elements.  Each element has 6 records, each preceded by 2
C    skipped template lines (a units/heading line and a rule), and the
C    whole set ends with 1 skipped terminator line ("*"):
C       1  J  NU  NR  NL  NC1  NC2  NPRINT                topology
C       2  XL  W  S  ZR  ZL  BW  R1  R2                   geometry
C       3  FMIN  G  POR  THI  THMX  ROC  RECS  DINT       infiltration
C       4  DEPNO RILLW RILLD ZRL RS RFR SIR               rills
C                  (for channels the record is skipped, not read)
C       5  COVER SHAPE PLANGLE PLANTBASE PLANTH DERO STONE   cover
C       6  D50 EROD SPLTEX COH RHOS PAVE SIGMAS MCODE     soil/erosion
C    Elements may appear in any order (J identifies the element).
C    Several checks here STOP the program or ask a question on the
C    keyboard; these would need handling in a GUI.
C-----------------------------------------------------------------------
      SUBROUTINE PAREAD
C
C              Reads data from file IREAD (formerly NAMELIST variables)
C
      LOGICAL LHEAD
      CHARACTER*1 HANS
      CHARACTER*80 HTITLE,HFIL(5)
      COMMON /HLESSL/ LHEAD
      COMMON /HLESSC/ HANS,HTITLE,HFIL
      COMMON /INFIL/ AL(60), THI(60), THMX(60), ROC(60), RECS(60),
     1        THR, THFC, TEMPW, stone(60), RVAL
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      COMMON /CNTRL/ nres, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /ESEM/ COVER(60), CUINT(60), ISHAPE(60), PLHGT(60),
     1        VPANG(60),PBASE(60)
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      COMMON /PONDS/ AP(10,10,3), TWP(10,10,3), EL(10,3), QSO(20,3),
     1   ELQ(20,3),ELZQ(3),ELVZ(3),NOPQ(3),MEL(3),NSEC(3),XSEC(10,3),
     2   ELZX(10), JOP(3), ELST(3), ALV, DIFUS(3), SIGMAS(60), NPART
      character*1 answ
      logical chan
C
C  Skip the 7 heading lines of the element-wise section.
      CALL READTM(IREAD,7)                                               READTM
C
      DO 10 I=1,NELE
      CALL READTM(IREAD,2)                                               READTM
C                                    NPRINT IS READ AND USED AS NP(J).
C                                    THE NAME NPRINT IS MAINTAINED IN
C                                    USER DIALOGUE ONLY, AN ARTIFACT FROM
C                                    NAMELIST DAYS.
C-----------------------------------------------------------------------
C  RECORD 1, topology (how flow is passed between elements):
C    NU        element feeding the upper end of this element
C    NR, NL    elements feeding laterally from the right / left
C    NC1, NC2  channels joining at the upper end of a channel
C    NPRINT    stored as NP(J):  > 0 writes this element's hydrograph;
C              7 also turns on detailed per-step diagnostics.
C    (the name NPRINT only survives in the user dialogue)
C-----------------------------------------------------------------------
       READ(IREAD,*) J,NU(J),NR(J),NL(J),NC1(J),NC2(J),np(j)
C             modified  Silsoe version (np(j) is needed for diagn.)
C     
C  NPNT (pond switch) is set to 0 and never changed in this routine, so
C  the pond reader below (PONDRD) is not reached and NPN(J) stays 0:
C  pond elements cannot be defined through this input file.
C  NRP(J) = 3 is the default output-printing level.
     	 npnt = 0
       nrp(j) = 3
C          
C                       Reading all Cards
        NPN(J) = 0
        CALL READTM(IREAD,2)                                         !READTM
C    ------------------------------------------------------------------
C  RECORD 2, geometry:  XL length (m), W width (m; W <= 0 marks a
C  channel), S slope, ZR/ZL right/left bank side slopes (channels),
C  BW bottom width (m, channels), R1 rill Manning n, R2 interrill n.
      READ(IREAD,*) XL(J),W(J),S(J),ZR(J),ZL(J),BW,R1(J),R2(J)         !4/93
C  R1 is rill Manning n, R2 (below) is interrill n
C  CHAN is true for a channel.  A(J) is derived from the bottom width BW
C  and the bank side slopes ZL, ZR.
        IF(W(J).LE.0.0 ) THEN
          A(J) = BW/( 1.0/ZL(J) + 1.0/ZR(J) )
          chan = .true.
        Else
          chan = .false.
        END IF
C  Checks: element length must be positive; a plane longer than the
C  characteristic length CLEN is reported.
      If(XL(J) .le. 0.) Then
        CALL FATALI('XL must be greater than 0 on element',J)
      End If
      If(CLEN .lt. XL(J)) Then
        Write(99,801) CLEN, XL(J)
  801   Format(' Parameter Error:'/
     &' Characteristic Length ',f6.1,' < Surface Length of'
     &,F6.1/)
       End If
       CALL READTM(IREAD,2)                                         !READTM
C  --------------------------------------------------------------------
C  RECORD 3, infiltration: FMIN final infiltration rate (mm/h), AL(J)
C  capillary drive G (mm), POR porosity, THI initial and THMX maximum
C  soil moisture, ROC rock-content reduction factor, RECS surface
C  microtopography storage (mm), DINTR interception depth (mm).
C                  READ INFILTRATION RELATED PARAMETERS
       READ(IREAD,*) FMIN(J),AL(J),POR,THI(J),THMX(J)
     +                         ,ROC(J),RECS(J),DINTR(J)                 !4/93
C             AL(J) = G 
C convert surface microtopography from mm to m
C  Microtopography: mm -> m.
       recs(j) = recs(j)/1000.
C
C  Consistency checks (each STOPs the program):
C    porosity >= THMX;  THI <= THMX;  0 <= ROC <= 1.
      If(POR .lt. THMX(J)) Then
        Write(*,803) J, POR, THMX(J)
 803  Format(' Element ',i2,' Parameter Error Stop:'/' Porosity ',F5.3,
     & ' must be at least THMX [ ',F5.3,']')
        CALL FATALI('Porosity less than THMX on element',J)
      Else If(THI(J) .gt. THMX(J)) Then
        Write(*,804) J, THI(J), THMX(J)
 804  Format(' Element ',I2,' Parameter Error Stop:'/
     & ' THR [',F5.3,'] must be at least THMX [ ',F5.3,']')
        CALL FATALI('THI greater than THMX on element',J)
      End If
      If(ROC(J) .gt. 1.0 .or. ROC(J) .lt. 0.0) Then
        Write(*,806) J, ROC(J)
  806 Format(' Element ',I2,' Parameter Error Stop:'/'  Value of ROC [ r
     &ead as',f8.4,'] must be between 0. and .999 ')
        CALL FATALI('ROC outside 0 to 0.999 on element',J)
      End If
C-----------------------------------------------------------------------
C  ROC > 0: ask the user whether to reduce soil infiltration storage by
C  that fraction.  RVAL = 1 - ROC if accepted, else 1.
C  NOTE: RVAL is one scalar (COMMON /INFIL/), not one value per element,
C  and is reset to 1.0 for every element read.  Only the value left by
C  the LAST element therefore survives to the run (used in eurozl when
C  the infiltration parameters are formed).  Check this is intended.
C-----------------------------------------------------------------------
C       Check with user to see if a reduction in soil storage due to
C         ROC is desired:
      RVAL = 1.0
      If(ROC(j) .gt. 0.0 ) Then
        change = roc(j)*100.
        IF(LHEAD) THEN
          answ = HANS
        ELSE
          Write(*,701) change,J
          Read(*,702) answ
        END IF
  701 Format(/'    Soil Infiltration Storage/Capillary Coefficient'/
     &'  will be decreased by ',f5.2,' percent, due to value of ROC ente
     &red'/'  for element ',I2,'.  Type Y to accept, or N to cancel: ')
  702 Format(a1)
        If(answ .eq. 'y' .or. answ .eq. 'Y') RVAL = (1.- ROC(J))
      End If
Continue with next card:
C  RECORD 4, rills (planes only):  DEPNO number of rills across the
C  width, RILLW rill width (m), RILLD rill depth (m), ZRL rill side
C  slope, RS, RFR surface roughness index, SIR interrill slope.
        CALL READTM(IREAD,2)                                            ! READTM
C   -------------------------------------------------------------------
      If(CHAN) Then
        Call READTM(IREAD,1)
      Else
        READ(IREAD,*) depno(j),rillw(J),rilld(j),ZRL(J),
C  ASR read removed to match Quinton;s version:
     1                    rs(j),RFR(J),SIR(J)                    !1/98
C  Roughness height RFH (mm) is estimated from RFR.  If positive, the
C  interrill Manning n (R2) is increased to include it, and the change
C  is reported on the screen.
C  Estimate roughness height (mm) from RFR:
            rfh = 60.*(1.-(29.-RFR(J))/29.)-6.                   !7/93
         if(rfh .gt. 0.) Then
C  Adjust R2 for RFH - (On Loan from Opus):                      !
           RRAN = 0.00087*RFh**(1.25)                                   !
C  write(*,'(6h Rran ,f7.4)') rran
           R2R = R2(j)
           R2(J) = (R2R**1.5 + RRAN**1.5)**(0.67)                       !
           rdif = r2(j) - R2R 
           write(*,*)
       write(*,"('  RFR adds ',f8.6,' to basic Manning n specified as '
     & ,F8.6/'  Total Interrill Manning n is ',F8.6)") rdif, r2r, r2(j)
         End If
         ASR(j) = ASR(j)/100.
C  Check that the rills fit within the plane width, else STOP.
C  Wid check:
         widr = DEPNO(J)*(zrl(J)*rilld(J)*2.+rillw(J))
         If(widr .ge. W(J) )Then
           Write(*,805) J, W(J), widr
 805  Format('  Element ',I2,' Parameter Error Stop:'/
     & ' Surface width ',F5.1,' is less than total rills: [',F5.1,']')
           CALL FATALI('Width less than total rills, element',J)
         End If
       End If
C  Check width, length, slope etc. (INSPEC code 4).
       Call INSPEC(4,J,0)         ! check width, x, s, etc
Continue with next Card:
         CALL READTM(IREAD,2)                                           ! READTM
C  RECORD 5, cover: COVER canopy cover fraction, ISHAPE plant shape,
C  VPANG plant angle, PBASE plant base, PLHGT plant height (cm),
C  DERO depth of the erodible soil layer (m; splash and flow
C  detachment stop once this much has been eroded), ISTONE (stone flag).
C-------------------------------------------------------------------
       READ(IREAD,*) cover(j),ishape(j),vpang(j),pbase(j), plhgt(j),
     &  DERO(J), istone
C  Plant height: cm -> m.
C  convert plant height in cm to internal m
       if (plhgt(j).gt.0.0) then
         plhgt(J)=plhgt(J)/100.
       end if
       CALL READTM(IREAD,2)                                             ! READTM
C  RECORD 6, soil and erosion: D50 median particle size (microns),
C  EROD (ERODGJ) soil erodibility, SPLTEX splash reduction parameter,
C  COH cohesion, RHOS particle density, PAVE pavement/stone fraction,
C  SIGMAS, MCODE (see the EUROSEM user manual for definitions).
C----------------------------------------------------------------------
       READ(IREAD,*) d50(j),erodgj(j),SPLTEX(J),COHE(J),
     1                       RHOS(J),PAVE(J),SIGMAS(J),mcode
C  Store the porosity in the channel or plane array, by element type.
       IF(chan) THEN
         PORCH(J) = POR
       ELSE
         PORPL(J) = POR
       END IF
C  ISTONE >= 1: stones are assumed to INCREASE the infiltration rate
C  (STONE = +1); otherwise to decrease it (STONE = -1).
       If(istone .ge. 1) Then
C  stones assumed to increase fmin    
         stone(j) = +1.
       Else
C  stones assumed to decrease fmin
         stone(j) = -1.
       End If
C  Warn (screen only) if D50 is outside 20-400 microns, the range for
C  which the transport-capacity formula is valid.
C Valid range of particle sizes warning
       If(d50(j) .lt. 20. .or. d50(j) .gt. 400.) Then
         Write(*,704) d50(j)
 704  Format(/'   Particle size of',F7.2,' microns is outside the range
     &(20 to 400) for which the transport capacity formula is valid.'/
     &'   Proceed with forewarning.')
       End If
C  PAVE >= 0.001: ask the user to confirm adjusting FMIN by the stone
C  fraction:  FMIN = FMIN * (1 + STONE * PAVE).
       If(PAVE(J) .ge. .001) Then
         change = PAVE(J)*100.*stone(j)
         FMN=FMIN(J)
         IF(LHEAD) THEN
           answ = HANS
         ELSE
           Write(*,703) change,J
           Read(*,702) answ
         END IF
  703 Format(/'    Soil final infiltration rate (K sat) '/
     &'  will be changed by ',f6.2,' percent, due to value of PAVE enter
     &ed'/'  for element ',I2,'.  Type Y to accept, or N to cancel: ')
         If(answ .eq. 'y' .or. answ .eq. 'Y')   FMIN(j) = FMN*(1. + 
     & STONE(J)*PAVE(J))
      End If
C  Plausibility test on SPLTEX against 0.3 / D50(mm): if it is outside
C  0.1 to 10 times that value, ask the user; any answer except Y STOPs.
       d5mm = d50(j)/1000.
       tests = 0.3/d5mm
       If(spltex(j) .gt. 10.*tests .or. spltex(j) .lt. 0.1*tests)
     & then
         IF(LHEAD) THEN
           answ = HANS
         ELSE
         Write(*,807) j,spltex(j)
  807 Format(/'  WARNING:  On Element ',I2,', value of splash reduction 
     &parameter'/
     &'  [',f9.4,'], is far out of range based on experimental data.'/
     &'  Consult Documentation, or  Enter Y to proceed regardless:')
         Read(*,'(a1)') answ
         END IF
         If(answ .ne. 'Y' .and. answ .ne. 'y')
     &     CALL FATALI('SPLTEX out of range, not accepted. Element',J)
       End If
C  Check the soil data (INSPEC code 6) and fill in defaults (DEFAUL).
       CALL INSPEC(6,J,0)                                               ! INSPEC
C
       CALL DEFAUL(J)
C
C  Skip the "*" line that ends this element's records.
       CALL READTM(IREAD,1)                                             ! READTM
       IF(NPNT.GE.1) THEN
         CALL PONDRD(J,NPNT,NP(J))
       ENDIF
   10 CONTINUE
C
C  Final input checks (INSPEC code 7).
      call inspec(7,0,0)
      RETURN
      END
C**************************************************************************
C
C  NOTE: READTM reads whole lines with format A1, so the content of
C  skipped lines is irrelevant, but their COUNT must be exact.
      SUBROUTINE READTM(INUN,N)
C
C       READS N RECORDS FROM INPUT FILE INUN
C       TO SKIP TEMPLATE PRETTY-PRINTING
C
      CHARACTER*1 DUM
      DO 10 I=1,N
        READ(INUN,100) DUM
100     FORMAT(A1)
 10   CONTINUE
      RETURN
      END
C**************************************************************************
C
C-----------------------------------------------------------------------
C  BLOCK DATA: initial values of the element arrays.
C    XL, W and A start at -1 (meaning "not yet read", tested by INSPEC);
C    the other arrays start at 0.  Also initial labels for the graph
C    (time, flow in mm/h, sediment in kg/min).
C-----------------------------------------------------------------------
      BLOCK DATA
C
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), rfr(60), DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      common /labels/ title, tlabl, qlabl, slabl, ifsed
C
      character*6 tlabl(2)
      Character*12 qlabl, slabl
      character*80 title
      logical ifsed
C
      DATA XL,W,S,R1,R2,FMIN,NL,NR,NU,NC1,NC2,ZL,ZR,A,
     +     NP,DINTR,NRP,SUMA
     +/60*-1.,60*-1.,60*0.,60*0.,60*0.,60*0.,60*0,
     +60*0,60*0,60*0,60*0,60*0.,60*0.,60*-1.,60*0.,
     +60*0.,60*0,60*0./
      DATA ifsed/.false./
      Data tlabl(1)/'Time, '/,qlabl/' Flow, mm/HR'/,
     + slabl/' KG/MIN     '/
      END
C**************************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE DEFAUL(I)
C    Fills in defaults for element I when values were left at zero or
C    negative:  links NU, NR, NL, NC1, NC2 = 0 (no link);  THI = 0.05;
C    NPART = 1;  AL = 1;  SIGMAS = 0;  RHOS = 2.65;  PORPL and PORCH =
C    0.4;  PAVE = 0.  A non-positive XL becomes -1 (flags "missing").
C-----------------------------------------------------------------------
      SUBROUTINE DEFAUL(I)
C
      COMMON /INFIL/ AL(60), THI(60), THMX(60), ROC(60), RECS(60),
     1        THR, THFC, tempw, stone(60), RVAL
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), rfr(60), DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /PONDS/ AP(10,10,3), TWP(10,10,3), EL(10,3), QSO(20,3),
     1   ELQ(20,3),ELZQ(3),ELVZ(3),NOPQ(3),MEL(3),NSEC(3),XSEC(10,3),
     2   ELZX(10), JOP(3), ELST(3), ALV, DIFUS(3), SIGMAS(60), NPART
C
      IF(NU(I).LE.0) NU(I)  = 0
      IF(NR(I).LE.0) NR(I)  = 0
      IF(NL(I).LE.0) NL(I)  = 0
      IF(NC1(I).LE.0)NC1(I) = 0
      IF(NC2(I).LE.0)NC2(I) = 0
      IF(XL(I).LE.0.)  XL(I)   = -1.
      IF(thi(I).LT.0.)  thi(I)   = 0.05
      IF(NPART.LE.0)    NPART = 1
      IF(AL(I).LE.0.)   AL(I) = 1.
C      IF(NP(I).LE.0)    NP(I)   = 1
      IF(SIGMAS(I).LE.0.) SIGMAS(I) = 0.
      IF(RHOS(I).LE.0.) RHOS(I) = 2.65
      IF(PORPL(I).LE.0.)PORPL(I) = 0.4
      IF(PORCH(I).LE.0.)PORCH(I) = 0.4
      IF(PAVE(I).LE.0.) PAVE(I)  = 0.0
      RETURN
      END
C************************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE PONDRD(J,NPND,NPRT)
C    Reads the definition of pond number NPND (element J) from the pond
C    input file (unit IPOND): cross-section positions, an outlet rating
C    table (elevation vs discharge), and for each cross-section a table
C    of elevation, area and width.  Checks the pond limits (at most 10
C    sections and 10 elevations) and that the top of the section tables
C    matches the top of the rating table (error 21).
C    NOTE: in this version PAREAD never calls PONDRD (see NPNT there).
C-----------------------------------------------------------------------
      SUBROUTINE PONDRD(J,NPND,NPRT)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60),ASR(60), rfr(60), DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      COMMON /CNTRL/ nres, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /PONDS/ AP(10,10,3), TWP(10,10,3), EL(10,3), QSO(20,3),
     1   ELQ(20,3),ELZQ(3),ELVZ(3),NOPQ(3),MEL(3),NSEC(3),XSEC(10,3),
     2   ELZX(10), JOP(3), ELST(3), ALV, DIFUS(3), SIGMAS(60), NPART
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
C
      CHARACTER*6 ISUBR
      data isubr/'PONDRD'/
C
C  Mark element J as pond number NPND: no plane length, print level
C  NPRT.
      JOP(NPND)  = J
      NPN(J)   = NPND
      XL(J)    = 0
      NRP(J)   = NPRT
      CALL READTM(IPOND,13)                                               READTM
      READ(IPOND,*,END=10000) NUMP,NSEC(NPND),MEL(NPND),ELST(NPND),
     +                        ELVZ(NPND),DIFUS(NPND)
      IF(NUMP .NE. NPND)     CALL ERRPO(ISUBR,18,0,0,' NPND ','  X   ')
      IF(NSEC(NPND) .GT. 10) CALL ERRPO(ISUBR,20,0,0,' NPND ','  X   ')
      IF(MEL(NPND) .GT. 10)  CALL ERRPO(ISUBR,21,0,0,' NPND ','  X   ')
10000 NS = NSEC(NPND)
      CALL READTM(IPOND,6)                                                READTM
      DO 20 L=1,NS
        READ(IPOND,*,END=10001) NXS,XSEC(NXS,NPND)
   20 CONTINUE
10001 CONTINUE
      CALL READTM(IPOND,8)                                                READTM
      READ(IPOND,*,END=9002) NOPQ(NPND),ELZQ(NPND)
 9002 NQ = NOPQ(NPND)
      CALL READTM(IPOND,6)
      DO 25 L=1,NQ
        READ(IPOND,*,END=9003) ELQ(L,NPND),QSO(L,NPND)
   25 CONTINUE
 9003 IF(ELQ(1,NPND).GT.ELZQ(NPND)) GO TO 30
      ELQ(1,NPND) = ELZQ(NPND) + 1.0E-5
      QSO(1,NPND) = 1.0E-6
C  The pond takes its particle density from the element that feeds it.
   30 NV = MAX0(NU(J),NC1(J),NC2(J))
      RHOS(J) = RHOS(NV)
      CALL READTM(IPOND,10)                                             ! READTM
      DO 40 K=1,NS
       CALL READTM(IPOND,4)                                             ! READTM
       NEL = MEL(NPND)
       DO 50 L=1,NEL
         READ(IPOND,*) NXS,EL(L,NPND),AP(NXS,L,NPND),
     -            TWP(NXS,L,NPND)
   50  CONTINUE
       DIFFEL = ABS(EL(NEL,NPND) - ELQ(NOPQ(NPND),NPND))
       IF(DIFFEL .GT. 0.00001) THEN
         CALL ERRPO(ISUBR,21,0,0,' NPND ',' NXS  ')
       ENDIF
       CALL READTM(IPOND,1)                                             ! READTM
   40 CONTINUE
      RETURN
      END
C***************************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE INSPEC(NREAD,J,JLST)
C    Checks the input data for obvious errors.  NREAD selects the group:
C       1  run control record           5  rainfall times
C       2  options record               6  infiltration / erosion values
C       3  computation order            7  flow links between elements
C       4  geometry of element J
C    Problems are reported through ERRPO, which writes to the auxiliary
C    file and normally STOPs the run.
C-----------------------------------------------------------------------
      SUBROUTINE INSPEC(NREAD,J,JLST)
C
C      THIS SUBROUTINE INSPECTS THE INPUT DATA FOR ANY OBVIOUS ERRORS
C
      COMMON /CNTRL/ nres, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /INFIL/ AL(60), THI(60), THMX(60), ROC(60), RECS(60),
     1        THR, THFC, tempw, stone(60), RVAL
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), rfr(60), DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
C
      DIMENSION ICONNECT(60),ICNT(60)
C                                                     Case 7 ADDED 1/00TOK
      CHARACTER*6 ISUBR
      data isubr/'INSPEC'/
C
      Select Case (NREAD)
        Case (1)
        IF(NELE.EQ.0) CALL ERRPO(ISUBR,1,0,0,' NELE ','  X   ')
        IF(CLEN.EQ.0.)CALL ERRPO(ISUBR,1,0,0,' CLEN ','  X   ')
        IF(TFIN.EQ.0.)CALL ERRPO(ISUBR,1,0,0,' TFIN ','  X   ')
C      IF(NRES.EQ.0) CALL ERRPO(ISUBR,1,0,0,' NRES ','  X   ')
        IF(DELT.EQ.0.)CALL ERRPO(ISUBR,1,0,0,' DELT ','  X   ')
        IF(NELE.GT.60)CALL ERRPO(ISUBR,9,0,0,' NELE ',' NELE ')
C        IF(NOPT.GT.1) CALL ERRPO(ISUBR,3,0,0,' NOPT ',' NOPT ')
        Case (2) 
        IF(NTIME.EQ.0) CALL ERRPO(ISUBR,1,0,0,'NTIME ','  X   ')
        IF(NTIME.GT.2) CALL ERRPO(ISUBR,3,0,0,'NTIME ',' NTIME')
       Case ( 3) 
        IF(NLOG(NELE).EQ.0) CALL ERRPO(ISUBR,4,NLOG(NELE),NELE,'  X   ',
     &  '  X   ')
       Case(4) 
        IF(XL(J).LT.0.) CALL ERRPO(ISUBR,2,0,J,'  XL  ','  X   ')
        IF(W(J) .LT.0.) CALL ERRPO(ISUBR,2,0,J,'  W   ','  X   ')
        IF(S(J) .EQ. 0. .and. rillw(j) .gt. 0.)
     &      CALL ERRPO(ISUBR,2,0,J,'  S   ','  X   ')
        IF(w(j) .gt. 0.001 .and. SIR(J) .EQ. 0.) 
     &       CALL ERRPO(ISUBR,2,0,J,' SIR  ','  X   ')
        IF(sir(j) .LE. S(j) .and. rillw(j) .gt. 0.) Then
          Write(*,901) J
  901 Format(/' Element ',I2,': Interrill slope toward rills must be gre
     &ater than slope along rills.'/'  Default Used: Interrill slope = 1
     &.4*rill_slope') 
          SIR(J) = 1.4*S(J)
        End If
        IF(RILLW(J) .gt. 0. .and. R1(J).EQ. 0.)
     &   CALL ERRPO(ISUBR,2,0,J,'  R1  ','  X   ')
        IF(W(J) .Eq. 0.) Then
          IF(XL(J).EQ.0..AND.J.NE.JLST) CALL ERRPO(ISUBR,6,J,0,'  X   ',
     & '  X   ')
          IF(ZL(J).EQ.0.) CALL ERRPO(ISUBR,2,0,J,'  ZL  ','  X   ')
          IF(ZR(J).EQ.0.) CALL ERRPO(ISUBR,2,0,J,'  ZR  ','  X   ')
          IF(A(J) .LT.0.) CALL ERRPO(ISUBR,2,0,J,'  A   ','  X   ')
        Else
          IF(XL(J).EQ.0.) CALL ERRPO(ISUBR,3,0,J,'  XL  ','  X   ')
        End If
       Case (5) 
        IF(ND.EQ.0) CALL ERRPO(ISUBR,1,0,0,'  ND  ','  X   ')
        DO 501 I=1,ND
         IF(TI(I).EQ.0..AND.I.NE.1) CALL ERRPO(ISUBR,8,0,I,'  TI  ',
     & '  I   ')
  501   CONTINUE
       Case(6) 
        If(FMIN(J).gt. 1.0E-08) Then
          IF(AL(J).LT.0.) CALL ERRPO(ISUBR,2,0,J,'  AL  ','  X   ')
          IF(THMX(J).EQ.0.) CALL ERRPO(ISUBR,2,0,J,'THMAX ','  X   ')
        End If
        IF(NEROS .LE. 0) RETURN
        IF(D50(J) .LE. 0.0) CALL ERRPO(ISUBR,22,0,J,' D50  ','  X   ')
        IF(RHOS(J) .LE. 1.0) CALL ERRPO(ISUBR,23,0,J,' RHOS ','  X   ')
C-----------------------------------------------------------------------
C  Case 7 - inter-element flow logic.  LEL is the final (outlet)
C  element, IEL the element being tested.  Checks that:
C    - the outlet element does not feed any other element (error 25);
C    - no element feeds a downstream element through more than one
C      link, and none feeds more than one downstream element
C      (errors 26, 28);
C    - every element except the outlet feeds a downstream element
C      (error 27: a message only, the run is not stopped).
C-----------------------------------------------------------------------
       Case (7)
C
C     CHECK INTER ELEMENT FLOW LOGIC  
C                                        (LEL IS FINAL DOWNSTREAM ELEMENT)
C                                        (IEL IS CURRENT UPSTREAM ELEMENT)
C
  700 LEL = NLOG(NELE)
      DO 999 K = 1,NELE
        ICONNECT(K) = 0
C
C            IS THE LAST DOWNSTREAM ELEMENT USED AS INPUT TO ANY ELEMENT?
C
        IEL = NLOG(K) 
        IF((NU(IEL) .EQ. LEL).OR.
     1     (NR(IEL) .EQ. LEL).OR.
     2     (NL(IEL) .EQ. LEL).OR.
     3     (NC1(IEL).EQ. LEL).OR.
     4     (NC2(IEL).EQ. LEL)) THEN
          CALL ERRPO(ISUBR,25,LEL,IEL,' LAST ',' PRIOR')
        ENDIF
C
C            IS AN UPSTREAM ELEMENT USED AS INPUT TO ONE AND ONLY ONE
C                DOWNSTREAM ELEMENT AND IS THE UPSTREAM ELEMENT ONE AND
C                ONLY ONE LATERAL OR UPSTREAM INPUT TO THAT DOWNSTREAM 
C                ELEMENT?
C
        DO 998 I = 1,NELE
          IELD = NLOG(I)
          ICNT(I) = 0
          IF(IEL .EQ. NU(IELD))  ICNT(I) = ICNT(I) + 1
          IF(IEL .EQ. NR(IELD))  ICNT(I) = ICNT(I) + 1
          IF(IEL .EQ. NL(IELD))  ICNT(I) = ICNT(I) + 1
          IF(IEL .EQ. NC1(IELD)) ICNT(I) = ICNT(I) + 1
          IF(IEL .EQ. NC2(IELD)) ICNT(I) = ICNT(I) + 1
          IF(ICNT(I) .GT. 1) CALL ERRPO(ISUBR,26,IEL,I,'     ','      ')
          IF(ICNT(I) .EQ. 1) ICONNECT(K) = ICONNECT(K) + 1
  998   CONTINUE
C
C          IF ICONNECT IS GREATER THAN 1 THEN UPSTREAM ELEMENT
C             IS CONNECTED TO MORE THAN ONE DOWNSTREAM ELEMENT
        IF(ICONNECT(K).GT.1) CALL ERRPO(ISUBR,28,IEL,0,'     ','     ')
C
C          IF ICONNECT IS ZERO (EXCEPT FOR FINAL DOWNSTREAM ELEMENT) THEN
C             UPSTREAM ELEMENT IS NOT CONNECTED TO ANY DOWNSTREAM ELEMENT
        IF(K.EQ.NELE) GOTO 999
        IF(ICONNECT(K).EQ.0) CALL ERRPO(ISUBR,27,IEL,0,'      ','     ')
  999 CONTINUE
       Case Default
        CALL GOTOER('inspec')                                           !GOTOER
      End Select
      RETURN
      END
C*************************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE ERRPO(ISUBR,I,IVAR,KVAR,CIVAR,CKVAR)
C    Writes diagnostic message number I, raised in routine ISUBR, to the
C    auxiliary file (unit IDIAGN).  IVAR/KVAR and CIVAR/CKVAR are the
C    numbers and names inserted in the message text.
C    The run is stopped afterwards (STOP) except for message 13 and
C    message 27, which return to the caller.
C-----------------------------------------------------------------------
      SUBROUTINE ERRPO(ISUBR,I,IVAR,KVAR,CIVAR,CKVAR)
C
      CHARACTER ISUBR*6,CIVAR*6,CKVAR*6
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
C
C      DATA IUP/29/
C
      WRITE(IDIAGN,1000) I,ISUBR
      Select Case (i)
        Case (:0, 29:)
         WRITE(IDIAGN,1010)
        Case (1)
          WRITE(IDIAGN,1200)
          WRITE(IDIAGN,1100) CIVAR
          WRITE(IDIAGN,1220)
        Case (2)
          WRITE(IDIAGN,1200)
          WRITE(IDIAGN,1210) CIVAR,KVAR
        Case (3)
          WRITE(IDIAGN,1200)
          WRITE(IDIAGN,1300) CIVAR,KVAR
        Case (4)
          WRITE(IDIAGN,1400) KVAR,KVAR
        Case (5)
          WRITE(IDIAGN,1500) IVAR,KVAR
        Case (6)
          WRITE(IDIAGN,1160) IVAR
        Case (7)
          WRITE(IDIAGN,1170) CIVAR,CKVAR
        Case (8)
          WRITE(IDIAGN,1180) CIVAR,KVAR
        Case (9)
          WRITE(IDIAGN,1190) IVAR
        Case (10)
          WRITE(IDIAGN,2000) IVAR
        Case (11)
          WRITE(IDIAGN,2100) IVAR
        Case (12)
          WRITE(IDIAGN,2200)
        Case (13)
          WRITE (IDIAGN,2300) IVAR
          RETURN
        Case (14)
          WRITE(IDIAGN,2400) IVAR
        Case (15)
          WRITE(IDIAGN,2500) IVAR
        Case (16)
          Continue
        Case (17)
          WRITE(IDIAGN,2700) IVAR,KVAR
        Case (18)
          WRITE(IDIAGN,2800) IVAR
        Case (19)
          WRITE(IDIAGN,2900) IVAR,KVAR
        Case (20)
          WRITE(IDIAGN,3000) IVAR
        Case (21)
          WRITE(IDIAGN,3100) IVAR
        Case (22)
          WRITE(IDIAGN,3200)
        Case (23)
          WRITE(IDIAGN,3300)
        Case (25)
          WRITE(IDIAGN,3500) IVAR,KVAR
        Case (26)
          WRITE(IDIAGN,3600) IVAR,KVAR
        Case (27)
          WRITE(*,3700)IVAR
          WRITE(IWRITE,3700) IVAR
          WRITE(IDIAGN,3700) IVAR
          RETURN
        Case (28)
          WRITE(IDIAGN,3800) IVAR
        Case Default
          write(*,*)' errpo'
          CALL GOTOER('errpo ')                                         !GOTOER
      End Select
      WRITE(IDIAGN,9900)
      CALL FATALI('Error stop: see auxiliary file, error no.',I)
 1000 FORMAT(' ',79('*')//' ERROR NO.',I3,' CALLED FROM ',A6)
 1010 FORMAT(' ERROR NUMBER OUT OF RANGE.  CALLED FROM ERRPO.')
 1100 FORMAT(' IT APPEARS THAT NO VALUE HAS BEEN INPUT FOR THE VARIABLE
     1',A8,' .')
 1200 FORMAT(' DATA CARD ERROR',/)
 1210 FORMAT(' IT APPEARS THAT A ZERO VALUE HAS BEEN INPUT FOR THE VARIA
     1BLE',A6,' FOR WATERSHED ELEMENT J=',I3)
 1220 FORMAT(' AT LEAST UNDER THE CONDITIONS SPECIFIED, INPUT FOR THIS
     1VARIABLE IS REQUIRED.')
 1300 FORMAT(' THE VALUE INPUT FOR THE VARIABLE ',A6,' IS ILLEGAL. THIS
     1VALUE OR ARRAY ELEMENT SUBSCRIPT(IF VARIABLE IS AN ARRAY) IS ',I4)
 1400 FORMAT(' APPARENTLY ALL WATERSHED ELEMENTS HAVE NOT BEEN ASSIGNED
     1TO THE ORDER OF PROCESSING ARRAY NLOG.'/' THERE ARE ',I3,' ELEMENT
     2S,BUT NLOG(',I3,')=0.')
 1500 FORMAT(' DATA CARDS FOR GROUPS 1ST AND 2ND ARE FOR TWO DIFFERENT W
     1.S. ELEMENTS. ON 1ST, J=',I3,' AND ON 2ND J=',I3)
 1160 FORMAT(' ELEMENT J=',I3,' HAS BEEN SPECIFIED AS AN ADDER CHANNEL
     1(XL=0.) BUT IS NOT THE LAST ELEMENT TO BE PROCESSED.')
 1170 FORMAT(' THE ARRAY ',A6,' IS REQUIRED INPUT ON GROUP CARD ',A6,' .
     1 APPARENTLY IT IS MISSING FROM THE DATA INPUT')
 1180 FORMAT(' FOR ARRAY ',A6,' EITHER ELEMENT ',I3,' IS ILLEGAL(=0.) OR
     1 THE ENTIRE ARRAY HAS NOT BEEN INPUT')
 1190 FORMAT(' THE NUMBER OF GEOMETRIC ELEMENTS EXCEEDS 60.  -',I10,
     1'-')
 2000 FORMAT(' NO CONVERGENCE AFTER ',I4,' ITERATION STEPS.')
 2100 FORMAT(' DERIVATIVE OF  X(',I4,') = 0.')
 2200 FORMAT(' BLANK CARD READ.')
 2300 FORMAT(' 5 CONSECUTIVE NEGATIVE VALUES WERE OBTAINED FOR NEW VALUE
     +S IN THE ITERATIVE SOLUTION OF  X(',I3,')')
 2400 FORMAT(' ELEMENT ',I2,' HAS A PLANE AND CHANNEL(S) AT UPPER END')
 2500 FORMAT(' ELEMENT ',I3,' HAS NOT BEEN PROCESSED YET')
 2700 FORMAT('  AREA OR DEPTH ON ELEMENT ',I2,' AT X(',I2,') IS CALCULAT
     1ED TO EXCEED MAXIMUM ALLOWED ')
 2800 FORMAT(' ',1X,'POND NUM. ',I2,' IN THE PARAMETER FILE DOES NOT COR
     1RESPOND TO THE ',/,' ',3X,' POND NUM. IN THE POND FILE ')
 2900 FORMAT(' ',1X,'THE AREA-TOP WIDTH DATA FOR POND ',I2,' X-S ',I2,
     1/,' ',    ' IS NOT DEFINED W/R AN ELEV. = TO THE TOP ELEV. OF THE
     2DEPTH-DISCHARGE RATING TABLE')
 3000 FORMAT(' ',1X,'THE NUM. OF X-S IN POND ',I2,' MUST BE LESS THAN OR
     1 EQUAL TO 10')
 3100 FORMAT(' ',1X,'THE NUM. OF X-S ELEV. IN POND ',I2, ' MUST BE LESS
     1THAN OR EQUAL TO 10')
 3200 FORMAT(' ',1X,'THE D50 PARTICLE SIZE MUST BE GREATER THAN ZERO')
 3300 FORMAT(' ',1X,        'THE DENSITY OF THE SEDIMENT MUST BE GREATER
     1 THAN 1.0')
 3500 FORMAT(' ERROR: LAST DOWNSTREAM ELEMENT',I5,  ' CANNOT BE INPUT TO
     1 ELEMENT',I5,'. CHECK PARAMETER FILE')
 3600 FORMAT(' ERROR: UPSTREAM ELEMENT',I5,' IS USED AS MULTIPLE INPUT T
     &O ELEMENT',I5,'. CHECK PARAMETER FILE')
 3700 FORMAT(' WARNING: ELEMENT',I5,' IS NOT CONNECTED TO ANY DOWNSTREAM
     & ELEMENTS')
 3800 FORMAT(' ERROR: ELEMENT',I5,' IS CONNECTED TO MORE THAN ONE DOWNST
     &REAM ELEMENT')
 9900 FORMAT(1H ,79('*'))
      END
C*************************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE CONVER(KEY)
C    Unit conversion.
C      KEY = 1  (called from READER, after reading) to INTERNAL units:
C          times to seconds, rainfall to rate in m/s, FMIN to m/s
C      KEY = 2  (called at the end of the run) rainfall data back to
C          input units.
C-----------------------------------------------------------------------
      SUBROUTINE CONVER(KEY)
C
      COMMON /CNTRL/ nres, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60),rfr(60),DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /CHAN/ A1(20),A2(20),QUB(1000),AUB(1000), CO1, CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /INFIL/ AL(60), THI(60), THMX(60), ROC(60), RECS(60),
     1        THR, THFC, tempw, stone(60), RVAL
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /PONDS/ AP(10,10,3), TWP(10,10,3), EL(10,3), QSO(20,3),
     1   ELQ(20,3),ELZQ(3),ELVZ(3),NOPQ(3),MEL(3),NSEC(3),XSEC(10,3),
     2   ELZX(10), JOP(3), ELST(3), ALV, DIFUS(3), SIGMAS(60), NPART
C
C  ALL surface flow INTERNAL CALCULATIONS ARE DONE IN M-SEC UNITS
C  All infiltration calculations are in mm/sec
C  KEY = 1: input -> internal units.
      IF(KEY.NE.2) Then
C  TFAC converts the input time unit to seconds.  (NTIME = 3 would need
C  3600 here; it is rejected by INSPEC.)  FFAC = 3.6E6 converts m/s to
C  mm/h; HFAC = 1000 converts m to mm.
        IF(NTIME.EQ. 1) THEN
C  TIME IS IN SECONDS.
          TFAC=1.
        ELSE
C  NTIME=2, TIME IS IN MINUTES.
          TFAC=60.
        END IF
C  METRIC UNITS, DEPTH IN MM. INFIL IN MM/HR
        FFAC=3.6e6
        HFAC=1000.
C
C  Rainfall: the cumulative-depth data (mm) of each gage become the
C  intensity in m/s over each interval:
C     QIDD(K) = (depth(K+1) - depth(K)) / (1000 * (t(K+1) - t(K)))
C  and the last interval is set to zero.  RNRMX = peak intensity.
C Convert gage depths to rates in m/sec
        RNRMX = 0.
        DO 11 I=1,NGAGES
          TIDD(I,1)=TIDD(I,1)*TFAC
          DO 10 K=1,ND-1
            TIDD(I,K+1)=TIDD(I,K+1)*TFAC
C  qidd becomes precip in m/sec
       QIDD(I,K)=(QIDD(I,K+1)-QIDD(I,K))/(HFAC*(TIDD(I,K+1)-TIDD(I,K)))
            If(qidd(I,K) .gt. rnrmx) rnrmx = QIDD(I,K)
 10       Continue
 11     QIDD(I,ND)=0.
C  Final infiltration rate FMIN (mm/h) -> FMINEW (m/s), every element.
        DO 20 IN=1,NELE
          I = NLOG(IN)
          FMINEW(I)=FMIN(I)/FFAC
C fminew() is in m/sec
C          RNDPTH(I) = RNDPTH(I)/HFAC  ! done elsewhere
 20     CONTINUE
C  Run length and time step to seconds.
        TFIN=TFIN*TFAC
        DELT=DELT*TFAC
        RETURN
      End If
C
CR!     KEY=2.. CONVERT UNITS BACK FOR CONVENIENT OUTPUT
C  KEY = 2: rainfall time back to input unit, intensity back to mm/h.
      If(NTIME .eq. 1) Then    ! Seconds in input
        DIV=1.
      Else If(NTIME .eq. 2) Then    ! MInutes in input
        DIV=60.
      Else If(NTIME .eq. 3) Then     ! Hours in input
        DIV=3600.
      Else
       CALL FATAL('NTIME value out of range')
      End If
      DO II=1,NGAGES
        DO I=1,ND
          TIDD(II,I)=TIDD(II,I)/DIV
          QIDD(II,I)=QIDD(II,I)*ffac
        End Do
      End Do
      RETURN
      END
C***************************************************************************
C
C  Called when a selector falls through a computed GO TO or SELECT CASE;
C  writes a message to the diagnostic file and STOPs.
	SUBROUTINE GOTOER(namef)
C***************************************************************************
C  THIS IS AN ESCAPE ROUTINE CALLED WHEN ARGUMENTS FALL THRU
C  COMPUTED "GO TO" STATEMENTS. THIS EMULATES THE RUN-TIME
C  TERMINATION THAT WOULD OCCUR IF SUCH AN EVENT OCCURRED
C  WHILE EXECUTING A FORTRAN-IV PROGRAM.
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      character*6 namef
C
	WRITE(IDIAGN,100)namef
100     FORMAT(/,5X,'Parameter out of range in module',a6)
	CALL FATAL('Selection parameter failed - see DIAGN output')
	END
c**************************************************************
C-----------------------------------------------------------------------
C  SUBROUTINE OUTPUT(J)
C    Writes the input data of element J (geometry, infiltration, rills,
C    cover, soil and erosion parameters) to the static summary file
C    (unit 99), under the run title for element 1.  Called from MAIN
C    before the element is computed.
C-----------------------------------------------------------------------
      subroutine output(J)
c     added by jq feb 93 to provide a static summary file
      COMMON /CNTRL/ nres, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /INFIL/ AL(60), THI(60), THMX(60), ROC(60), RECS(60),
     1        THR, THFC, tempw, stone(60), RVAL
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), rfr(60), DEPNO(60)
     4      ,rs(60), dero(60), sir(60)
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000),Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /ESEM/ COVER(60), CUINT(60), ISHAPE(60), PLHGT(60),
     1        VPANG(60),PBASE(60)
      COMMON /PONDS/ AP(10,10,3), TWP(10,10,3), EL(10,3), QSO(20,3),
     1   ELQ(20,3),ELZQ(3),ELVZ(3),NOPQ(3),MEL(3),NSEC(3),XSEC(10,3),
     2   ELZX(10), JOP(3), ELST(3), ALV, DIFUS(3), SIGMAS(60), NPART
C
      common /labels/ title, tlabl, qlabl, slabl, ifsed
      Character*12 qlabl, tlabl, slabl
      CHARACTER TITLE*80 
C
      if (j.eq.1) write(99,100)TITLE
C  STODJ is computed from the roughness index RFR with a formula the
C  original author "borrowed from PLANE"; it is only printed in the
C  summary below.
C  Borrow from PLANE for output purposes:
      STODJ = 0.001*exp(-6.66+0.27*rfr(j))                          !6/93
      write(99,101)J,NU(J),W(J),XL(J),S(J),R1(J),FMIN(J),al(j)
     - ,PORpl(j),thi(j),
     - THMX(J),ROC(J),RECS(J)*1000,DINTR(J),DEPNO(J),RS(J),RFR(J),
     - ZRL(J),RILLW(J),RILLD(J),COVER(J),ISHAPE(J),VPANG(J),
     - PBASE(J),PLHGT(J),D50(J), ERODGJ(J),SPLTEX(J),COHE(J),RHOS(J),
     - PAVE(J),SIGMAS(J), SIR(J), DERO(J), R2(J),stodj
 100  format(//,
     - 10X,'----------------------------------------------',/,
     - 10X,'|                                            |',/,
     - 10X,'*       EUROSEM 3 STATIC SUMMARY FILE        *',/,
     - 10X,'|                                            |',/,
     - 10X,'----------------------------------------------',///,
     - 10X,'RUN TITLE:',A80)
C
 101  format(/,
     - 6X,'INPUT DATA FOR ELEMENT ',I3, / ,
     - 6X,'=========================', / ,
     - ' ',2x,'   NU:',3X,I3,/,
     - ' ',2x,'    W:',3X,F9.2,' M',3X,'   XL:',3X,F9.2,' M',3X,'    S:'
     -,F9.2, / ,
     -' ',2x,' MANN:',3X,F9.3,'  ',3X,' FMIN:',3X,F9.2,' MM/HR'
     -,3X,'G:',F9.2,' MM',/,
     -' ',2x,'  POR:',3X,F9.2,'  ',3X,'  THI:',3X,F9.2,'  ',3X,' THMX:'
     -,F9.2,/,
     -' ',2x,'  ROC:',3X,F9.2,'  ',3X,' RECS:',3X,F9.2,' MM',2X,'DINTR:'
     -,F9.2,' MM',/,
     -' ',2x,'DEPNO:',3X,F9.2,'  ',3X,'   RS:',8X,f3.1,'  ',4X,'  RFR:',
     -F9.3,' % ',/
     - ' ',2x,'  ZLR:',3X,F9.2,'  ',3X,'RILLW:',3X,F9.2,' M',3X,'RILLD:'
     -,F9.2,' M',/,
     - ' ',2x,'COVER:',3X,F9.2,'  ',3X,'SHAPE:',6X,I3,'  ',6X,' PANG:'
     -,F9.2,' o',/,
     -' ',2x,'PBASE:',3X,F9.2,'  ',3X,'PHEIG:',3X,F9.2,' M',3X,'  D50:'
     -,F9.2,' um',/,' ',2x,' EROD:',3X,F9.2,' G/J',1X,'SPLTX:',3X,F9.2,
     -' ',3X,'   COH:',F9.2,' KPA',/,' ',2x,' RHOS:',3X,F9.2,'kgm3',1X,
     - ' PAVE:',3X,F9.2,'  ',3X,'SIGMA:',F9.2/
     -,' ',4X,'SIR:',3X,F9.3,6x,'DERO:   ',F9.2,' m'/
     -'  Derived parameters: ',3X,'MN(IR):   ',F9.3,'  SurfStor:',f11.6
     -,' mm'/)
      RETURN
      END

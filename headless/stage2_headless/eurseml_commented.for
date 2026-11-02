C-----------------------------------------------------------------------
C  EURSEML.FOR -- vegetation, interception and erosion physics (EUROSEM)
C
C    VEGWRI   writes the interception table for an element
C    RILLWI   changes rill width and depth with erosion / deposition
C    CAPAC    sediment transport capacity (Govers; Everaert for
C             interrill)
C    DRIPS    rainsplash and leaf-drip detachment
C    KINET    kinetic energy of rain and of leaf drip
C    INTERC   interception by the plant cover; throughfall, drip, stem
C             flow
C    SHIELDF  Shields-type critical value for particle movement
C-----------------------------------------------------------------------
C Old subroutine TEXTU removed since not used in latest version
C-----------------------------------------------------------------------
C  SUBROUTINE VEGWRI(J)
C    Writes the interception results for element J to unit 8 (written
C    as a literal, the same file as IDIAGN) if IVOUT is not 0: for each
C    rainfall interval the time, rain, throughfall, leaf drip, stem
C    flow, vegetation storage and the closure error, all in mm.
C-----------------------------------------------------------------------
      SUBROUTINE VEGWRI(J)
C
C     WRITES OUTUT FILES FOR THE ESEM MODEL
C
      COMMON /ESEM/ COVER(60), CUINT(60), ISHAPE(60), PLHGT(60),
     1        VPANG(60),PBASE(60)
      COMMON /CNTRL/ nres, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /EURO/ DINTRT(100), STEMD(100), DRIPD(100), TFALL(100),
     1        RAIN(100), RHRMM(100)
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
C
      dimension dintmm(100)
C
C
C     FOR EACH ELEMENT WHERE OUTPUT REQUESTED
C
      IF (ivout.EQ.0)THEN
            GOTO 10
      ELSE
C                          HEADINGS
           WRITE(8,1000)J
           WRITE(8,1003)
           WRITE(8,1001)
           WRITE(8,1002)
C
C     COMPUTE DEPTHS PER TIME DEPTH PAIRS
C
           DO 20 K=1,ND
C
C     WORK OUT PERCENT ERROR
C
C  Storage depth from m to mm.  VEGPER is the part of the interval rain
C  that is not accounted for by throughfall, drip, stem flow and
C  storage; ER expresses it as a percent of the rain.
      Dintmm(k)=dintrt(k)*1000.
      VEGPER=RAIN(K)-(TFALL(K)+DRIPD(K)+STEMD(K)
     &            +DINTmm(K))
      IF (abs(VEGPER).GT.0.0) THEN
      ER = 100.*vegper/RAIN(K)
      ELSE
         ER=0.0
      END IF
C
C     WRITE OUT DATA FOR EACH TIMEDEPTH PAIR
C
        WRITE(8,1004)TI(K+1),RAIN(K),TFALL(K),DRIPD(K),
     & STEMD(K),DINTmm(K),ER
C
 1000       FORMAT(//,25X,'INTERCEPTION DATA FOR ELEMENT',I3/2X)
 1003       FORMAT(20X,'ALL DATA EXPRESSED AS MM PER TIME DEPTH PAIR')
 1001       FORMAT(/,1X,'TIME',4X,' RAIN ',4X,'TFALL ',4X,' DRIP ',4X,
     & ' STEM ',4X,'VEGSTORE',4X,'ERROR')
 1002       FORMAT(1X,'MIN ',4X,'  MM  ',4X,'  MM  ',4X,'  MM  ',4X,
     & '  MM  ',4X,'  MM  ',4X,'   % ')
 1004       FORMAT(1X,F4.0,4(4X,F6.3),2(4X,F7.5))
   20       CONTINUE
        END IF
   10 CONTINUE
      RETURN
      END
C**********************************************************************
C-----------------------------------------------------------------------
C  SUBROUTINE RILLWI(J,K,DELA)
C    Changes the rill depth, width and width ratio at node K after
C    erosion or deposition.  DELA(K) is the area of sediment removed
C    from (+) or added to (-) the rill cross-section.
C      Erosion:    shared between the walls and the bottom until the
C                  rill depth reaches DERO (the depth of the erodible
C                  layer); after that only the walls move (widening).
C      Deposition: added to the bottom, making the rill shallower and
C                  narrower.  Filling a rill beyond its depth STOPs
C                  the run.
C    Also updates the interrill width WINTER and the rill width ratio
C    RWR (rill width / element width).
C-----------------------------------------------------------------------
      SUBROUTINE RILLWI(J,K,DELA)
C
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1, CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
C
      dimension dela(20)
C
C     ROUTINE CHANGES RILL DEPTH, WIDTH, AND RWR (RILL WIDTH RATIO=
C     RILLWIDTH/ ELEMENT WIDTH
C
C        DELA = LOCAL AREA OF SEDIMENT REMOVED FROM CIRCUMFERANCE
C                OF RILL
C
         IF(DELA(K) .GT. 0.) THEN
C  *NET EROSION AT SECTION. ASSUME EQUAL TAKEN FROM WALLS AND BOTTOM
C      UNTIL DERO IS REACHED
           If(RILLD1(K) .lt. DERO(J)) Then
             DD = DELA(K)/(RILLW1(K) + 2.*RILLD1(K))
             DW = 2.*DD
           Else
             DD = 0.
             DW = DELA(K)/(2.*RILLD1(K))
           End If
           RILLW1(K)=RILLW1(K) + DW
           RILLD1(K)=RILLD1(K) + DD
C           h2(k) = h2(k)*(1.-2.*dd/tw(k))                              !N8/91
         ELSE
C  *NET DEPOSITION. Add to bottom.
           dh = DELA(K)/RILLW1(K)
C  Note: In this case dh is negative, reduce depth of rill:
           RILLD1(K) = RILLD1(K) + DH                                   !N11/91
           rillw1(k) = rillw1(k) - 2.*zrl(j)*dh                         !12/92
           If(rilld1(k) .lt. 0.) Then
             write(*,99) K, dela(k), dh, rilld1(k)
 99   Format(/' dA for rill',i2,' is ',f8.4,', giving dh and rilld: ',
     &  2f8.5)
            write(*,*) dela(1),rilld1(k),tw(1),dh
            CALL FATAL('Excess rill fill. Check erosion parameters')
C            rilld1(k) = 0.
           End If
C           h2(k) = h2(k) + (1.-rillw1(k)/tw(k))*dh                     !N8/91
         END IF
         winter(k) = w(j)/depno(j)- RILLW1(k)
C
         RWR(k)=RILLW1(K)/W(J)
      RETURN
      END
C************************************************************************
C-----------------------------------------------------------------------
C  SUBROUTINE CAPAC(LRL,SPWR,AH,J,NZR,NM,CMX)
C    Sediment transport capacity CMX (volume concentration, m3/m3) at
C    nodes NZR to NM, from the stream power SPWR.
C      LRL   0 interrill area, 1 rills
C      SPWR  stream power at each node (cm/s for Govers, otherwise the
C            effective stream power used by Everaert)
C    Calculated in cgs units (G = 981 cm/s2).
C-----------------------------------------------------------------------
      SUBROUTINE CAPAC(LRL,spwr,AH,J,nzr,NM,cmx)
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1, CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
C      COMMON /EROSN/ V1(20), V2(20), CT(20), CMX(20), RDET(20), DLA(20),
C     1      DrlNET(20), Dirnet(20), suspir, susprl, SPOW(20)
      dimension ah(20),spwr(20),cmx(20)
      DATA g/981./
C
C LRL is 0 for interrill area, 1 for Rills
C
C-----------------------------------------------------------------------
C  RILLS, or MCODE = 0 (interrill too): Govers (1987) unit stream
C  power.  Below the critical value of 0.4 cm/s the capacity is zero;
C  above it (D50 in microns):
C      CMX  = ALFA * (SPWR - 0.4) ** BETA
C      ALFA = (2.1 / (D50 + 5)) ** 0.85
C      BETA = ((D50 + 5) / 280) ** 0.35
C-----------------------------------------------------------------------
      DO 77 K=nzr,nm
        If(LRL .ge. 1 .or. mcode .eq. 0) Then
C   RILLS USE EQUATIONS OF GOVERS unit stream power
          IF (SPWR(K) .GE. 0.4) THEN
C    New continuous method based on Govers (1987)
             dfifty = d50(j)
             alfagx = (2.1/(dfifty+5.))**(0.85)     !  revised 7/97
             betagx = ((dfifty+5.)/280.)**(.35)     !  revised 7/97
             cmx(k) = alfagx*(spwr(k)-0.4)**betagx
          ELSE
            CMX(K)=0.0
          END IF
C-----------------------------------------------------------------------
C  INTERRILL with MCODE >= 1: Everaert (1991).  The critical effective
C  stream power SPOWC comes from a Shields-type critical value (SHIELDF)
C  for the particle Reynolds number.  Above it, the sediment discharge
C  QSW (g/cm/s) follows Everaert's continuous relation, converted to a
C  volume concentration by dividing by the particle density RHOS and the
C  unit discharge WINQ, and limited to 0.3.
C-----------------------------------------------------------------------
        Else If(winq(k) .gt. 0.) Then
C This is Everaert's equations, used for interrill only:
 !         hcm = ah(k)*100./tw(k)
          hcm = winh(k)*100.
          vcm = 0.
          yc = 0.
          d50cm = d50(j)*1.e-4
          If(hcm .gt. 0.) Then
            vcm = winq(k)/hcm
C  Interrill flows Use Equations of Everaert(1991)
            spowc = 0.
C  Proposed critical value of Bagnold shear, based on Shields':
            ustar = sqrt(g*hcm*sir(j))
            reyp = ustar*d50cm/(xnu*1.e4)
            yc = shieldf(reyp)                                          !SHIELDF
            bspc = yc*(rhos(j)-1.)*g*d50cm*vcm
            Spowc = (0.5*bspc)**1.5/(hcm**.667)                         !res
          End If
          If(spwr(k) .gt. spowc .and. winq(k) .gt. 0.) Then
C gm/cm/sec:
C  Proposed continuous model for Everaert Data:
            cona = (19. - d50(j)/30.)*1.e-4
            qsw = cona*((spwr(k)-spowc)**(0.2/1.4) - 1.)**5.           
C winq is sq cm/sec: convert to conc
            cmx(k) = qsw/rhos(j)/winq(k)
            If(cmx(k) .gt. 0.3) cmx(K) = 0.3
          Else
            cmx(k) = 0.
          End If
C          If(k .ge. nm) write(99,'(9x,3g12.5)') hcm,hca
        Else
          cmx(k) = 0.
        End If
        IF (CMX(K).LE.0.0) THEN
          CMX(K)=0.0
        END IF
  77  CONTINUE
      RETURN
      END
C************************************************************************
C-----------------------------------------------------------------------
C  SUBROUTINE DRIPS(J,NM,LRL)
C    Rainsplash detachment RDET at each node, from the kinetic energy of
C    rain (RAINKE) and of leaf drip (DRIPKE):
C        detachment = KE * rain rate * exp(-SPLTEX * water depth in mm)
C                     * soil erodibility ERODGJ
C    (g/m2/s), converted to m of soil per second by dividing by a
C    particle density of 2.65 g/cm3.  The exponential reduces splash
C    where the water layer is deep.
C      LRL = 0  interrill area: depth WINH
C      LRL = 1  rill area: depth RILLH, weighted by the rill width
C               ratio RWR (interrill by 1 - RWR)
C    IZ is the index of the current rainfall interval.
C    NOTE: RMMSEC is the GROSS rainfall rate (from KINET), so the rain
C    term applies to all the rain, not only the direct throughfall.
C    Worth checking against the EUROSEM documentation.
C-----------------------------------------------------------------------
      Subroutine DRIPS(J,nm,LRL)
C
      COMMON /KINETI/DRIPKE(100),RAINKE(100)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000),Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      COMMON /EROSN/ V1(20), V2(20), CT(20), CMX(20), RDET(20), DLA(20),
     1      DrlNET(20), Dirnet(20), suspir, susprl, SPOW(20)
      common /intrill/ intrt, kint, d1(20), d2(20), vir(20), ql3(20), 
     1                 cbal, dcx, slrat, izr
C
      COMMON /EURO/ DINTRT(100), STEMD(100), DRIPD(100), TFALL(100),    !5/95
     1        RAIN(100), RHRMM(100)                                     !
C
      dimension rillh2(20), xinh2(20), RDETRL(100),RDETIN(100),
     1      DRDTRL(100), DRDTIN(100), CUDTDR(20), 
     2      XINDET(20), CUDTRN(20),spro(20),cmu(20)
      SAVE IZ
C
C  LRL is 0 for interrill, 1 for rill area
C
C  Find the rainfall interval IZ for the current time (restarted at
C  each element, L <= 2).
       IG = IGAGES(J)                                                   !CORRS
       IF (L.LE.2) THEN
          IZ=1
       END IF
       IF (T(L).GT.TIDD(IG,IZ+1)) THEN                                  !corrs
          IZ=IZ+1
       END IF
C  RMMSEC: rain rate in mm/s.  DMMSEC: the part of it that falls as
C  leaf drip.
       rmmsec = rhrmm(IZ)/3600.                                         !5/95
       dmmsec = 0.                                                      !
       If(Rain(IZ) .gt. 0.) DMMSEC = RMMSEC*DRIPD(IZ)/RAIN(IZ)          !
       DO 77 K=1,nm
C
         RDET(k) = 0.
C      FOR THE RILL AREA ONLY
C  Rill area: rain and drip detachment, damped by the rill flow depth,
C  and weighted by the rill width ratio RWR.  Converted from g/m2/s to
C  m/s of soil.
       If(LRL .eq. 1) Then
C                                    Convert to mm
         RILLH2(K)=RILLH(K)*1.e3
         Dampf = exp(-SPLTEX(J)*RILLH2(K))                              !5/95
         RDETRL(K)=RAINKE(IZ)*rmmsec*DAMPF*RWR(K)*ERODGJ(J)
         DRDTRL(K)=DRIPKE(IZ)*DMMSEC*DAMPF*RWR(K)*ERODGJ(J)
         CUDTRN(K)=RDETRL(K)+RDETIN(K)+CUDTRN(K)
C
C Convert g/M**2/sec into M/sec
         RDET(k) = rwr(k)*(RDETRL(k)+DRDTRL(k))/2.65*1.e-6
       End If
C
C  Interrill area: the same, using the interrill flow depth WINH
C  (cm to mm).
C    FOR THE Interrill CASE :
C
C get interrill depth from cm. to mm:
         XINH2(K)=WINH(K)*10.
         Dampf = exp(-SPLTEX(J)*XINH2(K))
         RDETIN(K)=(RAINKE(IZ)*rmmsec*DAMPF)*(1.-RWR(K))*erodgj(j)      !5/95
         DRDTIN(K)=DRIPKE(IZ)*DMMSEC*DAMPF*(1.-RWR(K))*erodgj(j)        !5/95
         CUDTDR(K)=DRDTRL(K)+DRDTIN(K)+CUDTDR(K)
         XINDET(K) = (RDETIN(k)+DRDTIN(k))*1.e-6/2.65    ! in m/sec
C this is redone now to take balanced input from interrill area (3/92)
C-----------------------------------------------------------------------
C  Rills fed by interrill areas too short to route explicitly
C  (INTRT = 0): a steady approximation of the sediment concentration
C  CBAL that the interrill area delivers into the rill.  It compares
C  the transport capacity of the interrill flow (CAPAC) with the splash
C  supply, allowing for deposition at the settling velocity VS reduced
C  by the cohesion factor FBETA.
C-----------------------------------------------------------------------
         If(lrl .ge. 1. .and. INTRT .LE. 0) Then                      !12/92
C Short unrouted interrills indicated:
C    first find transport capacity of interrill stuff entering rill:
          If(mcode .ge. 1) Then
C                     Bagnold effective stream power (Everaert):
            spowb = grav*100.*winq(k)*sir(j)
            spro(1) = spowb**(1.5)/(winh(k))**.667
          Else
C                     Govers (rills) metric stream power (cm/sec):
            virr = winq(k)/winh(k)
            SPRO(k) = virr*sir(j)*100.
          End If
          call capac(0,spro,rillh2,j,1,1,cmu)
!           
           CBAL = xindet(k)/vs   ! this is rain splash transport
!                 This is a steady approx.
           fq = 1.
           If(ql2(k) .gt. 0.) Then
           ratl = 1./(sin(atan(sir(j)/s(j))))
           wqx = ql2(k)/ratl
           crat = (cmu(1)-cbal)/cmu(1)
           cgu = vs*fbeta(crat,cohe(j))
           splsh = cbal*wqx
           ap = splsh + cgu*cmu(1)
           bp = wqx + cgu
           fq = exp(-bp*2./wqx)
           cout = (ap + (bp*cbal - ap)*fq)/bp
           cbal = cout
           end if
C           write(99,'(6g12.4)') cbal,cmu(1),winq(k),winh(k),fq,cout
         Else If(LRL .le. 0) Then                                       !
           Rdet(k) = xindet(k)  + rdet(k)                               !
         End If                                                         !
C
  77  CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C  SUBROUTINE KINET(J)
C    Kinetic energy of rain and of leaf drip for each rainfall
C    interval K:
C      RHRMM(K)   gross rainfall intensity (mm/h)
C      RAINKE(K)  rain energy per mm of rain (J/m2/mm):
C                    8.95 + 8.44 * LOG(intensity)
C                 zero if the interval's rain depth is 0.087 mm or less
C      DRIPKE(K)  leaf-drip energy: (15.8 * sqrt(PH) - 5.87) times
C                 the drip depth (mm), PH = plant height in m; zero
C                 if PH is zero
C    CHANGED (2026-10-29): the code originally used LOG, the NATURAL
C    logarithm.  Brandt (1990) uses log to base 10, as confirmed from
C    the EUROSEM documentation, so the line now uses LOG10.  The cut-off
C    0.08701 is where 8.95 + 8.44*log10(I) = 0, which supports this.
C    Effect: lower rain energy (at 100 mm/h, 25.8 instead of 47.8
C    J/m2/mm) and so less splash detachment.  Sediment results will fall
C    against the original baseline; water results should not change.
C    NOTE: the 0.08701 test is applied to RAIN(K), the depth in the
C    interval, not to the intensity RHRMM(K).  The clamp to zero below
C    covers negative values, so the result is the same.
C-----------------------------------------------------------------------
      SUBROUTINE KINET(J)
C
      COMMON /ESEM/ COVER(60), CUINT(60), ISHAPE(60), PLHGT(60),
     1        VPANG(60),PBASE(60)
      COMMON /KINETI/DRIPKE(100),RAINKE(100)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /EURO/ DINTRT(100), STEMD(100), DRIPD(100), TFALL(100),
     1        RAIN(100), RHRMM(100)
C
      DO 77 K=1,ND
         KP = K + 1
         ramimm = RAIN(K)/(TI(KP)-TI(K))  ! mm/min
C
         RHRMM(K)=ramimm*60.0      !  mm/h
C
         IF (PLHGT(J).LE.0.0) THEN
            DRIPKE(K)=0.0
         ELSE
           IF (RAIN(K).GT.0.0) THEN
             DRIPKE(K)=((15.8*SQRT(PLHGT(J)))-5.87)*DRIPD(K)
     &                /RAIN(K)
           ELSE
             DRIPKE(K)=0.0
           END IF
         END IF
         IF (RHRMM(K).LE.0.0) THEN
            RAINKE(K)=0.0
         ELSE
           IF (RAIN(K).GT.0.08701) THEN
             RAINKE(K)=(8.95+(8.44*LOG10(RHRMM(K))))
     &                *TFALL(K)/RAIN(K)
           ELSE
              RAINKE(K)=0.0
           END IF
         END IF
         IF (DRIPKE(K).LT.0.0)  DRIPKE(K)=0.0
         If (rainke(k) .lt. 0.0) rainke(k) = 0.
77    CONTINUE
      RETURN
      END
C*************************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE INTERC(J)
C    Interception by the plant cover on element J (written 19/4/90 by
C    J.N. Quinton).  From the maximum interception depth DINTR (mm) and
C    the cover fraction COVER it finds, for each rainfall interval:
C      - the water held on the vegetation (storage), which grows with
C        the cumulative rain CUMD (mm):
C            storage = COVER * DINTR * (1 - exp(-cumulative rain))
C      - DINTRT(K)  interception in the interval (returned in m),
C        CUINT(J)  total interception (mm)
C      - the rain intercepted but not stored is split into stem flow
C        (by plant shape ISHAPE and angle VPANG) and leaf drip DRIPD;
C        TFALL is the direct throughfall (rain not hitting the cover)
C    Finally VEGWRI writes the table.
C-----------------------------------------------------------------------
      SUBROUTINE INTERC(J)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), rilld(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /ESEM/ COVER(60), CUINT(60), ISHAPE(60), PLHGT(60),
     1        VPANG(60),PBASE(60)
      COMMON /EURO/ DINTRT(100), STEMD(100), DRIPD(100), TFALL(100),
     1        RAIN(100), RHRMM(100)
C
      DIMENSION vegsto(100),dinter(100),dripst(100)
C
C     THIS SUBROUTINE TAKES A SINGLE VALUE OF MAXIMUM INTERCEPTION DEPTH
C     FOR EACH PLANE (DINTR(J); mm), AND THE % COVER (COVER(J)). IT THEN
C     COMPUTES THE AMOUNT OF VEGETATIVE STORAGE AS A FUNCTION OF
C     CUMALATIVE RAIN. IT RETURNS A VALUE FOR THE DEPTH OF INTERCEPTION
C     (DINTRT(K))FOR EACH TIME DEPTH PAIR
C     N.B ALL CALCULATIONS WILL BE DONE IN FEET. WHEN WRITTING NEED TO
C     CONVERT
C     WRITTEN 19/4/90 BY J.N.QUINTON
C
C  Plant angle from degrees to radians (VPANG(J) is changed in place, so
C  this must only be done once per element).
      VPANG(J)=(VPANG(J))*(3.142/180.0)    ! Convert to radians
C
      IGE = IGAGES(J)
      DICMAX=DINTR(J)
      CUINT(J)=0.0
      VEGSTO(1) = 0.
C------------------------------------------------------------------------
      DO 10 K=2,ND
       KM = K - 1
      RAIN(KM) = cumd(ige,K) - cumd(ige,KM)   !Interval raindepth in mm
C     QIDD = rate thru interval K in m/sec
C
       dinter(KM)=COVER(J)*RAIN(KM)
C
c     DINTRT = RETURNED INTERCEPTION DEPTH
c     dinter= ACTUAL AMOUNT INTERCEPTED in INTERVAL K
C
       vegsto(K)=COVER(J)*DICMAX*(1-EXP(-cumd(ige,k)))
       DINTRT(KM)=vegsto(K)-vegsto(KM)
       CUINT(J)=CUINT(J)+DINTRT(KM)
C
C     THIS PART OF THE ROUTINE SPLITS THE RAIN INTO STEMFLOW,DRIP AND
C     DIRECT TFALL.JNQ 3.5.90
C
       dripst(KM) = dinter(KM)-DINTRT(KM)
C  Stem flow by plant shape: 0 none; 1 and 2 a fraction of the
C  intercepted water that depends on the plant angle.  Any other value
C  is replaced by 1 with a message.
  20   IF (ISHAPE(J).EQ.1) THEN
         STEMD(KM) = dripst(KM)*(COS(VPANG(J)))*
     &   ((SIN(VPANG(J)))**2)
        ELSE IF (ISHAPE(J).EQ.2) THEN
          STEMD(KM)=dripst(KM)*(COS(VPANG(J)))
        Else If (ISHAPE(J) .eq. 0) Then
          STEMD(KM) = 0.
       ELSE
         WRITE(*,*)'BAD PLANT SHAPE FACTOR. 1 ASSUMED'
         ISHAPE(J)=1
         GOTO 20
       END IF
       DRIPD(KM)=dripst(KM)-STEMD(KM)
       IF (RAIN(KM).LE.0.0) THEN
         TFALL(KM)=0.0
       ELSE
         TFALL(KM)= RAIN(KM)-DINTER(KM)
       END IF
C     CONVERT MM TO M
C
       DINTRT(KM)=DINTRT(KM)/1000.
  10  CONTINUE
C     --------------------------------------------------------------------
      CALL VEGWRI(J)       
      RETURN
      END
C*************************************************************************
C-----------------------------------------------------------------------
C  FUNCTION SHIELDF(REYN)
C    Dimensionless critical value for particle movement as a function
C    of the particle Reynolds number (a piecewise approximation of the
C    Shields curve):
C        REYN <= 10          0.08 / REYN**0.40
C        10 < REYN <= 500    0.022 * REYN**0.16
C        REYN > 500          0.060
C-----------------------------------------------------------------------
      FUNCTION SHIELDF(REYN)
C* DIMENSIONLESS CRITICAL TRACTIVE FORCE FOR PARTICLE MOVEMENT:
      IF(REYN.GT.10.) GO TO 61
      SHIELDF = 0.08/REYN**0.40
      GO TO 63
   61 IF(REYN.GT.500.) GO TO 62
      SHIELDF = 0.022*REYN**0.16
      GO TO 63
   62 SHIELDF = 0.060
   63 CONTINUE
      RETURN
      END

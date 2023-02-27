C Old subroutine TEXTU removed since not used in latest version
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
            stop' excess rill fill. Check erosion parameters'
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
       IG = IGAGES(J)                                                   !CORRS
       IF (L.LE.2) THEN
          IZ=1
       END IF
       IF (T(L).GT.TIDD(IG,IZ+1)) THEN                                  !corrs
          IZ=IZ+1
       END IF
       rmmsec = rhrmm(IZ)/3600.                                         !5/95
       dmmsec = 0.                                                      !
       If(Rain(IZ) .gt. 0.) DMMSEC = RMMSEC*DRIPD(IZ)/RAIN(IZ)          !
       DO 77 K=1,nm
C
         RDET(k) = 0.
C      FOR THE RILL AREA ONLY
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
             RAINKE(K)=(8.95+(8.44*LOG(RHRMM(K))))
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



C************************************************************************
C
      SUBROUTINE ADD (J)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT, IRDATA
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1,CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(2000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(2000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
C
C     THIS ROUTINE DOES TWO THINGS  1) IT CALCULATES THE UPPER BOUNDARY
C     AREA FOR INPUT INTO THE CHANNEL ROUTINE FOR CASES WHERE THERE IS A
C     EITHER A CONVERGING PLANE OR CHANNELS CONTRIBUTING TO THE UPPER
C     BOUNDARY OF CHANNEL J...2)IT ADDS TOGETHER ALL LATERAL INFLOW FROM
C     ANY CONTRIBUTING LATERAL PLANES.
C     NOTE THAT THE INDEX (N) IN THIS ROUTINE IS ON TIME.
      EXTERNAL IMPAUB
      NOQL=1
      NPA=0
      NCHN=0
        SUMA(J) = 0.
      NUT=NU(J)
      NC1T=NC1(J)
      NC2T=NC2(J)
C
C      ALSO RESET QUB, QL
      DO 10 N=1,NI
      CUB(N)=0.0
        QUB(N)=0.0
        QL(N) =0.0
   10 AUB(N)=0.0
C
C     CHECK FOR CONTRIB. ELEMENTS AT UPPER END AND SEND TO APPROP. LOOP
   20 IF (NUT.NE.0) GO TO 30
      IF (NC1T.NE.0) THEN
          NCHN=NCHN+1
          SUMA(J) = SUMA(J) + SUMA(NC1T)
        END IF
      IF (NC2T.NE.0) THEN
          NCHN=NCHN+1
          SUMA(J) = SUMA(J) + SUMA(NC2T)
        END IF
      IF (NCHN.EQ.0) GO TO 230
      IF (NCHN-1) 50,50,70
C
C     CONVERGING PLANE AT UPPER END OF CHANNEL
   30 IF (NC1T.NE.0.OR.NC2T.NE.0) CALL errpo ('  ADD ',14,J,0,'   X  ',
     & '   X  ')
      IF (NB(NUT).EQ.0) CALL errpo ('  ADD ',15,NUT,0,'   X  ','   X  ')
        SUMA(J) = SUMA(J) + SUMA(NUT)
      MU=NB(NUT)
      DO 40 N=1,NI
      MM=N-1
      QUB(N)=QS(MM+MU)
      CUB(N)=SCON(MM+MU)
   40 CONTINUE
      NB(NUT)=0
      NUT=0
      GO TO 90
C
C     ONE CHANNEL AT UPPER END
   50 M=NB(NC1T+NC2T)
      IF(NB(NC1T+NC2T).EQ.0) CALL errpo('   ADD  ',15,NC1T+NC2T,0,
     &'   X  ','   X  ')
C
      DO 60 N=1,NI
      MM=N-1
      QUB(N)=QS(MM+M)
      CUB(N)=SCON(MM+M)
   60 CONTINUE
C
      NB(NC1T+NC2T)=0
      NC1T=0
      NC2T=0
      NCHN=0
      GO TO 90
C
C     TWO CHANNELS AT UPPER END
   70 M1=NB(NC1T)
      M2=NB(NC2T)
      IF(NB(NC1T).EQ.0) CALL errpo('  ADD ',15,NC1T,0,'   X  ','   X  ')
      IF(NB(NC2T).EQ.0) CALL errpo('  ADD ',15,NC2T,0,'   X  ','   X  ')
      DO 80 N=1,NI
      MM=N-1
      MX=M1+MM
      MY=M2+MM
      QUB(N)=QS(MX)+QS(MY)
      IF (QUB(N).LE.0.) GO TO 80
      CUB(N)=(QS(MX)*SCON(MX)+QS(MY)*SCON(MY))/QUB(N)
   80 CONTINUE
      NB(NC1T)=0
      NB(NC2T)=0
      NC1T=0
      NC2T=0
      NCHN=0
C
   90 Continue
C
C     TRAPEZOIDAL CASE -- UPPER BOUND AREA (AUB)
      IF (XL(J).EQ.0) RETURN
      AUB(1)=0.
      IEND=50
      XST=(QUB(2)*B*CO1/ALPHA(1))**(1./POWER(1))
C
C  ************************************************* N=2,NI LOOP *******
      DO 130 N=2,NI
      INDEX=N
      IF (QUB(N).EQ.0.) Then
        AUB(N)=0.
      Else
        IF (B .le. .000001) Then
          AUB(N)=((QUB(N)/ALPHA(1))*(2.*CO2*CO2/CO1)**
     &    (.5*(POWER(1)-1.)))**(1./(.5*(POWER(1)+1.)))
        Else !GO TO 130
  110     IF (XST.EQ.0.) XST=0.001
          CALL ITER (AUB(N),FA,DERFA,IMPAUB,XST,0.0001,IEND,IER,10000.)
          IER=IER+1
        
          IF (IER.EQ.3) Then
            AUB(N)=0.
            IER=1
          Else IF (IER.EQ.4) Then
            AUB(N)=0.
            IER=1
          End If
          IF (AUB(N).LT.0.0001) AUB(N)=0.
          XST=AUB(N)
          If(IER .eq. 2) Then
            CALL errpo ('  ADD ',10,IEND,0,'   X  ','   X  ')
          Else If(IER .eq. 3) Then
            CALL errpo ('  ADD ',11,N,0,'   X  ','   X  ')     
          Else If(IER .eq. 4) Then
            CALL errpo ('  ADD ',13,N,0,'   X  ','   X  ')
          Else If(IER .eq. 5) Then
            CALL errpo ('  ADD ',17,J,1,'   X  ','   X  ')
          Else If(IER .gt. 5) Then
            CALL GOTOER ('IER')
          End If
C           
C          GO TO (130,170,180,190,200), IER
        End If
      End If
  130 CONTINUE
C        ********************* END LOOP
      GO TO 20
C
C     CHECK FOR NUMBER OF CONTRIB. PLANES ON SIDES AND SEND TO CORRECT
C     LOOP.  AFTER THIS IS DONE, RETURN.
  230 NRT=NR(J)
      NLT=NL(J)
      IF (NRT.NE.0) THEN
          NPA=NPA+1
          SUMA(J) = SUMA(J) + SUMA(NRT)
        END IF
      IF (NLT.NE.0) THEN
          NPA=NPA+1
          SUMA(J) = SUMA(J) + SUMA(NLT)
        END IF
      IF (NPA.EQ.0) GO TO 300
      IF (NPA-1) 240,240,260
C
C     ONE SIDE LATERAL INFLOW ONLY
  240 M=NB(NLT+NRT)
      IF(NB(NRT+NLT).EQ.0) CALL errpo(' ADD ',15,NRT+NLT,0,'  X  ',
     &'   X  ')
C
      DO 250 N=1,NI
      MM=N-1
      QL(N)=(QS(MM+M))/XL(J)
      CONL(N)=SCON(MM+M)
  250 CONTINUE
C
      NPA=0
      NB(NLT+NRT)=0
      NRT=0
      NLT=0
      GO TO 320
C
C     BOTH LEFT AND RIGHT LATERAL INFLOW
  260 MR=NB(NRT)
      ML=NB(NLT)
      IF (NB(NRT).EQ.0) CALL errpo ('  ADD ',15,NRT,0,'   X  ','   X  ')
      IF (NB(NLT).EQ.0) CALL errpo ('  ADD ',15,NLT,0,'   X  ','   X  ')
      DO 290 N=1,NI
      MM=N-1
      MX=MM+MR
      MY=MM+ML
      QL(N)=(QS(MX)+QS(MY))/XL(J)
      IF (QL(N)) 280,280,270
  270 CONL(N)=(SCON(MX)*QS(MX)+SCON(MY)*QS(MY))/QL(N)/XL(J)
      GO TO 290
  280 CONL(N)=0.
  290 CONTINUE
      NPA=0
      NB(NRT)=0
      NB(NLT)=0
      NRT=0
      NLT=0
      RETURN
C
C     NO LATERAL INFLOW
  300 DO 310 N=1,NI
      QL(N)=0.0
      CONL(N)=0.
  310 CONTINUE
      NOQL=2
  320 RETURN
C
  330 FORMAT (//,1X,120(1H*),/,5X,'PIPE NO. ',I2,' WITH DIAMETER=',F7.3,
     1 ' HAS EXCEEDED ITS FLOW CAPACITY OF QMAX=',1PE12.5,/,5X,'OVERFLOW
     2  OCCURRED AT TIME STEP ',I3,' WITH QUB=',1PE12.5,/1X,120(1H*),1X)
      END
C*******************************************************************
      SUBROUTINE BISECT(A,B,FCT,EP,XOP,FOP)
C
C  THIS SUBROUTINE WILL LOCATE THE MIN. OF AN EXTERNAL ONE-DIMEN-
C  SIONAL FUNCTION BY USING THE BISECTION METHOD.
C
C  CALL SEQUENCE: CALL BISECT(A,B,FCT,EP,XOP,FOP)
C
C  CALLED BY: SUBROUTINE ITER
C
C  VARIABLES         TYPE USAGE DEFINITION
C  ---------         ---- ----- ----------
C  A                 R*4    I   LEFT ENDPOINT OF STARTING INTERVAL
C  B                 R*4    I   RIGHT   "     "     "         "
C  FCT               C*     I   EXTERNAL FUNCTION SUB. NAME
C  EP                R*4    I   STOPPING CRITERIA EPSILON
C  XOP               R*4    O   OPTIMUM OF VARIABLE X
C  FOP               R*4    O   OPT. FUNCTION VALUE AT VARIABLE XOP
C  XL                R*4   INT  LEFT END INTERVAL X
C  XR                R*4   INT  RIGHT END INTERVAL X
C  XM                R*4   INT  MIDDLE OF INTERVAL
C  FXL               R*4   INT  FUNCTION VALUE AT XL
C  FXR               R*4   INT     "       "    " XR
C  FXM               R*4   INT     "       "    " XM
C  ITER              I*4    O   ITERATION COUNTER
C /IO/
C  IDIAGN            I*4   COM  UNIT NUMBER FOR DIAG. ERROR WRITE
C
C
C  EXTERNAL REFERENCES: CALL FCT(X,F,DERF)
C-------------------------------------------------------------------
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT, IRDATA
C
      REAL*4 A,B,EP,XOP,FOP,XL,XR,XM,FXL,FXR,FXM
      REAL*4 DERF
C
      ITER = 0
      XL = A
      XR = B
C     BISECT THE INTERVAL
      XM = XL + ( (XR - XL)/2. )
C     COMPUTE THE INITIAL FUNCTION VALUES
      CALL FCT(XL,FXL,DERF)
      CALL FCT(XR,FXR,DERF)
      CALL FCT(XM,FXM,DERF)
C
C     * LOOP UNTIL THE INTERVAL IS LESS THAN EPSILON (EP)
   10   CONTINUE
C       * IF ZERO CROSSING IS BETWEEN LEFT AND MIDDLE VALUE
          IF( (FXL*FXM) .LT. 0.) THEN
C            THEN UPDATE VALUES AND SET XR = XM
             XR = XM
             FXR = FXM
          ELSE
C         ELSE SET XL = XM AND FXL = FXM
             XL = XM
             FXL = FXM
          ENDIF
C       * ENDIF
C
C       BISECT AND COMPUTE NEW FUNCTION VALUE
        XM = XL + ( (XR - XL)/2. )
        CALL FCT(XM,FXM,DERF)
C
C       * IF STOPPING CRITERIA MET
          IF( (ABS(XR-XL)) .LT. EP) THEN
C            THEN EXIT THE LOOP
             GO TO 20
          ELSE
C         ELSE CONTINUE LOOPING
             ITER = ITER + 1
             IF(ITER .GT. 40) THEN
                WRITE(IDIAGN,900)
  900           FORMAT(' ','STOP - MORE THAN 40 ITER. IN',
     -                     ' BISECTION AS CALLED FROM ITER',/)
                STOP
             ELSE
                GO TO 10
             ENDIF
          ENDIF
C       * ENDIF
C     * ENDLOOP
C
   20 CONTINUE
C     COMPUTE OPTIMUN VALUES
      XOP = XM
      FOP = FXM
C
      RETURN
      END
C****************************************************************************
C
      SUBROUTINE CHAINF(F0,F1,DELT,FMIN,AL)
C
C  THIS SUBROUTINE CALCULATES ACCUMULATED INFILTRATION AT EACH CHANNEL
C  NODE.
C
C    F1 = ACCUM. INFILT. AT TIME T+DELT IN FT. (OUTPUT)
C    F0 = ACCUM.    "    AT TIEM T      IN FT. (INPUT)
C---------------------------------------------------------------------------
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT, IRDATA
C
      DATA TOL /0.0005/
C
      FX = F0 + FMIN*DELT
C
      DO 16 K=1,30
         ERF = FX + AL*EXP(-FX/AL) - FMIN*DELT - F0 - AL*EXP(-F0/AL)
         DERFDX = 1.0 - EXP(-FX/AL)
         CORR = ERF/DERFDX
         FX = FX - CORR
CD!      STOPPING CRITERIA OF WALRUS*5
         IF(ABS(CORR/FX) .LT. TOL) GO TO 20
C
   16 CONTINUE
C     NO CONVERGENCE AFTER 30 ITERATIONS
      WRITE(IDIAGN,600)
      WRITE(IWRITE,600)
  600 FORMAT(//,' ','NO CONVERGENCE AFTER 30 ITER. IN  SUB. CHAINF',/)
C
   20 CONTINUE
      F1 = FX
      RETURN
      END
C**************************************************************************
C
      SUBROUTINE CHANNL (J)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT, IRDATA
      COMMON /CNTRL/ NRES, NOPT, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), QOB(100), TOB(100),
     1        NO, WTRAIN(60), QIDD(20,100), TIDD(20,100), MAXND,
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(2000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1,CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
C      COMMON /LAWS/ ATURB, PTURB, ALAM, PLAM, HTRANS, QTRANS
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(2000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /INFIL/ AL(60), SI(60), SMAX(60), ROC(60), RECS(60),
     1        thr, thfc, tempw, stone(60)
      COMMON /EROSN/ V1(20), V2(20), CT(20), CMX(20), RDET(20), DLA(20),
     1      DrlNET(20), Dirnet(20), suspir, susprl, SPOW(20)
      COMMON /BALANC/ STORA(60),EINF(60),CINF2(20), STODEP(60)
C
      CHARACTER*7 TYPE
      DATA NTELE /60/
C
C       WRITE OUT HEADER
        WRITE(IDIAGN,700) J
      JJCHAN=J
      NK=MAX1((15.*XL(J)/CLEN),5.)
      IF (NK.LE.LENH1) GO TO 10
      WRITE (IDIAGN,440) NK
      NK=15
   10 DX=XL(J)/(FLOAT(NK)-1.)
      NI=IFIX(TFIN/DELT)+1
      IF (NI.GT.LENQS) GO TO 360
C       INITIALIZE VOLUME BALANCE VARIABLES
C         SPLIT INFLOW INTO UPSTREAM INFLOW AND LAT. INFLOW
C              - SPLIT OUTFLOW INTO INFIL. AND RUN OUT
        SUMIN = 0.0
        SUMLAT = 0.0
        SOUT = 0.0
        STOR = 0.0
        SINF = 0.0
        EINF(J) = 0.0
        DO 15 KK=1,20
           CINF2(KK) = 0.0
   15   CONTINUE
C
      IF (NEROS.GT.1) VS=VSETL(RHOS(J),D50(J),XNU)
      POWER(1) = 5./3.
      ALPHA(1) = SQRT(S(J))/R1(J)
      If(alpha(1) .le. 0. .or. power(1) .le. 0.) 
     &            stop 'zero alpha or power in CHANNL'
      SEDTOT=0.
C
   30 CO1=1./ZR(J)+1./ZL(J)
      CO2=(1.+1./(ZR(J)*ZR(J)))**0.5+(1.+1./(ZL(J)*ZL(J)))**0.5
      B=A(J)
      BWID=B*CO1
C
   50 CALL ADD (J)
      IF (XL(J)) 400,400,60
   60 IF (NEROS.LE.1) GO TO 90
      DO 70 K=1,NK
         DRLNET(K)=0.
   70 Continue 
C
C     TRAPEZOIDAL SHAPED CROSS-SECTION
   90 Q(1)=0.0
      T(1)=0.0
      CONC(1)=0.
      DO 100 K=1,NK
      V1(K)=0.
C      INIT. FC1 FOR ALL NODES
        IF(FMINEW(J) .GE. 1.0E-08) THEN
C          INIT. CHAN. CUMM. INFIL. ARRAY TO APPROX. 0.0012 INCHES (A SMALL
C          NUMBER) TO AVOID DIVIDE BY ZERO WHEN COMPUTING INFIL. RATES
           FC1(K) = 0.0001
           FC2(K) = 0.0001
        ELSE
C       ELSE SET IT = 0.0 SO IMPREV. CASE CAN BE DETECTED IN IMPCHA
           FC1(K) = 0.0
           FC2(K) = 0.0
        ENDIF
      A1(K)=0.0
      A2(K)=0.0
  100   CONTINUE
C                                             BEGIN TRAPEZ. CHANNEL LOOP
      DO 160 L=2,NI
      T(L)=T(L-1)+DELT
  110 A2(1)=AUB(L)
      DT=DELT
      CALL IMPLCT (NK,J)
      DO 120 K=1,NK
      DUM=(B*B+2.*A2(K)/CO1)**.5
      WPR=(DUM-B)*CO2+B*CO1
  120 V2(K)=ALPHA(1)*(A2(K)/WPR)**(POWER(1)-1.)
      Q(L)=V2(NK)*A2(NK)
      IF (NP(J).GE.2) THEN
         WRITE (IDIAGN,450) L,T(L),L,Q(L),L,QL(L)
        ENDIF
      IF (NEROS.LE.1) GO TO 130
      CT(1)=CUB(L)
      DY = DX
      CALL SEDCOM (J,QL(L),dy,A2,A1,CONL(L),NK,2)                       !SEDCOM
      CONC(L)=CT(NK)
C       SUMIN = CUMM. UPPER BOUND. INFLOW VOLUME
C       SUMLAT= CUMM. LATERAL        "      "
  130   SUMIN = SUMIN + ( (QUB(L) + QUB(L-1)) * DT )/2.0
        SUMLAT= SUMLAT + ( (QL(L) + QL(L-1)) * DT * XL(J) )/2.0
C       COMPUTE LOOP STOP FOR TRAP. RULE
        NKM2 = NK - 2
        SINF2 = 0.0
        STOR = 0.0
C       * IF ONE OR MORE MIDDLE NODES
          IF(NKM2 .GE. 1) THEN
C            THEN LOOP SUMMING THE STORAGE AND INFIL. TERMS
             DO 96 K=1,NKM2
                STOR = STOR + A2(K+1)
                SINF2 = SINF2 + CINF2(K+1)
   96        CONTINUE
          ENDIF
C       * ENDIF
C       ADD IN END NODE CONTRIBUTIONS
        STOR = STOR + ( A2(1) + A2(NK) )/2.0
        SINF2 = SINF2 + ( CINF2(1) + CINF2(NK) )/2.0
C       CONVERT TO VOLUME
        STOR = STOR * DX
        SINF = SINF2 * DT * DX
C       SAVE THE CUMM. INFIL. OVER TIME FOR THIS ELEMENT
        EINF(J) = EINF(J) + SINF
C       SOUT = CUMM. OUTFLOW VOLUME, INCREASED AT EACH TIME STEP
        SOUT = SOUT + (Q(L) + Q(L-1))*DT/2.0
C       COMPUTE PERCENT ERROR
C       * IF NONZERO INFLOW
          IF( (SUMIN+SUMLAT) .NE. 0.0) THEN
C            THEN COMPUTE THE ERROR
             ER = ((SUMIN+SUMLAT) - (SOUT+EINF(J)) - STOR)/
     -             (SUMIN+SUMLAT) * 100.0
          ELSE
             ER = 0.0
          ENDIF
C
      IF (NP(J).LT.2) GO TO 140
      WRITE (IDIAGN,460) (A2(K),K=1,NK)
        IF (FMINEW(J) .GE. 1.0E-08) THEN
          WRITE (IDIAGN,464) (CINF2(K),K=1,NK)
        ENDIF
      IF (NEROS.GE.2) THEN
          WRITE (IDIAGN,465) (CMX(K),K=1,NK)
          WRITE (IDIAGN,470) (CT(K),K=1,NK)
        WRITE (IDIAGN,480) (DRLNET(K),K=1,NK)
        END IF
        WRITE(IDIAGN,625) SUMIN,SUMLAT,EINF(J),SOUT,STOR,ER
  140   DO 150 K=1,NK
        V1(K)=V2(K)
        A1(K)=A2(K)
  150   CONTINUE
      IF (NEROS.LE.1) GO TO 160
      SEDTOT=SEDTOT+(Q(L)*CONC(L)+Q(L-1)*CONC(L-1))*.5*RHOS(J)*1000.*DT
  160 CONTINUE
C
  290 CONTINUE
      IF (NP(J).EQ.0) GO TO 340
        BWT=A(J)*(1./ZL(J)+1./ZR(J))
      WRITE (IDIAGN,520) XL(J),S(J)
      WRITE (IDIAGN,530) ZL(J),ZR(J),BWT
C  310 IF (NREST.GT.1) GO TO 320
      WRITE (IDIAGN,550) R1(J)
C      GO TO 330
C  320 WRITE (IDIAGN,560) R1(J)
  330   CONTINUE
        IF(FMINEW(J) .GT. 1.0E-08) THEN
C          CONVERT FMINEW TO MM/HR FROM M/S
           FM = FMINEW(J) * 3.6e6
           WRITE (IDIAGN,565) FM
           G = AL(J)/PORCH(J)
           WRITE (IDIAGN,567) G,PORCH(J)
        ELSE
           WRITE (IDIAGN,568)
        ENDIF
      WRITE (IDIAGN,570) NC1(J),NC2(J),NL(J),NR(J),NU(J)
  340 LASTNB=NB(1)
      DO 350 NE=2,NTELE
         IF (NB(NE).GT.LASTNB) LASTNB=NB(NE)
  350 CONTINUE
      MBT=LASTNB+NI
      IF (LASTNB.EQ.0) MBT=1
      IF (MBT.LE.LENQS) GO TO 370
  360 WRITE (IDIAGN,580)
  370 DO 380 L=1,NI
      MM=L-1
      QS(MM+MBT)=Q(L)
      SCON(MM+MBT)=CONC(L)
  380 CONTINUE
      IF (NEROS.LE.1) GO TO 390
      WRITE (IDIAGN,590) D50(J),RHOS(J),PORCH(J
     1 ),PAVE(J)
      WRITE (IDIAGN,600) (DRLNET(K),K=1,NK)
  390 CONTINUE
      NB(J)=MBT
C       SAVE THE STORAGE FOR GLOBAL WATER BALANCE
        STORA(J) = STOR
        TYPE = 'CHANNEL'
C          FORM INFLOW AND OUTFLOW TOTALS
        SUMINT = SUMIN + SUMLAT
        SOUTT  = SOUT + EINF(J)
        WRITE(IDIAGN,605) SUMINT,SOUTT,STOR,ER
        WRITE (IWRITE,615) J,TYPE,ER,SEDTOT
        GO TO 420
C
C     IF XL=0 WE MERELY OUTPUT THE ADDED UPPER BOUND DISCHARGE (QUB
C     WHICH WAS CALCULATED IN ADD.
  400 T(1)=0.0
      Q(1)=QUB(1)
      CONC(1)=CUB(1)
      DO 410 L=2,NI
      T(L)=T(L-1)+DELT
      Q(L)=QUB(L)
      CONC(L)=CUB(L)
  410   CONTINUE
C     CODE FOR PROPER STORAGE OF AN INTERMED. ADDER CHAN.
      LASTNB=NB(1)
      DO 800 NE=2,NTELE
         IF (NB(NE).GT.LASTNB) LASTNB=NB(NE)
  800 CONTINUE
      MBT=LASTNB+NI
      IF (LASTNB.EQ.0) MBT=1
      IF (MBT.LE.LENQS) GO TO 820
  810 WRITE (IWRITE,580)
  820 DO 830 L=1,NI
      MM=L-1
      QS(MM+MBT)=Q(L)
      SCON(MM+MBT)=CONC(L)
  830 CONTINUE
  420 CONTINUE
      RETURN
C  430 STOP 6666
C
  440 FORMAT (/,' ',10X,'NK CALCULATED TO BE',I3,', BUT SET AT 15')
  450 FORMAT (/,' ','T(',I3,')=',E11.4,1X,'  Q(',I3,')=',E11.4,
     1                     1X,'  QL(',I3,')=',E11.4)
  460 FORMAT (' ','A2(K=1,NK)=',8F8.4/(10F8.4))
  464 FORMAT (' ','CINF2(K=1,NK)=',6(E10.3,1X)/(7(E10.3,1X)))
  465   FORMAT (' ','CMX(K=1,NK)=',7F9.6/(8F9.6))
  470 FORMAT (' ','CT(K=1,NK)=',7F9.6/(8F9.6))
  480 FORMAT (' ','DANET(K=1,NK)=',7F9.5/(8F9.5))
  495 FORMAT (' ','A1(K=1,NK)',8F8.4/(10F8.4))
  500 FORMAT (' ','A2(K=1,NK)',8F8.4/(10F8.4))
  520 FORMAT (/,' ',5X,'GEOM. PARAMETERS ARE  L=',F7.1,'  S=',F7.4)
  530 FORMAT (' ',5X,'TRAP. X-S  LEFT SLOPE=',F6.3,' RIGHT SLOPE='
     1 ,F6.3,' BOTT. WID.=',F7.2)
  550 FORMAT (' ',5X,'ROUGHNESS COEF. IS MANNINGS N= ',F5.3)
C  560 FORMAT (' ',5X,'ROUGHNESS COEF. IS CHEZY    C= ',F6.1)
  565 FORMAT (' ',5X,'INFILT. PARAMETER IS FMINEW=',F8.5)
  567 FORMAT (' ',5X,'INFILT. PARA. G=',F7.4,' POR=',F7.4)
  568 FORMAT (' ',5X,'IMPERVIOUS CHANNEL')
  570 FORMAT(' ',5X,'CONTRIB. CHANNEL NUMBERS: NC1=',I4,' NC2=',I4,/,
     1' ',5X,'CONTRIB. PLANE NUMBERS:  LEFT=',I4,'  RIGHT=',I4,
     2       '  UPPER=',I4)
  580 FORMAT (' QS NEEDS TO BE DIMENSIONED LARGER')
  590 FORMAT (' ',5X,'EROSION PARAMETERS ARE ---',/,
     3    ' ',6X,' D50=',G9.3,' RHOS=',F6.2,' SURF.POR.=',F6.2,
     4    ' PAVE. FAC.=',F6.3)
  600 FORMAT (' ',5X,'ACCUMUL. CHAN. BOTTOM DEPOSIT. OR EROSION
     1(NEG.)AT EACH NODE ',
     2                  '(SQ.FT.)',/,' ',(8F10.6))
  605 FORMAT (' ',/,16X,'*** WATER BALANCE AT END OF CHANNEL ***',/,
     1          ' ',14X,'<INFLOW BASED ON RUNIN + LATERAL INFLOW>',/,
     2  ' ',2X,'INFLOW=',E10.3,' OUTFLOW=',E10.3,' STOR.=',E10.3,
     3      ' ERROR=',E10.3,' %')
  615   FORMAT (' ',9X,I3,6X,A7,5X,E10.3,5X,F11.3)
  625   FORMAT(' ','  UPPER INFLOW= ',E10.3,'  LAT. INFLOW=',E10.3,
     1             '  INFIL. OUTFLOW=',E10.3,/,
     2         ' ','  CHAN. OUTFLOW=',E10.3,'      STORAGE=',E10.3,
     3             '           ERROR=',E10.3 ,' %')
  700 FORMAT(' ',8X,'** CHANNEL NO ',I4,' DIAGNOSTIC INFORMATION ** ',/)
      END
C****************************************************************************
      FUNCTION FCP(QC,R,APAR,FKO)
C  FINDS FLUX CAPACITY CORRESP. TO GIVEN PROFILE WATER
C
      IF (R.LE.FKO) THEN
        WRITE(6,100)
  100     FORMAT(' CALL FCP FOR R<KS. ERROR ')
          STOP 66
      ELSE
        RP=1.001*R
          IF(QC .GT. 0.0) THEN
C       CHECK FOR OVERFLOW FOR SATURATION OR NEAR SATURATION
             TEST = QC/APAR
             IF(TEST .GE. 80.0) THEN
                FCM = FKO
             ELSE
              EF=EXP(QC/APAR)
              FCM=FKO*EF/(EF-1.)
             ENDIF
          ELSE
             FCM=10000.
          ENDIF
        FCP=AMIN1(FCM,RP)
      ENDIF
      RETURN
      END
C*************************************************************************
C
      SUBROUTINE FINTG (DT,QSL,RF,FCL,FMN,APAR,QSN)
C
C INTEGRATES ALONG TIME SCALE UNTIL CHANGE IN SOIL STORAGE IS FOUND.
C *  REVISED VERSION 8/88 BY RES
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT, IRDATA
C
      DQR=RF*DT
      IF (RF.LT.FCL) THEN
C * ADJUST TO INTEGRATE FROM Qp:
        QP = QCP(RF,APAR,FMN)
        CHEQ = QSL + DQR
        IF(CHEQ .LE. QP) THEN
C*  NO RUNOFF:
          QSN = CHEQ
          GO TO 10
        ELSE
C * RUNOFF AFTER QP:
          QSO = QP
          DTU = DT - (QP-QSL)/RF
          dqr = rf*dtu
        END IF
      ELSE
C  RUNOFF FROM RAIN CONSISTENTLY GREATER THAN FCL:
        DTU = DT
        QSO = QSL
      END IF
      QST = QSO + DQR
      FCT = FCP(QST,RF,APAR,FMN)
      If(FCT/fmn .ge. 50.) then  ! use imbibition approx:
        tz = 0.5*QSO*QSO/APAR/FMN
        QSN = SQRT((tz+dtu)*2.*APAR*FMN)
      Else
       ITR = 0
        QBM = FMN*DTU
        QF2 = (FCL + FCT)*DTU + QSO
        QStst = qf2                                                     !10/91
        fnz = FOF(QSO,APAR)
    1   ITR = ITR + 1
C*  NEWTON RAPHSON ITERATION TO SOLVE FOR QF2:
        fnc = FOF(QF2,APAR)
        OFN = fnc - fnz - QBM
        DOFN = 1. - EXP(-QF2/APAR)
        ERFN = OFN/DOFN
        QFT = QF2 - ERFN
        If(QFT .lt. QSO) Then
          QF2 = 0.2*QF2 + 0.8*QSO
        Else
          QF2 = QFT
        End If
        If(itr .gt. 10) 
     &  WRITE(IDIAGN,191) QSo,qf2,fnc,FCL,DOFN,OFN,ERFN,dtu
        IF(ITR .GT. 20) GO TO 92
        IF(ABS(ERFN) .GT. 0.0005*QStst) GO TO 1                         !10/91
        QSN = QF2
C
      End If
   10 CONTINUE
      RETURN
C
   92   WRITE(IDIAGN,192) QSo,qf2,fnc,FCL,DOFN,OFN,ERFN,dtu
     &                   ,rf,dt,qf2,fmn,apar
  191 Format(7g13.6)
  192  FORMAT(' 20 ITRS IN FINTG:QStst,RF,FCL,DOFN,OFN,ERFN:'/,(7G13.6))
      STOP 'Iteration limit in FINTG'
      END
C-------------------------------------
       FUNCTION FOF(Q,AP)
       FOF = Q + AP*EXP(-Q/AP)
       RETURN
       END
C *****************************************
      FUNCTION FSRT(RS,SMX,SR)
C* GET STAEDY SATURATION AS FUNCTION OF FLOW, F .LT. KS:
C         ASSUME EPSL = 10.
      FSRT = SR + (SMX - SR)*RS**0.1
      RETURN
      END
C*************************************************************************
C*************************************************************************
C
      SUBROUTINE IMPAUB (X,FX,DERF)
C
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1,CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(2000), LENQ, LENQS, LENH1, L, QL2(20), JU
C
      dqdact(AL,A,D1,D2,P,PM1,C)=AL*(A/D2)**PM1*(P-(PM1*A*C/(D1*D2)))
C
      N=0
   10 KFLAG=0
      J=INDEX
      IF (X.EQ.0..AND.QUB(J).EQ.0.) QUB(J)=0.0000001
      P=POWER(1)
      PM1=P-1.
        C=CO2/CO1
      DUM1=(B*B+2.*X/CO1)**0.5
      DUM2=(DUM1-B)*CO2+B*CO1
      FX=ALPHA(1)*X**P/DUM2**PM1-QUB(J)
C
C     CALCULATE DERIVATIVE OF ERROR FUNCTION (DERF)
      DERF=dqdact(ALPHA(1),X,DUM1,DUM2,P,PM1,C)
C     DUE TO MISBEHAVIOR OF THE ERROR FUNCTION IN SOME CASES, THE FOL -
C     LOWING CORRECTION OF X MAY BE NECESSARY FOR CONVERGENCE
      IF (QUB(J-1).EQ.0.) GO TO 20
      IF (QUB(J)/5..GT.QUB(J-1)) KFLAG=1
      IF (KFLAG.EQ.0) GO TO 20
      IF (DERF.LT.0..AND.FX.LT.0.) GO TO 30
   20 RETURN
   30 X=X*10.
      N=N+1
      IF (N.GT.2) CALL errpo ('IMPAUB',10,10,0,'   X  ','   X  ')
      GO TO 10
      END
C*************************************************************************
C
      SUBROUTINE IMPCHA (X,FX,DERF)
C
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1,CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), QOB(100), TOB(100),
     1        NO, WTRAIN(60), QIDD(20,100), TIDD(20,100), MAXND,
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(2000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(2000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /INFIL/ AL(60), SI(60), SMAX(60), ROC(60), RECS(60),
     1        thr, thfc, tempw, stone(60)
      COMMON /BALANC/ STORA(60),EINF(60),CINF2(20), STODEP(60)
C
      DIMENSION DQDA(2)
C
      dqdact(AP,A,D1,D2,P,PM1,C)=AP*(A/D2)**PM1*(P-(PM1*A*C/(D1*D2)))
      ZFUNC(B,AREA,CO1)=(B*B+2.*AREA/CO1)**0.5
C       NEW FUNCTIONS FOR CHANNEL INFLITRATION
        PRFUNC(A,B,BW,CO1)=AMIN1( (-B+SQRT(B*B+2.*A/CO1))/
     -  (0.15*BW**0.5), 1.0 )
        WPFUNC(B,A,CO1,CO2)=(-B+SQRT(B*B+2.*A/CO1))*CO2 + B*CO1
C
C       NEW FUNCTIONS FOR CHANNEL INFLITRATION
        BW = BWID
C
      IF (B.EQ.0.) B=.000001
      J=INDEX
        JP1 = J+1
      C=CO2/CO1
      XK=B*CO1
      P=POWER(1)
      PM1=P-1.
      PM2=P-2.
      THETA1=1.-THETA
      AREA1=(A1(J)+A1(J+1))/2.
      DUM11=ZFUNC(B,AREA1,CO1)
      TW(J)=DUM11*CO1
      DUM21=(DUM11-B)*CO2+XK
      DQDA(1)=dqdact(ALPHA(1),AREA1,DUM11,DUM21,P,PM1,C)
      AREA2=(X+A2(J))/2.
      DUM12=ZFUNC(B,AREA2,CO1)
C  TW IS THE TOP WIDTH
      TW(J+1)=DUM12*CO1
C  DUM22 IS THE WETTED PERIMETER
      DUM22=(DUM12-B)*CO2+XK
C      LOOP AROUND THIS IF NO INFILTRATION - NOTE: DETECT
C              IMPERV. CASE USING FC1 AS SET IN SUBR. CHANNL FOR IMPREV.
        IF(FC1(J) .GT. 0.0) THEN
c        GAL=AL(JJCHAN)*(SMAX(JJCHAN)-SI(JJCHAN))/12.
C          COMPUTE INFILT EST. AT JP1
C          CHECK FOR NEG. AREA
           IF(X .GT. 0.0) THEN
              xfljp1 = (FC2(JP1) - FC1(JP1))/DT * (WPFUNC
     -        (B,X,CO1,CO2) * PRFUNC(X,B,BW,CO1) + WPFUNC
     -        (B,A1(JP1),CO1,CO2) * PRFUNC(A1(JP1),B,BW,CO1))
     -        * 0.5
           ELSE
              xfljp1 = (FC2(JP1) - FC1(JP1))/DT * (WPFUNC
     -        (B,A1(JP1),CO1,CO2) * PRFUNC(A1(JP1),B,BW,CO1))
     -        * 0.5
           ENDIF
        ENDIF
CD!
      DQDA(2)=dqdact(ALPHA(1),AREA2,DUM12,DUM22,P,PM1,C)
C
C     CALCULATE SECOND DERIVATIVE OF DISCHARGE EQN WITH RESPECT TO X
C     (D2Q/DADX)
      IF (AREA2.EQ.0.) GO TO 10
      FAC=(AREA2/DUM22)
      TERM1=FAC**PM2*PM1*(1.-X*C/(DUM12*DUM22))/(2.*DUM22)
      TERM2=FAC**PM1*C*PM1*(1.-AREA2*(C+DUM22/(CO1*DUM12))/(DUM12*DUM22)
     1 )/(2.*DUM12*DUM22)
      D2QDAX=ALPHA(1)*(TERM1-TERM2)
      GO TO 20
   10 D2QDAX=0.
   20 CONTINUE
        IF(FC1(J) .GT. 0.0) THEN
           AOUT = (XINFLJ + xfljp1) * 0.5
        ELSE
           AOUT = 0.0
        ENDIF
C
      FX=(AREA2-AREA1)/DT+(THETA/DX)*DQDA(2)*(X-A2(J))+(THETA1/DX)*DQDA(
     1 1)*(A1(J+1)-A1(J))-0.5*(QL(L-1)+QL(L))+AOUT
      DERF=1./(2.*DT)+(THETA/DX)*(DQDA(2)+D2QDAX*(X-A2(J)))
      RETURN
      END


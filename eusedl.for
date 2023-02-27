      SUBROUTINE POND (J)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1,CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), rfr(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /PONDS/ AP(10,10,3), TWP(10,10,3), EL(10,3), QSO(20,3),
     1   ELQ(20,3),ELZQ(3),ELVZ(3),NOPQ(3),MEL(3),NSEC(3),XSEC(10,3),
     2   ELZX(10), JOP(3), ELST(3), ALV, DIFUS(3), SIGMAS(60), NPART
      COMMON /CNTRL/ NRES, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /BALANC/ STORA(60),EINF(60),CINF2(20), STODEP(60)
      Common /pldat/ ifad,ndat,tdat(500),qdat(500),sdat(500),sumrmm,
     1                sedtot
C
      DIMENSION VOL(10,3),DELX(10,3),SAR(10,3),STOR(1000),SUMIN(1000),
     1 SUMOUT(1000), DVS(7), CPTS(7), VSTL(7)
C
        CHARACTER*5 TYPE
      Logical ifad
C
C     TOLERANCE CHECK FOR EQUALITY
      DATA TOLEL,TOLEQ,TOLB/1.0E-04,1.0E-06,1.0E-04/
      DATA IBMAX/40/
C     UNITS CONVERSION CONSTANTS
      DATA NTELE /60/
C
C ELST = ELEV OF W.S. AT BEGINNING OF SIMULATION
C ELVZ = ELEV. OF POND BOTTOM
C ELZQ = ELEV. BELOW WHICH OUTFLOW DOES NOT OCCUR
C NPN(J) IS POND NO. WHICH IS ELEMENT J
C JOP(N) IS ELEMENT NO. CORRESPONDING TO POND N
C NOPQ IS NO. OF POINTS ON QOUT RATING TABLE
C QSO(L) IS OUTFLOW RATING, L=1,NOPQ
C ELQ(L) IS OUTFLOW RATING DEPTH, L=1,NOPQ
C VOL(J) IS VOL. OF POND AT ELEVATION J, J = 1 TO MEL
C MEL IS NO OF ELEVATIONS AT WHICH X-SEC AREAS ARE GIVEN, AND VOL CALCUL
C MS IS INDEX OF FORST X-SEC AREA IN ARRAY K
C AP(K,J) IS X-SEC AREA AT STATION K AND ELEV. J
C EL(J) IS ELEVATION J, SAME FOR ALL STATIONS
C TWP(K,J) IS TOP WIDTH AT STA. K, ELEV. J
C NSEC IS NO. OF STATIONS, K=1,NSEC
C
C
      NI = IFIX(TFIN/DELT) + 1
      CALL ADD(J)
C
      NPND = NPN(J)
      ME = MEL(NPND)
      NS = NSEC(NPND)
      MQ = NOPQ(NPND)
      CALL VOLCAL(NPND,VOL,SAR,DELX)
C  WRITE OUT THE RATING TABLE
C     CHECK FOR PROPER UNITS
         WRITE(IDIAGN,350)
  350    FORMAT(/,' ',7X,' POND OUTFLOW RATING TABLE ',/,
     -            ' ',7X,'   ELEV. (M)      Q (CMS)    ',/,
     -            ' ',7X,'   ----------     -------    ')
      DO 370 K=1,MQ
         QR = QSO(K,NPND)
         ELR= ELQ(K,NPND)
         WRITE(IDIAGN,360) ELR,QR
  360    FORMAT(' ',9X,F10.3,3X,F10.3)
  370 CONTINUE
C
C  FIND INITIAL POND VOLUME
      CALL VOFH(NPND,ME,VOL1,VOL,ELST(NPND),EL,1,ELVZ(NPND),IS)
      VSTR = VOL1
C  INITIALIZE OUTFLOW FOR THE STARTING ELEVATION
      IF(ELST(NPND) .LE. ELZQ(NPND)) THEN
         Q(1) = 0.
      ELSE
         CALL VOFH(NPND,MQ,QINT,QSO,ELST(NPND),ELQ,1,ELZQ(NPND),IS)
         Q(1) = QINT
      ENDIF
C  CALCULATE SEDIMENT PARTICLE CLASSES
      IF(NEROS.LE.1) GO TO 50
      JU = MAX0(NU(J),NC1(J),NC2(J))
      ICS = 2
      LR = 1
      CALL SEDIV(D50(JU),SIGMAS(JU),NPART,DVS)
      DO 4 K=1,NPART
      CPTS(K) = 0.
    4 VSTL(K) = VSETL(RHOS(J),DVS(K),XNU)
      PART = FLOAT(NPART)
      IF(NP(J).LE.1) GO TO 50
      WRITE(IDIAGN,211) NPART,(DVS(I),I=1,NPART)
  211 FORMAT(' SEDIMENT SIMULATED IN',I2,' PARTS, SIZES (m) = ',
     17F8.6)
      WRITE(IDIAGN,215) (VSTL(K),K=1,NPART)
  215 FORMAT(' SETTLING VEL. IN M/S. =',7(1X,G11.5))
   50 CONTINUE
      SUMOUT(1) = 0.
      SUMIN(1) = 0.
      STOR(1)  = VOL1
      ZE1 = ELST(NPND)
      IT = 1
      CONL(1) = DT
      SEDTOT = 0.
C       VNOTCH = VOLUME AT NOTCH ELEVATION
C       DTT    = TIME (0 < DTT < DT) AT WHICH POND DEPTH = NOTCH ELEV.
C       FLOWIN = VOLUME OF INFLOW DURING DT
C       floout=    "   "  OUTFLOW  "    DT
C       flonew= VOLUME OF FLOW ABOVE WEIR NOTCH
C       VOL1   = VOLUME OF POND AT TIME T (LM IN DO LOOP)
C       VOL2   = VOLUME OF POND AT TIME T + DT (L IN DO LOOP)
C       VOLL   = VOLUME OF POND AT ZE2L
C       VOLU   = VOLUME OF POND AT ZE2U
C       ZE1    = POND ELEVATION AT TIME T (LM IN DO LOOP)
C       ZE2L   = LOWER BOUND ON POND ELEV. AT TIME T + DT (L IN DO LOOP)
C       ZE2U   = UPPER BOUND ON POND ELEV. AT TIME T + DT (L IN DO LOOP)
C       ZE2    = POND ELEVATION AT TIME T + DT (L IN DO LOOP)
C       ZT     = PERTURBED ZE2 DEPTH
C       ZMID   = MIDPOINT OF Z IN BISECTION
C       EF     = ERROR FUNCTION
C       EFL    = LOWER INTITAL ERROR FUNCTION
C       EFU    = UPPER INITIAL ERROR FUNCTION
C       DZ     = CHANGE TO ZE2 IN BISECTION SCHEME
C**C       AREA1  = POND AREA AT ZE1
C**C       AREA2  = POND AREA AT ZE2
C       storno= POND STORAGE BETWEEN NOTCH AND CURRENT POND ELEV. (ZE1)
C       DIFFQ  = (FLOWIN - 0.0)
C       DIFFEL = (ZE1 MINUS THE NOTCH ELEVATION)
C
      CALL VOFH(NPND,ME,VNOTCH,VOL,ELZQ(NPND),EL,1,ELVZ(NPND),IS)
C**C     CALCULATE THE AREA AT THE START ELEVATION ELST=ZE1
C**      CALL VOFH(NPND,ME,AREA1,SAR,ZE1,EL,1,ELVZ(NPND),IS)
C     INITIALIZE DTT
      DTT = DT
C
C MAIN LOOP THRU TIME INCREMENTS*********************
      DO 100 L=2,NI
      LM = L - 1
      SUMIN(L) = SUMIN(LM) + 0.5*DT*(QUB(L)+QUB(LM))
      FLOWIN = 0.5 * DT * (QUB(L) + QUB(LM))
C     FORM DIFFERENCES FOR REAL VAR. EQUALITY COMPARISONS
      DIFFEL = ABS(ZE1 - ELZQ(NPND))
      DIFFQ  = ABS(FLOWIN - 0.0)
C
C     FIVE CASES FOR A GENERAL SOLUTION MUST BE EXAMINED
C
      IF(ZE1 .LT. ELZQ(NPND) .AND. FLOWIN .GT. 0.0) THEN
C        THEN EXAMINE CASE 1 OR 2: WATER LEVEL BELOW NOTCH ELEV. (ELZQ)
C        THEN CHECK TO SEE IF INFLOW WILL CAUSE OUTFLOW
         storno = VNOTCH - VOL1
         IF(FLOWIN .LE. storno) THEN
C CASE 1:   THEN CASE 1: INFLOW WILL NOT CAUSE POND TO FILL; NO OUTFLOW
            VOL2 = VOL1 + FLOWIN
            CALL HOFV(NPND,ME,ZE2,EL,VOL2,VOL,1,ELVZ(NPND),IS)
            Q(L) = 0.0
C           UPDATE ELEV. AND VOLUME
            ZE1 = ZE2
            VOL1 = VOL2
            GO TO 90
         ELSE
C CASE 2:ELSE CASE 2: START ELEV < NOTCH ELEV. AND INFLOW WILL CAUSE
C                     OUTFLOW IN THIS CASE THE TIME INCREMENT MUST BE
C                     PARTIIONED SOLVE FOR DTT USING QUADRATIC EQUATION
            BB = QUB(LM)
            AA = 0.5 * (QUB(L) - QUB(LM))/DT
            CC = -storno
            DTP = (-BB + SQRT(BB*BB - 4.0 * AA * CC)) * 0.5/AA
            flonew = FLOWIN - storno
C           COMPUTE UPPER BOUND ON ELEV. (ZE2U) ASSUMING NO OUTFLOW
            VOL2 = VOL1 + FLOWIN
            VOLU = VOL2
C           CHECK UPPER BOUND ON VOLUME
            IF(VOLU .GT. VOL(MEL(NPND),NPND) ) THEN
               VOLU = VOL(MEL(NPND),NPND)
            ENDIF
            CALL HOFV(NPND,ME,ZE2U,EL,VOLU,VOL,1,ELVZ(NPND),IS)
C           SINCE WE ARE JUST CROSSING THE NOTCH LOWER BOUND (ZE2L)
C             IS THE NOTCH
            ZE2L = ELZQ(NPND)
            VOL1 = VNOTCH
            VOLL = VOL1
            FLOWIN = flonew
            DTT = DT - DTP
            GO TO 5
         ENDIF
      ELSE IF(ZE1 .LT. ELZQ(NPND) .AND. FLOWIN .LE. 0.0) THEN
C CASE 3:START ELEV. < NOTCH ELEV. AND NO INFLOW => NO OUTFLOW
         VOL2 = VOL1
         ZE2  = ZE1
         Q(L) = 0.0
         GO TO 90
      ELSE IF(DIFFEL .LE. TOLEL .AND. DIFFQ .LE. TOLEQ) THEN
C CASE 4:START ELEV. = NOTCH AND NO INFLOW => NO OUTFLOW
         VOL2 = VOL1
         ZE2  = ZE1
         Q(L) = 0.0
         GO TO 90
      ELSE IF(ZE1 .GE. ELZQ(NPND) .AND. FLOWIN .GE. 0.0) THEN
C CASE 5:START ELEV. .GE. NOTCH ELEV. WITH ZERO OR POSITIVE INFLOW
C        THEN SET UP FOR THE ITERATIVE SOLUTION
C        COMPUTE LOWER BOUND ASSUMING CONSTANT OUTFLOW OF Q(LM)
C          AND NO INFLOW
         floout = Q(LM) * DTT
         VOLL = VOL1 - floout
         CALL HOFV(NPND,ME,ZE2L,EL,VOLL,VOL,1,ELVZ(NPND),IS)
C        IF NEW ZE2L IS BELOW THE NOTCH
         IF(ZE2L .LE. ELZQ(NPND)) THEN
            ZE2L = ELZQ(NPND)
            VOLL = VNOTCH
         ENDIF
C        COMPUTE THE UPPER BOUND ASSUMING NO OUTFLOW
         VOLU= VOL1 + FLOWIN
C        CHECK UPPER BOUND ON VOLUME
         IF(VOLU .GT. VOL(MEL(NPND),NPND) ) THEN
C           THEN RESET TO MAXIMUM POND VOLUME - DELTA
            VOLU = (VOL(MEL(NPND),NPND) - 0.01)
         ENDIF
         CALL HOFV(NPND,ME,ZE2U,EL,VOLU,VOL,1,ELVZ(NPND),IS)
C        WITH ZE2L,ZE2U DEFINED GO TO BISECTION (LABEL 5)
      ENDIF
C
C     WITH THE ZE2 BRACKETED BY ZE2L AND ZE2U USE BISECTION TO
C     COMPUTE FINAL ZE2 WITHIN TOLERANCE (TOLB)
C
    5 CONTINUE
C     FOR THE ERROR FUNCTION EF, Q2 AND VOL2 WILL VARY WITH
C       CHANGING STAGE ZE2
C     COMPUTE THE INITIAL ERROR FUNCTION VALUES FOR THE LOWER AND
C       UPPER BRACKET ON ZE2
C     FIRST GET THE ASSOCIATED DISCHARGES
      CALL VOFH(NPND,MQ,Q2L,QSO,ZE2L,ELQ,IT,ELZQ(NPND),IQ)
      CALL VOFH(NPND,MQ,Q2U,QSO,ZE2U,ELQ,IT,ELZQ(NPND),IQ)
      EFL = FLOWIN - 0.5 * (Q(LM) + Q2L) * DTT - (VOLL - VOL1)
      EFU = FLOWIN - 0.5 * (Q(LM) + Q2U) * DTT - (VOLU - VOL1)
C     CHECK TO SEE THAT ZERO IS BRACKETED
      IF( (EFL*EFU) .GE. 0.0 ) THEN
         ZE2L = ELZQ(NPND)
         ZE2U = EL(MEL(NPND),NPND) - 0.001
         VOLL = VNOTCH
         VOLU = VOL(MEL(NPND),NPND) - 0.001
         CALL VOFH(NPND,MQ,Q2L,QSO,ZE2L,ELQ,IT,ELZQ(NPND),IQ)
         CALL VOFH(NPND,MQ,Q2U,QSO,ZE2U,ELQ,IT,ELZQ(NPND),IQ)
         EFL = FLOWIN - 0.5 * (Q(LM) + Q2L) * DTT - (VOLL - VOL1)
         EFU = FLOWIN - 0.5 * (Q(LM) + Q2U) * DTT - (VOLU - VOL1)
      ENDIF
      IF(EFL .LT. 0.0) THEN
         ZE2 = ZE2L
         DZ  = ZE2U - ZE2L
      ELSE
         ZE2 = ZE2U
         DZ  = ZE2L - ZE2U
      ENDIF
C     BEGIN BISECTION LOOP
      IB = 1
   7  CONTINUE
      DZ = DZ * 0.5
      ZMID = ZE2 + DZ
C     COMPUTE ASSOC. ERROR FUNCT.
      CALL VOFH(NPND,ME,VOL2,VOL,ZMID,EL,IS,ELVZ(NPND),I)
      IS = I
C     GET DISCHARGE AS A FUNCTION OF ZE2
      CALL VOFH(NPND,MQ,Q2,QSO,ZMID,ELQ,IT,ELZQ(NPND),IQ)
      IT = IQ
      EF = FLOWIN - 0.5 * (Q(LM) + Q2) * DTT - (VOL2 - VOL1)
      IF(EF .LE. 0.0) THEN
         ZE2 = ZMID
      ENDIF
      IF(ABS(DZ) .LT. TOLB .OR. EF .EQ. 0.0) THEN
C        THEN EXIT THE DO LOOP
         GO TO 10
      ELSE
C        ITERATE
         IB = IB + 1
         IF(IB .LT. IBMAX) THEN
            GO TO 7
         ELSE
            WRITE(IDIAGN,310)
            WRITE(*,310)
  310       FORMAT(//,' ',10X,' TOO MANY BISECTIONS IN SUBR. POND',//)
            STOP
         ENDIF
      ENDIF
C
   10 CONTINUE
C     UPDATE AND RESETS
      ZE1 = ZE2
      VOL1 = VOL2
      Q(L) = Q2
      STOR(L) = VOL2
      DTT = DT
C
   90 IF(NP(J).NE.7) GO TO 20
      IF(L .EQ. 2) THEN
C        THEN WRITE OUT A HEADER
         WRITE(IDIAGN,320)
  320    FORMAT(//,' ',10X,'DETAILED PRINTOUT OF POND COMPUTATIONS',/,
     -             ' ',10X,'--------------------------------------',/)
      ENDIF
      WRITE(IDIAGN,201) T(L),QUB(L),Q(L),ZE2,VOL2
  201 FORMAT(' ','TIME =',F9.2,3X,'Q - IN =',F8.3,3X,
     -'Q - OUT =',F8.3,/,
     -' ',5X,'STAGE =',F10.2,3X,'POND VOLUME =',F10.3)
   20 CONTINUE
      QOUT = 0.5*DT*(Q(L)+Q(LM))
      STOR(L) = VOL2
      SUMOUT(L) = SUMOUT(LM) + QOUT
      IF(NEROS.LE.1) GO TO 100
      IF(SUMOUT(L).GT.VSTR) GO TO 21
C ** OUTFLOW LESS THAN ORIG. STORAGE VOLUME
      CONC(L) = 0.
      TELAP = T(L)
      GO TO 60
   21 NETOUT = SUMOUT(L)-VSTR
   22 IF(SUMIN(LR).GE.NETOUT) GO TO 23
      LR = LR + 1
      GO TO 22
   23 IF(SUMIN(LR-1).LT.NETOUT) GO TO 24
      LR = LR-1
      GO TO 23
   24 LRM = LR-1
      FRACT = (NETOUT-SUMIN(LRM))/(SUMIN(LR)-SUMIN(LRM))
      TIN = T(LRM) + DT*FRACT
C* TELAP IS FLOW-THRU TIME FOR WATER LEASVING POND AT TIME L:
      TELAP = T(L)-TIN
      CONCIN = CUB(LRM) + (CUB(LR)-CUB(LRM))*FRACT
      STAV = (STOR(LR)+STOR(L))*.5
      CALL HOFV(NPND,ME,HAV,EL,STAV,VOL,ICS,ELVZ(NPND),IS)
      ICS = IS
      CALL VOFH(NPND,ME,SAV,SAR,HAV,EL,ICS,ELVZ(NPND),IS)
      DEPAV = STAV/SAV
      CONC(L) = CONCIN
      DO 30 K=1,NPART
      DSETL = VSTL(K)*TELAP
      RAT = DSETL/DEPAV
      F1 = 1.05 + RAT
      PCTS = 0.5*(F1 - SQRT(F1*F1 - 4.*RAT))
      CONC(L) = CONC(L) - CONCIN/PART*PCTS
  30  CPTS(K) = CONCIN*(1.-PCTS)/PART
      SEDTOT =SEDTOT + DT*500.*RHOS(J)*(Q(L)*CONC(L)+Q(LM)*CONC(LM))
      LR = LRM
   60 IF(NP(J).NE.7) GO TO 35
C* CONL(L) IS BORROWED AS AN ARRAY TO STORE VALUES OF TELAP
      CONL(L) = TELAP
      WRITE(IDIAGN,212) CUB(L),CONC(L)
  212 FORMAT(' ',5X,'INPUT CONC. =',F8.6,1X,' TOTAL OUT CONC. =',F8.6)
      WRITE(IDIAGN,213) (CPTS(K),K=1,NPART)
  213 FORMAT(' ',5X,'CONC. OF PART. CLASSES =',7(F8.6,1X))
   35 CONTINUE
  100 CONTINUE
      IF(NEROS.LE.1) GO TO 110
C ******** APPROXIMATE DIFFUSION IN POND
      NM = NI - 1
      DO 70 L=2,NM
        TDI = 2.*SQRT(DIFUS(NPND)*CONL(L))*CONL(L)/XSEC(NS,NPND)
        DLT = DT*10.
        TDI = AMIN1(TDI,DLT)
        IF(TDI .ge. DT) Then
          NAV = TDI*2./DT
          NAV = (NAV/2)*2 + 1
          NHAF = NAV/2
          IF(L .gt. NHAF) Then
            IF((L+NHAF).GT.NI) Then
              LK = NI
              LS = L - NI + L
              LT = NI - LS + 1
            Else
              LS = L - NHAF
              LK = L + NHAF
              LT = NAV
            End If
          Else
            LK = 2*L - 1
            LS = 1
            LT = LK
          End If
          STOR(L) = 0.
          DO K = LS,LK
            STOR(L) = STOR(L) + CONC(K)
          End Do
          STOR(L) = STOR(L)/FLOAT(LT)
        Else
          STOR(L) = CONC(L)
        End If
   70 CONTINUE
      DO 79 L=2,NM
   79 CONC(L) = STOR(L)
      GO TO 110
C
C     VOLUME BALANCE COMP.
C
  110 CONTINUE
      STPOND = VOL2 - VSTR
C     SAVE STORAGE FOR GLOBAL VOL BAL.
      STORA(J) = STPOND
      IF(SUMIN(NI) .LT. 1.0E-06) THEN
         ERP = 0.0
      ELSE
         ERP = (SUMIN(NI) - SUMOUT(NI) - STPOND)/SUMIN(NI) * 100.0
      ENDIF
      TYPE = 'POND'
      WRITE(IWRITE,210) J,TYPE,ERP,SEDTOT
  210 FORMAT(' ',9X,I3,6X,A5,7X,E10.3,5X,F11.3)
C**      NPQ = NOPQ(NPND)
C
         elstrt = ELST(NPND)
         VSTART = VSTR
         VEND = VOL2
         WRITE(IDIAGN,410) elstrt,VSTART,VEND
  410    FORMAT(/,' ',10X,'POND INIT. SURF. ELEV.(M)=',F8.3,
     -   2X,'INIT. AND FINAL VOL. (M**3) = ',2(F11.3,2X),/)
      LASTNB=NB(1)
C         DO OVER ALL ELEMENTS
      DO 106 NE=2,NTELE
      IF(NB(NE).GT.LASTNB) LASTNB=NB(NE)
  106 CONTINUE
      MBT=LASTNB+NI
      IF(LASTNB.EQ.0) MBT=1
      IF(MBT.LE.LENQS) GO TO 107
   17 WRITE(IDIAGN,108)
  108 FORMAT(' QS NEEDS TO BE DIMENSIONED LARGER')
  107 DO 111 L=1,NI
      MM=L-1
      QS(MM+MBT)=Q(L)
      SCON(MM+MBT)=CONC(L)
  111 CONTINUE
      NB(J)=MBT
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE HOFV (N,M,HOUT,H,VIN,VS,IS,HZR,I)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
C
      DIMENSION VS(10,3), H(10,3)
      I=IS
C THIS SUBROUTINE DOES LOG. INTERPOLATION FOR TABULATED FUNCTION H(VIN)
C   WHERE H = HZR AT VIN = 0.
      IF (VIN.GE.VS(1,N)) GO TO 10
      HOUT=HZR
      GO TO 40
   10 IF (VS(I,N).GT.VIN) GO TO 20
      I=I+1
      IF (I.GT.M) GO TO 50
      GO TO 10
   20 IF (VIN.GE.VS(I-1,N)) GO TO 30
      I=I-1
      IF (I.LE.1) GO TO 50
      GO TO 20
   30 IL=I-1
      RAT=ALOG(VIN/VS(IL,N))/ALOG(VS(I,N)/VS(IL,N))
      HOUT=HZR+(H(IL,N)-HZR)*EXP(RAT*ALOG((H(I,N)-HZR)/(H(IL,N)-HZR)))
   40 RETURN
   50 WRITE (IDIAGN,60) VIN
      STOP 9107
C
   60 FORMAT (' VOL OUT OF RANGE IN HOFV =',E12.5)
      END
C-----------------------------------------------------------------------
      SUBROUTINE SEDCOM(J,qil,DZ,AH,AHL,CR,NM,LCH)
C
C  THIS VERSION IS REVISED FOR EUROSEM BY R.SMITH, 4/91.
C
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1, CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), rfr(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
C
      COMMON /EROSN/ V1(20), V2(20), CT(20), CMX(20), RDET(20), DLA(20),
     1      DrlNET(20), Dirnet(20), suspir, susprl, SPOW(20)
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      common /intrill/ intrt, kint, d1(20), d2(20), vir(20), ql3(20), 
     1                 cbal, dcx, slrat, izr
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
C
      DIMENSION AH(20), AHL(20), SUP(20), supl(20,2), rdetl(20,2),
     1          cl(20,2), cml(20,2), eronet(20), vn(20)
      REAL IPAV
      data supl,rdetl,cl,cml/40*0.,40*0.,40*0.,40*0./
C
C  LCH   is 0 for interrill area
C           1 for rills
C           2 for channels
C *** *** *** DT IS IN SECONDS, AH IS m**2
      IPAV=1.-PAVE(J)
        WI = THETA
      WC=1.-WI
      susprl = 0.
      LI = LCH+1
      If(LI .gt. 2) LI=2
C GET TRANSPORT CAPACITY FOR LOCAL CONDITIONS
      CALL CAPAC(lch,spow,ah,J,1,NM,cmx)                                    !CAPAC
      qirc = 0.
      qrc = 0.
      IF (LCH .LE. 0) THEN
C ESTIMATE SPLASH EROSION FOR UPLAND FLOW:
        CALL drips (J,NM,LCH)                                           !DRIPS
C The next lines match initial conditions to splash erosion 
C    RDET from DRIPS in cu. m per unit of bed area per sec.
C
        do 5 i = 1,nm
          If( -dirnet(i) .ge. dero(J)) rdet(i) = 0.
          IF (AHL(i).LE.0.0 .and. ql2(i) .gt. 1.e-6) then
              cl(i,LI) = rdet(i)/(ql2(i)+vs)                            !New
          end if
          IF (i .eq. 1 .and. AH(1).LE.0.0 .and. ct(1) .le. 0.) ! **
     +       ct(1) = rdet(1)/(ql2(1)+vs)                                !New
          eronet(i) = dirnet(i)*(1.-porpl(j))
          if(kint .ge. 1) then
            vn(i) = vir(i)
          else
            vn(i) = v2(i)
          end if
    5   continue
        bw = tw(i)            !  new
      ELSE
C OR GET LATERAL INFLOW FOR CHANNELS or RILLS:
       Do 9 i=1,NM
        eronet(i) = drlnet(i)*(1.-porpl(j))
        rdet(I) = 0.
        vn(i) = v2(i)
    9  continue
        If(LCH .le. 1) Then
C                     Rainsplash into rills:
          Call DRIPS(J,nm,lch)                                          !DRIPS
          qlm = qil
          qlml = qll
          If(INTRT .le. 0) cr = cbal
          qls = rdet(1)*bwid + qil*cr
          qrc = qil*cr
          qirc = qll*crl
          qlsl = rdetl(1,LI)
          BW = rillw1(1)
          CRL = CR
          qll = qil
        Else
C   True Channels
          qlm = ql(l)
          qlml = ql(l-1)
          qls = conl(l)*qlm
          qlsl = conl(l-1)*qlml
         do 4 i=1,nm
           eronet(i) = drlnet(i)*(1.-porch(j))
    4    Continue
          BW = B
        End If
C   For rills or channels, qil is lateral inflow in sq m/sec
        DO I=1,NM
          RDET(I)=QIL*CR + rdet(i)*bw     ! rdet units now m2/sec
        End Do
C   For rills or channels, RDET is in sq m/sec
        IF (AHL(1) .LE. 0.0)
     1        CL(1,LI)=(QLSL)/(QLML+VS*BW)
        IF (AH(1) .LE. 0.0)
     2        CT(1)=(QLS)/(QLM+VS*BW)
      END IF
C*********************************************************  I=1,NM  LOOP
      DO 80 I=1,NM
 !?      If(I .le. 1) supl(i,LI) = 0.
        DTRF = 0.
        DTRC = 0.
        AI = AH(I)
        AL = AHL(I)
        V  = VN(I)
        VL = V1(I)
        IL = I - 1
        CP = CL(I,LI)
        CGU = 0.
        PAVF=1.0
C* FIND MEAN TRANSPORT CAPACITY TO ANTICIPATE DEPOSIT OR EROS. CASE:
C
        VSU = VS
        fldet = SUPL(I,LI)
        IF(AI .LT. 1.E-8) GO TO 50
C
C      if(danet(i)*(1.-porch(j)) .le. fldet*dt) pavf = ipav
      If(eronet(i) .le. fldet*dt) pavf = ipav
      PRO = 1.
      if(cl(i,li) .lt. cml(i,LI)) then                                     !cml
C Erosion indicated: Was this node depositing at the last step?
C        if(cml(I,LI) .le. cl(i,LI)) supl(i,LI) = 0.
        IF(eronet(i) .le. 0. .and. cmx(i) .gt. 0.) Then
          If(LCH .le. 0 .and. -dirnet(i) .gt. dero(J)) Then
            pavf = 0.
          Else
C Erodibility reduction for cohesive soil:
           cdef = (cml(i,LI)-cl(i,LI))/cmx(i)                               !cml
           pro =  fbeta(cdef,cohe(j))                                    !FBETA
          End If
C        If(lch.eq.0 .and. i.eq.5)write(idiagn,797)cdef,pro
C 797  Format('  cdef(5),pro: ',2g12.5)
        End If
      else
C Deposition indicated. Was this node formerly eroding?
        if(cml(i,LI) .gt. cl(i,LI)) supl(i,LI) = 0.
CC        If(ahl(i).gt. 0.1e-6) vst = ahl(i)/tw(i)/dt
CC        vsu = amin1(vst,vs)
      end if
C
      if(ah(i) .le. 0.) then
        write(*,979) i,cl(i,LI),ai,cmx(i)
 979  Format(' i,cl,a,cmx:',i3,3g11.4)
       stop' neg. or zero depth with deposition in sedcom'
      end if
C
   18   CGU = VSU*PRO
C*
        CEX = CMX(I)
        DTRF = PAVF*CGU*TW(I)
C**                      TO DISTINGUISH PLANES/CHANNELS
        IF(LCH .LE. 0) THEN
            DTRC = DTRF*CEX + PAVF*RDET(I)
        ELSE
            DTRC = DTRF*CEX + RDET(I)
        ENDIF
   20 CONTINUE
C
        IF(I .LE. 1) GO TO 60
C        IL = I - 1
  40  Continue
      supc = wi*dtrc + wc*supl(i,li)                                     try
C      AF=(CP*AL+CPM*ALM-CT(IL)*AM)*.5/DT
      AF = (CP*AL)/DT
      BF=(CT(IL)*AM*VM*WI - (CP*AL*VL-CPM*ALM*VLM)*WC)/DZ
C      DF=(0.5/DT+V*WI/DZ)
      DF = (1./DT + V*WI/DZ)
      CT(I) = (SUPC+AF+BF)/(AI*DF +  wi*DTRF)                            try
C
      IF (CT(I).LT.0.) CT(I)=0.
      CT(I)=AMIN1(CT(I),0.5)
      GO TO 70
C*      ZERO DEPTH CASE:
   50   IF (I .GT. 1) then
          CT(I) = 0.0
        end if
        IF(AI .GT. AL) CT(I) = CL(I,LI)
      GO TO 70
C*      FIRST NODE ONLY:
   60 CONTINUE
      IF (AI.LE.0.) THEN
C*      FIRST NODE IS ZERO DEPTH:
        DTRC = PAVF*RDET(I)
          IF(AH(2) .LE. 0.) DTRC = 0.
        CT(1)=CL(1,LI)
      ENDIF
   70 CONTINUE
C
      AM = AI
      ALM=AL
      VM=V
      VLM=VL
      CPM=CP
C   Sup and Supl are net supply to flow from bed in cu m/m length/ sec.
      SUP(I) = DTRC - CT(I)*DTRF
C
      IF(I .EQ. 1) THEN 
  ! calculation of sup(1) for AI .le. 0. moved from here
        volsed = 0.
      Else
        volsed = 0.5*(AI*CT(I)+AM*CT(IL))*DZ
      ENDIF
      IF (LCH.LE.1) THEN
        if(AI .gt. 0.) THEN
          prop = 1.                                                     *N11/91
          If(LCH .eq. 1) Then     ! RILLS:
            If(tw(i) .gt. rillw1(i) .and. rilld1(i) .gt. 0.) Then     
              prop = (rillw1(i) + 2.*(zrl(j)+1.)*rilld1(i))/(tw(i)+        "
     &                2.*rilld1(i))                                        "
              if(prop. gt. 1.) prop = 1.                                   "
            End If
          End If
          DLA(I) = prop*(WC*(SUPL(I,LI)-qirc)+WI*(SUP(I)-qrc))*DT/
     &            (1.-PORPL(J))                                            "
        ELSE
          DLA(I) = 0.
        ENDIF
        If(LCH .le. 0) Then
          dirnet(i) = dirnet(i) - dla(i)
          suspir = suspir + VOLSED
C          if(i .ge. nm) write(*,102) suspir,volsed
  102  Format(30x,2g12.5)
        Else
          drlnet(i) = drlnet(i) - dla(i)
          susprl = susprl + VOLSED
        End If
      ELSE                                ! Channels only
C* DELETE AMOUNT THAT COMES FROM OUTSIDE SOURCE FOR DELH IN CHANN:
        IF (AI .GT. 0.) THEN
          DLA(I) = (WC*(SUPL(I,LI)-RDETL(I,LI))+
     *              WI*(SUP(I)-RDET(I)))*DT/(1.-PORCH(J))
        ELSE
          DLA(I) = 0.
        ENDIF
        IF(PAVE(J) .gt. 0.99) THEN
          IF(DLA(I) .LT. 0.0) DLA(I)=0.0
        ENDIF
        drlnet(i) = drlnet(i) - dla(i)
        susprl = susprl + VOLSED
      END IF
      IF(AH(1).LE.0. .and. ct(1) .le. 0.) SUP(1) = 
     &           RDET(1) + VSU*tw(1)*(CMX(1)-CT(1))
C
   80 CONTINUE
C ************************************************************* END LOOP
      DO 90 I=1,NM
        SUPL(I,LI)=SUP(I)
        IF(I .GT. 1 .AND. CT(I) .LE. 0.) SUPL(I,LI) = 0.
        RDETL(I,LI) = RDET(I)
        CML(I,LI)=CMX(I)
   90   CL(I,LI)=CT(I)
        IF(NP(J) .EQ. 7 ) THEN
           If(lch.le.0) Then
           WRITE(IDIAGN,96) (d2(I),I=1,NM)
C           WRITE(IDIAGN,92) (pavf*rdet(I),I=1,NM)
           End If
           ptim = t(L)/60.
           qlat = ql2(nm)*1000.*3600.
           write(IDIAGN,95) ptim,qlat,sup(nm),pavf*rdet(nm),cbal
   95 Format(5(2x,g13.4))
           write(idiagn,93) (cmx(i),i=1,nm)
           WRITE(IDIAGN,94) (ct(I),I=1,NM)
 !          write(idiagn,94) qil,cr
C           End If
        ENDIF
   91   FORMAT(' SUP(I=1,NK)=',5(G11.4)/(6G11.4))
   92   FORMAT(' rdt(I=1,NK)=',5(G11.4)/(6G11.4))
   93   FORMAT(' cmx(I=1,NK)=',6(F10.6)/(7F10.6))
   94   FORMAT(' ct(I=1,NK)=', 6(F10.6)/(7F10.6))
   96   FORMAT(' d2(I=1,NK)=', 6(F10.6)/(7F10.6))
        RETURN
        END
C-----------------------------------------------------------------------
      Function FBETA(CDEL, COH)
C This subroutine calculates a reduction factor for the erosion
C   resistance due to soil cohesion, based on experimental data
C   of Govers and others, as designed by the Eurosem Committee
C      CDEL is Scaled Concentration deficit, (Cmx-C)/Cmx
C      COH is soil cohesion in kPa
C
      az = 0.79
      bz = 8.50e-1
      If(COH .le. 1.) Then
        Bo = .335
      Else
        Bo = az*exp(-bz*COH)
      End If
      If(Bo .gt. 1.0) Bo = 1.0
C      Fbeta = Bo + (1.-Bo)*exp(-25.*CDEL*cdel)
      Fbeta = Bo
C      write(*,10) COH,CDEL,fbeta
   10 Format(20x' coh, cs, Beta:',3g13.5)
      Return
      End
C-----------------------------------------------------------------------
      SUBROUTINE SEDIV (MEDS,SIGS,NPART,DVS)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      DIMENSION DVS(7), V(200)
C  THIS SUBROUTINE COMPUTES FRACTIONS OF A SEDIMENT SAMPLE IN
C  N EQUAL PARTS GIVEN A MEDIAN AND STD. DEV. AND ASSUMING LOG
C  NORMAL DISTRIBUTION
      REAL MUS, MUSQ, MUL, MEDS, MEDSQ
        FV(Y)=EXP(MUL+SDL*Y)
      VARS=SIGS*SIGS
      MEDSQ=MEDS*MEDS
      MUSQ=0.5*(MEDSQ+SQRT(MEDSQ*MEDSQ+4.*MEDSQ*VARS))
      MUS=SQRT(MUSQ)
      VARL=ALOG(VARS/MUSQ+1.)
      SDL=SQRT(VARL)
      MUL=ALOG(MUS)-.5*VARL
      ND=(NPART/2)*2+1
      IF (ND.NE.NPART) NPART=ND
C ** GENERATE DISTRIBUTION SAMPLE SIZE N
      N=25*NPART
      NU=N+1
      DO 10 I=2,NU,2
      CALL URAN (R1)
      CALL URAN (R2)
      F1=SQRT(-2.*ALOG(R1))
      F2=6.28318531*R2
      Y1=F1*COS(F2)
      Y2=F1*SIN(F2)
      V(I-1)=FV(Y1)
      V(I)=FV(Y2)
   10 CONTINUE
C ORDER SAMPLE DIAMETERS
      NM=N-1
      DO 30 I=1,NM
      VMIN=V(I)
      IP=I+1
      DO 20 J=IP,N
      IF (V(J).GE.VMIN) GO TO 20
      V(I)=V(J)
      V(J)=VMIN
      VMIN=V(I)
   20 CONTINUE
   30 CONTINUE
      ICHECK=25
C ** GET WEIGHTED DIAM. ACCORDING TO APPROX FALL VELOCITY
      FSUM=0.
      K=1
      DO 40 I=1,N
      FSUM=FSUM+V(I)*V(I)
      IF (I.LT.ICHECK) GO TO 40
      DVS(K)=SQRT(FSUM/25.)
      K=K+1
      ICHECK=ICHECK+25
      FSUM=0.
   40 CONTINUE
      IM=(ND+1)/2
      IF (ABS(DVS(IM)-MUS).GT.1.E-04) WRITE (IDIAGN,50) DVS(IM),MUS
      RETURN
C
   50 FORMAT (//,' ',10X,'RELEVANT SEDIMENT CHARACTERISTICS',/,
     -             ' ',10X,'---------------------------------',/,
     -             ' ','CALCULATED MODE AND MEAN = ',2F12.6,/)
      END
C**************************************************************************
C
      SUBROUTINE  URAN(V)
C
      DIMENSION K(4)
C
      DATA K(1),K(2),K(3),K(4)/37,2283,3,9911/
C
      K(4)    = 3*K(4)+K(2)
      K(3)    = 3*K(3)+K(1)
      K(2)    = 3*K(2)
      K(1)    = 3*K(1)
      I       = K(1)/1000
      K(1)    = K(1)-I*1000
      K(2)    = K(2)+I
      I       = K(2)/100
      K(2)    = K(2)-100*I
      K(3)    = K(3)+I
      I       = K(3)/1000
      K(3)    = K(3)-I*1000
      K(4)    = K(4)+I
      I       = K(4)/100
      K(4)    = K(4)-100*I
      V       = (((FLOAT(K(1))* .001+FLOAT(K(2)))*.01
     1           + FLOAT(K(3)))*.001+FLOAT(K(4)))*.01
      RETURN
      END
C*********************************************************************
      SUBROUTINE VOFH (N,M,VOUT,VS,XIN,X,IS,XZR,I)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
C
      DIMENSION VS(20,3), X(20,3)
      I=IS
C THIS SUBROUTINE DOES LOG. INTERPOLATION FOR TABULATED FUNCTION VS(XIN)
C    WHERE VS = 0 AT X = XZR
      IF (XIN.GE.X(1,N)) GO TO 10
      VOUT=0.
      GO TO 40
   10 IF (X(I,N).GT.XIN) GO TO 20
      I=I+1
      IF (I.GT.M) GO TO 50
      GO TO 10
   20 IF (XIN.GE.X(I-1,N)) GO TO 30
      I=I-1
      IF (I.LE.1) GO TO 50
      GO TO 20
   30 IL=I-1
      RAT=ALOG((XIN-XZR)/(X(IL,N)-XZR))/ALOG((X(I,N)-XZR)/(X(IL,N)-XZR))
      VOUT=VS(IL,N)*EXP(RAT*ALOG(VS(I,N)/VS(IL,N)))
   40 RETURN
   50 WRITE (IDIAGN,60) XIN
      STOP 9106
C
   60 FORMAT ('  ELV OUT OF RANGE IN VOFH =',E12.5)
      END
C*********************************************************************
      SUBROUTINE VOLCAL (N,VOL,SAR,DX)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /CNTRL/ NRES, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /PONDS/ AP(10,10,3), TWP(10,10,3), EL(10,3), QSO(20,3),
     1   ELQ(20,3),ELZQ(3),ELVZ(3),NOPQ(3),MEL(3),NSEC(3),XSEC(10,3),
     2   ELZX(10), JOP(3), ELST(3), ALV, DIFUS(3), SIGMAS(60), NPART
C
      DIMENSION VOL(10,3), DX(10,3), SAR(10,3)
C
      NS=NSEC(N)
      DO 10 K=2,NS
        KL=K-1
   10 DX(K,N)=XSEC(K,N)-XSEC(KL,N)
      VOL(1,N)=0.01
      SAR(1,N)=.05
      M=MEL(N)
      DO 80 J=1,M
       VOL(J,N)=0.
       SAR(J,N)=0.
       DO 70 K=2,NS
        KL=K-1
        AH=AP(K,J,N)
        AB=AP(KL,J,N)
        IF (AH.GT.0.) GO TO 20
        IF (AB.LE.0.) GO TO 70
        KLL=KL-1
        IF (KLL.LE.0) GO TO 30
        ABB=AP(KLL,J,N)
        IF (ABB.LE.AB) GO TO 30
        DEL=AB*DX(KL,N)/(ABB-AB)
        GO TO 50
   20   IF (AB.GT.0.) GO TO 30
        KP=K+1
        IF (KP.GT.NS) GO TO 30
        AF=AP(KP,J,N)
        IF (AF.GT.AH) GO TO 40
   30   DEL=DX(K,N)
        GO TO 60
   40   DEL=AH*DX(KP,N)/(AF-AH)
   50   DEL=AMIN1(DEL,DX(K,N))
   60   VOL(J,N)=VOL(J,N)+DEL/3.*(AH+AB+SQRT(AH*AB))
        SAR(J,N)=SAR(J,N)+DEL*0.5*(TWP(K,J,N)+TWP(K-1,J,N))
   70  CONTINUE
   80 CONTINUE
      WRITE (IDIAGN,90) N,N
      WRITE(IDIAGN,110)
      DO 95 K=1,M
        ELP = EL(K,N)
        VP  = VOL(K,N)
        SP  = SAR(K,N)
        WRITE (IDIAGN,120) ELP,VP,SP
   95 CONTINUE
      RETURN
C
   90 FORMAT(' ',15X,'*** POND NO. ',I2,' DIAGNOSTIC INFORMATION ***',
     -//,' ',2X,'ELEV., VOLUMES AND SURFACE AREAS FOR POND',I2,/)
  110 FORMAT  (' ',5X,' (M)      (M**3)           (M**2)        ',/,
     -         ' ',5X,'-----     -------      -------------')
  120 FORMAT(' ',3X,F8.3,1X,F11.3,4X,F11.3)
      END
C-----------------------------------------------------------------------
      FUNCTION VSETL (SS,D,VNU)
C  Drag coefficient calculation of fall velocity
C     asuming spherical particle
C  Input vnu is in sq m/sec   Input D in microns.
      grav = 9.81
      dm = D*1.e-6  ! diameter changed to meters
C
      ITER=0
      CA=24*VNU/dm
C      CK=42.91*(SS-1.)*D
      CK = grav*1.3333*(SS-1.)*dm
      CB=3.*SQRT(VNU/dm)
      CC=.34
      VS=8.E05*dm*dm
   10 F = CA*VS+CB*VS**1.5+CC*VS*VS-CK
      FP = CA+1.5*CB*SQRT(VS)+2*CC*VS
      DVS = F/FP
      IF (ABS(DVS) .gt. 0.0001*VS) Then
        VS = VS-DVS
        ITER = ITER+1
        IF (ITER .GT. 40) GO TO 40
        IF (VS .LE. 0.) VS=0.5*(VS+DVS)
        GO TO 10
      Else
        VSETL = VS
        RETURN
      End If
   40 WRITE (6,50) SS,D,VNU,DVS
      STOP ' Convergence failure in finding settling velocity'
C
   50 FORMAT (' 40 ITERS WITH SS, D, VNU, AND DVS =',4G13.5)
      END



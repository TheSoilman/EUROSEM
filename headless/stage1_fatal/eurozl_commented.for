C-----------------------------------------------------------------------
C  EUROZL.FOR -- infiltration, implicit flow solver and helpers
C
C    UNIF    resample an outflow series onto equal time steps
C    XPLINF  infiltration at every node (Smith-Parlange method)
C    IMPLCT  implicit solution of the kinematic-wave flow equation, for
C            planes (depth) and channels (area)
C    IMPOCF  residual and derivative for the plane Newton solution
C    ITER    Newton-Raphson root finder (with bisection fallback)
C    PERIM   flow geometry of a plane, rill or rill with overflow
C    QCP     soil-profile storage at the start of ponding
C    HFUN    depth of a trapezoid from its area
C-----------------------------------------------------------------------
C**************************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE UNIF(Q,T,NIN,QO,NI,DELT,CONC,SCON)
C    The flow solution uses uneven time steps (the step is cut at
C    rainfall breakpoints and at TFIN).  The next element needs the
C    outflow on EVEN steps.  UNIF takes the NIN values Q(T) and
C    CONC(T) and, by linear interpolation, returns NI values at
C    intervals DELT in QO and SCON.  The first and last values are
C    copied as they are.
C-----------------------------------------------------------------------
      SUBROUTINE UNIF (Q,T,NIN,QO,NI,DELT,CONC,SCON)
C
      DIMENSION Q(NIN), T(NIN), QO(NI), CONC(NIN), SCON(NI)
C
C     THIS SUBROUTINE TAKES THE VALUES OF Q AT THEIR CORRESPONDING TIME
C     CONVERTS THEM INTO VALUES WITH EQUAL TIME INTERVALS (DELT).  THESE
C     VALUES ARE STORED IN QO.
      I=2
      K=1
      QO(1)=Q(1)
      SCON(1)=CONC(1)
      TO=0.
   10 TO=TO+DELT
   20 IF (T(K).GE.TO) GO TO 30
      K=K+1
      GO TO 20
   30 IF (ABS(T(K)-TO).LE.1.E-5) GO TO 50
      QO(I)=Q(K-1)+(TO-T(K-1))/(T(K)-T(K-1))*(Q(K)-Q(K-1))
      SCON(I)=CONC(K-1)+(TO-T(K-1))/(T(K)-T(K-1))*(CONC(K)-CONC(K-1))
   40 I=I+1
      IF (I.GE.NI) GO TO 60
      GO TO 10
   50 QO(I)=Q(K)
      SCON(I)=CONC(K)
      GO TO 40
   60 QO(NI)=Q(NIN)
      SCON(NI)=CONC(NIN)
      RETURN
      END
C***********************************************************************
C
C-----------------------------------------------------------------------
C  SUBROUTINE XPLINF(J,RA,QF,NS,NIR)
C    Infiltration at every node of element J for one time step.
C    Smith-Parlange infiltration ("decay function by the
C    Parlange-Smith method 2"), with redistribution of soil moisture
C    during breaks in the rain.  Version for spatially variable
C    interrill flow: two-dimensional variation (NS nodes down the slope
C    by NIR interrill sub-nodes).
C
C    ARGUMENTS
C      J    element
C      RA   rainfall rate at the node (m/s), after interception
C      QF   (out) net lateral inflow at each node: rain minus
C           infiltration (m/s).  Negative while surface water is still
C           draining into the soil.
C      NS   number of nodes down the slope
C      NIR  number of interrill sub-nodes (1 = main flow path only)
C
C    INTERNAL UNITS: mm and minutes (RFM mm/min, QF mm/min until the
C    last loop, which converts back to m/s).
C
C    SOIL PARAMETERS (per element J)
C      FMN = FMIN(J)   final (saturated) infiltration rate Ks
C      AL(J)           capillary drive G
C      THI, THMX       initial and maximum soil moisture
C      APAR = G * (THMX - THIN) * RVAL   storage parameter; THIN is the
C                      current initial moisture of the node
C      TOTF            cumulative infiltration at the node
C      FCL             infiltration capacity at the node
C      NRO             state of the node:  0 no ponding yet,
C                      1 supply infiltrated fully this step after
C                      ponding began, 2 runoff (excess rain),
C                      -1 recently in redistribution (rain below Ks)
C      THIN            current soil moisture, rising or falling with
C                      the rain history
C
C    STATE KEPT BETWEEN CALLS
C      TOTF, FCL, NRO, APAR (the SAVE statement) hold the history of
C      every node.  They are reset on the first time step of each
C      element (T(L) <= 1.01*DT).
C
C    NOTE: RVAL, the rock-content factor for the profile, is ONE value
C    for the whole run (see PAREAD in emainl): it is not set per
C    element here.  The line that would do so (RVAL = 1 - ROC(J)) is
C    commented out further down.
C    NOTE: stops the run ('bad thin') if the moisture reaches THMX.
C-----------------------------------------------------------------------
      SUBROUTINE XPLINF (J,RA,QF,NS,NIR)
C
C      OPERATING DIMENSIONS ARE MM AND MINUTES
C
C  Version for spatially variable interrill flow: 2 dimensional variation
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /INFIL/ AL(60), thi(60), thmx(60), ROC(60), RECS(60),
     1        THR, thfc, tempw, stone(60), RVAL
C
      common /intrill/ intrt, kint, d1(20), d2(20), vir(20), ql3(20), 
     1                 cbal, dcx, slrat, izr
C
C      common /splice/ thinp(1000)
      DIMENSION QF(20), THIN(20,20), FCL(20,20), NRO(20,20), 
     1 TOTF(20,20), S1(20), APAR(20,20), DRATE(20,20), TIN(20,20)
C      Logical FIRST
C!!
      SAVE TOTF,FCL,NRO,AVF,apar,FMN
C!!
C  DPAR is the redistribution-rate parameter, 0.95, based on Walnut
C  Gulch (Lucky Hills 106) data.
C     SET DPAR = 0.95 AS BASED ON WALNUT GULCH LUCKY HILLS 106 DATA
      DATA DPAR/0.95/
C
C  DECAY FUNCTION BY PARLANGE-SMITH METHOD 2
C  CONVERT sec to MIN and m/sec TO MM/MIN:
C  Convert the rain rate from m/s to mm/min, and the time step from
C  seconds to minutes.
      RFM=RA*60000.
      DMT=DT/60.
C      TM=T(L)/60.
      fmn = fmin(j)
C
C  First time step of this element: set every node to its starting
C  state (capacity large, no runoff, no cumulative infiltration, soil
C  moisture THI).
      IF (T(L).LE.1.01*DT) THEN
C  FIRST TIME THROUGH:
        JRP=0
C  Available Porosity correction for soil profile rocks:
C        RVAL = 1. - ROC(J)
        DO 10 K=1,NS
        DO 10 M=1,nir
          FCL(K,M) = AMAX1(RFM,1000.)
          APAR(K,M) = AL(J)*(thmx(J)-thi(J))*RVAL
          NRO(K,M)=0
          TOTF(K,M)=0.
          TIN(K,M) =0.
          THIN(K,M)=thi(J)
   10   Continue
      ENDIF
C
C  Surface water depth that supplies infiltration: SU is the depth at
C  the first node and NN the number of nodes in use.  With interrill
C  routing (NIR > 1) the interrill depths D1 are used instead of H1.
      If(NIR .le. 1) Then
        su = h1(1)
        nn = ns
        miz = 1
        mzp = 2
C  Can't properly use RECS(j) (as read) for rill flow:
C       If(KINT .gt. 0) recsf = 1.e-5
      Else
        su = d1(1)
        nn = nir
        miz = 2
        mzp = 3
      End If
C  Frozen-soil surrogate: if the water temperature is below zero the
C  profile can store only (THMX - THIN) * DERO.  Once the cumulative
C  infiltration reaches that, the capillary drive is set to zero, so
C  the infiltration rate falls to the constant rate FMN.
C  Surrogate modification for frozen soil at depth DERO(J):             !4/93
      If(TEMPW .lt. 0.) Then 
        CAPC = (THMX(J) - Thin(1,MIZ))*DERO(J)
        If(TOTF(1,miz) .ge. CAPC)   AL(J) = 0.
      End If
C
C  No capillary drive (AL = 0), or negligible Ks: constant infiltration
C  rate.  Infiltration is the smaller of Ks and the rain rate.
      IF (AL(J) .LE. 0. .or. FMN .LE. 1.E-7) THEN
C ****  CONSTANT INFIL RATE INDICATED BY  AL = 0.
        AVF = AMIN1(FMN,RFM)
        FCL(1,miz) = AVF
        QF(miz) = RFM - AVF
        TOTF(1,miz)=TOTF(1,miz) + AVF*DMT
        I = amin0(NS,2)
        M = amin0(mzp,nir)
        NRO(1,miz) = 1
        GO TO 90
      ENDIF
C
C OBTAIN LOCAL WATER DEPTHS IN INCHES
C      SU=THETA*H2(1)+(1.-THETA)*H1(2)
C  RECS is in m
C       recsf = recs(j)*1000.
C  Local surface water depth at each node, in mm (S1).
      S1(miz) = 1000.*SU
      AVF = RFM
      DO 20 I=mzp,NN
C      Do 20 k=2,NN
      If(nir .le. 1) Then
        Hef = H1(I)
      Else
        Hef = d1(i-1)
      End If
      S1(I)=1000.*Hef
C       AVF(I)=RFM
   20 CONTINUE
C      If(NIR .gt.1) WRite(*,108)(s1(i),i=1,NN)
  108 Format(10G12.4)
C
C######## BEGIN LOOP THROUGH ALL NODES
C
C-----------------------------------------------------------------------
C  LOOP OVER THE NODES (I down the slope, M interrill sub-nodes).
C  For each node:
C    SURFIN  surface water available to infiltrate, in mm/min (S1 / DT
C            / RECES; RECES is set to 1 here, so the full depth counts)
C    RFT     total supply = rain rate + SURFIN
C  Then one of three cases:
C    RFT > FMN   supply exceeds Ks: infiltration from the capacity
C                (FCP, FINTG), with the excess leaving as runoff
C    RFT <= FMN  supply below Ks: no ponding; the soil redistributes
C                its moisture (below)
C    nodes with identical history are not recomputed: the result is
C    copied down the slope at label 90.
C-----------------------------------------------------------------------
      DO 80 I = 1,NS
      DO 80 M = miz,NIR
        If(ns .le. 1) MI = M
        If(NIR .le. 1) MI = I
        reces = 0.0
        SURFIN = 0.
        IF(NRO(I,M) .EQ. 0) THEN
          QF(MI) = 0.
C NRO = 0.: DON'T REPEAT CALCULATIONS FOR IDENTICAL CASES:
          IF (JRP.LE.0) THEN
            IF (MI.GT.miz) THEN
              IF (S1(MI-1).LE.0.) GO TO 90
            ENDIF
          ENDIF
        ELSE
C* EFFECTIVE SURFACE INFLOW SUPPLY TO NODE IN mm/MIN.:
C         RATH = S1(MI)/RECSF
C         IF(RATH .LT. 0.2) RATH = 0.2
          IF(S1(MI) .GT. 0.) THEN
C           RECES = AMIN1(1.,RATH)
            reces = 1.
            SURFIN=S1(MI)/DMT/RECES                                       *N1091
C           If(m.eq.4.and.i.eq.1)
C     +   write(*,109) nro(i,m),mi,s1(mi),reces,rfm
 109  Format('nro, mi',2i3,4g13.5)
          END IF
        END IF
C * TOTAL EFFECTIVE AVAILABLE INFLUX SUPPLY:
C  Total supply to the node.
        RFT = RFM + SURFIN
C LOGICAL CHECK OF CURRENT RUNOFF CONDITION
C  Supply exceeds Ks.  FCP gives the infiltration capacity for the
C  current cumulative infiltration TOTF; FINTG integrates it over the
C  step to give the new TOTF (QSN) and the mean infiltration rate AVF.
C  QF = rain - infiltration; if this is not positive everything
C  infiltrates (NRO = 1) and THIN is reset in preparation for a break
C  in the rain.
        IF(RFT .GT. FMN) THEN
          DTE=DMT
          FCL(I,M) = FCP(TOTF(I,M),RFT,APAR(I,M),FMN)
          CALL FINTG (DTE,TOTF(I,M),RFT,FCL(I,M),FMN,APAR(I,M),QSN)
          AVF = (QSN-TOTF(I,M))/DTE
          TOTF(I,M) = QSN
          FCL(I,M) = FCP(QSN,RFT,APAR(I,M),FMN)
          QF(MI)=RFM - AVF
          IF(QF(MI) .LE. 0.) THEN
            NRO(I,M) = 1
C * RESET THIN(I) IN PREP. FOR HIATUS IF NECESS:
            If(RFM .gt. FMN) Then
            vqc = qcp(rfm,apar(I,M),fmn)
            qsu = min(vqc,qsn)
            THIN(I,M) = 0.99*(thi(j) +(thmx(J)-thi(j))
     &                    *(QSU/vqc)**.15)
C            If(qsn .gt. vqc) Then
C              write(IDIAGN,945)L,qsn,vqc,thmx(j),fcl(i,m),rfm,
C     & avf,dte, surfin, fmn
C  945 Format(3x,I3,(5g13.5))
C            End If
            End If
          ELSE
C  JRP = 1 indicates redistribution has begun lower on the slope; the
C  initial saturations after the break are written to the diagnostic
C  file once.  Otherwise the node is in runoff (NRO = 2); if the first
C  node of an interrill column shows uniform rain excess the remaining
C  nodes copy it.
            IF(JRP .EQ. 1) THEN
              IF(NRO(I,M) .LT. 2. AND. I .EQ. NS) THEN
                WRITE(IDIAGN,130) (THIN(K,M),K=1,NS)
                JRP = 2
              END IF
              NRO(I,M) = 2
            ELSE
              NRO(I,M) = 2
C* TEST FOR UNIFORM RF EXCESS AND SHORTCUT:
              IF(MI .GT. MIZ) THEN
                IF(QF(MI) .GT. 0.  .AND. s1(MIZ) .LE. 0.) GO TO 90
              END IF
            END IF
          END IF
C-----------------------------------------------------------------------
C  RFT <= FMN: rain (plus surface water) is below the saturated rate Ks,
C  so there is no ponding.
C    THULT   the moisture the soil tends towards at this supply rate:
C            0.99 * (THR + (THMX - THR) * (RFT / FMN) ** 0.1)
C    DRATE   rate at which the moisture THIN approaches THULT.  After
C            runoff it depends on the cumulative infiltration TOTF;
C            before runoff on the rain accumulated since (TIN).
C    THIN is moved towards THULT by  (THULT - THIN)*(1 - exp(-ADR*DT)).
C    QF is negative: surface water of depth S1 drains into the soil at
C    1.2 * S1 / DT.
C    TOTF is reset to zero and APAR is recomputed from the new THIN, so
C    that the next burst of runoff starts from the redistributed state.
C-----------------------------------------------------------------------
        ELSE
C * EFFECTIVE R .LT. KS:
C              * SOME REDISTRIBUTION AND CHANGE IN INIT. SAT.:
          thult = .99*(THr+(THMX(J)-thr)*(RFT/FMN)**0.1)
          IF(NRO(I,M) .GT. 0) THEN
C    MAKE REDUCTION RATE A FUNCTION OF KS (FMN) AS WELL NOTE:
C   DPAR SET = 0.95 (DATA Stmt. ABOVE) AS BASED ON WALNUT GULCH LH106 DATA
            IF(I .GT. 1) JRP = 1
            TIN(I,M)  = 0.
C      USE THR AND THFC FROM /SOIL/ AS COMP. FROM REGRESSION
C            SM = AMAX1(THR,thi(J))
            THULT = AMAX1(THI(j),THULT)
C  Following empirical expression converted to mm>
C            DRATE(I,M) = SQRT(DPAR * FMN)*EXP(-2.*TOTF(I,M))
            DRATE(I,M) = SQRT(DPAR*FMN/25.)*EXP(-TOTF(I,M)/12.)
            QF(MI) = -1.2*S1(MI)/DMT
          Else
C No runoff yet, add to TIN
            TIN(I,M) = TIN(I,M) + RFT*DMT
            If(THULT .gt. THI(J)) Then
              DRATE(I,M) = SQRT(DPAR*FMN/25.)*EXP(-TIN(I,M)/3.)
            Else
              DRATE(I,M) = SQRT(DPAR*FMN/25.)*.1
            End If
C * NO RUNOFF YET
C            QF(MI) = 0.
            QF(MI) = -1.2*S1(MI)/DMT
          END IF
C            NEW EXPONENTIAL FORM OF thINT DECAY
          ADR = DRATE(I,M) * (1. - RECES)
          othi = thin(i,m)
          THIN(I,M) = THIN(I,M) + (THULT - THIN(I,M))*
     @                       (1.0-EXP(-ADR*DMT))
C  Safety check: the moisture must stay below THMX, otherwise the run
C  stops with 'bad thin' and a diagnostic line.
          If(THIN(i,m) .ge. thmx(j)) then
       write(idiagn,946) i,nro(i,m),L,thult,drate(i,m),adr,rft,fmn
     & ,othi,thin(i,m)
  946 Format(2x,3i3,/2x,9g13.5)
           CALL FATAL('Bad THIN value in infiltration routine')
          End IF
          APAR(I,M) = AL(J)*(thmx(J)-THIN(I,M))*RVAL
          NRO(I,M) = -1
          TOTF(I,M) = 0.
        END IF
      IF(I .eq. 1 .and. M .eq. 1) Then
C        write(*,"(' R, k, thi, THIN, thult: ',5F12.6)") 
C     &    RFM, fmn,THI(j),THIN(1,1), thult
C        read(*,*) 
C        thinp(L) = thin(1,1)
      End If
   80 CONTINUE
C
C################END LOOP THROUGH SURFACE NODES #########
      GO TO 110
C  Nodes skipped above (identical cases) copy the infiltration state
C  of the node above, or of the previous sub-node.
   90 Continue
      DO 100 K=I,NS
      DO 100 N=M,NIR
      If(ns .le. 1) Then
        MI = N
        FCL(K,N)=FCL(K,N-1)
        NRO(K,N)=NRO(K,N-1)
        IF(NRO(K,N) .GE. 2) THIN(K,N) = THIN(K,N-1)
        TOTF(K,N)=TOTF(K,N-1)
C        
      Else If(NIR .le. 1) Then 
        MI = K
        FCL(K,N)=FCL(K-1,N)
        NRO(K,N)=NRO(K-1,N)
        IF(NRO(K,N) .GE. 2) THIN(K,N) = THIN(K-1,N)
        TOTF(K,N)=TOTF(K-1,N)
      End If
      QF(MI)=QF(MI-1)
  100 CONTINUE
C
C  Convert QF from mm/min to m/s before returning to PLANE.
  110 CONTINUE
C CONVERT UNITS BEFORE RETURNING TO SUBROUTINE PLANE:
      DO 120 I=1,NN
  120 QF(I)=QF(I)/60000.
      RETURN
  130 FORMAT (5X,' INITAL SATS. AFTER HIATUS: '/10(2X,F10.8))
  140 FORMAT (' ',9X,'TOTAL ACCUM INFIL(mm.)=',F8.4,' AT NODE',I3)
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C  SUBROUTINE IMPLCT(NM,J)
C    Implicit finite-difference solution of the kinematic-wave flow
C    equation for one time step, for element J with NM nodes:
C          dA/dt + dQ/dx = ql          Q = ALPHA * FV(depth)
C    with ALPHA = sqrt(S) / n (Manning) and FV from PERIM.  Weighting
C    between old and new time levels is THETA.
C
C    PLANE (W(J) <> 0): the unknown is the depth H2(K+1) at node K+1,
C      found from the node above, node by node down the slope.  Each
C      depth is a root of IMPOCF, found by ITER (Newton-Raphson).  Where
C      the flow front has not yet arrived the depth is found by tracing
C      the upstream characteristic in 5-second sub-steps.  With
C      interrill routing (KINT >= 1) the same code solves the interrill
C      strips (D1, D2, QL3 and slope DCX instead of H1, H2, QL2, DX).
C      RESULTS: H2, V2, A2, VR (rill velocity), WIDW (flow width).
C
C    CHANNEL (W(J) = 0): the unknown is the flow area A2(K+1) of a
C      trapezoidal channel (bottom width B, side slopes CO1, CO2), with
C      infiltration loss through the bed and banks (CHAINF).  The
C      residual is IMPCHA (in Euroa.for).  RESULT: A2 and the
C      infiltration CINF2.
C
C    ITER flag IER (after IER = IER + 1):
C          1 converged    2 no convergence in IEND passes (error 10)
C          3 zero derivative (error 11 for channels; for planes the
C            flag is accepted and the depth kept)
C          4 negative trial values (treated as zero depth, recession)
C          5 depth above the limit XMAX (error 17)
C-----------------------------------------------------------------------
      SUBROUTINE IMPLCT (NM,J)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /CNTRL/ nres, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1,CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /BALANC/ STORA(60),EINF(60),CINF2(20), STODEP(60)
      COMMON /INFIL/ AL(60), thi(60), thmx(60), ROC(60), RECS(60),
     1        THR, thfc, tempw, stone(60), RVAL
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2         SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /EROSN/ V1(20), V2(20), CT(20), CMX(20), RDET(20), DLA(20),
     1      DrlNET(20), Dirnet(20), suspir, susprl, SPOW(20)
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      common /intrill/ intrt, kint, d1(20), d2(20), vir(20), ql3(20), 
     1                 cbal, dcx, slrat, izr
      dimension y2(20), y1(20), qlt(20)
C
      EXTERNAL IMPOCF, IMPCHA
      DATA EPS, IEND,ch,cd /0.00001,50,0.,0./
C
C       DEFINE CHANNEL FUNCTIONS
      dqdact(AP,AAA,b1,b2,P,PM1,C)=
     -  AP*(AAA/b2)**PM1*(P-(PM1*AAA*C/(b1*b2)))
      ZFUNC(B,AREA,CO1)=(B*B+2.*AREA/CO1)**0.5
C       NEW FUNCTIONS FOR CHANNEL INFLITRATION
        PRFUNC(AAA,B,BW,CO1)=AMIN1( (-B+SQRT(B*B+2.*AAA/CO1))/
     -  (0.15*BW**0.5), 1.0 )
        WPFUNC(B,AAA,CO1,CO2)=(-B+SQRT(B*B+2.*AAA/CO1))*CO2 + B*CO1
C
C       NEW CODE FOR CHANNEL INFLITRATION
        BW = BWID
C
      nmM1=nm-1
C
C  W(J) = 0 identifies a channel (label 110); otherwise a plane.
      IF (W(J).EQ.0) GO TO 110
C====================================================
C         ***  PLANE CASE  ***
C -- SOLUTION FOR ADVANCE TIME DEPTH(H2). --
C====================================================
C Alternate use for interrill routing:
C  With interrill routing (KINT >= 1) work on the interrill arrays
C  (Y = D, QLT = QL3, step DCX, characteristic speed CD).  Otherwise
C  use the rill / plane arrays (Y = H, QLT = QL2, step DX, speed CH).
      Do 5 k=1,nm
        If(KINT .ge. 1) Then
          y1(k) = d1(k)
          y2(k) = d2(k)
          qlt(k) = ql3(k)
          dy = dcx
          csf = slrat
          cx = cd
        Else
          y1(k) = h1(k)
          y2(k) = h2(k)
          qlt(k) = ql2(k)
          dy = dx
          csf = 1.
          cx = ch
        End If
  5   Continue
      IF (Y1(nm).LE.0.) CX=0.
C           +++++++++++++++++++++++ Begin Loop ++++++++++
C  Loop down the plane, node K to K+1.  XST is the starting value of the
C  depth for the iteration.  QTEST is the volume of lateral inflow over
C  the step; if it is negative (recession) and larger than the depth
C  already there, the node is dry.
      DO 100 K=1,nmM1
      INDEX=K
      KP1=K+1
      XST = Y2(K)
      IF (XST.EQ.0.) XST=Y1(KP1)
      If(kint .le. 0) Then
        QTEST=0.5*(QLt(K)+QLt(KP1))*DT
      Else
        QTEST = ql3(nm)*dt
      End If
      IF (QTEST.GT.0.) GO TO 10
      xst = 0.5*(Y2(k) + Y1(kp1))
      QTEST=-QTEST
      IF (QTEST.LT.Y1(KP1)) GO TO 10
        IF (Y1(K) .GT. 0.0 .OR. Y2(K) .GT. 0.0) GO TO 10
        Y2(KP1)=0.
        GO TO 95
   10 continue
C ****characteristic method thru first node:
C  First node, flow not yet established (Y2(1) = 0): the upstream
C  characteristic is traced from the top in 5 s sub-steps (CX = its
C  distance down the slope).  Until it reaches the first node
C  spacing (DY) the depth at node 2 is built directly from the lateral
C  inflow (area A = A + DDT * QLT * width, depth from HFUN for
C  rills with side slopes).
      if(K .le. 1) then
        ZERO = 0.0
        IF (CX .lt. Dy .and. Y2(1) .le. ZERO) then
C trace upstream shock carefully, using 5 second subdivisions:
          ddt = 5.
          udt = 0.
          yu = y1(2)
          IER=1
C          if(yu .le. 0.) yu = qlt(2)*0.5*ddt
   15     Continue
          call perim(1,j,2,power(2),Yu,tw2,au,dum,dAVdh,vrl)
C          if(yu .le. 0.) Then
C            write(*,*) kint,yu,tw2,davdh
C            read(*,*)
C          End If
          cv = csf*alpha(2)*dAVdh/tw2
          CX=CX+CV*ddT
          udt = udt + ddt
C     Trace Characteristic:  cx is position of upstream char. at time t
          twr = t(L) - dt + udt
          If(cx .gt. dy) go to 40
          if(rillw1(2) .le. 0.) then
            y2(2)=yu + ddT*QLt(2)
          else
            If(kint .ge. 1) Then
              wd2 = 1.
              wdfl = 1.
            Else
              wd2 = winter(2) + rillw1(2)
              wdfl = tw2
            End If
            a2L = au
            a2t = A2L + ddT*QLT(2)*WD2
            if(zrl(j) .gt. 0.01 .and. kint .le. 0) then
              if(a2t .lt. 0.) a2t = 0.
              y2(2) = hfun(rillw1(2),zrl(j),a2t)                        !new
            else
              y2(2) = a2t/wdfl
            end if
          End If
          If(udt .lt. dt) Then
            Yu = amax1(Y2(2),0.)
            If((udt + ddt) .ge. dt) ddt = dt - udt + .01
            go to 15
          End if
          go to 40
        end if
      end if
C  Normal case: solve the implicit equation for the new depth Y2(K+1)
C  by Newton-Raphson (ITER with IMPOCF).  Convergence tolerance 1.E-5,
C  at most IEND = 50 passes, depth limit 3 m.
      CALL ITER (y2(KP1),FH2,DERFH2,IMPOCF,XST,.00001,IEND,IER,3.)
C diagnos: message indicating that neg depth found during recession
C      If(ier .eq. 3) Then
C        If(y2(kp1) .lt. 0.) Then
C          write(idiagn,39) kp1,y2(kp1)
C        Else
C          write(idiagn,38) kp1,y2(kp1),rilld1(kp1)
C        End If
C      End If
C  38  Format(' At I=',I2,', Depth ',f10.5,', Rilldepth ',f10.5)
C  39  Format(' At I=',i2,', depth ',f10.5)
C
C  Shift IER to the code used below (see the list in the header).
      IER=IER+1
   40 y2(KP1)=AMAX1(0.,y2(KP1))
C
C     INTERPRET ERROR FLAG FOR TOO MANY NEG. TRIAL VALUES AS CONVERGENCE
C     AT ZERO.  ONLY POSSIBLE DURING RECESSION PERIOD OF HYDROGRAPH.
C  Negative trial values during a recession mean the depth is
C  converging to zero: set it to zero (if no later node has water) and
C  go on.  A failure to converge (IER = 2) or a depth above the limit
C  (IER = 5) leaves through the error exits at labels 310 and 360.
      If(IER .eq. 4 .and. (QLT(K)+QLT(KP1)) .le. 0.) Then
        If(KP1 .lt. NM) Then
          DO 50 JK=KP1,NM
            IF(y2(JK).GT.0.) Then
              y2(KP1) = 0.
              GO TO 95
            End If
   50     CONTINUE
        Else
          y2(KP1) = 0.
          GO TO 95
        End If
      End If
C      GO TO (95,310,95,330,360), IER
      If(IER .eq. 2) Then
        go to 310
      Else If(IER .eq. 5) Then
        Go to 360
      Else If(IER .ne. 1 .and. IER .ne. 3 .and. IER .ne. 4) Then
         write(*,*)' ierzp'
        CALL GOTOER('IER_zP')
      End If
C
   95 Continue
C  Flow area, velocity and flow width from the new depth (PERIM).
C  Velocity = ALPHA * FV / A.
      Call PERIM(0,J,KP1,power(kp1),y2(KP1),twp,ah,FV,dum2,vrl)
      a2(KP1) = ah
      vnew = 0.
      coefq = csf*alpha(kp1)
      if(ah .gt. 0.) vnew = coefq*FV/ah
      If(KINT .le. 0) Then
        h2(KP1) = y2(KP1)
        V2(KP1) = vnew
        vr(kp1) = vrl*coefq
        widw(kp1) = twp
        If(k .le. 1) ch = cx
      Else
        d2(KP1) = y2(KP1)
        vir(kp1) = vnew
        if(k .le. 1) cd = cx
      End If
C
  100 CONTINUE
C          ++++++++++++++ End Loop +++++++++++++
      RETURN
C
C====================================================
C    ***  TRAPEZOIDAL CHANNEL CASE  ***
C -- SOLUTION FOR ADVANCE TIME AREA(A2). --
C====================================================
C-----------------------------------------------------------------------
C  CHANNEL CASE.  Trapezoidal channel: wetted perimeter, area-width
C  relations and the Manning wave speed DQDA = dQ/dA (statement
C  function DQDACT) are written as functions of the area.  XC is
C  the distance the front has travelled during the step while the
C  channel is still dry above node 1.
C-----------------------------------------------------------------------
  110 Continue
      IF (B.EQ.0.) B=.000001
      C=CO2/CO1
      XK=B*CO1
      P=POWER(1)
        ALP=ALPHA(1)
      PM1=P-1.
C
C       INITIALIZE XC
      IF(A2(nm) .LE. 0.0) XC = 0.0
C                       ************************** K=1,NMK1 LOOP *******
      DO 190 K=1,nmM1
        INDEX=K
        KP1=K+1
        IF(K .EQ. 1) THEN
C          THEN COMPUTE NEW STARTING VALUE FOR UPPER BOUNDARY
           DUMJP1 = ZFUNC(B,A1(KP1),CO1)
           DUMJ   = ZFUNC(B,A1(K),CO1)
           WPJP1  = (DUMJP1-B)*CO2 + B*CO1
           WPJ    = (DUMJ-B)*CO2 + B*CO1
           WPPOW1 = WPJP1**(P-1.)
           WPPOW  = WPJ**(P-1.)
           DADX = ( ALP*(A1(KP1)**P)/WPPOW1 - ALP*(A1(K)**P)/WPPOW )/DX
           XST = (QL(L-1)+QL(L))*DT - 2.0*DT*DADX - A2(K)+A1(KP1)+A1(K)
        ELSE
C          NOTE: A2(KP1) IS USED FOR IN THE FOLLOWING EQUATION BECAUSE IN
C              THE FIRST STEP IT = 0 AND AFTERWARD IT HOLDS THE VALUE FROM
C              THE PREVIOUS TIME STEP
           XST=(A2(K)+A2(KP1))/2.0
        ENDIF
C  Starting value XST for the area iteration: at node 1 from the
C  upstream boundary and the lateral inflow, elsewhere the average of
C  the old areas at nodes K and K+1.
        IF(XST .LE. 0.0) THEN
           XST = 0.0
        ENDIF
        IF(NOQL.EQ.1) THEN
C          THEN THERE IS A CONTRIBUTING PLANE SO SHOULD HAVE LATERAL INFLOW
           CONTINUE
        ELSE
C       ELSE IF NO LATERAL INFLOW
           IF (A1(KP1) .LE. 0.) XST=0.
        ENDIF
C       UNITS FOR NEW VARAIBLES
C  QLAT = mean lateral inflow over the step.  GAL = G * (THMX - THI) /
C  12, the capillary parameter used by the channel infiltration.
        QLAT =  (QL(L-1) + QL(L)) * 0.5
        IF(FMINEW(J) .GE. 1.0E-08) THEN
           GAL = AL(J) * (thmx(J)-thi(J))/12.0
        ENDIF
C  Node 1: infiltration at the upper boundary is computed first, then
C  the dry-bed start-up ("zone A": area from the lateral inflow
C  alone, iterated up to 20 times, label 116).  If the front has not yet
C  moved one node spacing, the zone A value is kept and the iteration
C  skipped.
        IF( K .EQ. 1) THEN
          IF(FMINEW(J) .GE. 1.0E-08) THEN
C             THEN COMPUTE UPPER BOUNDARY INFILTRATION CAP.
C             FIRST COMPUTE WETTED PERIMETERS
            DUML1 = ZFUNC(B,AUB(L-1),CO1)
            DUML  = ZFUNC(B,AUB(L),CO1)
            WPL1  = (DUML1-B)*CO2 + B*CO1
            WPL   = (DUML-B)*CO2 + B*CO1
C             COMPUTE INFIL CAP [L**2/T] AT THE UPPER BOUNDARY
            IF( (A1(K)+A2(K)+QLAT) .LE. 0.0) THEN
               FC2(K) = FC1(K)
               AOUT = 0.0
C  CHECK (not changed): in this branch and in the one at label 115 the
C  channel infiltration parameter is read as FMINEW(I), but the index I
C  is not set anywhere in this routine, so it holds whatever value it
C  has (undefined: in practice often zero, which is outside the array).
C    Other places in the same
C  routine use FMINEW(J), which looks like what was meant.  This matters
C  only for channel elements with infiltration; it does not affect the
C  test case, which has only planes.  Compiling with -fcheck=bounds
C  should flag it if a channel is run.
            ELSE
               CALL CHAINF( FC1(1),FC2(1),DELT,FMINEW(I),GAL )
               AOUT = ((FC2(1) - FC1(1))/DELT) * 0.5 * (WPL1 *
     -         PRFUNC(AUB(L-1),B,BW,CO1)+WPL*PRFUNC(AUB(L),B,BW,CO1))
            ENDIF
            CINF2(K) = AOUT
          ELSE
C             NO INFILTRATION
            CINF2(K) = 0.0
          ENDIF
  115     CONTINUE
C          ZONE A SOLUTION
          IF( (A1(K)+A2(K)) .LE. 0.0) THEN
            IF(FMINEW(I) .GE. 1.0E-08) THEN
              IF( (A1(K)+A2(K)+QLAT) .LE. 0.0) THEN
                FC2(K) = FC1(K)
                AOUT = 0.0
              ELSE
                CALL CHAINF( FC1(2),FC2(2),DELT,FMINEW(I),GAL )
                AOUT = (FC2(1) - FC1(1))/DT
              ENDIF
            ELSE
              AOUT = 0.0
            ENDIF
            WPL=WPFUNC(B,A1(2),CO1,CO2)*PRFUNC(A1(2),B,BW,CO1)
            AX = A1(2) + (QLAT*DT)
            DO 116 NCOUNT=1,20
              FAX = AX - A1(2) - DT*(QLAT - AOUT*(WPL +
     -              WPFUNC(B,AX,CO1,CO2)*PRFUNC(AX,B,BW,CO1))*0.5)
              IF( ABS(FAX) .LT. 0.00001 ) GO TO 117
              AX = AX - FAX
  116       CONTINUE
            If(NP(J) .eq. 7) Then
              WRITE(IDIAGN,600)
              WRITE(IWRITE,600)
            End If
  600    FORMAT(/,' ','MORE THAN 20 ITER. FROM ZONE A SOLUT.',/,
     -                 ' ','FOR TRAP. CHAN. IN SUB. IMPLCT ',/)
  117       A2(2) = AX
            DUM11=ZFUNC(B,A2(2),CO1)
            DUM21=(DUM11-B)*CO2+XK
            DQDA=dqdact(ALPHA(1),A2(2),DUM11,DUM21,P,PM1,C)
            XC = XC + DQDA*DT
            IF(XC .LT. DX .AND. A2(1) .LE. 0.0) THEN
C                THEN DO A ZONE A SOLUTION (I.E. SKIP ITER AND USE A2(2)
C                AS COMPUTED ABOVE
              IER = 1
              GO TO 125
            ENDIF
          ENDIF
        ENDIF
  120   CONTINUE
C  Infiltration loss for the reach from node K to K+1 (XINFLJ, kept
C  fixed for the iteration in IMPCHA), then the Newton solution for the
C  new area A2(K+1).
        IF(FMINEW(J) .GE. 1.0E-08) THEN
           IF( (A1(K)+A2(K)+QLAT) .LE. 0.0) THEN
                FC2(KP1) = FC1(KP1)
           ENDIF
C          COMPUTE XINFLJ SO IT IS NOT RECOMPUTED FOR ALL ITER. IN IMPCHA
           XINFLJ = (FC2(K) - FC1(K))/DT * (WPFUNC(B,A1(K),CO1,CO2)
     -            * PRFUNC(A1(K),B,BW,CO1) + WPFUNC(B,A2(K),CO1,CO2)
     -            * PRFUNC(A2(K),B,BW,CO1)) * 0.5
           CALL CHAINF( FC1(KP1),FC2(KP1),DELT,FMINEW(J),GAL )
        ENDIF
        CALL ITER (A2(KP1),FA2,DERFA2,IMPCHA,XST,EPS,IEND,IER,10000.)
        IER=IER+1
  125   CONTINUE
C       SAVE THE INFIL. AT EACH NODE (UNITS OF AOUT ARE [L**2/T])
C       CHECK FOR NEG. A2(KP1)
C  Infiltration at each node, CINF2 (area units per time): the change
C  in cumulative infiltration times the wetted perimeter.
        IF(A2(KP1) .GT. 0.0) THEN
           CINF2(KP1) = (FC2(KP1) - FC1(KP1))/DT *
     -     ((WPFUNC(B,A1(KP1),CO1,CO2) * PRFUNC(A1(KP1),B,BW,CO1) +
     -     WPFUNC(B,A2(KP1),CO1,CO2)*PRFUNC(A2(KP1),B,BW,CO1))*0.5)
        ELSE
           CINF2(KP1) = (FC2(KP1) - FC1(KP1))/DT *
     -     (WPFUNC(B,A1(KP1),CO1,CO2) * PRFUNC(A1(KP1),B,BW,CO1)
     -     *0.5)
        ENDIF
C  Negative or zero area at the node: decide whether it is the
C  pre-runoff period or a recession, and set the remaining nodes to
C  zero or continue.  At the end of the node loop the cumulative
C  infiltration FC2 becomes FC1 for the next step.
        IF (IER .EQ. 4) GO TO 130
        IF (A2(KP1).EQ.0.) GO TO 130
        IF (A2(KP1).LT.EPS) A2(KP1)=0.
        GO TO 170
C
C     CHECK TO SEE IF NEG. VALUE IS DUE TO PRE- RUNOFF PERIOD OR TO
C     RECESSION
C  130   KP1=K+1
  130  Continue
C        If(IER .eq. 4) Then
C          write(*,"('  kp1,a2m,a2: ',i3,3g13.5)")kp1,a2(K),a2(KP1)
C          write(IDIAGN,"('  kp1,a2m,a2: ',i3,3g13.5)")kp1,a2(K),a2(KP1)
C        read(*,*)
C        End If
        DO 140 JK=KP1,NM
          IF (A1(JK).GT.0.) GO TO 160
  140   CONTINUE
        IF (NOQL.LE.1.AND.FMINEW(JJCHAN).LE.0.0) GO TO 170
        IF(QLAT .LE. 0.0) THEN
          DO 150 KJ=KP1,NM
            A2(KJ)=0.0
              FC2(KJ) = FC1(KJ)
              CINF2(KJ) = 0.0
  150     CONTINUE
        ELSE
C********** SET NEG. A2(KP1) TO 0.0 IN THE CASE : QLAT >0.0
           IF (A2(KP1) .LT. 0.0) A2(KP1) = 0.0
           GO TO 190
        ENDIF
C  ************RESET CUMM. INFILTRATION ARRAY and RETURN
        DO 155 KK=1,nm
           FC1(KK) = FC2(KK)
  155   CONTINUE
        RETURN
C                   
  160   A2(KP1)=0.
        GO TO 190
C  170   GO TO (190,310,320,340,360), IER
  170 Continue
        If(IER .eq. 4) Then
C  Moved into Loop:
          IF ((A1(K) .NE. A1(K+1)) .or. (A2(K) .LE. 
     1      ((QL(L-1)+QL(L))*.5*DT))) GO TO 350
C
  180     A2(K+1)=(QL(L)+QL(L-1))*.5*DT
          WRITE (IDIAGN,400) K
        Else If(IER .eq. 2) Then
          Go To 310
        Else If(IER .eq. 3) Then
          Go TO 320
        Else If(IER .eq. 5) Then
          Go To 360
        Else If(IER .ne. 1) Then
          write(*,*) ier,' ierzch'
          CALL GOTOER('IER_CH')
        End If
  190 CONTINUE
C **************End K Loop ******************
C       RESET CUMM. INFILTRATION ARRAY
        DO 195 KK=1,nm
           FC1(KK) = FC2(KK)
  195   CONTINUE
C            ************************************************** END LOOP
      RETURN
C
C  Error exits, from ITER failures:  10 no convergence, 11 zero
C  derivative, 17 depth above the limit.  ERRPO writes the message and
C  stops the run.
  310 CALL ERRPO('IMPLCT',10,IEND,0,'  X   ','  X   ')
  320 IP1=INDEX+1
      CALL ERRPO('IMPLCT',11,IP1,0,'  X   ','  X   ')
  350 CALL ERRPO('IMPLCT',13,IP1,0,'  X   ','  X   ')
      RETURN
  360 IP1=INDEX+1
      CALL ERRPO('IMPLCT',17,I,IP1,'  X   ','  X   ')
  370 IP1=INDEX+1
      CALL ERRPO('IMPLCT',16,I,IP1,'  X   ','  X   ')
      RETURN
C
  380 FORMAT (2X,'  NO POSITIVE ROOT AT J= ',I3,/,10X,'H2(K+1)=H1(K+1)+(
     1QL)(DT)')
  390 FORMAT (2X,' Zero DERF or zero depth AT I= ',I2,', at ',F10.2,' Se
     1c. ON PLANE ',I2)
  400 FORMAT (2X,' NO POSITIVE ROOT AT J=',I3,/,10X,'A2(K+1)=(QL(1)+QL(2
     1))(.5)(DT)')
      END
C***************************************************************************
C-----------------------------------------------------------------------
C  SUBROUTINE IMPOCF(X,FX,DERFX)
C    Residual FX and derivative DERFX of the implicit plane equation,
C    as a function of the unknown depth X at node I+1 (I = INDEX, set
C    by IMPLCT):
C        FX = (A2_i + A(X) - A1_i - A1_i+1)
C             + (2 DT / DY) * (THETA * (ALPHA_i+1 * FV(X)
C                                        - ALPHA_i * FV2_i)
C                              + (1 - THETA) * (...old time level...))
C             - PC
C    where A is the flow area from PERIM, 1 = old and 2 = new time
C    level, and PC = DT * (QL_i + QL_i+1) * flow width is the lateral
C    inflow, already averaged over the step.
C    DERFX = TW + (2 DT / DY) * ALPHA * THETA * dFV/dh.
C    Modified 1/91 to use hydraulic radii as variable and the wetted
C    perimeter WP(h).  With interrill routing (KINT >= 1) the unit-width
C    form with DCX, D1, D2, QL3 is used.
C-----------------------------------------------------------------------
      SUBROUTINE IMPOCF (X,FX,DERFX)
C
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      common /intrill/ intrt, kint, d1(20), d2(20), vir(20), ql3(20), 
     1                 cbal, dcx, slrat, izr
C
C Modified 1/91 to use hydraulic radii, as variable, and use
C  WP(h) as wetted perimeter.  X is still h
C
      I=INDEX
      J = JU
      IP1=I+1
      P=POWER(I)
      PP=POWER(I+1)
C      PPM1=PP-1.
      TH=THETA
      T1=1.-THETA
      If(KINT .ge. 1) Then
        DY = DCX
        HLI = D1(I)
        HLP = D1(IP1)
        HNI = D2(I)
        WD = 1.
        PC = dt*(ql3(i)+ql3(ip1))
        fs = slrat
      ELSE
        DY = DX
        HLI = H1(I)
        HLP = H1(IP1)
        HNI = H2(I)
        wd = rillw1(ip1) + winter(ip1)
        PC = DT*(QL2(I)+QL2(IP1))*wd
        fs = 1.0
      END IF
      Call Perim(0,J,I,P,HLI,tw1,AL1,FVL1,dum2,vrl)
      Call Perim(0,J,I,P,HNI,tw2,aL2,FVL2,dum2,vrl)
      Call Perim(0,J,IP1,PP,HLP,twp,AP1,FVP1,dum2,vrl)
      Call Perim(1,J,IP1,PP,X,twx,APX,FVAX,dFVdh,vrl)
      PA = AL2 + APX - AL1 - AP1
      PB = (2.*DT/DY)*(TH*(ALPHA(IP1)*FVAX - ALPHA(I)*FVL2)
     &    +  T1*(ALPHA(IP1)*FVP1 - ALPHA(I)*FVL1))*fs
C* LATERAL FLOWS ARE ALREADY AVERAGED OVER TIME:
      FX=PA+PB-PC
C
      DERFX = TWX + 2.*DT/Dy*fs*ALPHA(IP1)*TH*(DFVDH)
      RETURN
      END
C****************************************************************************
C-----------------------------------------------------------------------
C  SUBROUTINE ITER(X,F,DERF,FCT,XST,EPS,IEND,IER,XMAX)
C    Newton-Raphson solution of F(X) = 0, where FCT(X,F,DERF) returns
C    the function and its derivative (here IMPOCF or IMPCHA).
C      XST   starting value        IEND  maximum number of passes
C      XMAX  upper limit for X     IER   result flag (before the
C                                        caller adds 1):
C            0 converged   1 no convergence in IEND passes
C            2 zero derivative   3 negative-X limit or change of sign
C            4 X above XMAX five times
C    Features:
C      - a negative X is replaced by a fraction of the old value
C        (up to 5 times)
C      - after 15 passes, if the function has changed sign since the
C        start, BISECT (in Euroa.for) finds the root within 1.E-6
C      - normal exit when |F| < 1.E-8, or when both the relative change
C        in X and the relative F/X are below 0.0005
C    NOTE: the argument EPS is set to TOL but never used: the 0.0005
C    tolerance is written into the code.
C    NOTE: writes to the screen (unit *) when it takes more than 15
C    passes, and to unit 8 (written as a literal) for a zero derivative.
C-----------------------------------------------------------------------
       SUBROUTINE ITER (X,F,DERF,FCT,XST,EPS,IEND,IER,XMAX)
C
C      DCG 11/29/88 - MODIFIED TO CALL BISECTION AFTER 15 NEWTON-
C                     RAPHSON ITERATION AND WHEN A CROSS
C                     OVER OF THE FUNCTION (NEG. TO POS. OR VICE VERSA)
C                     IS ENCOUNTERED
C      DCG 1/25/89  - MODIFY TO EXIT WHEN F .LT. 1.0E-08
C                     AND PREVENT DIVIDE BY ZERO FOR RELATIVE ERROR
C                     CHECK
C
C--------------------------------------------------------------------------
C     PREPARE ITERATION
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /BALANC/ STORA(60),EINF(60),CINF2(20), STODEP(60)
C
      EXTERNAL FCT
      DATA BISTOP /1.0E-06/
C
      IER=0
      NC=0
      LC=0
      X=XST
      TOL=X
      CALL FCT (TOL,F,DERF)
C       KEEP ORIGINAL FUNCTION AND START VALUE TO CHECK ZERO CROSS OVER
        FOLD = F
        XOLD = XST
      DX=F/DERF
      X=X-DX
      NSIGN=0
      IF (DERF.LT.0.) NSIGN=1
C     START ITERATION LOOP
      DO 100 I=1,IEND
      rlxf = 1.                     ! -0.5*float(i-1)/float(iend)
      IF (F) 10,120,10
   10 IF (DERF) 20,130,20
C     IF X TAKES A NEGATIVE VALUE, CORRECT X TO BE HALF ITS OLD VALUE.
C     IF NOT, MAKE SURE NC=0 AND CONTINUE.
   20 IF (X) 30,40,40
   30 NC=NC+1
        IF ((NC-6) .LT. 0) THEN
           X=(X+DX)/(1.+FLOAT(NC))
           GO TO 70
        ELSE
          If(X .lt. -0.001) Then
            Write(idiagn,652)
  652 Format  (' Iteration limit in subroutine ITER reached, call for 
     1assistance')
          End If
          GO TO 150
        ENDIF
C**     IF (NC-5) 70,150,150
   40 NC=0
      IF (X-XMAX) 60,60,50
   50 LC=LC+1
      X=0.9*XMAX
      IF (LC-5) 70,70,160
   60 LC=0
   70 TOL=X
      CALL FCT (TOL,F,DERF)
      NCK=0
      IF (DERF.LT.0.) NCK=1
      IF (NSIGN-NCK.NE.0) GO TO 140
      DX=F/DERF
      X=X-DX*rlxf
      TOL=EPS
C       FOR VERY SMALL FUNCTION VALUE - NORMAL EXIT
        IF(ABS(F) .LT. 1.0E-08) GO TO 120
       IF (I.LT.15) GO TO 80
C       CHECK FOR ZERO CROSSING AND IF IT EXISTS CALL BISECTION
        IF( (FOLD*F) .LE. 0.) THEN
           IF(XOLD .LT. X) THEN
              CALL BISECT(XOLD,X,FCT,BISTOP,XF,FF)
           ELSE
              CALL BISECT(X,XOLD,FCT,BISTOP,XF,FF)
           ENDIF
C          SET OPT VALUE AND EXIT
           X = XF
           F = FF
           GO TO 120
        ENDIF
C
       WRITE (*,170) I,F,X,derf
C       USE STOPPING CRITERIA OF WALRUS (RELATIVE) REDUCE BY FACTOR OF 5
   80   CONTINUE
C       if(DERF .LT. 0.) WRITE(6,170) I,F,X,DERF
C       PREVENT DIVIDE BY ZERO
        IF(X .LE. 1.0E-15) GO TO 90
        ERRH = ABS(DX/X)
        ERRF = ABS(F/X)
        IF(ERRH .LT. 0.0005 .AND. ERRF .LT. 0.0005) GO TO 120
   90   CONTINUE
C       RESET XOLD AND FOLD
        XOLD = X
        FOLD = F
C       AFTER 49 ITERATIONS AND VERY SMALL AREA THEN EXIT NORMALLY
        IF(I .GT. 48 .AND. X .LE. 1.0E-12) GO TO 120
C
  100 CONTINUE
C     END OF ITERATION LOOP
      GO TO 110
  110 IER=1
  120 RETURN
  130 IER=2
      Write(8,131) derf, f, x, dx
  131 Format(' derf,f,x,dx:',4g11.4)
      RETURN
  140 WRITE (IDIAGN,180)
Czzzz
      write(*,189)x,tol,xst,f,derf
  189 Format(' x,tol,xst,f,derf',5g12.5)
Czzzz
  150 IER=3
      RETURN
C ***  FLAGGED RETURN AFTER CONVERGENCE TO GT. MAX VALUE OF X
  160 IER=4
      RETURN
C
  170 FORMAT (' I',i2,'; F(x) ',G12.5,';  x ',G12.5,' dF',g12.5)
  180 FORMAT (11X,'ITER/ERROR: Change of sign in the objective function
     1implies a multiple solution -  check geometry')
      END
C***************************************************************************
C-----------------------------------------------------------------------
C  SUBROUTINE PERIM(IFDER,J,I,PN,H,WID,AH,FV,DFVDH,VRL)
C    Flow geometry at node I of element J for a flow depth H.  Returns
C      AH     flow area, summed over the rill and any overflow
C      WID    total flow width at the surface
C      FV     discharge function, so that Q = ALPHA * FV, with
C             ALPHA = sqrt(S)/n and FV = A * R**(PN-1) = WP * R**PN
C             (PN = 5/3 for Manning)
C      DFVDH  dFV/dh (only if IFDER >= 1)
C      VRL    velocity in the rill, FV / AH
C    Three cases:
C      1  interrill routing (KINT >= 1): unit width, FV = H**PN
C      2  no rill (depth RILLD1 <= 1.E-5): a flat surface of width
C         WINTER, FV = WINTER * H**PN
C      3  rill with bottom width BW = RILLW1 and side slope Z = ZRL(J):
C           H <= HC (rill depth)  trapezoid of top width BW + 2*Z*H
C           H > HC                the rill is full and the water spreads
C             onto the interrill surface on each side.  The width of the
C             overflow, (H - HC) / XS with XS = SIR(J) (the interrill
C             slope across the rill), is limited to half of what is left
C             of the rill spacing W / DEPNO after the rill.
C    NOTE: Z is used as the horizontal run per unit depth (top width =
C    BW + 2*Z*H); the "h/v" wording of the original header is the
C    reciprocal of that.
C-----------------------------------------------------------------------
      subroutine PERIM(IFDER,J,I,PN,H,wid,AH,FV,dFVdh,vrl)
C
C  Finds effective discharge function (FV as function of
C  flow depth, and derivative when flag IFDER is positive.  Uses two
C  rill models, one for vertical sides, and another for furrow type rills
C  Rill side slope is defined as tangent h/v = ZRL(j)
C  PN is hydraulic exponent
C  HC is depth of rill, rilld1(i)
C  WIDC is width at top of rill
C  WOVF is width of overflow, if H > HC
C  WID is total flow width, = 2.*WOVF + WIDC for overflow case
C  hwov = width of overflow on one side
C  vrl = rill velocity Q/A
C
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), DERO(60), SIR(60)
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      common /intrill/ intrt, kint, d1(20), d2(20), vir(20), ql3(20), 
     1                 cbal, dcx, slrat, izr
C
       xs = sir(j)
C       data xs/.04/
      If(KINT .ge. 1) then
        ah = h
        wid = 1.
        FV = ah**PN
        dFVdh = PN*ah**(PN-1.)
        RETURN
      End If
      HC = RILLD1(I)
      If(HC .le. 1.e-5) then
C **      No Rill depth
        wid = winter(i)
        ah = h*wid
        FV = wid*h**PN
        if(h .gt. 0.) vrl = FV/ah
        dFVdh = wid*PN*h**(PN-1.)
      Else 
        BW = RILLW1(I)
        z = zrl(j)
        cf = SQRT(1.+ z*z)
        If(H .le. HC) Then
C ** simple trapezoid case
          wid = BW + 2.*z*H
          ah = H*(wid + BW)/2.
          WP = BW + H*2.*cf
          rh = ah/WP
          FV =  WP*rh**PN
          vrl = 0.
          If(ah .gt. 0.) vrl = FV/ah
          drdh = (wp*wid - ah*2.*cf)/wp/wp
          DFVDH = (PN*WP*rh**(pn-1.) + 2.*cf*rh**pn)*drdh
        Else
          WC = BW + 2.*z*HC
          hwmax = 0.5*(w(j)/depno(j)-wc)
          WPC = BW + H*2.*cf
          AC = HC*(BW+WC)/2.
          A1 = AC + (H-HC)*WC
          rh1 = A1/WPC
          fv1 = WPC*rh1**PN
          vrl = fv1/A1
          hwov = (H-HC)/xs
          If(hwov .le. hwmax) Then
            a2 = hwov*(H-HC)
            rh2 = (H-HC)/2.
            If(IFDER.ge.1) dfv2 = (PN+1.)*(H-HC)**PN/xs/(2.**(PN-1.))
          Else
            hwov = hwmax
            HL = hwmax*xs
            a2 = hwmax*0.5*(H-HC+ H-HL)
            rh2 = a2/hwmax
            If(IFDER.ge.1)dfv2 = 2.*hwov*PN*(rh2)**(PN-1.)
          End If
          fv2 = 2.*hwov*rh2**pn
          AH = A1+A2
          WID = WC + 2.*hwov
          FV = fv1 + fv2
          If(IFDER .ge. 1) Then
            dfv1 = pn*wpc*rh1**(PN-1.)
            DFVDH = dfv1 + dfv2
          End If
        End If
      End if
      RETURN
      END
C-------------------------------------------------
C-----------------------------------------------------------------------
C  FUNCTION QCP(FLUX,APAR,FKO)
C    Cumulative infiltration (the profile storage) at which ponding
C    begins under a steady supply FLUX larger than Ks = FKO:
C          QCP = APAR * ln( FLUX / (FLUX - FKO) )
C    A flux not above Ks never ponds: returns zero with a message.
C-----------------------------------------------------------------------
      FUNCTION QCP(FLUX,APAR,FKO)
C
C  FINDS PROFILE STORAGE AS A FUNCTION OF FLUX
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      IF (FLUX.LE.FKO) THEN
        QCP=0.
        WRITE (IDIAGN,10)
      ELSE
C     FLUX IS GREATER THAN FMIN:
        QCP=APAR*ALOG(FLUX/(FLUX-FKO))
      ENDIF
      RETURN
C
   10 FORMAT (' NO DEPTH ASSOCIATED WITH FLUX .LT. FMIN')
      END
C------------------------------------------------------------
C-----------------------------------------------------------------------
C  FUNCTION HFUN(BW,Z,AR)
C    Flow depth of a trapezoid with bottom width BW and side slope Z
C    (horizontal run per unit depth) for a flow area AR:
C          h = ( sqrt(BW**2 + 4*Z*AR) - BW ) / (2*Z)
C    Zero if the area is not positive.
C-----------------------------------------------------------------------
      Function hfun(bw,z,ar)
C  finds depth of trapezoidal area given bottom width, bw,
C    side slope, z, and area a
      If(ar .gt. 0.) Then
        hfun = (sqrt(bw*bw + 4.*z*ar)-bw)/2./z
      Else
        hfun = 0.
      End If
      return
      end

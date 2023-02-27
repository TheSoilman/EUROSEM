C*************************************************************************
      SUBROUTINE PLANE (J)
C
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      COMMON /CNTRL/ NRES, NTIME, NELE, CLEN, DELT,
     1       NLOG(60), NEROS, NPN(60)
      COMMON /GEOM/ XL(60), W(60), S(60), R1(60), R2(60), FMIN(60),
     1       NL(60), NR(60), NU(60), NC1(60), ZL(60), ZR(60),
     2       A(60), NC2(60), NP(60), DINTR(60), NRP(60),
     3       SUMA(60), RILLW(60), RILLD(60), ASR(60), RFR(60), DEPNO(60)
     4      ,RS(60), dero(60), sir(60)
      COMMON /EVENT/ TFIN, ND, QI(100), TI(100), NO, WTRAIN(60), 
     1        QIDD(20,100), TIDD(20,100), cumd(20,100), MAXND, 
     2        IGAGES(60), NGAGES, RNDPTH(60), RNRMX
      COMMON /PLANE1/ H1(20), H2(20), QL(1000), ALPHA(20), POWER(20),
     1    T(1000), Q(1000), HUB(1000), DX, DT, INDEX, THETA, XNU, GRAV,
     2    NB(60), QS(4000), LENQ, LENQS, LENH1, L, QL2(20), JU
      COMMON /CHAN/ A1(20),A2(20), QUB(1000),AUB(1000), CO1,CO2, B, NI,
     1       NOQL, TW(20), FMINEW(60), JJCHAN, FC1(20), FC2(20), XINFLJ
      COMMON /INFIL/ AL(60), SI(60), SMAX(60), ROC(60), RECS(60),
     1        THR, THFC, tempw, stone(60), RVAL
      COMMON /SED/ D50(60), RHOS(60), PORPL(60), PORCH(60), PAVE(60),
     2        SCON(4000),CONC(1000), CUB(1000), BWID, CONL(1000), VS,
     3         SPLTEX(60), ERODGJ(60),ITEX(60), COHE(60)
      COMMON /BALANC/ STORA(60),EINF(60),CINF2(20), STODEP(60)
      COMMON /RILL1/ RILLW1(20),rilld1(20),winter(20),ZRL(60),RWR(20),
     1        rillh(20),winh(20),winq(20),widw(20),vr(20),mcode,nk,
     2        arrild(20),arrilw(20)
      COMMON /ESEM/ COVER(60), CUINT(60), ISHAPE(60), PLHGT(60),
     1         VPANG(60),PBASE(60)
      COMMON /EROSN/ V1(20), V2(20), CT(20), CMX(20), RDET(20), DLA(20),
     1      DrlNET(20), Dirnet(20), suspir, susprl, SPOW(20)
      COMMON /EURO/ DINTRT(100), STEMD(100), DRIPD(100), TFALL(100),
     1        RAIN(100), RHRMM(100) 
      COMMON /KINETI/DRIPKE(100),RAINKE(100)
C
      common /intrill/intrt,kint,d1(20),d2(20),vir(20),ql3(20),cbal,
     1        dcx,slrat,izr
      Common /pldat/ ifad,ndat,tdat(500),qdat(500),sdat(500),sumrmm,
     1                sedtot
C
      DIMENSION XIN(50), RF(1000), ar1(20), ar2(20), qmph(20), qmmph(20)
     1         , smda(20), arill(20), qinrl(20)
C
      Character*1 line(4)   !, answer
      CHARACTER*5 TYPE
      logical ifad, rills
      Real*4 inrlbal
C
       DATA  NTELE /60/
      line(1) = char(179)
      line(2) = char(47)
      line(3) = char(196)
      line(4) = char(92)
C
C              WRITE OUT HEADER
      WRITE(IDIAGN,800) J
C
C              CALCULATE STORAGE DEPTH BASED ON DOWNSLOPE ROUGHNESS
C              STORAGE IS IN (METERS)
C
      STODEP(J) = 0.001*exp(-6.66+0.27*rfr(j))                          !6/93
C
C              NEED TO SAVE THIS FOR LATER USE IN VOLUME BALANCE CALCS
C
      stdpm=0.
      STOCAP = STODEP(J)
      write(idiagn,'(12h stocap(m): ,f10.6)') stocap
C
C              RECALCULATE KSAT USING HOLTAN'S MODIFICATION
C
      If (pbase(j).gt.0.0) then
C          write(*,701) J, fmin(j), FMIN(J)*(1./(1.-pbase(J)))
C          read(*,702) answer
C          if (answer .eq.'y'. or.answer.eq. 'Y') then
              FMIN(J)=FMIN(J)*(1./(1.-pbase(J)))
C          end if 
      end  if
 701  Format(/'  For element ',I2,', FMIN (mm/hr) will be changed'/
     &'  from ',F7.3,' to ',F7.3,', due to the value of PBASE entered. '
     &/ '  Type Y to continue with change, or N for no change: ')
 702  Format (A1)
C
      JU = J
      SUM1=0.
      SUM2=0.
      srain=0.0
      SUMQIN=0.0
      sedin = 0.
      sqout=0.0
      sinfil=0.0
        sumdir = 0.
        sumqirr = 0.
      SEDTOT=0.
C             CALCULATE NUMBER OF DISTANCE INCREMENTS (NK) 
C             AND THE DISTANCE INCREMENT SIZE (DX).
C
      ndx = 0
      NK = MAX(INT(15.*XL(J)/CLEN),5)
      IF (NK .gE. LENH1) Then
        WRITE (IDIAGN,470) NK
        NK=15
      End If
      DX=XL(J)/(FLOAT(NK)-1.)
      Write(IDIAGN,770)NK, dx
C
C     INITIALIZE
      IGAGED=IGAGES(J)
        IF(IGAGED .EQ. 0) THEN
C          THEN WRITE WARNING AND GAGE FOR PLANE ASSIGNED TO GAGE 1
           IGAGED = 1
           WRITE(*,15) J
           WRITE(IDIAGN,15) J
   15      FORMAT(//,6X,' *** WARNING ***',/,
     -     ' A RAINGAGE FOR ELEMENT ',I3,' WAS NOT ASSIGNED IN THE ',
     -     'INPUT FILE',/,' THE RAINGAGE ASSOC. WITH THE ELEMENT',
     -     'HAS BEEN ASSIGNED TO RAINGAGE 1',//)
        ENDIF
      CLAT=0.0
C
C       COMPUTE RELATIVE RESIDUAL SOIL MOISTURE (SR) AND RELATIVE
C       "FIELD CAPACITY" SOIL MOISTURE (THFC) BASED ON REGRESSION
C       RELATIONSHIPS AS DERIVED FROM SOIL CLASS MEAN VALUES FROM
C       RAWLS AND BRAKENSIEK, THESE ARE USED IN XPLINF ONLY TO ESTI-
C       MATE REDISTRIBUTION.
C          SINCE THE REGRESSION EQUATIONS ARE BASED ON
C       THE TEXTURAL SAT. HYD. COND. THEREFORE OBTAIN THE TEXTURAL
C       BASED SAT. HYD. COND. AS FOLLOWS
C
C                                   Ks(input)
C                     Ks(text.) = ------------
C                                 2.0 * (1-Vr)
C
C       NOTE: 2.0 IS A CORRECTION FOR AIR ENTRAPMENT
C             (1-Vr) IS A CORRECTION FOR VOLUME OF ROCK
C             FMIN, FMTEXT ARE IN (mm/hr)
C
        IF (FMIN(J).GE. 1.0E-08) THEN
           FMTEXT = FMIN(J) / ( 2.0 * (1. - ROC(J)) )
           THR = porpl(j)*(0.08122 - (0.05458 * ALOG10(FMTEXT)))
           THFC = porpl(j)*(0.45157 - (0.25361 * ALOG10(FMTEXT)))
           FMIN(J)=FMIN(J)/60.
        ENDIF
C      
        Call INTERC(J)
C      DO 310 JKL =1,NELE
        IGE = IGAGES(J)
C Inserted
C        DCG-2/25/88 SUBTRACT INTERCEPTION SO GLOBAL VOLUME BALANCE IS CORRECT
C        FOR TFIN LESS THAN RAINFALL TIME WE MUST INTERPOLATE QIDD TO GET
C        THE PROPER RAINFALL DEPTH FOR THE GLOBAL WATER BAL.
        QRAT = cumd(IGE,ND)
        DINADJ = 0.0
C        IF(IFLAGT .EQ. 0) THEN
        If(tfin .lt. tidd(ige,nd)) then
C               THEN TIDD(IGE,ND) IS .GE. TFIN SO CHECK INTERPOLATION
          DINADJ = cuint(J)
          DO 19 II=2,ND
            TDEL = TIDD(IGE,II) - TIDD(IGE,II-1)
            QDEL = CUMD(IGE,II)-CUMD(IGE,II-1)
            IF(TFIN.GE.TIDD(IGE,II-1).AND.TFIN.LT.TIDD(IGE,II)) THEN
C                 THEN INTERPOLATE PROPER RAINFALL DEPTH
              TDIFF = TFIN - TIDD(IGE,II-1)
              TRAT = TDIFF/TDEL
              QRAT=CUMD(IGE,II-1)+TRAT*(QDEL)
              IF(CUINT(J) .GE. 1.0E-07) THEN
                DINADJ = DINADJ - (DINADJ * TRAT)
C              ELSE
C                GO TO 305
              ENDIF
            ELSE
              IF(CUINT(J) .GE. 1.0E-07) THEN
                QTEM = WTRAIN(J)*CUMD(IGE,II+1)
                IF(QTEM .GE. DINADJ) THEN
                  DINADJ = 0.0
                ELSE
                  DINADJ = DINADJ - QTEM
                ENDIF
              ENDIF
            ENDIF
   19     CONTINUE
        ENDIF
C  rain depth on element J in m
        RNDPTH(J) = (QRAT*WTRAIN(J) - (CUINT(J) - DINADJ))/1000.
C  310 CONTINUE
C
         DINTR(J) = CUINT(J)
        CALL KINET(J)
        IF(NRP(J).EQ.1 .OR. NRP(J).EQ.3) THEN
            TYPE = 'MM/HR'
            WRITE(IDIAGN,730) J,TYPE
        END IF
        DO 20 K=1,ND
            QI(K)=WTRAIN(J)*QIDD(IGAGED,K)
            IF(DINTRT(k).LT.1.0E-07) GO TO 25
            DELTAT = TIDD(IGAGED,K+1) - TIDD(IGAGED,K)
            IF(QI(K)*DELTAT .GE. DINTRT(K)) THEN
               QI(K) = QI(K) - DINTRT(K)/DELTAT
            ELSE
               QI(K) = 0.0
            END IF
   25       CONTINUE
            IF(NRP(J).EQ.1 .OR. NRP(J).EQ.3) THEN
                TIME = TIDD(IGAGED,K)/60.
                QIT = QI(K)*3.6e6
                WRITE(IDIAGN,740) TIME,QIT,rainke(k),dripke(k)
            END IF
C
   20    CONTINUE
C                       
        TYPE = '(MM) '
        WRITE (IDIAGN,480) J,IGAGED,WTRAIN(J),DINTR(J),TYPE
C
      rills = .false.
      If(depno(j) .gt. 0. .and. rillw(j) .gt. 0.) Then
        rills = .true.
C  wpr is width of flow per rill .     depno is no of rills total
        wpr = w(j)/depno(j)                                               N11/91
C rillw() is width of each 'rill', not total         : 1/93
        rillwj = rillw(j)
C uwinj is interrill width between rill (bottom) edges:
        uwinj = wpr - rillwj
        If(rilld(j) .le. 0. .and. zrl(j) .le. 0.) zrl(j) = 1.
      Else
C  No rills in simulation:
        LCH = 0
        wpr = w(j)
        rillwj = 0.
        R1(J) = R2(J)
        s(j) = sir(j)
        uwinj =  w(j)
C set to 1 for use as multiplier later on:
        depno(j) = 1.
      End If
      rilldj = 0.
      if(rilld(j) .le. 0.) then   !  calculate rill depths from ASR:
C        if(rillw(j) .gt. 0.) rilldj = ((1.-asr(j))/2.)*w(j)/depno(j)
  !  This assumes the depressions have zrl(j) side slopes
        if(rillw(j) .gt. 0.) then
           asu = asr(j)/2.
           fac = sqrt(zrl(j)**2 + 1.) - zrl(j)
           rilldj = (asu/fac)*w(j)/depno(j)
        end if
      else
        rilldj = rilld(j)
      end if
C
       POWERJ = 5./3.
       ALPHAJ = SQRT(S(J))/R1(J)
       If(alphaj .le. 0. .or. powerj .le. 0.) stop 'zero alpha or power'
C
C  Max potential runoff in m^2/sec:   Check for dt too large:
      RXM = (RNRMX - FMINEW(J))*(XL(J)-DX)
      If(RXM .gt. 0.) Then
        hmt = (RXM/ALPHAj)**(1./powerj)
        vel = RXM/hmt
        deltl = DX/vel*2.
        If(deltl .lt. delt) Then
          Write(*,891) delt/60.,deltl/60.
 891  Format('  Specified DT of',f6.2,' min. is larger than'/
     &' Courant Limit of ',F6.2,':  Results may be seriously'/
     &' affected for large runoff coefficients'/)
          Write(*,*)' RETURN to continue:'
          Read(*,*)
        End If
C  No potential for runoff:
      Else
        Write(*,916) J
  916 Format(' Rainfall rates all less than minimum infiltration rate.'/
     &'  There will be no runoff on element ',I2 )!,'    RETURN to continue
C     &: ' ) 
C        Read(*,*)
      End If
C
      INTRT = 0
      LCH = 0
      contrib = 1.
C Check to see if interrill routing is indicated:
      If(rillwj .gt. 0.) Then
        slrat = R1(j)/R2(j)*SQRT(sir(j)/s(j))
C        ratl = 2.
        ratl = 1./(sin(atan(sir(j)/s(j))))
        ofl = uwinj/2.*ratl
		write(99,890) ratl, slrat, ofl
  890   Format(' Interrill Flow Length Ratio ',f6.2/
     &         ' Interrill Velocity Ratio    ',F6.4/
     &         ' Overland flow length (m)    ',F6.2/)
 !       if(uwinj .ge. 3. ) Then
        If(ofl .gt. 1.5) Then
		   	
C  Interrill flow long enough for routing:
          nkir = int(sqrt(4.* ofl)) + 1      !  Revised 1/98
          if(nkir .lt. 5) nkir = 5
          if(nkir .gt. 19) nkir = 19
          nkirp = nkir + 1
          dcx = ofl/(nkir-1)
          contrib = 1.-rillwj/wpr
          INTRT = 1
          Write(IDIAGN,771) nkir,ofl
        Else
C  No Interrill Routing
          cbal = 0.
          nkir = 0
          Write(idiagn,772) ofl
        End If
        LCH = 1
C        
      Else
        Write(IDIAGN,773)
      End If
C  
      KINT = 0
      nkm = max(nk,nkir)
C                         
      DO 70 K=1,NKM
C
        If(K .le. NK) Then
          xld = (k-1)*dx
          pxl = 0.25*xl(j)
          RILLW1(K) = rillwj
          winter(k) = uwinj
C     CALCULATE RILL DEPTHS FROM ACROSS SLOPE ROUGHNESS
C
          RILLD1(k) = 0.
          if(rillw(j) .gt. 0) then
C   CALCULATE THE INITIAL RILL WIDTH RATIO (changes with rill overtopping)
            RWR(K)=rillw(j)/WPR
            tw(k) = rillwj
            widw(k) = rillwj
C
            RILLD1(K) = rilldj
C Proposed variable rill width for upper plane:
C  rs() added 3/92
            if(rs(j) .le. 0 .and. nu(j) .le. 0) 
     &               rilld1(k) = rilldj*sqrt((xld+pxl)/(xl(j)+pxl))
C
C          Call PROFIL(j,k)
          Else
C  No Rills:
           rwr(k) = 0.
            tw(k) = wpr
            widw(k) = wpr                    ! Add 7/96
            depno(j) = 1.
C
          end if
          QL2(K)=QI(1)
          H1(K)=0.0
          V1(K)=0.0
          a1(k) = 0.
          a2(k) = 0.
          arrild(k)=rilld1(k)
          arrilw(k)=rillw1(K)
          smda(k) = 0.
          arill(k) = rilld1(k)*(rillw1(k) + zrl(k)*rilld1(k))
        End If
        ar2(k) = 0.
c                        archive rill d and w for use in hydwri 
        ar1(k) = 0.
C      widql(k) = w(j)
        h2(k) = 0.
        d1(k) = 0.
        vir(k) = 0.
        d2(k) = 0.
        power(k) = powerj
        alpha(k) = alphaj
        IF (NEROS.LE.0) GO TO 70
        DRLNET(K) = 0.
        DIRNET(K) = 0.
   70 CONTINUE
C
      If(rillw(j) .gt. 0.) Then                                         !4/92
        Write (IDIAGN,601) wpr,zrl(j),rillw1(1),rilld1(1),rillw1(nk),
     +    rilld1(nk)   
      End If                                                            !4/92
      IF (NEROS.Ge.1) THEN
        VS=VSETL(RHOS(J),D50(J),XNU)
        IF(NP(J) .GE. 2) Then
          write(idiagn,"(/' xnu: ',g13.4)") xnu
          WRITE(IDIAGN,631) VS
        End If
      END IF
      KB=2
      L=2
      qupt = 0.
      qsup = 0.
      T(1)=0.0
      Q(1)=0.0
      T(L)=T(1)+DELT
      RF(L)=QI(1)
      DO 80 I=1,LENQ
   80 CONC(I)=0.0
C
C      CHECK FOR OTHER PLANES CONTRIBUTING.  NU(J)=0 - NO CONTRIBUTORS.
      Write(*,*)'      '
      If(NU(J) .ne. 0) Then
C
C     OTHER PLANES ARE CONTRIBUTING -- CHECK TO SEE IF CORRECT LAW IS IN
C     EFFECT -- CALCULATE UPPER BOUND DEPTHS (HUB) FOR ALL TIME INCREMNT
        JJ=NU(J)
        SUMA(J) = SUMA(JJ) + ( W(J)*XL(J) )
        MM=NB(JJ)-1
        IF (NB(JJ).EQ.0) CALL ERRPO ('PLANE',15,JJ,0,'X','X')
        DO 100 M=1,NI
           IF (NEROS.GE.1) THEN
             CUB(M)=SCON(M+MM)
           ENDIF
C Eurosem: redefine hub as discharge.  HUB not used elsewhere
           hub(M) = QS(M+MM)
  100    CONTINUE
        HUB(MB)=HUB(NI)
        IF (NEROS.GE.1) CUB(MB)=CUB(NI)
        NB(JJ)=0
      Else
C
C     NO PLANES ARE CONTRIB. -- DEFINE NI (NO. OF TIME INCREMENTS)-- SET
C     ALL UPPER BOUND DEPTHS (HUB) TO ZERO.
        MB=IFIX(TFIN/DELT)+2
        SUMA(J) = XL(J)*W(J)
        IF (MB.GT.LENQ) GO TO 370
        NI=MB-1
        DO 120 M=1,MB
          HUB(M)=0.
          CUB(M)=0.
  120   CONTINUE
C
      End If
C  BEGIN TIME STEP LOOP:::::::::::::::::::::::::::::::::::::::::::::::::
C
  130   CONTINUE
C          Clock
      ndx = ndx + 1
      If(ndx.gt.4) ndx = 1
      Write(*,917)j,line(ndx)
  917 Format('+','    Simulating Plane ',I2,'...',a1)
C       INITIALIZE Cumulative ELEMENT INFILTRATION
        EINF(J) = 0.0
      RF(L)=QI(KB-1)
      RF(L-1)=RF(L)
      IF (T(L).LE.TFIN) GO TO 140
      IF (T(L-1)+.0001.GE.TFIN) GO TO 340
      T(L)=TFIN
  140 IF (T(L).LT.TIDD(IGAGED,KB)-.0001) GO TO 150
      T(L)=TIDD(IGAGED,KB)
C
C     CALCULATE INITIAL ADVANCED TIME DEPTH (H2(1)) USING HUB AND A TIME
C     INTERPOLATION.
  150 MT=IFIX(T(L)/DELT)+1
C
      qupt= HUB(MT)+(HUB(MT+1)-HUB(MT))*((T(L)-(DELT*FLOAT(MT-1)))/DELT)
CC  Find initial h corresp. to discharge "qupt":
C
      qup = qupt/depno(j)
      if(qup .le. 1.e-6) then
        h2(1) = 0.
        a2(1) = 0.
        v2(1) = 0.
      else if(RILLW1(1) .gt. 0.) then
        cz = sqrt(1. + zrl(j)*zrl(j))
        wr = RILLW1(1)
        dr = RILLD1(1)
        pw = power(1)
        alp = alpha(1)
        wpril = wr + 2.*cz*dr
        aril = dr*(wr + zrl(j)*dr)
        hril = aril/wpril
        qrm = alpha(1)*wpril*hril**pw
        if(qrm .ge. qup) then
          ht = dr*(qup/qrm)**(1./pw)
        else
          ht = ((qup-qrm)/alp/(winter(1)+wr))**(1./pw) + dr
        end if
        iter = 1
  151   continue
          call perim(1,j,1,pw,ht,twd,at,VQ,dVVdh,vrf)
          F = alp*VQ - qup
          DFDH = alp*dVVdh
          dlH = F/DFDH
          ht = ht - dlH
          if(abs(dlH) .le. 0.001*ht) go to 152
          if(ht .le. 0.) ht = ht + 0.5*dlH
          iter = iter + 1
          if(iter .gt. 20) stop ' GT. 20 iters on UB/ plane'
          go to 151
  152   continue
        h2(1) = ht
        vr(1) = vrf*alp
        widw(1) = twd
        v2(1) = qup/at
        a2(1) = at
      else
C No rills:
        pwi = 1./power(1)
        widw(1) = wpr
        h2(1) = (qupt/(w(j)*alpha(1)))**pwi
        a2(1) = h2(1)*w(j)
        v2(1) = qupt/a2(1)
      end if
C
      IF (NEROS.GE.1) Then
C********************************INITIALIZE CT(1) 
        ct(1) = 0.
        IF (H2(1) .GT. 0.0) THEN
          CT(1)=CUB(MT)+(CUB(MT+1)-CUB(MT))*((T(L)
     &          -(DELT*FLOAT(MT-1)))/DELT)
          qsup = ct(1)*qupt
        END IF
      End If
C
  160 DT=T(L)-T(L-1)
C
      suspir = 0
C  use KINT as code for options in XPLINF:
      If(INTRT .gt. 0 .or. rillwj .gt. 0.) kint = 1
      CALL XPLINF (J,RF(L),QL2,NK,1)                 ! XPLINF rill or main
C
      IF (NP(J) .EQ. 7) THEN
        PPT=QI(KB-1)*3.6e6
        WRITE (IDIAGN,490) T(L),PPT
      END IF
C
C      STORE DEP HAS UNITS OF M, QL M/SEC SO NEED TO DIVIDE
C      STODEP BY DT, THE TIME IN SECS BETWEEN TIME DEPTH PAIRS
        QIO = ql2(nk-1)
        DQL2 = 0.
        STOF = 0.
        stov = 0.
        IF(QIO .GT. 0.) THEN
          IF(STOCAP .GT. 1.E-6) THEN
            QIOL = QIO*contrib
            DQL2 = STOCAP/DT
            IF (DQL2 .GT. QIOL) DQL2 = QIOL
            DSTO = DQL2*DT
            STOCAP = STOCAP - DSTO
            STDPM = STDPM + DSTO
          END IF
        ELSE IF(QIO .LT. 0.) THEN
          IF(STDPM .GT. 1.E-6) THEN
C            DSTO = -QIO*DT*0.15*STDPM/STODEP(J)
            DSTO = -QIO*DT*0.15*0.5
            IF(DSTO .GT. STDPM) DSTO = STDPM
            STDPM = STDPM - DSTO
C Infil depth from storage
            STOF = DSTO
            STOCAP = STOCAP + DSTO
          END IF
        END IF
C
      If(INTRT .ge. 1) Then
        KINT = nkir
        Call XPLINF(J,RF(L),qinrl,1,NKIRP)                   ! XPLINF Interrill
        do 171 k=1,nkir
          kp = k+1
          ql3(k) = qinrl(k+1)
          v1(k) = vir(k)
          fr = 1.
          If(ql3(k) .lt. -1.e-6) Then
            If(recs(j) .gt. 0.) fr = sqrt(2.*d2(k)/recs(j))
            fr = amin1(1., fr)
            ql3(k) = ql3(k)*fr
            qinrl(k+1) = ql3(k)
          Else If(dql2 .gt. 0.) Then
            ql3(k) = ql3(k) - dql2
            If(ql3(k) .le. 0.) ql3(k) = 0.
C            If(k .ge. nkir) dql2 = 0.
          End If
          qmph(k) = ql3(k)*3.6e6
  171   continue
C
        izr = i
        if(rf(l) .gt. 0. .or. d1(nkir) .gt. 0.)
     +    Call IMPLCT(nkir,j)                                           !IMPLCT
C  IMPLCT gets new interrill values d2() and vir()
        xinrl = 0.
        stol = 0.
        Do 172 k=1,nkir
          ar2(k) = d2(k)
          ar1(k) = d1(k)
          tw(k) = 1.
          d1(k) = d2(k)
          If(k.lt.nkir) Then
            kp = k+1
            qia = 0.5*(ql3(k)+ql3(kp))
            stol = stol + d2(k) + d2(kp)
            xinrl = xinrl + (rf(L)-qia- dql2)*dcx/ofl
          end IF
          If(NEROS .ge. 1) Then
            winh(k) = d2(k)*100.
            winqm = vir(k)*d2(k)
C  winq in cm**2/sec for use in CAPAC
            winq(k) = winqm*1.e4
            spow(k) = 0.
            If(d2(k) .gt. 1.e-6) Then
             If(mcode .ge. 1) Then
C Bagnold effective stream power (Everaert):
               spowb = grav*100.*winq(k)*sir(j)
               spow(k) = spowb**(1.5)/(d2(k)*100.)**.667
             Else
C Govers (rills) metric stream power (cm/sec):
               SPOW(k) = vir(k)*sir(j)*100.        ! corrected 1/98
             End If
            End If
          End If
  172   Continue
C  storage on interrill area:
        stov = stov + stol*dcx*xl(j)/ratl
C
	
        If(NEROS .ge. 1) Then
          Call SEDCOM(J,RF(L),dcx,AR2,AR1,0.,NKIR,0)                    !SEDCOM
          cll = clir
C interrill supply conc for call to SEDCOM for rills:
          clir = ct(nkir)
          clat = 0.5*(cll+clir)
        End If
C disch per unit length along rill, both sides:
        qlo = 2.*winqm
        qltol = qlt
        qlt = qlo/ratl
        qlrav = 0.5*(qltol + qlt)
C net sed from interrill
        qsirr = xl(j)*(qlt*clir + qltol*cll)*0.5
      End If
C
      DO 180 I=1,NK
C   CALCULATE RECESSION FACTOR AND REDUCE INFILTRATION BY THIS
C   FACTOR IF LATERAL INFLOW IS NEGATIVE (I.E. RAINFALL CONCLUDED)
C   NOTE THAT RECS IS A USER INPUT PARAMETER  REFLECTING PLANE
C   MICROTOPOGRAPHY.  IT IS ASSUMED THAT THE AREA CONTRIBUTING TO
C   INFILTRATION VARIES LINEARLY WITH DEPTH DURING RECESSION.
C
      RECES=AMIN1(1.,sqrt(2.*H1(I)/RECS(J)))              ! 3/95
      IF (QL2(I).LT.-1.0E-6) THEN
C
        If(INTRT .ge. 1) Then
          rwf = 1.
        Else If(LCH .gt. 0) Then
C This should be (width of flow) / (rill spacing):
          rwf = RWR(i)
        Else
          rwf = RECES
        End If
C
C        qlt = rf(l)
        QL2(I)=QL2(I)*rwf
      ELSE
        QL2(I) = QL2(I) - DQL2
      ENDIF
C
      v1(i) = v2(i)
      XIN(I) = (RF(L)-QL2(I)-dql2)*dt
      QmmPH(I) = QL2(I)*3.6e6
      If(INTRT .ge. 1) Then
        eftw = rillw1(i)
C         
C *Net input redefined: total from interrill and rill areas:
C    discharge per unit area along rill area
        ql2(i) = (qlrav + ql2(i)*eftw)/wpr
        XIN(I) = (xin(i)*eftw +  xinrl*dt*winter(i))/wpr ! infil dpth in m
      End If
C
  180 CONTINUE
      kint = 0
      IF (NP(J) .eq. 7) Then
  !        xinrl = xinrl*3.6e6
        If(INTRT .ge.1) write (idiagn,511) (qmph(i),i=1,nkir)
C temp change for comparison:
        WRITE (IDIAGN,511) (Qmmph(I),I=1,NK)
C       WRITE INFIL ave. rate
        WRITE (IDIAGN,500) (XIN(I)/dt*3.6e6,I=1,NK)
      End If
C
  190 CONTINUE
C                MASS BALANCE AND ACTIVITY TEST
       QLS = 0.
       H1S = 0.
      DO 200 K=1,NK
      IF(H1(K) .LE. 0.) THEN
        IF(QL2(K) .LE. 1.E-8) THEN
          QL2(K) = 0.
        ELSE
          QLS = QLS + QL2(K)
          IF(NEROS .Ge. 1) THEN
            IF(Ct(K) .GT. 0.2) Ct(K) = 0.2
          END IF
        END IF
      ELSE
       H1S = H1S + H1(K)
      END IF
C
  200 CONTINUE
      H1S = H1S + H2(1)
      IF(H1S .GT. 0. .OR. QLS .GT. 0.) GO TO 220
      IF(T(L).GE.TIDD(IGAGED,KB)-.0001) KB = KB+1
      GO TO 235
  220 IF (T(L).LT.TIDD(IGAGED,KB)-.0001) GO TO 230
      KB=KB+1
  230   CONTINUE
C
C     IMPLCT RETURNS WITH NEW RILL DEPTHS (H2 and V2 (I=2,NK)) --
C
      CALL IMPLCT (NK,J)                                                !implct
CR!
  235   CONTINUE
      DO 260 K=1,NK
C default (min.) value:
       rwr(k) = rillw(j)/wpr
       if(rilldj .le. 0.) rwr(k) = 0.
      IF (H2(K) .le. 0.) Then
        h2(k) = 0.
        v2(k) = 0.
        a2(k) = 0.
      End If
C
      If(rilldj .gt. 0.)   rwr(k) = tw(k)*depno(j)/w(j)
      ar1(k) = a1(k)
      ar2(k) = a2(k)
C
  260 CONTINUE
C
C           CALCULATE TOTAL DISCHARGE (CMS) FROM PLANE J FOR TIME INCREMENT L.
      Q(L)=V2(NK)*ar2(nk)*depno(j)
C
      IF (NEROS.GE.1) THEN
C
      DO 789 K=1,NK
C
C  Estimate rill and overrill depths, q, and v for splash
C
C  Smooth velocity to reduce sawtooth cmx:
        if(k .gt. 1 .and. k .lt. nk) then
          vu = 0.25*(v2(k-1) + 2.*v2(k) + v2(k+1))
        else if(k .eq. nk) then
          vu = 0.5*(v2(nk) + v2(k-1))
        else
          vu = v2(k)
        end if
C
        If(LCH .le. 0) Then

C  No Rills Case
          winh(k) = h2(k)*100.
          if(winh(k) .gt. 1.e-4) Then
            winq(k) = vu*h2(k)*1.e4
            If(mcode .ge. 1) Then
C Bagnold effective stream power (Everaert):
              spowb = grav*100.*winq(k)*s(j)
              spow(k) = spowb**1.5/(h2(k)*100.)**.667
            Else
C Govers (rills) metric stream power (cm/sec):
              SPOW(k) = vu*s(j)*100.
            End If
          Else      
            spow(k) = 0.
          End If
C
        Else
          tw(k) = widw(k)
C
          If(LCH .eq. 1 .and. INTRT .eq. 0) then
C use steady flow estimate for interrill flow depth:                    !3/92
C q into rill:
            If(ql2(k) .gt. 0.) Then
              qlir = ql2(k)*ofl*10                   
              winh(k) = (qlir/alpha(k)/slrat)**(1./power(k))*100. !IN CM.
C              qlrav = qlir*2./ratl
C              stov = 0.
            Else If(ar2(k) .gt. 0.) Then
              winh(k) = winh(k) + .5*ql2(k)*dt
              if(winh(k) .gt. 0.) then
                qlir = slrat*alpha(k)*(winh(k)/100.)**power(k)  !corrected 1/98 
C                qlrav =qlir*2./ratl
              else
                qlir = 0.
C                qlrav = ql2(k)*2.
                winh(k) = 0.
              end if
            Else
C              qlrav = 0.
              qlir = 0.
              winh(k) = 0.
            End If
            winq(k) = amax1(qlir,0.)*10000.   !  use winq in capac in cm**2/s
  !          
          End if
C
C rills metric stream power (cm/sec):
          SPOW(k) = vr(k)*s(j)*100.
          rillh(k) = h2(k)
        End If
C
789   CONTINUE
C
      If(LCH .eq. 1 .and. INTRT .le. 0) Then
            If(ql2(nk) .gt. 0.) Then
              qlir = ql2(nk)*ofl*10                   
              qlrav = qlir*2./ratl
              stov = 0.
            Else If(ar2(nk) .gt. 0.) Then
C              winh(nk) = winh(nk) + .5*ql2(nk)*dt
              if(winh(nk) .gt. 0.) then
                qlir = slrat*alpha(k)*(winh(nk)/100.)**powerj  !corrected 1/98 
                qlrav =qlir*2./ratl
              else
                qlrav = ql2(k)*2.
              end if
            Else
              qlrav = 0.
            End If
      End If
C
C  Sediment transport calculation
C
        dv = dx
        CALL SEDCOM (J,qlrav,dv,ar2,ar1,CLAT,NK,LCH)                    !SEDCOM
        CONC(L)=CT(NK)
        qsav = .5*(Q(L)*CONC(L)+Q(L-1)*CONC(L-1))
        If(LCH .le. 0) qsirr = qsav
C                        SEDTOT IN Kg
       SEDTOT=SEDTOT+ qsav*RHOS(J)*1000.*DT
       sedin = sedin + qsup*dt
C
      END IF
C
      srain=srain+RF(L)*XL(J)*W(J)*DT
      SUMQIN=SUMQIN+(V2(1)*ar2(1) + V1(1)*ar1(1))*0.5*DT*depno(j)
      dvout = (Q(L)+Q(L-1))*0.5*DT
      sqout = sqout+dvout
      SUM1 = srain+SUMQIN
      SUM2 = dvout+SUM2
C       ARRAY XIN(K) FOR VOL. BAL.
  280 CONTINUE
      XINS=(XIN(1)+XIN(NK))/2.
      STOR=(ar2(1)+ar2(NK))/2.                                          Modfd
      NKM2=NK-2
      IF (NKM2 .ge. 1) Then
       DO 290 K=1,NKM2
         XINS=XINS+XIN(K+1)
  290  STOR=STOR+ar2(K+1)
      End If
      STOR=(STOR*DX + stov)*depno(j) + STDPM*w(j)*xl(j)                 Modfd
C      write(idiagn,774) stdpM*1000.
  774 Format('   stdpm (mm):',f10.6)
      XINS=XINS*DX*W(J) + stof*w(j)*xl(j)                               Modfd
      sinfil=sinfil+XINS                                                Modfd
      SUM2=SUM2+XINS
      If(neros .ge. 1) then
C
C     EUROSEM ROUTINES
C                                 -------------RILLWIDEN
        if(rillw(j) .gt. 1.e-4) then
          sx = 1./sir(j)
          do 305 k = 1,nk
            smda(k) = smda(k) + dla(k)
            if(abs(smda(k)) .gt. 0.05*arill(k) .or. T(L) .ge. 
     &        (TFIN-.1) ) then
C
              CALL RILLWI(J,k,smda)                                     *10/92
              arill(k) = RILLD1(k)*(RILLW1(k) + RILLD1(k)*zrl(j))
              twrill = rillw1(k) + 2.*rilld1(k)*zrl(j)
C **Update h2 based on new rill profile:                                *n8/91
              acr = arill(k)
              if(ar2(k) .gt. acr) then                                  *n8/91
                apm = ar2(k) - acr
                hpm = hfun(twrill,sx,apm)
                h2(k) = rilld1(k) + hpm
              else                                                      *n8/91
                h2(k) = hfun(rillw1(k),zrl(j),ar2(k))                   *n8/91
              end if
              smda(k) = 0.
            end if
  305     Continue
        End If
C             ---------------  Check sed balance / Volume Basis
        stout = sedtot/rhos(j)/1000.
        nv = nk
        dxu = dx
        If(INTRT .ge. 1) then
          nv = nkir
          dxu = dcx
        End If
        If (INTRT .ge. 1 .or. LCH .le. 0) Then
C Routed interrill areas or No Rills:
          sumdir = 0.
          Do 306 k=2,nv
            sumdir = sumdir + dirnet(k)
  306     Continue
          sumdir = sumdir*(1.-porpl(j))*dxu*depno(j)
          sumqirr = sumqirr + qsirr*dt*max(depno(j),1.0)
        Else If(INTRT .le. 0 .and. LCH .eq. 1) Then
C  Unrouted, short interrill areas:
C            vinrl = winh(2)*ofl*dx/ratl
          sumqirr = sumqirr + qlrav*cbal*xl(j)*depno(j)*dt
          sumdir = - sumqirr
          suspir = 0.
          if(NP(j) .eq. 7) write(IDIAGN,'(" cbal ",f10.6)') cbal
C            write(*,*) suspir
C            read(*,*)
        End If

        sumdrl = 0.
        If(rillw(j) .gt. 1.e-4) Then
          If(INTRT .ge. 1) Then
            sumdir = sumdir*2.*xl(j)/ratl
            suspir = suspir*2.*xl(j)/ratl*depno(j)
          End If
C          sumdrl = 0.
          Do 307 k=2,nk
C  drlnet is net eros (-), depos (+) sectional area at location k
C      along rills.  [For unrilled areas, the alanog is dirnet(), 
C      in loop 306, above.  Both cases need to be accounted for ]
            sumdrl = sumdrl + drlnet(k)*dx
  307     Continue
C   Multiply by (1.-porpl(j))*depno(j) to get solid volume,
C   and multiply cumulative value by sediment density to get weight 
C   lost/gained up to that distance, for output purposes:
          susprl = susprl*depno(j)
          sumdrl = sumdrl*(1.-porpl(j))*depno(j)
        End If
        rillbal = stout + sumdrl + susprl - sedin
        inrlbal = sumqirr + sumdir + suspir 
        If(NP(j) .eq. 7 .or. T(L) .ge. (TFIN-0.1)) Then
          Write(idiagn,780) sumdir, suspir, sumqirr, inrlbal
          If(rills) Write(idiagn,781) sumdrl, susprl, stout, rillbal
        End If
      End if
C
      IF (SUM1.NE.0.) THEN
        er = (SUM1-SUM2-STOR)
        PER=ER/SUM1*100.0
      ELSE
        PER=0.
      ENDIF
      If(NP(J) .eq. 7) Then
C        If(intrt .gt. 0) Write(IDIAGN,518) (d2(k),k=1,nkir)
        WRITE (IDIAGN,520) (h2(K),K=1,NK)
        Write (IDIAGN,521) (widw(k),k=1,nk)
        write(idiagn,519) (v2(k),k=1,nk)
        IF (NEROS.GE.1) THEN
C          WRITE (IDIAGN,550) (CMX(K),K=1,NK)   ! These are output in SEDCOM
C          WRITE (IDIAGN,530) (CT(K),K=1,NK)
          write(idiagn,'(10x,g12.4)')spow(nk)
          If(LCH .ge. 1) Then                                ! change 7/96
            WRITE (IDIAGN,540) (DRLNET(K),K=1,NK)            !
          Else                                               !
            WRITE (IDIAGN,540) (DIRNET(K),K=1,nkir)          !  2/98
          End If                                             !
        ENDIF
C
C          VOL BAL. WRITE
C        write(idiagn,571) stov*depno(j),stof
  571  Format(' interrill and micro stor:',2g12.4)
        WRITE (IDIAGN,560) srain,SUMQIN,sinfil,sqout
        WRITE (IDIAGN,570) SUM1,SUM2,STOR,pER
      End If
C
C   INCREMENT L AND INTERNAL TIME  -- REDEFINE LATERAL INFLOW (QL),AND
C   SET ADVANCE TIME DEPTHS (H2) EQUAL TO THE KNOWN DEPTHS (H1)
      L=L+1
      IF (L .GT. LENQ) Then
        WRITE (IDIAGN,580)
        STOP 'Too Many Q Time Steps: Exceeds Dimensions'
      End If
      T(L)=T(L-1)+DELT
      DO 330 K=1,NK
        A1(k) = A2(K)
        H1(K)=H2(K)
  330 continue
      GO TO 130
C
CR!     END TIME LOOP ::::::::::::::::::::::::::::::::::::::::::::::::::
C
  340 CONTINUE
      L=L-1
      LASTNB=NB(1)
      IF (NELE.GT.1) Then
C       DO OVER ALL ELEMENT
        DO 350 NE=2,NTELE
          IF (NB(NE).GT.LASTNB) LASTNB=NB(NE)
  350   CONTINUE
      End If
      MBT=LASTNB+NI
      IF (LASTNB.EQ.0) MBT=1
      MBE=MBT+NI
      IF (MBE .gt. LENQS) Go To 370
      CALL UNIF (Q,T,L,QS(MBT),NI,DELT,CONC,SCON(MBT))
C                 CODE READ THE UNIFORM TIME INCREMENTED OUTFLOWS
C                    BACK INTO ARRAY Q AND ALSO EVEN TIME INCREMENTS INTO
C                    ARRAY T
        IEND = NI + 1
        DO 385 II=1,IEND
           ID = MBT + (II-1)
           Q(II) = QS(ID)
           T(II) = FLOAT(II-1) * DELT
  385   CONTINUE
C       SAVE TOTAL INFILTRATION AND STORAGE FOR GLOBAL VOL BAL
        EINF(J) =sinfil
C Total storage on plane at end of time in cu. m.:
        STORA(J)=stor
C        METRIC CONVERSION
      NB(J)=MBT
      IF (NP(J) .gt. 0) Then
        wpr = w(j)/depno(j)
        XLT=XL(J)
        WT=W(J)
        RECST=RECS(J)                                                   !11/92
        ALT=AL(J)
        FMINT=FMIN(J)*60.
        SUM1T=SUM1
        SUM2T=SUM2
        STORT=STOR
        WRITE (IDIAGN,600) XLT,WT,S(J)
        If(Rillw(j) .gt. 0.)Write(IDIAGN,601) wpr,zrl(j),rillw1(1),
     +         rilld1(1),rillw1(nk),rilld1(nk)
        WRITE (IDIAGN,660) R1(J)
        IF (FMIN(J).GT.1.0E-08) THEN
C        
          WRITE (IDIAGN,700) FMINT,alt,SI(J),PORPL(J),
     1                        SMAX(J),ROC(J),RECST
        ELSE
          WRITE (IDIAGN,710)
        ENDIF
        IF (NEROS .ge. 1) Then
         WRITE (IDIAGN,630) D50(J),RHOS(J),PORPL(J),PAVE(J)
         If(LCH .ge. 1) Then
           WRITE (IDIAGN,640) (DRLNET(K),K=1,NK)
         Else
           Write(IDIAGN,640) (DIRNET(K),K=1,NK)
         End If
        End If
        WRITE (IDIAGN,610) SUM1T,SUM2T,STORT,PER
      End IF
C
C          WRITE OUT SUMMARY INFOR
      TYPE = 'PLANE'
      WRITE (IWRITE,605) J,TYPE,PER,SEDTOT
      stodmm = stodep(j)*1000.
      write(iwrite,606) stodmm
      cmtkg = -rhos(j)*1000.
      sednet = sedtot + sedin*cmtkg
      areaco=10./(w(j)*xl(j))
C      areanet = 10./suma(j)
      Write(99,1020)
      If(rillw(j) .gt. 1.e-4) Then
        write(99,1021)sumdrl*cmtkg, sumdrl*cmtkg*areaco
      End If
      write(99,1022) sumdir*cmtkg,
     - sumdir*cmtkg*areaco, sednet, sednet*areaco
      IF (NP(J).EQ.0) THEN
           WRITE (IDIAGN,650)
           GO TO 460
        ENDIF
  460 CONTINUE
      RETURN
C
  370 WRITE (IDIAGN,580)
      STOP 'Too many time steps for storage dimensions'
 1020 FORMAT(//,3X,'EROSION SUMMARY',/,3X,'---------------')
 1021 Format(3X,'GROSS RILL EROSION',9X,F9.3,' kg',f9.3,' t/ha')
 1022 Format(3x,'GROSS INTERRILL EROSION',4X,F9.3,' kg',f9.3,' t/ha'//,
     - 3X,'NET EROSION/DEPOSITION',5X,F9.3,' kg',f9.3,' t/ha',/,
     - 3x,'(a minus denotes deposition)',/)     
  470 FORMAT (/,' ',10X,'NK CALCULATED TO BE',I3,', BUT SET AT 15')
  480 FORMAT(/,' ',10X,'THE RAIN GAGE FOR PLANE',I3,' IS GAGE NO.',I3,/,
     1' ',7X,'PPCT. WEIGHT IS',F5.2,'   INTERCEPTION IS',F5.2,2X,A5)
  490 FORMAT (5X,'TIME IS',F11.5,' SEC., RAIN RATE IS ',F15.5,' MMPH',
     1      18('*'))
  500 FORMAT (1X,'INFILmmh=',6F9.6/(8F9.6))
  511 FORMAT (1X,'QL(mmPH)=',6F9.5/(8F9.5))
  518 FORMAT (1X,'d2(I=1,Ni)=',6F10.5/(8F9.5))
  519 Format(' v2(k=1,nk)=', 6f9.7/(8f9.7))
  520 FORMAT (1X,'H2(I=1,NK)=',6F10.5/(8F9.5))
  521 FORMAT (1X,'WD(I=1,NK)=',6F10.4/(8F9.4))
  530 FORMAT (1X,'CT(I=1,NK)=',6F9.6/(8F9.6))
  540 FORMAT (1X,'DANET(I=1,NK)=',5(G11.5,1X)/(' ',6(G11.5,1X)))
  550 FORMAT (1X,'CMX(I=1,NK)=',6F9.6/(8F9.6))
  560 FORMAT(' RAIN=',E10.3,' UPSTREAM IN=',E10.3,' INFIL=',E10.3
     1          ,' RO OUT=',E10.3)
  570 FORMAT (' INFLOW=',E10.3,' OUTFLOW=',E10.3,' STOR. =',E10.3,' ERRO
     1R=',E10.3,' %',/)
  580 FORMAT ("   You have exceeded EUROSEM's storage limit. Try reducin
     1g the storm length,"/' increasing the timestep, or reducing the ',
     2' complexity of the storm.')
  600 FORMAT (/,' ',5X,'GEOM. PARAMETERS ARE  L=',F7.1,
     1'  W=',F7.1,'  S=',F7.4)
  601 FORMAT(' Every ',f5.2,' m there is a rill with sideslope 'f5.2,
     +/'    Width(m) and Depth(m) at Top of slope:   ',2f10.5,/
     +'    Width and Depth at Bottom:               ',2f10.5/)
  605 FORMAT (' ',9X,I3,6X,A5,7X,E10.3,5X,F11.3)
  606 Format (9x,f8.4,' mm Inactive Storage Capacity on plane')
  610 FORMAT (' ',/,16X,'**** WATER BALANCE AT END OF PLANE ****',/,
     1  ' ',11X,'<INFLOW BASED ON (PPT*GAGE WT) - INTER. + RUNON>',/,
     2  ' ',2X,'INFLOW=',E10.3,' OUTFLOW=',E10.3,' STOR.=',E10.3,
     3      ' ERROR=',E10.3,' %',/)
  620 FORMAT (' TOTAL SOIL LOSS IN KGS. =',F11.3/)
  630 FORMAT (' ',5X,'EROSION PARAMETERS ARE ---',/,
     3    ' ',6X,' D50=',G9.3,' RHOS=',F6.2,' POR=',F6.2,
     4    ' PAVE.FAC.=',F6.3)
  631  FORMAT(/'   Particle settling velocity is ',F8.6,' m/s'/)
  640 FORMAT (' ',5X,'ACCUMUL. SURFACE DEPOSIT. OR EROSION (NEG.)
     1AT EACH NODE (m.)'/,
     2          ' ',6(G12.5,1X))
  650 FORMAT (30X)
  660 FORMAT (' ',5X,'ROUGHNESS: MANNINGS N=',F5.3)
  700 FORMAT (' ',5X,'INFILT. PARAMETERS ARE FMIN=',F10.5,' mm/h;  G=',
     1F8.3,' mm'/5x,'POR=',F7.4,'  SMAX=',F7.4,'   SI=',F7.4,'   ROC='
     2,F6.3,'  RECS=',F7.2,' mm')
  710 FORMAT (' ',5X,'IMPERVIOUS PLANE')
  730   FORMAT (//,' ',20X,'RAINFALL HYETOGRAPH FOR PLANE NO.',I4,/,
     +    ' ',10X,'(AFTER INTERCEPTION REMOVED)'10X,'Kinetic Energy (J/m
     +2/mm)', /,' ',8X,'TIME (MIN)',5X,'INTENSITY(',A5,')',6x,
     +           'Rain',12X,'Leaf Drip' )
  740 FORMAT (8X,F8.1,13X,F8.2,2(10x,f8.3))
  770 Format(/' Downslope flow distance divided into ',I2,' nodes',
     & F5.2,' m apart.')
  771 Format(/' Interrill Flow routed explicitly, using',I2,' nodes'/
     &'  with total interrill distance of ',F6.2,' m'/)
  772 Format(/' Short Interrill flow length (',F5.2,' m): not dynamicall
     &y routed'/)
  773 Format(/' Surface contains no explicit rills'/)
  780 Format('  INtRill eros, susp,    sedout, and Bal. (m*3):'/,
     + 2x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)
  781 Format('  Rill eros,   susp,     sedout, and Bal. (m*3):'/,
     + 2x,f9.5,2x,f9.5,2x,f9.5,2x,f9.5)
  800 FORMAT(' ',8X,'*** PLANE NO. ',I4,' DIAGNOSTIC INFORMATION ***')
      END
C-------------------------------------------------



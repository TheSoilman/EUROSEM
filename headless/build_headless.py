#!/usr/bin/env python3
"""
build_headless.py  --  makes the three staged, headless-capable copies of
the commented EUROSEM sources.

    python3 build_headless.py            (run from this folder)

Reads  ../*_commented.for  (comments-only versions + the LOG10 change)
Writes stage1_fatal/  stage2_headless/  stage3_mbfix/  stage4_rval/  (each complete)

Each stage includes the one before it, so a problem can be bisected:
  stage 1  every STOP becomes CALL FATAL(...): the message goes to the
           screen, the auxiliary file and a status file, and the exit
           code is 2.  A normal run ends through ENDOK (status "OK").
           No change to results, still interactive.
  stage 2  adds the control file eurosem.ctl.  When it exists, every
           keyboard prompt is skipped and file names (now up to 80
           characters) come from the file.  Without it the program asks
           questions exactly as before.
  stage 3  PLANE sets MB and NI itself when an element receives upstream
           inflow (they used to be left over from the previous element).
  stage 4  RVAL (rock-content factor from ROC) becomes one value per
           element; before, only the last element's ROC took effect.
Every edit must match exactly once, or the script stops.
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, ".."))
FILES = ["Euroa", "emainl", "eplanel", "eurozl", "eurseml", "eusedl"]


def load(name):
    return open(os.path.join(SRC, name + "_commented.for"), encoding="latin-1").read()


def rep(text, old, new, label):
    n = text.count(old)
    if n != 1:
        sys.exit(f"EDIT FAILED ({label}): found {n} matches")
    for ln in new.split("\n"):
        body = ln.expandtabs(6)
        if not ln.startswith(("C", "c", "*")) and len(body.rstrip()) > 72:
            sys.exit(f"LINE TOO LONG ({label}): {len(body.rstrip())} cols: {ln}")
    return text.replace(old, new)


# ---------------------------------------------------------------- stage 1
def stage1(t):
    t["Euroa"] = rep(t["Euroa"], "                STOP\n             ELSE",
        "                CALL FATAL('Over 40 iterations in bisection (ITER)')\n             ELSE", "E400")
    t["Euroa"] = rep(t["Euroa"], "     &            stop 'zero alpha or power in CHANNL'",
        "     &            CALL FATAL('Zero alpha or power in CHANNL')", "E568")
    t["Euroa"] = rep(t["Euroa"], "          STOP 66",
        "          CALL FATAL('FCP called with R less than KS')", "E838")
    t["Euroa"] = rep(t["Euroa"], "      STOP 'Iteration limit in FINTG'",
        "      CALL FATAL('Iteration limit in FINTG')", "E942")

    m = t["emainl"]
    m = rep(m, "      STOP 'NORMAL COMPLETION'", "      CALL ENDOK", "m450")
    m = rep(m, "        Write(*,*)' dimension of output array exceeded. '\n        Stop",
        "        Write(*,*)' dimension of output array exceeded. '\n"
        "        CALL FATAL('Dimension of output array exceeded (ND>100)')", "m501")
    m = rep(m, "        Stop ' Parameter Error:  TFIN/DELT must be less than 1000'",
        "        CALL FATAL('TFIN/DELT must be less than 1000')", "m941")
    m = rep(m, "        Stop ' Parameter Fatality: XL must be .gt. 0.'",
        "        CALL FATALI('XL must be greater than 0 on element',J)", "m1172")
    m = rep(m, "     & ' must be at least THMX [ ',F5.3,']')\n        Stop\n",
        "     & ' must be at least THMX [ ',F5.3,']')\n"
        "        CALL FATALI('Porosity less than THMX on element',J)\n", "m1200")
    m = rep(m, "     & ' THR [',F5.3,'] must be at least THMX [ ',F5.3,']')\n        Stop\n",
        "     & ' THR [',F5.3,'] must be at least THMX [ ',F5.3,']')\n"
        "        CALL FATALI('THI greater than THMX on element',J)\n", "m1205")
    m = rep(m, "     &ead as',f8.4,'] must be between 0. and .999 ')\n        Stop\n",
        "     &ead as',f8.4,'] must be between 0. and .999 ')\n"
        "        CALL FATALI('ROC outside 0 to 0.999 on element',J)\n", "m1211")
    m = rep(m, "is less than total rills: [',F5.1,']')\n           Stop\n",
        "is less than total rills: [',F5.1,']')\n"
        "           CALL FATALI('Width less than total rills, element',J)\n", "m1270")
    m = rep(m, "         If(answ .ne. 'Y' .and. answ .ne. 'y') stop\n",
        "         If(answ .ne. 'Y' .and. answ .ne. 'y')\n"
        "     &     CALL FATALI('SPLTEX out of range, not accepted. Element',J)\n", "m1346")
    m = rep(m, "      STOP' Error Stop.  Read Diagnosis in Auxiliary Output File'",
        "      CALL FATALI('Error stop: see auxiliary file, error no.',I)", "m1774")
    m = rep(m, "       Stop ' NTIME value out of range'",
        "       CALL FATAL('NTIME value out of range')", "m1921")
    m = rep(m, "STOP ' Selection parameter failed - see DIAGN output'",
        "CALL FATAL('Selection parameter failed - see DIAGN output')", "m1947")
    t["emainl"] = m

    p = t["eplanel"]
    p = rep(p, "       If(alphaj .le. 0. .or. powerj .le. 0.) stop 'zero alpha or power'",
        "       If(alphaj .le. 0. .or. powerj .le. 0.)\n"
        "     &   CALL FATALI('Zero alpha or power on element',J)", "p351")
    p = rep(p, "          if(iter .gt. 20) stop ' GT. 20 iters on UB/ plane'",
        "          if(iter .gt. 20)\n"
        "     &      CALL FATAL('More than 20 iterations on UB, plane')", "p637")
    p = rep(p, "        STOP 'Too Many Q Time Steps: Exceeds Dimensions'",
        "        CALL FATAL('Too many Q time steps: exceeds dimensions')", "p1184")
    p = rep(p, "      STOP 'Too many time steps for storage dimensions'",
        "      CALL FATAL('Too many time steps for storage dimensions')", "p1299")
    t["eplanel"] = p

    t["eurozl"] = rep(t["eurozl"], "           stop ' bad thin'",
        "           CALL FATAL('Bad THIN value in infiltration routine')", "z372")
    t["eurseml"] = rep(t["eurseml"], "            stop' excess rill fill. Check erosion parameters'",
        "            CALL FATAL('Excess rill fill. Check erosion parameters')", "s149")

    u = t["eusedl"]
    u = rep(u, "            STOP\n         ENDIF",
        "            CALL FATAL('Too many bisections in POND')\n         ENDIF", "u342")
    u = rep(u, "      STOP 9107", "      CALL FATAL('Volume out of range in HOFV (pond rating)')", "u534")
    u = rep(u, "       stop' neg. or zero depth with deposition in sedcom'",
        "       CALL FATAL('Zero or negative depth with deposition in SEDCOM')", "u777")
    u = rep(u, "      STOP 9106", "      CALL FATAL('Elevation out of range in VOFH (pond rating)')", "u1101")
    u = rep(u, "      STOP ' Convergence failure in finding settling velocity'",
        "      CALL FATAL('Convergence failure finding settling velocity')", "u1216")
    t["eusedl"] = u


EHLESS_1 = r"""C=======================================================================
C  EHLESS.FOR  -  support routines for running EUROSEM without a
C  keyboard (for a web service or batch job).
C
C  FATAL(MSG)     stops the run with a message.  The message goes to the
C                 screen, to the auxiliary file (if open) and to the file
C                 eurosem.status as "ERROR: ...".  Exit code is 2.
C  FATALI(MSG,I)  the same, with an integer (element or error number)
C                 added to the end of the message.
C  ENDOK          normal end of run: writes "OK" to eurosem.status.
C  HEADLS         (stage 2) looks for the control file eurosem.ctl.
C
C  A wrapper should run each case in its own empty folder, so that an old
C  eurosem.status can never be mistaken for the status of a new run.
C=======================================================================
      SUBROUTINE FATAL(MSG)
      CHARACTER*(*) MSG
      COMMON /IO/ IREAD, IWRITE, IDIAGN, JREAD, IPOND, IVOUT
      LOGICAL OP
      WRITE(*,'(2A)') ' FATAL ERROR: ',MSG(1:LEN_TRIM(MSG))
      INQUIRE(UNIT=IDIAGN,OPENED=OP)
      IF(OP) WRITE(IDIAGN,'(2A)') ' FATAL ERROR: ',MSG(1:LEN_TRIM(MSG))
      OPEN(97,FILE='eurosem.status',STATUS='UNKNOWN')
      WRITE(97,'(2A)') 'ERROR: ',MSG(1:LEN_TRIM(MSG))
      CLOSE(97)
      STOP 2
      END
C-----------------------------------------------------------------------
      SUBROUTINE FATALI(MSG,I)
      CHARACTER*(*) MSG
      CHARACTER*12 NUM
      WRITE(NUM,'(I12)') I
      CALL FATAL(MSG(1:LEN_TRIM(MSG))//' '//ADJUSTL(NUM))
      END
C-----------------------------------------------------------------------
      SUBROUTINE ENDOK
      OPEN(97,FILE='eurosem.status',STATUS='UNKNOWN')
      WRITE(97,'(A)') 'OK'
      CLOSE(97)
      STOP 'NORMAL COMPLETION'
      END
"""

EHLESS_2 = r"""C-----------------------------------------------------------------------
C  HEADLS  -  headless-run switch.
C  If the file eurosem.ctl exists in the working folder the run is
C  headless: LHEAD is set true, and these 7 lines are read:
C     1  run title (up to 80 characters)
C     2  parameter file name        (up to 80 characters)
C     3  rainfall file name
C     4  dynamic output file name
C     5  auxiliary (diagnostic) output file name
C     6  static output file name
C     7  Y or N: answer given to the questions about ROC, PAVE and the
C        SPLTEX range check.  Y accepts the adjustment, N declines it
C        (a declined SPLTEX warning stops the run).
C  If the file does not exist LHEAD is false and the program behaves as
C  it always has, asking its questions at the keyboard.
C-----------------------------------------------------------------------
      SUBROUTINE HEADLS
      LOGICAL LHEAD, EX
      CHARACTER*1 HANS
      CHARACTER*80 HTITLE, HFIL(5)
      COMMON /HLESSL/ LHEAD
      COMMON /HLESSC/ HANS, HTITLE, HFIL
      LHEAD = .FALSE.
      INQUIRE(FILE='eurosem.ctl',EXIST=EX)
      IF(.NOT.EX) RETURN
      OPEN(10,FILE='eurosem.ctl',STATUS='OLD')
      READ(10,'(A)',ERR=90,END=90) HTITLE
      DO 10 K=1,5
        READ(10,'(A)',ERR=90,END=90) HFIL(K)
   10 CONTINUE
      READ(10,'(A)',ERR=90,END=90) HANS
      CLOSE(10)
      LHEAD = .TRUE.
      RETURN
   90 CALL FATAL('eurosem.ctl is incomplete: it needs 7 lines')
      END
"""

DECL = ("      LOGICAL LHEAD\n      CHARACTER*1 HANS\n"
        "      CHARACTER*80 HTITLE,HFIL(5)\n"
        "      COMMON /HLESSL/ LHEAD\n      COMMON /HLESSC/ HANS,HTITLE,HFIL\n")


# ---------------------------------------------------------------- stage 2
def stage2(t):
    m = t["emainl"]
    # MAIN: declarations, then skip banner / title / graph question
    m = rep(m, "      LOGICAL IFSED, IFAD, ifgrf\n", "      LOGICAL IFSED, IFAD, ifgrf\n" + DECL, "main decl")
    m = rep(m, "      WRITE(*,136)\n136   FORMAT(//,6X,",
        "C  Headless run: title comes from eurosem.ctl, no questions asked.\n"
        "      CALL HEADLS\n      IF(LHEAD) THEN\n        TITLE = HTITLE\n"
        "        IFGRF = .FALSE.\n      ELSE\n      WRITE(*,136)\n136   FORMAT(//,6X,", "main skip")
    m = rep(m, "      if(ans.eq. 'Y' .or. ans .eq. 'y') Ifgrf = .true.\n",
        "      if(ans.eq. 'Y' .or. ans .eq. 'y') Ifgrf = .true.\n      END IF\n", "main end if")
    # INTSEL: 80-character names, headless block, skip EUROM.FIL write
    m = rep(m, "      CHARACTER*12 inpt,FIL1,FIL2,FIL3,FIL4,FIL5,fil6,FIL7,dtitl,BLANK\n",
        "      CHARACTER*80 inpt,FIL1,FIL2,FIL3,FIL4,FIL5,fil6,FIL7,dtitl,BLANK\n" + DECL, "intsel decl")
    m = rep(m, "C  Read the remembered names from the last run, if EUROM.FIL has any.\n",
        "C  Headless run: names come from eurosem.ctl (see HEADLS); no prompts\n"
        "C  and no EUROM.FIL.\n      IF(LHEAD) THEN\n        FIL4 = HFIL(1)\n"
        "        FIL5 = HFIL(2)\n        FIL1 = HFIL(3)\n        FIL2 = HFIL(4)\n"
        "        FIL6 = HFIL(5)\n        GO TO 3\n      END IF\n"
        "C  Read the remembered names from the last run, if EUROM.FIL has any.\n", "intsel headless")
    m = rep(m, "c                               FOR NEXT TIME USE:\n      REWIND 9\n",
        "c                               FOR NEXT TIME USE:\n      IF(LHEAD) RETURN\n      REWIND 9\n", "intsel return")
    # PAREAD: three questions
    m = rep(m, "      SUBROUTINE PAREAD\nC\nC              Reads data from file IREAD (formerly NAMELIST variables)\nC\n",
        "      SUBROUTINE PAREAD\nC\nC              Reads data from file IREAD (formerly NAMELIST variables)\nC\n" + DECL, "pareaad decl")
    m = rep(m, "        Write(*,701) change,J\n        Read(*,702) answ\n",
        "        IF(LHEAD) THEN\n          answ = HANS\n        ELSE\n"
        "          Write(*,701) change,J\n          Read(*,702) answ\n        END IF\n", "roc q")
    m = rep(m, "         Write(*,703) change,J\n         Read(*,702) answ\n",
        "         IF(LHEAD) THEN\n           answ = HANS\n         ELSE\n"
        "           Write(*,703) change,J\n           Read(*,702) answ\n         END IF\n", "pave q")
    m = rep(m, "         Write(*,807) j,spltex(j)\n",
        "         IF(LHEAD) THEN\n           answ = HANS\n         ELSE\n"
        "         Write(*,807) j,spltex(j)\n", "spltex q1")
    m = rep(m, "         Read(*,'(a1)') answ\n",
        "         Read(*,'(a1)') answ\n         END IF\n", "spltex q2")
    t["emainl"] = m

    p = t["eplanel"]
    p = rep(p, "      Character*1 line(4)   !, answer\n",
        "      Character*1 line(4)   !, answer\n" + DECL, "plane decl")
    p = rep(p, "          Write(*,*)' RETURN to continue:'\n          Read(*,*)\n",
        "          If(.not.LHEAD) Then\n            Write(*,*)' RETURN to continue:'\n"
        "            Read(*,*)\n          End If\n", "courant")
    t["eplanel"] = p


# ---------------------------------------------------------------- stage 3
def stage3(t):
    p = t["eplanel"]
    p = rep(p, "        JJ=NU(J)\n        SUMA(J) = SUMA(JJ) + ( W(J)*XL(J) )\n",
        "C  MB and NI are set here, as in the ELSE branch, so they do not\n"
        "C  depend on a previous element having been computed.\n"
        "        MB=IFIX(TFIN/DELT)+2\n        NI=MB-1\n"
        "        IF (MB.GT.LENQ) GO TO 370\n"
        "        JJ=NU(J)\n        SUMA(J) = SUMA(JJ) + ( W(J)*XL(J) )\n", "mb fix")
    # update the explanatory note (comment lines only)
    a = p.index("C  NOTE: in this branch MB is used")
    b = p.index("C  be safer to set MB = NI + 1 here.") + len("C  be safer to set MB = NI + 1 here.")
    p = p[:a] + ("C  NOTE: MB and NI are now set at the top of this branch (they\n"
                 "C  used to be left over from the previous element, which only worked\n"
                 "C  because elements are computed in upstream-first order).") + p[b:]
    t["eplanel"] = p


# ---------------------------------------------------------------- stage 4
def stage4(t):
    """RVAL (rock-content factor on soil storage) becomes one value per
    element, RVAL(60), at the end of COMMON /INFIL/.  Before, a single
    scalar was reset for every element read, so only the last element's
    ROC took effect.  Runs where ROC = 0 on every element are unchanged."""
    for f, want in (("emainl", 6), ("eplanel", 1), ("eurozl", 2)):
        n = len(re.findall(r"stone\(60\), RVAL\n", t[f]))
        if n != want:
            sys.exit(f"RVAL declarations in {f}: found {n}, expected {want}")
        t[f] = re.sub(r"stone\(60\), RVAL\n", "stone(60), RVAL(60)\n", t[f])
    m = t["emainl"]
    m = rep(m, "      RVAL = 1.0\n      If(ROC(j) .gt. 0.0 ) Then",
            "      RVAL(J) = 1.0\n      If(ROC(j) .gt. 0.0 ) Then", "rval init")
    m = rep(m, "Y') RVAL = (1.- ROC(J))", "Y') RVAL(J) = (1.- ROC(J))", "rval set")
    a = m.index("C  that fraction.  RVAL = 1 - ROC if accepted, else 1.")
    b = m.index("C  the infiltration parameters are formed).  Check this is intended.") \
        + len("C  the infiltration parameters are formed).  Check this is intended.")
    m = m[:a] + ("C  that fraction.  RVAL(J) = 1 - ROC(J) if accepted, else 1.\n"
                 "C  NOTE: RVAL is now one value per element (it used to be a single\n"
                 "C  scalar, so only the last element's ROC took effect).  It is used\n"
                 "C  in eurozl (XPLINF) when the infiltration parameters are formed.") + m[b:]
    t["emainl"] = m
    z = t["eurozl"]
    z = rep(z, "          APAR(K,M) = AL(J)*(thmx(J)-thi(J))*RVAL\n",
               "          APAR(K,M) = AL(J)*(thmx(J)-thi(J))*RVAL(J)\n", "z156")
    z = rep(z, "          APAR(I,M) = AL(J)*(thmx(J)-THIN(I,M))*RVAL\n",
               "          APAR(I,M) = AL(J)*(thmx(J)-THIN(I,M))*RVAL(J)\n", "z374")
    a = z.index("C    NOTE: RVAL, the rock-content factor for the profile, is ONE value")
    b = z.index("C    commented out further down.") + len("C    commented out further down.")
    z = z[:a] + ("C    NOTE: RVAL(J) = 1 - ROC(J) is the rock-content factor for element\n"
                 "C    J.  It is set per element in PAREAD (emainl).  Before stage 4 it\n"
                 "C    was one value for the whole run, so only the last element's ROC\n"
                 "C    took effect.") + z[b:]
    t["eurozl"] = z


def write(stage_dir, t, ehless):
    os.makedirs(os.path.join(HERE, stage_dir), exist_ok=True)
    for k, v in t.items():
        open(os.path.join(HERE, stage_dir, k + "_commented.for"), "w", encoding="latin-1", newline="\n").write(v)
    open(os.path.join(HERE, stage_dir, "ehless.for"), "w", encoding="latin-1", newline="\n").write(ehless)


t = {k: load(k) for k in FILES}
stage1(t);  write("stage1_fatal", t, EHLESS_1)
stage2(t);  write("stage2_headless", t, EHLESS_1 + EHLESS_2)
stage3(t);  write("stage3_mbfix", t, EHLESS_1 + EHLESS_2)
stage4(t);  write("stage4_rval", t, EHLESS_1 + EHLESS_2)
print("built stage1_fatal, stage2_headless, stage3_mbfix, stage4_rval")

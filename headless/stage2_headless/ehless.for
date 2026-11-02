C=======================================================================
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
C-----------------------------------------------------------------------
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

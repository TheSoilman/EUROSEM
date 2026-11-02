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

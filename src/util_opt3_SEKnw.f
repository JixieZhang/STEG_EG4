********************************************************************************
* This file contains utility routines (mostly integration related) to be
* used in conjunction with the polarized radiative correction program called
* RCSLACPOL, written for experiments E142, E143, E154, and E155.
* 2/94, LMS.
********************************************************************************
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,
     >                  DERIVS,RKQC)
*****************************************************************************
* Runge-Kutta driver program with adaptive stepsize control. Integrate 
* the NVAR starting values YSTART from X1 to X2 with accuracy EPS, storing
* immediate results in the common block /PATH/. H1 should be set as a guessed
* first stepsize, and HMIN as the minimum allowed stepsize (can be zero).
* On output NOK and NBAD are the number of good and bad (but retried and
* fixed) steps taken, and YSTART is replaced by values at the end of the
* integration interval. DERIVS is the user-supplied subroutine for 
* calculating the right-hand side derivative, while RKQC is the name of the
* stepper routineto be used. PATH contains its own information about how
* often an intermediate value is to be stored.
*****************************************************************************
      IMPLICIT NONE

      INTEGER MAXSTP, NMAX
      REAL*8 TWO, ZERO, TINY
      PARAMETER (MAXSTP = 10000, NMAX = 10,TWO = 2.D0,
     >           ZERO = 0.D0, TINY = 1.D-30)

      INTEGER KMAX, NVAR, NOK, NBAD, KOUNT, I, NSTP
      REAL*8 YSTART(NVAR), X1, X2, EPS, H1, HMIN, DXSAV, XP(200),
     >       YP(10,200), YSCAL(NMAX), Y(NMAX), DYDX(NMAX), H,
     >       X, XSAV, HDID, HNEXT
      EXTERNAL DERIVS
      COMMON /PATH/ KMAX, KOUNT, DXSAV, XP, YP
     
      X = X1
      H = SIGN(H1,X2-X1)
      NOK = 0
      NBAD = 0
      KOUNT = 0
      DO I = 1,NVAR
        Y(I) = YSTART(I)
      ENDDO
      XSAV = X - DXSAV*TWO
      DO NSTP = 1,MAXSTP
         CALL OPT_DERIVS(X,Y,DYDX,DERIVS) !kp: Defined below in this file itself.
        DO I = 1,NVAR
          YSCAL(I) = ABS(Y(I)) + ABS(H*DYDX(I)) + TINY
        ENDDO
        IF(KMAX.GT.0)THEN
          IF(ABS(X - XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX - 1)THEN
              KOUNT = KOUNT + 1
              XP(KOUNT) = X
              DO I = 1,NVAR
                YP(I,KOUNT) = Y(I)
              ENDDO
              XSAV = X
            ENDIF
          ENDIF
        ENDIF
        IF((X + H - X2)*(X + H - X1).GT.ZERO) H = X2 - X
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK = NOK + 1
        ELSE
          NBAD = NBAD + 1
        ENDIF
        IF((X - X2)*(X2 - X1).GE.ZERO)THEN
          DO 14 I = 1,NVAR
            YSTART(I) = Y(I)
14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT = KOUNT+1
            XP(KOUNT) = X
            DO 15 I = 1,NVAR
              YP(I,KOUNT) = Y(I) 
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) PRINT*, 'Stepsize < min in ODEINT'
        H = HNEXT
      ENDDO
      PRINT*, 'Too many steps.'
      RETURN
      END
********************************************************************************
      SUBROUTINE ODEINT_INT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,
     >                  DERIVS,RKQC_INT)
*****************************************************************************
* This routine is the same as ODEINT, but in can be called internally
* to ODEINT for multiple integration purposes.
*
* Runge-Kutta driver program with adaptive stepsize control. Integrate 
* the NVAR starting values YSTART from X1 to X2 with accuracy EPS, storing
* immediate results in the common block /PATH/. H1 should be set as a guessed
* first stepsize, and HMIN as the minimum allowed stepsize (can be zero).
* On output NOK and NBAD are the number of good and bad (but retried and
* fixed) steps taken, and YSTART is replaced by values at the end of the
* integration interval. DERIVS is the user-supplied subroutine for 
* calculating the right-hand side derivative, while RKQC_INT is the name of the
* stepper routineto be used. PATH contains its own information about how
* often an intermediate value is to be stored.
*****************************************************************************
      IMPLICIT NONE

      INTEGER MAXSTP, NMAX
      REAL*8 TWO, ZERO, TINY
      PARAMETER (MAXSTP = 10000, NMAX = 10,TWO = 2.D0,
     >           ZERO = 0.D0, TINY = 1.D-30)

      INTEGER KMAX, NVAR, NOK, NBAD, KOUNT, I, NSTP
      REAL*8 YSTART(NVAR), X1, X2, EPS, H1, HMIN, DXSAV, XP(200),
     >       YP(10,200), YSCAL(NMAX), Y(NMAX), DYDX(NMAX), H,
     >       X, XSAV, HDID, HNEXT
      EXTERNAL DERIVS
      COMMON /PATH_INT/ KMAX, KOUNT, DXSAV, XP, YP
     
      X = X1
      H = SIGN(H1,X2-X1)
      NOK = 0
      NBAD = 0
      KOUNT = 0
      DO I = 1,NVAR
        Y(I) = YSTART(I)
      ENDDO
      XSAV = X - DXSAV*TWO
      DO NSTP = 1,MAXSTP
        CALL OPT_DERIVS_INT(X,Y,DYDX,DERIVS)
        DO I = 1,NVAR
          YSCAL(I) = ABS(Y(I)) + ABS(H*DYDX(I)) + TINY
        ENDDO
        IF(KMAX.GT.0)THEN
          IF(ABS(X - XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX - 1)THEN
              KOUNT = KOUNT + 1
              XP(KOUNT) = X
              DO I = 1,NVAR
                YP(I,KOUNT) = Y(I)
              ENDDO
              XSAV = X
            ENDIF
          ENDIF
        ENDIF
        IF((X + H - X2)*(X + H - X1).GT.ZERO) H = X2 - X
        CALL RKQC_INT(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK = NOK + 1
        ELSE
          NBAD = NBAD + 1
        ENDIF
        IF((X - X2)*(X2 - X1).GE.ZERO)THEN
          DO 14 I = 1,NVAR
            YSTART(I) = Y(I)
14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT = KOUNT+1
            XP(KOUNT) = X
            DO 15 I = 1,NVAR
              YP(I,KOUNT) = Y(I) 
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) PRINT*, 'Stepsize < min in ODEINT_INT'
        H = HNEXT
      ENDDO
      PRINT*, 'Too many steps.'
      RETURN
      END
********************************************************************************
      SUBROUTINE ODEINT_EX1(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,
     >                  DERIVS,RKQC_EX1)
*****************************************************************************
* This routine is the same as ODEINT, but in can be called internally
* to ODEINT for multiple integration purposes.
*
* Runge-Kutta driver program with adaptive stepsize control. Integrate 
* the NVAR starting values YSTART from X1 to X2 with accuracy EPS, storing
* immediate results in the common block /PATH/. H1 should be set as a guessed
* first stepsize, and HMIN as the minimum allowed stepsize (can be zero).
* On output NOK and NBAD are the number of good and bad (but retried and
* fixed) steps taken, and YSTART is replaced by values at the end of the
* integration interval. DERIVS is the user-supplied subroutine for 
* calculating the right-hand side derivative, while RKQC_EX1 is the name of the
* stepper routineto be used. PATH contains its own information about how
* often an intermediate value is to be stored.
*****************************************************************************
      IMPLICIT NONE

      INTEGER MAXSTP, NMAX
      REAL*8 TWO, ZERO, TINY
      PARAMETER (MAXSTP = 10000, NMAX = 10,TWO = 2.D0,
     >           ZERO = 0.D0, TINY = 1.D-30)

      INTEGER KMAX, NVAR, NOK, NBAD, KOUNT, I, NSTP
      REAL*8 YSTART(NVAR), X1, X2, EPS, H1, HMIN, DXSAV, XP(200),
     >       YP(10,200), YSCAL(NMAX), Y(NMAX), DYDX(NMAX), H,
     >       X, XSAV, HDID, HNEXT
      EXTERNAL DERIVS
      COMMON /PATH_EX1/ KMAX, KOUNT, DXSAV, XP, YP
     
      X = X1
      H = SIGN(H1,X2-X1)
      NOK = 0
      NBAD = 0
      KOUNT = 0
      DO I = 1,NVAR
        Y(I) = YSTART(I)
      ENDDO
      XSAV = X - DXSAV*TWO
      DO NSTP = 1,MAXSTP
        CALL OPT_DERIVS_EXT1(X,Y,DYDX,DERIVS)
        DO I = 1,NVAR
          YSCAL(I) = ABS(Y(I)) + ABS(H*DYDX(I)) + TINY
        ENDDO
        IF(KMAX.GT.0)THEN
          IF(ABS(X - XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX - 1)THEN
              KOUNT = KOUNT + 1
              XP(KOUNT) = X
              DO I = 1,NVAR
                YP(I,KOUNT) = Y(I)
              ENDDO
              XSAV = X
            ENDIF
          ENDIF
        ENDIF
        IF((X + H - X2)*(X + H - X1).GT.ZERO) H = X2 - X
        CALL RKQC_EX1(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK = NOK + 1
        ELSE
          NBAD = NBAD + 1
        ENDIF
        IF((X - X2)*(X2 - X1).GE.ZERO)THEN
          DO 14 I = 1,NVAR
            YSTART(I) = Y(I)
14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT = KOUNT+1
            XP(KOUNT) = X
            DO 15 I = 1,NVAR
              YP(I,KOUNT) = Y(I) 
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) PRINT*, 'Stepsize < min in ODEINT_EX1'
        H = HNEXT
      ENDDO
      PRINT*, 'Too many steps.'
      RETURN
      END
********************************************************************************
      SUBROUTINE ODEINT_EX2(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,
     >                  DERIVS,RKQC_EX2)
*****************************************************************************
* This routine is the same as ODEINT, but in can be called internally
* to ODEINT for multiple integration purposes.
*
* Runge-Kutta driver program with adaptive stepsize control. Integrate 
* the NVAR starting values YSTART from X1 to X2 with accuracy EPS, storing
* immediate results in the common block /PATH/. H1 should be set as a guessed
* first stepsize, and HMIN as the minimum allowed stepsize (can be zero).
* On output NOK and NBAD are the number of good and bad (but retried and
* fixed) steps taken, and YSTART is replaced by values at the end of the
* integration interval. DERIVS is the user-supplied subroutine for 
* calculating the right-hand side derivative, while RKQC_EX2 is the name of the
* stepper routineto be used. PATH contains its own information about how
* often an intermediate value is to be stored.
*****************************************************************************
      IMPLICIT NONE

      INTEGER MAXSTP, NMAX
      REAL*8 TWO, ZERO, TINY
      PARAMETER (MAXSTP = 10000, NMAX = 10,TWO = 2.D0,
     >           ZERO = 0.D0, TINY = 1.D-30)

      INTEGER KMAX, NVAR, NOK, NBAD, KOUNT, I, NSTP
      REAL*8 YSTART(NVAR), X1, X2, EPS, H1, HMIN, DXSAV, XP(200),
     >       YP(10,200), YSCAL(NMAX), Y(NMAX), DYDX(NMAX), H,
     >       X, XSAV, HDID, HNEXT
      EXTERNAL DERIVS
      COMMON /PATH_EX2/ KMAX, KOUNT, DXSAV, XP, YP
     
      X = X1
      H = SIGN(H1,X2-X1)
      NOK = 0
      NBAD = 0
      KOUNT = 0
      DO I = 1,NVAR
        Y(I) = YSTART(I)
      ENDDO
      XSAV = X - DXSAV*TWO
      DO NSTP = 1,MAXSTP
        CALL OPT_DERIVS_EXT2(X,Y,DYDX,DERIVS)
        DO I = 1,NVAR
          YSCAL(I) = ABS(Y(I)) + ABS(H*DYDX(I)) + TINY
        ENDDO
        IF(KMAX.GT.0)THEN
          IF(ABS(X - XSAV).GT.ABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX - 1)THEN
              KOUNT = KOUNT + 1
              XP(KOUNT) = X
              DO I = 1,NVAR
                YP(I,KOUNT) = Y(I)
              ENDDO
              XSAV = X
            ENDIF
          ENDIF
        ENDIF
        IF((X + H - X2)*(X + H - X1).GT.ZERO) H = X2 - X
        CALL RKQC_EX2(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK = NOK + 1
        ELSE
          NBAD = NBAD + 1
        ENDIF
        IF((X - X2)*(X2 - X1).GE.ZERO)THEN
          DO 14 I = 1,NVAR
            YSTART(I) = Y(I)
14        CONTINUE
          IF(KMAX.NE.0)THEN
            KOUNT = KOUNT+1
            XP(KOUNT) = X
            DO 15 I = 1,NVAR
              YP(I,KOUNT) = Y(I) 
15          CONTINUE
          ENDIF
          RETURN
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) PRINT*, 'Stepsize < min in ODEINT_EX2'
        H = HNEXT
      ENDDO
      PRINT*, 'Too many steps.'
      RETURN
      END
********************************************************************************
      SUBROUTINE RKQC(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
****************************************************************************
* Fifth-order Runge-Kutta step with monitoring of local truncation error to
* ensure accuracy and adjust stepsize. Input are the dependent variable Y of
* length N and its derivative DYDX at the starting value of the independent
* variable X. Also input are the stepsize to be attempted HTRY, the required
* accuracy EPS, and the vector YSCAL against which the error is scaled.
* On output, Y and X are replaced by their new values, HDID is the stepsize
* which was actually accomlished, and HNEXT is the estimated next stepsize.
* DERIVS is the user-supplied subroutine that computes the right-hand
* side derivatives.
****************************************************************************
      IMPLICIT NONE
      INTEGER NMAX
      REAL*8 FCOR, ONE, SAFETY, ERRCON
      PARAMETER (NMAX = 10,FCOR = 0.0666666667D0,
     >           ONE = 1.D0, SAFETY = 0.9D0, ERRCON = 6.D-4)
      EXTERNAL DERIVS

      INTEGER N, I
      REAL*8 Y(N), DYDX(N), X, HTRY, EPS, YSCAL(N), HDID, HNEXT,
     >       YTEMP(NMAX), YSAV(NMAX), DYSAV(NMAX), PGROW, PSHRNK,
     >       XSAV, H, HH, ERRMAX

      PGROW = -0.2D0
      PSHRNK = -0.25D0
      XSAV = X
      DO I = 1,N
        YSAV(I) = Y(I)
        DYSAV(I) = DYDX(I)
      ENDDO
      H = HTRY
1     HH = 0.5D0*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X = XSAV + HH
      CALL OPT_DERIVS(X,YTEMP,DYDX,DERIVS)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X = XSAV + H
      IF(X.EQ.XSAV) PRINT*, 'Stepsize not significant in RKQC'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX = 0.D0
      DO I = 1,N
        YTEMP(I) = Y(I) - YTEMP(I)
        ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
      ENDDO
      ERRMAX = ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H = SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID = H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT = SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT = 4.D0*H
        ENDIF
      ENDIF
      DO I = 1,N
        Y(I) = Y(I) + YTEMP(I)*FCOR
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE RKQC_INT(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
****************************************************************************
* This routine is the same as RKQC, but it can be called internally in
* order to do multiple integration.
*
* Fifth-order Runge-Kutta step with monitoring of local truncation error to
* ensure accuracy and adjust stepsize. Input are the dependent variable Y of
* length N and its derivative DYDX at the starting value of the independent
* variable X. Also input are the stepsize to be attempted HTRY, the required
* accuracy EPS, and the vector YSCAL against which the error is scaled.
* On output, Y and X are replaced by their new values, HDID is the stepsize
* which was actually accomlished, and HNEXT is the estimated next stepsize.
* DERIVS is the user-supplied subroutine that computes the right-hand
* side derivatives.
****************************************************************************
      IMPLICIT NONE
      INTEGER NMAX
      REAL*8 FCOR, ONE, SAFETY, ERRCON
      PARAMETER (NMAX = 10,FCOR = 0.0666666667D0,
     >           ONE = 1.D0, SAFETY = 0.9D0, ERRCON = 6.D-4)
      EXTERNAL DERIVS

      INTEGER N, I
      REAL*8 Y(N), DYDX(N), X, HTRY, EPS, YSCAL(N), HDID, HNEXT,
     >       YTEMP(NMAX), YSAV(NMAX), DYSAV(NMAX), PGROW, PSHRNK,
     >       XSAV, H, HH, ERRMAX

      PGROW = -0.2D0
      PSHRNK = -0.25D0
      XSAV = X
      DO I = 1,N
        YSAV(I) = Y(I)
        DYSAV(I) = DYDX(I)
      ENDDO
      H = HTRY
1     HH = 0.5D0*H
      CALL RK4_INT(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X = XSAV + HH
      CALL OPT_DERIVS_INT(X,YTEMP,DYDX,DERIVS)
      CALL RK4_INT(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X = XSAV + H
      IF(X.EQ.XSAV) PRINT*, 'Stepsize not significant in RKQC_INT'
      CALL RK4_INT(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX = 0.D0
      DO I = 1,N
        YTEMP(I) = Y(I) - YTEMP(I)
        ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
      ENDDO
      ERRMAX = ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H = SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID = H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT = SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT = 4.D0*H
        ENDIF
      ENDIF
      DO I = 1,N
        Y(I) = Y(I) + YTEMP(I)*FCOR
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE RKQC_EX1(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
****************************************************************************
* This routine is the same as RKQC, but it can be called internally in
* order to do multiple integration.
*
* Fifth-order Runge-Kutta step with monitoring of local truncation error to
* ensure accuracy and adjust stepsize. Input are the dependent variable Y of
* length N and its derivative DYDX at the starting value of the independent
* variable X. Also input are the stepsize to be attempted HTRY, the required
* accuracy EPS, and the vector YSCAL against which the error is scaled.
* On output, Y and X are replaced by their new values, HDID is the stepsize
* which was actually accomlished, and HNEXT is the estimated next stepsize.
* DERIVS is the user-supplied subroutine that computes the right-hand
* side derivatives.
****************************************************************************
      IMPLICIT NONE
      INTEGER NMAX
      REAL*8 FCOR, ONE, SAFETY, ERRCON
      PARAMETER (NMAX = 10,FCOR = 0.0666666667D0,
     >           ONE = 1.D0, SAFETY = 0.9D0, ERRCON = 6.D-4)
      EXTERNAL DERIVS

      INTEGER N, I
      REAL*8 Y(N), DYDX(N), X, HTRY, EPS, YSCAL(N), HDID, HNEXT,
     >       YTEMP(NMAX), YSAV(NMAX), DYSAV(NMAX), PGROW, PSHRNK,
     >       XSAV, H, HH, ERRMAX

      PGROW = -0.2D0
      PSHRNK = -0.25D0
      XSAV = X
      DO I = 1,N
        YSAV(I) = Y(I)
        DYSAV(I) = DYDX(I)
      ENDDO
      H = HTRY
1     HH = 0.5D0*H
      CALL RK4_EX1(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X = XSAV + HH
      CALL OPT_DERIVS_EXT1(X,YTEMP,DYDX,DERIVS)
      CALL RK4_EX1(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X = XSAV + H
      IF(X.EQ.XSAV) PRINT*, 'Stepsize not significant in RKQC_EX1'
      CALL RK4_EX1(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX = 0.D0
      DO I = 1,N
        YTEMP(I) = Y(I) - YTEMP(I)
        ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
      ENDDO
      ERRMAX = ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H = SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID = H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT = SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT = 4.D0*H
        ENDIF
      ENDIF
      DO I = 1,N
        Y(I) = Y(I) + YTEMP(I)*FCOR
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE RKQC_EX2(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
****************************************************************************
* This routine is the same as RKQC, but it can be called internally in
* order to do multiple integration.
*
* Fifth-order Runge-Kutta step with monitoring of local truncation error to
* ensure accuracy and adjust stepsize. Input are the dependent variable Y of
* length N and its derivative DYDX at the starting value of the independent
* variable X. Also input are the stepsize to be attempted HTRY, the required
* accuracy EPS, and the vector YSCAL against which the error is scaled.
* On output, Y and X are replaced by their new values, HDID is the stepsize
* which was actually accomlished, and HNEXT is the estimated next stepsize.
* DERIVS is the user-supplied subroutine that computes the right-hand
* side derivatives.
*
* NOTE 4/2000 frw  function OPT_DERIVS_EXT2 was called without argument
*                  DERIVS which is used everywhere else   WHY?
*                  fixed, but no answer (old code?)
*
****************************************************************************
      IMPLICIT NONE
      INTEGER NMAX
      REAL*8 FCOR, ONE, SAFETY, ERRCON
      PARAMETER (NMAX = 10,FCOR = 0.0666666667D0,
     >           ONE = 1.D0, SAFETY = 0.9D0, ERRCON = 6.D-4)
      EXTERNAL DERIVS

      INTEGER N, I
      REAL*8 Y(N), DYDX(N), X, HTRY, EPS, YSCAL(N), HDID, HNEXT,
     >       YTEMP(NMAX), YSAV(NMAX), DYSAV(NMAX), PGROW, PSHRNK,
     >       XSAV, H, HH, ERRMAX

      PGROW = -0.2D0
      PSHRNK = -0.25D0
      XSAV = X
      DO I = 1,N
        YSAV(I) = Y(I)
        DYSAV(I) = DYDX(I)
      ENDDO
      H = HTRY
1     HH = 0.5D0*H
      CALL RK4_EX2(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X = XSAV + HH
      CALL OPT_DERIVS_EXT2(X,YTEMP,DYDX,DERIVS)   ! note applies HERE
      CALL RK4_EX2(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X = XSAV + H
      IF(X.EQ.XSAV) PRINT*, 'Stepsize not significant in RKQC_EX2'
      CALL RK4_EX2(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX = 0.D0
      DO I = 1,N
        YTEMP(I) = Y(I) - YTEMP(I)
        ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
      ENDDO
      ERRMAX = ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H = SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID = H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT = SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT = 4.D0*H
        ENDIF
      ENDIF
      DO I = 1,N
        Y(I) = Y(I) + YTEMP(I)*FCOR
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
****************************************************************************
* Given values for N variables Y and their derivatives DYDX known at X, use 
* the fourth order Runge-Kutta method to advance the solution over an
* interval H and return the incremented variables as YOUT, which need not 
* be a distinct array from Y. The user supplies the subroutine 
* DERIVS(X,Y,DYDX) which returns derivatives DYDX at X.
****************************************************************************
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX = 10)

      INTEGER N, I
      REAL*8 Y(N), DYDX(N), YOUT(N), YT(NMAX), DYT(NMAX), DYM(NMAX),
     >       X, H, HH, H6, XH
      EXTERNAL DERIVS      

      HH = H*0.5D0
      H6 = H/6.D0
      XH = X + HH
      DO I = 1,N
        YT(I) = Y(I) + HH*DYDX(I)
      ENDDO
      CALL OPT_DERIVS(XH,YT,DYT,DERIVS)
      DO I = 1,N
        YT(I) = Y(I) + HH*DYT(I)
      ENDDO
      CALL OPT_DERIVS(XH,YT,DYM,DERIVS)
      DO I = 1,N
        YT(I) = Y(I) + H*DYM(I)
        DYM(I) = DYT(I) + DYM(I)
      ENDDO
      CALL OPT_DERIVS(X + H,YT,DYT,DERIVS)
      DO I = 1,N
        YOUT(I) = Y(I) + H6*(DYDX(I) + DYT(I) + 2.D0*DYM(I))
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE RK4_INT(Y,DYDX,N,X,H,YOUT,DERIVS)
****************************************************************************
* This subroutine is the same as RK4, but it can be called internally for
* multiple integration purposes.
*
* Given values for N variables Y and their derivatives DYDX known at X, use 
* the fourth order Runge-Kutta method to advance the solution over an
* interval H and return the incremented variables as YOUT, which need not 
* be a distinct array from Y. The user supplies the subroutine 
* DERIVS(X,Y,DYDX) which returns derivatives DYDX at X.
****************************************************************************
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX = 10)

      INTEGER N, I
      REAL*8 Y(N), DYDX(N), YOUT(N), YT(NMAX), DYT(NMAX), DYM(NMAX),
     >       X, H, HH, H6, XH
      EXTERNAL DERIVS      


      HH = H*0.5D0
      H6 = H/6.D0
      XH = X + HH
      DO I = 1,N
        YT(I) = Y(I) + HH*DYDX(I)
      ENDDO
      CALL OPT_DERIVS_INT(XH,YT,DYT,DERIVS)
      DO I = 1,N
        YT(I) = Y(I) + HH*DYT(I)
      ENDDO
      CALL OPT_DERIVS_INT(XH,YT,DYM,DERIVS)
      DO I = 1,N
        YT(I) = Y(I) + H*DYM(I)
        DYM(I) = DYT(I) + DYM(I)
      ENDDO
      CALL OPT_DERIVS_INT(X + H,YT,DYT,DERIVS)
      DO I = 1,N
        YOUT(I) = Y(I) + H6*(DYDX(I) + DYT(I) + 2.D0*DYM(I))
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE RK4_EX1(Y,DYDX,N,X,H,YOUT,DERIVS)
****************************************************************************
* This subroutine is the same as RK4, but it can be called internally for
* multiple integration purposes.
*
* Given values for N variables Y and their derivatives DYDX known at X, use 
* the fourth order Runge-Kutta method to advance the solution over an
* interval H and return the incremented variables as YOUT, which need not 
* be a distinct array from Y. The user supplies the subroutine 
* DERIVS(X,Y,DYDX) which returns derivatives DYDX at X.
****************************************************************************
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX = 10)

      INTEGER N, I
      REAL*8 Y(N), DYDX(N), YOUT(N), YT(NMAX), DYT(NMAX), DYM(NMAX),
     >       X, H, HH, H6, XH
      EXTERNAL DERIVS      


      HH = H*0.5D0
      H6 = H/6.D0
      XH = X + HH
      DO I = 1,N
        YT(I) = Y(I) + HH*DYDX(I)
      ENDDO
      CALL OPT_DERIVS_EXT1(XH,YT,DYT,DERIVS)
      DO I = 1,N
        YT(I) = Y(I) + HH*DYT(I)
      ENDDO
      CALL OPT_DERIVS_EXT1(XH,YT,DYM,DERIVS)
      DO I = 1,N
        YT(I) = Y(I) + H*DYM(I)
        DYM(I) = DYT(I) + DYM(I)
      ENDDO
      CALL OPT_DERIVS_EXT1(X + H,YT,DYT,DERIVS)
      DO I = 1,N
        YOUT(I) = Y(I) + H6*(DYDX(I) + DYT(I) + 2.D0*DYM(I))
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE RK4_EX2(Y,DYDX,N,X,H,YOUT,DERIVS)
****************************************************************************
* This subroutine is the same as RK4, but it can be called internally for
* multiple integration purposes.
*
* Given values for N variables Y and their derivatives DYDX known at X, use 
* the fourth order Runge-Kutta method to advance the solution over an
* interval H and return the incremented variables as YOUT, which need not 
* be a distinct array from Y. The user supplies the subroutine 
* DERIVS(X,Y,DYDX) which returns derivatives DYDX at X.
****************************************************************************
      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX = 10)

      INTEGER N, I
      REAL*8 Y(N), DYDX(N), YOUT(N), YT(NMAX), DYT(NMAX), DYM(NMAX),
     >       X, H, HH, H6, XH
      EXTERNAL DERIVS      

      HH = H*0.5D0
      H6 = H/6.D0
      XH = X + HH
      DO I = 1,N
        YT(I) = Y(I) + HH*DYDX(I)
      ENDDO
      CALL OPT_DERIVS_EXT2(XH,YT,DYT,DERIVS)
      DO I = 1,N
        YT(I) = Y(I) + HH*DYT(I)
      ENDDO
      CALL OPT_DERIVS_EXT2(XH,YT,DYM,DERIVS)
      DO I = 1,N
        YT(I) = Y(I) + H*DYM(I)
        DYM(I) = DYT(I) + DYM(I)
      ENDDO
      CALL OPT_DERIVS_EXT2(X + H,YT,DYT,DERIVS)
      DO I = 1,N
        YOUT(I) = Y(I) + H6*(DYDX(I) + DYT(I) + 2.D0*DYM(I))
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE OPT_DERIVS(X,Y,DYDX,DERIVS)
********************************************************************************
* This routine is designed to speed up the integration routine by checking
* the argument X used in the last ten calls to the user supplied routine
* DERIVS. If the current argument is identical to one of the past ten calls
* to DERIVS return the saved derivatives DYDX at X.
* - JNF, 14Apr94
******************************************************************************* 
      IMPLICIT NONE
      INTEGER*4 NMAX
      PARAMETER (NMAX = 10)

      INTEGER*4 I, ILAST/0/
      REAL*8 X, Y, DYDX, XLAST(NMAX), DYDXLAST(NMAX)
      EXTERNAL DERIVS
      COMMON/ LATEST/ XLAST ,DYDXLAST
                
      DO I = 1,NMAX 
        IF(X.EQ.XLAST(I)) THEN
          DYDX = DYDXLAST(I)
          RETURN
        ENDIF
      ENDDO
      CALL DERIVS(X,Y,DYDX)
      ILAST = ILAST + 1
      IF(ILAST.EQ.11) ILAST = 1
      XLAST(ILAST) = X
      DYDXLAST(ILAST) = DYDX
      RETURN
      END	
********************************************************************************
      SUBROUTINE OPT_DERIVS_INT(X,Y,DYDX,DERIVS)
********************************************************************************
* This routine is designed to speed up the integration routine by checking
* the argument X used in the last ten calls to the user supplied routine
* DERIVS. If the current argument is identical to one of the past ten calls
* to DERIVS return the saved derivatives DYDX at X.
* - JNF, 14Apr94
******************************************************************************* 
      IMPLICIT NONE
      INTEGER*4 NMAX
      PARAMETER (NMAX = 10)
      INTEGER*4 I, ILAST/0/
      REAL*8 X, Y, DYDX, XLAST(NMAX), DYDXLAST(NMAX)
      EXTERNAL DERIVS
      COMMON/ LATEST_INT/ XLAST ,DYDXLAST
                
      DO I = 1,NMAX 
        IF(X.EQ.XLAST(I)) THEN
          DYDX = DYDXLAST(I)
          RETURN
        ENDIF
      ENDDO
      CALL DERIVS(X,Y,DYDX)
      ILAST = ILAST + 1
      IF(ILAST.EQ.11) ILAST = 1
      XLAST(ILAST) = X
      DYDXLAST(ILAST) = DYDX
      RETURN
      END	
********************************************************************************
      SUBROUTINE OPT_DERIVS_EXT1(X,Y,DYDX,DERIVS)
********************************************************************************
* This routine is designed to speed up the integration routine by checking
* the argument X used in the last ten calls to the user supplied routine
* DERIVS. If the current argument is identical to one of the past ten calls
* to DERIVS return the saved derivatives DYDX at X.
* - JNF, 14Apr94
******************************************************************************* 
      IMPLICIT NONE
      INTEGER*4 NMAX
      PARAMETER (NMAX = 10)

      INTEGER*4 I, ILAST/0/
      REAL*8 X, Y, DYDX, XLAST(NMAX), DYDXLAST(NMAX)
      EXTERNAL DERIVS
      COMMON/ LATEST_EX1/ XLAST ,DYDXLAST
                
      DO I = 1,NMAX 
        IF(X.EQ.XLAST(I)) THEN
          DYDX = DYDXLAST(I)
          RETURN
        ENDIF
      ENDDO
      CALL DERIVS(X,Y,DYDX)
      ILAST = ILAST + 1
      IF(ILAST.EQ.11) ILAST = 1
      XLAST(ILAST) = X
      DYDXLAST(ILAST) = DYDX
      RETURN
      END	
********************************************************************************
      SUBROUTINE OPT_DERIVS_EXT2(X,Y,DYDX,DERIVS)
********************************************************************************
* This routine is designed to speed up the integration routine by checking
* the argument X used in the last ten calls to the user supplied routine
* DERIVS. If the current argument is identical to one of the past ten calls
* to DERIVS return the saved derivatives DYDX at X.
* - JNF, 14Apr94
******************************************************************************* 
      IMPLICIT NONE
      INTEGER*4 NMAX
      PARAMETER (NMAX = 10)

      INTEGER*4 I, ILAST/0/
      REAL*8 X, Y, DYDX, XLAST(NMAX), DYDXLAST(NMAX)
      EXTERNAL DERIVS
      COMMON/ LATEST_EX2/ XLAST ,DYDXLAST
                
      DO I = 1,NMAX 
        IF(X.EQ.XLAST(I)) THEN
          DYDX = DYDXLAST(I)
          RETURN
        ENDIF
      ENDDO
      CALL DERIVS(X,Y,DYDX)
      ILAST = ILAST + 1
      IF(ILAST.EQ.11) ILAST = 1
      XLAST(ILAST) = X
      DYDXLAST(ILAST) = DYDX
      RETURN
      END	
******************************************************************************

      SUBROUTINE QUADMO1(ANSWER,FUNCT,LOWER,UPPER,EPSLON,NLVL, mtarg)

      INTEGER NLVL, LEVEL, MINLVL/3/, MAXLVL/24/, RETRN(50), I
      REAL*8 ANSWER,LOWER,UPPER,EPSLON
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),
     >   FMRX(50), ESTRX(50), EPSX(50)
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,
     >   AREA, ABAREA,   M, COEF, ROMBRG,   EPS, DUM
      EXTERNAL FUNCT

         LEVEL = 0
         NLVL = 0
         ABAREA = 0.0
         L = LOWER
         R = UPPER
         CALL FUNCT(L,DUM,FL)
         CALL FUNCT(0.5*(L+R),DUM,FM)
         CALL FUNCT(R,DUM,FR)
         EST = 0.0
         EPS = EPSLON
  100 LEVEL = LEVEL+1
      M = 0.5*(L+R)
      COEF = R-L
      IF(COEF.NE.0) GO TO 150
         ROMBRG = EST
         GO TO 300
  150 CALL FUNCT(0.5*(L+M),DUM,FML)
      CALL FUNCT(0.5*(M+R),DUM,FMR)
      ESTL = (FL+4.0*FML+FM)*COEF
      ESTR = (FM+4.0*FMR+FR)*COEF
      ESTINT = ESTL+ESTR
      AREA=DABS(ESTL)+DABS(ESTR)
      ABAREA=AREA+ABAREA-DABS(EST)
      IF(LEVEL.NE.MAXLVL) GO TO 200
         NLVL = NLVL+1
         ROMBRG = ESTINT
         GO TO 300
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.
     1         (LEVEL.LT.MINLVL))  GO TO 400
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0
  300    LEVEL = LEVEL-1
         I = RETRN(LEVEL)
         VALINT(LEVEL, I) = ROMBRG
         GO TO (500, 600), I
  400    RETRN(LEVEL) = 1
         MX(LEVEL) = M
         RX(LEVEL) = R
         FMX(LEVEL) = FM
         FMRX(LEVEL) = FMR
         FRX(LEVEL) = FR
         ESTRX(LEVEL) = ESTR
         EPSX(LEVEL) = EPS
         EPS = EPS/1.4
         R = M
         FR = FM
         FM = FML
         EST = ESTL
         GO TO 100
  500    RETRN(LEVEL) = 2
         L = MX(LEVEL)
         R = RX(LEVEL)
         FL = FMX(LEVEL)
         FM = FMRX(LEVEL)
         FR = FRX(LEVEL)
         EST = ESTRX(LEVEL)
         EPS = EPSX(LEVEL)
         GO TO 100
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)
      IF(LEVEL.GT.1) GO TO 300
      ANSWER = ROMBRG /12.0D0
      RETURN
      END

******************************************************************************

      SUBROUTINE QUADMO2(ANSWER,FUNCT,LOWER,UPPER,EPSLON,NLVL)

      INTEGER NLVL, LEVEL, MINLVL/3/, MAXLVL/24/, RETRN(50), I
      REAL*8 ANSWER, LOWER,UPPER,EPSLON
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),
     >   FMRX(50), ESTRX(50), EPSX(50)
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,
     >   AREA, ABAREA,   M, COEF, ROMBRG,   EPS, DUM
      EXTERNAL FUNCT

         LEVEL = 0
         NLVL = 0
         ABAREA = 0.0
         L = LOWER
         R = UPPER
         CALL FUNCT(L,DUM,FL)
         CALL FUNCT(0.5*(L+R),DUM,FM)
         CALL FUNCT(R,DUM,FR)
         EST = 0.0
         EPS = EPSLON
  100 LEVEL = LEVEL+1
      M = 0.5*(L+R)
      COEF = R-L
      IF(COEF.NE.0) GO TO 150
         ROMBRG = EST
         GO TO 300
  150 CALL FUNCT(0.5*(L+M),DUM,FML)
      CALL FUNCT(0.5*(M+R),DUM,FMR)
      ESTL = (FL+4.0*FML+FM)*COEF
      ESTR = (FM+4.0*FMR+FR)*COEF
      ESTINT = ESTL+ESTR
      AREA=DABS(ESTL)+DABS(ESTR)
      ABAREA=AREA+ABAREA-DABS(EST)
      IF(LEVEL.NE.MAXLVL) GO TO 200
         NLVL = NLVL+1
         ROMBRG = ESTINT
         GO TO 300
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.
     1         (LEVEL.LT.MINLVL))  GO TO 400
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0
  300    LEVEL = LEVEL-1
         I = RETRN(LEVEL)
         VALINT(LEVEL, I) = ROMBRG
         GO TO (500, 600), I
  400    RETRN(LEVEL) = 1
         MX(LEVEL) = M
         RX(LEVEL) = R
         FMX(LEVEL) = FM
         FMRX(LEVEL) = FMR
         FRX(LEVEL) = FR
         ESTRX(LEVEL) = ESTR
         EPSX(LEVEL) = EPS
         EPS = EPS/1.4
         R = M
         FR = FM
         FM = FML
         EST = ESTL
         GO TO 100
  500    RETRN(LEVEL) = 2
         L = MX(LEVEL)
         R = RX(LEVEL)
         FL = FMX(LEVEL)
         FM = FMRX(LEVEL)
         FR = FRX(LEVEL)
         EST = ESTRX(LEVEL)
         EPS = EPSX(LEVEL)
         GO TO 100
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)
      IF(LEVEL.GT.1) GO TO 300
      ANSWER = ROMBRG /12.0D0
      RETURN
      END
******************************************************************************

      SUBROUTINE QUADMO3(ANSWER,FUNCT,LOWER,UPPER,EPSLON,NLVL, mtarg)

      INTEGER NLVL, LEVEL, MINLVL/3/, MAXLVL/24/, RETRN(50), I
      REAL*8 ANSWER, LOWER,UPPER,EPSLON
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),
     >   FMRX(50), ESTRX(50), EPSX(50)
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,
     >   AREA, ABAREA,   M, COEF, ROMBRG,   EPS, DUM
	  
      real*8 mtarg              ! csk
	  
      EXTERNAL FUNCT

         LEVEL = 0
         NLVL = 0
         ABAREA = 0.0
         L = LOWER
         R = UPPER
         CALL FUNCT(L,DUM,FL,mtarg)
         CALL FUNCT(0.5*(L+R),DUM,FM,mtarg)
         CALL FUNCT(R,DUM,FR,mtarg)
         EST = 0.0
         EPS = EPSLON
  100 LEVEL = LEVEL+1
      M = 0.5*(L+R)
      COEF = R-L
      IF(COEF.NE.0) GO TO 150
         ROMBRG = EST
         GO TO 300
  150 CALL FUNCT(0.5*(L+M),DUM,FML,mtarg)
      CALL FUNCT(0.5*(M+R),DUM,FMR,mtarg)
      ESTL = (FL+4.0*FML+FM)*COEF
      ESTR = (FM+4.0*FMR+FR)*COEF
      ESTINT = ESTL+ESTR
      AREA=DABS(ESTL)+DABS(ESTR)
      ABAREA=AREA+ABAREA-DABS(EST)
      IF(LEVEL.NE.MAXLVL) GO TO 200
         NLVL = NLVL+1
         ROMBRG = ESTINT
         GO TO 300
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.
     1         (LEVEL.LT.MINLVL))  GO TO 400
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0
  300    LEVEL = LEVEL-1
         I = RETRN(LEVEL)
         VALINT(LEVEL, I) = ROMBRG
         GO TO (500, 600), I
  400    RETRN(LEVEL) = 1
         MX(LEVEL) = M
         RX(LEVEL) = R
         FMX(LEVEL) = FM
         FMRX(LEVEL) = FMR
         FRX(LEVEL) = FR
         ESTRX(LEVEL) = ESTR
         EPSX(LEVEL) = EPS
         EPS = EPS/1.4
         R = M
         FR = FM
         FM = FML
         EST = ESTL
         GO TO 100
  500    RETRN(LEVEL) = 2
         L = MX(LEVEL)
         R = RX(LEVEL)
         FL = FMX(LEVEL)
         FM = FMRX(LEVEL)
         FR = FRX(LEVEL)
         EST = ESTRX(LEVEL)
         EPS = EPSX(LEVEL)
         GO TO 100
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)
      IF(LEVEL.GT.1) GO TO 300
      ANSWER = ROMBRG /12.0D0
      RETURN
      END
******************************************************************************

      SUBROUTINE QUADMO4(ANSWER,FUNCT,LOWER,UPPER,EPSLON,NLVL,mtarg)

      INTEGER NLVL, LEVEL, MINLVL/3/, MAXLVL/24/, RETRN(50), I
      REAL*8 ANSWER, LOWER,UPPER,EPSLON
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),
     >   FMRX(50), ESTRX(50), EPSX(50)
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,
     >   AREA, ABAREA,   M, COEF, ROMBRG,   EPS, DUM

	  real*8 mtarg ! csk
	  
      EXTERNAL FUNCT

         LEVEL = 0
         NLVL = 0
         ABAREA = 0.0
         L = LOWER
         R = UPPER
         CALL FUNCT(L,DUM,FL,mtarg)
         CALL FUNCT(0.5*(L+R),DUM,FM,mtarg)
         CALL FUNCT(R,DUM,FR,mtarg)
         EST = 0.0
         EPS = EPSLON
  100 LEVEL = LEVEL+1
      M = 0.5*(L+R)
      COEF = R-L
      IF(COEF.NE.0) GO TO 150
         ROMBRG = EST
         GO TO 300
  150 CALL FUNCT(0.5*(L+M),DUM,FML,mtarg)
      CALL FUNCT(0.5*(M+R),DUM,FMR,mtarg)
      ESTL = (FL+4.0*FML+FM)*COEF
      ESTR = (FM+4.0*FMR+FR)*COEF
      ESTINT = ESTL+ESTR
      AREA=DABS(ESTL)+DABS(ESTR)
      ABAREA=AREA+ABAREA-DABS(EST)
      IF(LEVEL.NE.MAXLVL) GO TO 200
         NLVL = NLVL+1
         ROMBRG = ESTINT
         GO TO 300
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.
     1         (LEVEL.LT.MINLVL))  GO TO 400
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0
  300    LEVEL = LEVEL-1
         I = RETRN(LEVEL)
         VALINT(LEVEL, I) = ROMBRG
         GO TO (500, 600), I
  400    RETRN(LEVEL) = 1
         MX(LEVEL) = M
         RX(LEVEL) = R
         FMX(LEVEL) = FM
         FMRX(LEVEL) = FMR
         FRX(LEVEL) = FR
         ESTRX(LEVEL) = ESTR
         EPSX(LEVEL) = EPS
         EPS = EPS/1.4
         R = M
         FR = FM
         FM = FML
         EST = ESTL
         GO TO 100
  500    RETRN(LEVEL) = 2
         L = MX(LEVEL)
         R = RX(LEVEL)
         FL = FMX(LEVEL)
         FM = FMRX(LEVEL)
         FR = FRX(LEVEL)
         EST = ESTRX(LEVEL)
         EPS = EPSX(LEVEL)
         GO TO 100
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)
      IF(LEVEL.GT.1) GO TO 300
      ANSWER = ROMBRG /12.0D0
      RETURN
      END
******************************************************************************

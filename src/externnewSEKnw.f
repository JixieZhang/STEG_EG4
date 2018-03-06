*******************************************************************************
* File: extern.f, contains two external RC calculation programs.
* EXTERN1 is for special case where TA=0.0
* EXTERN2 is for full integration for given TA, and TB
*******************************************************************************


       SUBROUTINE EXTERN1NEW(TA,ANSWER, mtarg)
*******************************************************************************
* This subroutine calculates externally radiated cross sections for the case 
* where the target thickness before scattering (TB) is negligible.
* 
* This is done by calculating a convolution integral over
* dEP where the integrand is:
*    SIGMAR(E0,EP)*I(EP,EPFINAL,TA)
* TA is the total target radiation length (units of rl) after scattering.  
* I is the energy loss probability distribution given by ELOSS_PROB.
* ANSWER is the result of the integration.
* E0 is the beam energy.
* EPFINAL is the final(measured) scattered momentum.
* SIGMAR = internally radiated cross section function.
*
* The limits of integration for E' are EPFINAL to Ep_elastic(E)
* However:
* The probabilty function cannot be evaluated at the point I(E,E,T) because
* it diverges here. So, to avoid integration problems, the integral is
* divided into two parts. The first part is the integral evaluated from
* EPFINAL + DELTA to Ep_elastic(E). The second part is the analytically
* integrated part from EPFINAL to EPFINAL + DELTA.
*
* !!!!! It is assumed that the common KINLIST is filled before this routine 
* !!!!! is called. 
* 4/94, LMS.
*******************************************************************************
       IMPLICIT NONE
       INCLUDE 'instruct.inc'

       INTEGER NOK, NBAD, NLVL
       REAL*8 TA, ANSWER, E0, EPFINAL, THETA, E, Z, TBB, TAA,
     >        DELTA, W2, FZ, BZ, SIN2, PI/3.141592654D0/, IB,
     >        BT, GAMMA, EPS1, PROB_EP, NORM, EPMAX,
     >        EPMIN, PROBTOT_EP, SIGNOM, SIGMAR, I1, I2,
     >        EPSS/5.D-5/, H1/0.01D0/, HMIN/1.0D-30/,
csk     >        MP/0.939D0/, ELOSS_PROB, XCUT, W2CUT ! Changed MP SEK apr12
     >        MP/0.938272D0/, ELOSS_PROB, XCUT, W2CUT
     
     	real*8 mtarg ! csk
     	
       LOGICAL FL_ION/.FALSE./ 
       COMMON /KINLIST/E0, EPFINAL, THETA, E, Z,XCUT
       COMMON /KINLIST2/ DELTA, W2, TBB, TAA
       COMMON /B_ELOSS/ FZ, BZ
       EXTERNAL RKQC_EX1, BINTE, EXTERNI2NEW, EXTERNI2

       TBB = 0.D0
       TAA = TA
       DELTA = 0.01D0         !Nominal
       IF(DELTEST) DELTA = 0.05D0

       SIN2 = (SIN(THETA*PI/360.D0))**2

! Call once to define Z dependnet variables
       IB = ELOSS_PROB(E0,EPFINAL,TA,Z,FL_ION)


! PROB_EP is the analytically calculated integral of I(E',EPFINAL,TA) 
! for E':EPFINAL to EPFINAL+DELTA
       BT = BZ*TA
       EPS1 = DELTA/EPFINAL
       NORM = ( (1.D0 - BT/(1.D0 + BT))*(1.D0 + FZ)
     >        + 3.D0*BT/(4.D0*(2.D0 + BT)))
       if (EPS1 .gt. 0.0D0) then
         PROB_EP = (EPS1**(BT))*( (1.D0 - BT*EPS1/
     >       (1.D0 + BT))*(1.D0 + FZ)  + 3.D0*BT*EPS1*EPS1/
     >       (4.D0*(2.D0 + BT)))
         PROB_EP = PROB_EP/NORM
      else
        PROB_EP = 0.0D0
	  endif

       EPMAX = E0/(1.D0 + (2.D0*E0/mtarg)*SIN2) ! csk
       EPMIN = EPFINAL
       EPS1 = (EPMAX - EPMIN)/EPFINAL
       if (EPS1 .gt. 0.0D0) then
         PROBTOT_EP = (EPS1**(BT))*( (1.D0 - BT*EPS1/
     >       (1.D0 + BT))*(1.D0 + FZ)  + 3.D0*BT*EPS1*EPS1/
     >       (4.D0*(2.D0 + BT)))
         PROBTOT_EP = PROBTOT_EP/NORM
      else
        PROBTOT_EP = 0.0D0
	  endif

! I1 is the analytically evaluated integral between the limits 
! E':EPFINAL to EPFINAL+DELTA. It is estimated for this bin assuming
! that SIGMAR can be evaluated at the nominal kinematics and pulled
! outside the integral.
       SIGNOM = SIGMAR(E0,EPFINAL,THETA)
       I1 = SIGNOM*PROB_EP

! I2 is the integral evaluated between the limits EPFINAL+DELTA to 
! Ep_elastic. 
       EPMAX = E0/(1.D0 + (2.D0*E0/mtarg)*SIN2) !csk
       W2CUT = MP*MP + (1.D0/XCUT - 1.D0)*
     >     (4.D0*E0*E0*SIN2/(1.D0 + 2.D0*E0*SIN2/(MP*XCUT)))
       IF(.NOT.TAIL_ON) EPMAX = (MP*MP + 2.D0*MP*E0 - 
     >        W2CUT)/(2.D0*MP + 4.D0*E0*SIN2)
       EPMIN = EPFINAL + DELTA
       I2 = 0.D0
       IF(EPMAX.GT.EPMIN) CALL QUADMO3(I2,EXTERNI2,EPMAX,EPMIN,
     >                     1.D-3,NLVL, mtarg)
       I2 = -I2

       ANSWER = I1 + I2
c      WRITE(6,'(4(1PE13.5))') I1,I2,SIGNOM,I1+I2

       RETURN
       END
*******************************************************************************


       SUBROUTINE EXTERN2NEW(TB,TA,ANSWER, mtarg)
*******************************************************************************
* This subroutine calculates externally radiated cross sections.
* This is done by calculating a double convolution integral over
* dE and dEP where the integrand is:
*    I(E0,E,TB)*SIGMAR(E,EP)*I(EP,EPFINAL,TA)
* TB is the total target radiation length (units of rl) before scattering.  
* TA is the total target radiation length (units of rl) after scattering.  
* I is the energy loss probability distribution given by ELOSS_PROB.
* ANSWER is the result of the double integration.
* E0 is the initial energy.
* EPFINAL is the final(measured) scattered momentum.
* SIGMAR = internally radiated cross section function.
* This routine is called by EXT1 which decides how to integrate over target.
*
* The limits of integration for E are E_elastic to E0.
* The limits of integration for E' are EPFINAL to Ep_elastic(E)
* However:
* The probabilty function cannot be evaluated at the point I(E,E,T) because
* it diverges here. So, to avoid integration problems, the integral is
* divided into parts. The double integral has integration limits given
* by E: E_elastic to E0 - DELTA and E': EFINAL + DELTA to Ep_elastic(E).
* The remaining pieces are reduced to single integrations, and the total
* probability for the bin of size DELTA is calculated analytically. 
*
* !!!!! It is assumed that the common KINLIST is filled before this routine 
* !!!!! is called. 
* 4/94, LMS.
*******************************************************************************
       IMPLICIT NONE
       INCLUDE 'instruct.inc'

       INTEGER NOK, NBAD, NLVL
       REAL*8 TA, TB, ANSWER, E0, EPFINAL, THETA, 
     >        ALITTLE/1.0D-14/, EPS/3.D-6/, EPSS/1.D-3/, H1/0.01D0/, 
csk     >        HMIN/1.0D-30/, MP/0.939D0/, TBB, TAA, ! changed MP sek apr12
     >        HMIN/1.0D-30/, MP/0.938272D0/, TBB, TAA,
     >        E, Z, DELTA, FZ, EPS1, GAMMA, DE, INTP,
     >        SIGMAR, SIGMARNEW, BT, SIGNOM, TEST, DUM, BZ,
     >        PI/3.141592654D0/, MID, MPITHR/1.1518568D0/,SIN2,
     >        IB, ELOSS_PROB, I1, I2, I3, I4, I3A, I2A, YLIM,
     >        I44, I444, PROB_E, PROB_EP, EMIN, EMAX, EPMIN, EPMAX, 
     >        W2MIN, W2MAX, W2, W2CUT1/1.0D0/, W2CUT2, PROBTOT_E,
     >        PROBTOT_EP, NORM, DELTA2, XCUT, W2CUT, VARMIN, VARMAX
     
     	real*8 mtarg ! csk
     	     
       LOGICAL FL_ION/.FALSE./ 
       COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT
       COMMON /KINLIST2/ DELTA, W2, TBB, TAA
       COMMON /B_ELOSS/ FZ, BZ
       EXTERNAL RKQC_EX1, BINTE, EXTERNI2NEW, EXTERNI3NEW, EXTERNI4NEW,
     >          EXTERNI2

! Define functions for integrating probabilities:
       INTP(YLIM,FZ,BT) = (YLIM**BT)*((1.D0 - 
     >       YLIM*BT/(1.D0 + BT))*(1.D0+FZ) + 
     >       3.D0*YLIM*YLIM*BT/(4.D0*(2.D0 + BT)))

!test
c       TB = 0.00001D0
c       TA = 0.00001D0

       TBB = TB
       TAA = TA
       DELTA = 0.01D0         !Nominal
       IF(DELTEST) DELTA = 0.05D0
!Delta2 used to avoid problems when super close to elastic peak
       DELTA2 = 1.D-15
       DELTA2 = 1.D-5

       SIN2 = (SIN(THETA*PI/360.D0))**2

! Call once to define Z dependnet variables
       IB = ELOSS_PROB(E0,E0/2.D0,TB,Z,FL_ION)

! PROB_E is the analytically calculated integral of I(E0,E,TB) 
! for E:E0-DELTA to E0
       BT = BZ*TB
       EPS1 = DELTA/E0 
       NORM = INTP(1.D0,FZ,BT)
	   if (EPS1 .gt. 0.0D0) then
         PROB_E = INTP(EPS1,FZ,BT)/NORM
       else
         PROB_E = 0.0D0
       endif

       EMAX = E0
       EMIN = EPFINAL/(1.D0 - (2.D0*EPFINAL/mtarg)*SIN2) ! csk
       EPS1 = (EMAX - EMIN)/E0 
	   if (EPS1 .gt. 0.0D0) then
         PROBTOT_E = INTP(EPS1,FZ,BT)/NORM
       else
         PROBTOT_E = 0.0D0
       endif

! PROB_EP is the analytically calculated integral of I(E',EPFINAL,TB) 
! for E':EPFINAL to EPFINAL+DELTA
       BT = BZ*TA
       EPS1 = DELTA/EPFINAL
       NORM = INTP(1.D0,FZ,BT)
	   if (EPS1 .gt. 0.0D0) then
         PROB_EP = INTP(EPS1,FZ,BT)/NORM
       else
         PROB_EP = 0.0D0
       endif

       EPMAX = E0/(1.D0 + (2.D0*E0/mtarg)*SIN2) ! csk
       EPMIN = EPFINAL
       EPS1 = (EPMAX - EPMIN)/EPFINAL
       if (EPS1 .gt. 0.0D0) then
          PROBTOT_EP = INTP(EPS1,FZ,BT)/NORM
       else
          PROBTOT_EP = 0.0D0
       endif

c       print*,'externew2:EPMAX,EPMIN: ',EPMAX,EPMIN
c       print*,'externnew2:EPS1,FZ,BT,NORM: ',EPS1,FZ,BT,NORM
c       print*,'externew2:EPS1**BT: ',EPS1**BT
c       print*,'externew2: INTP(EPS1,FZ,BT) = ',(EPS1**BT)*((1.D0 - 
c     >       EPS1*BT/(1.D0 + BT))*(1.D0+FZ) + 
c     >       3.D0*EPS1*EPS1*BT/(4.D0*(2.D0 + BT)))
       
c     print*,'externnewSED: theta=',theta           !kp: 6/1/12

! I1 is the double integral evaluated between the limits E:E0-DELTA to E0
! and E':EPFINAL to EPFINAL+DELTA. It is estimated for this bin assuming
! that SIGMAR can be evaluated at the nominal kinematics and pulled
! outside the integral.
        SIGNOM = SIGMAR(E0,EPFINAL,THETA)
        I1 = SIGNOM*PROB_E*PROB_EP

! I2 is the double integral evaluated between the limits E:E0-DELTA to E0
! and E':EPFINAL+DELTA to Ep_elastic(E0). By treating the E integration as 
! for the I1 calculation, the double integral is reduced to a single integral.
       EPMAX = E0/(1.D0 + (2.D0*E0/mtarg)*SIN2) - DELTA2 !csk
       W2CUT = MP*MP + (1.D0/XCUT - 1.D0)*
     >     (4.D0*E0*E0*SIN2/(1.D0 + 2.D0*E0*SIN2/(MP*XCUT)))
csk      if (W2CUT .le. 0.0D0) W2CUT = 0.001 ! sek apr12
       IF(.NOT.TAIL_ON) EPMAX = (MP*MP + 2.D0*MP*E0 - 
     >        W2CUT)/(2.D0*MP + 4.D0*E0*SIN2)
       EPMIN = EPFINAL + DELTA
       I2 = 0.D0
! CSK
 	   if(EPMAX .le. EPFINAL) goto 7708 ! sek apr12
 	   BT = BZ*TA
       VARMIN = ((EPMIN - EPFINAL)/EPMIN)**BT
       VARMAX = ((EPMAX - EPFINAL)/EPMAX)**BT
C Note: EXTERNI2NEW yields slightly different results than EXTERNI2. I could
C not determine that one result was "more correct" than the other so I am 
C leaving both in for now. 3/97, LMS.
C       IF(EPMAX.GT.EPMIN) CALL QUADMO3(I2,EXTERNI2NEW,VARMAX,
C     >                     VARMIN,EPSS,NLVL)
       EPMIN = EPFINAL/(1.D0 - VARMIN**(1.D0/BT))
       EPMAX = EPFINAL/(1.D0 - VARMAX**(1.D0/BT))
       IF(EPMAX.GT.EPMIN) CALL QUADMO3(I2,EXTERNI2,EPMAX,EPMIN,
     >                     1.D-3,NLVL, mtarg)
       I2 = -(I2*PROB_E)
 7708	continue !SEK apr12      

! I3 is the double integral evaluated between the limits 
! E:E_elastic to E0-DELTA and E':EPFINAL to EPFINAL+DELTA. By treating 
! the E' integration as for the I1 calculation, the double integral is 
! reduced to a single integral.
       EMAX = E0 - DELTA
       EMIN = EPFINAL/(1.D0 - (2.D0*EPFINAL/mtarg)*SIN2) + DELTA2 !csk

       W2CUT = MP*MP + (1.D0/XCUT - 1.D0)*
     >     (4.D0*EPFINAL*EPFINAL*SIN2)/
     >     (1.D0 - 2.D0*EPFINAL*SIN2/(MP*XCUT))
       if (W2CUT .le. 0.0D0) W2CUT = 0.001 ! sek apr12

       IF(.NOT.TAIL_ON) EMIN = (W2CUT - MP*MP + 
     >       2.D0*MP*EPFINAL)/(2.D0*MP - 4.D0*EPFINAL*SIN2)
       I3 = 0.D0
       IF(EMAX.GT.EMIN.AND.PROB_EP.GT.0.001D0.and.E0.gt.EMIN) THEN
           DE = (EMAX - EMIN)/10.D0
           BT = BZ*TB
           VARMIN = ((E0 - EMAX)/E0)**BT
           VARMAX = ((E0 - EMIN)/E0)**BT
           CALL QUADMO3(I3,EXTERNI3NEW,VARMIN,VARMAX,EPSS,NLVL, mtarg)
c           CALL ODEINT_EX1(I3,1,VARMIN,VARMAX,EPSS,
c     >              H1,HMIN,NOK,NBAD,EXTERNI3NEW,RKQC_EX1)
       ENDIF
 
       I3 = I3*PROB_EP

! I4 is the double integral evaluated between the limits 
! E:E_elastic to E0-DELTA and E':EPFINAL+DELTA to Ep_elastic(E0).
! This is the cpu eating integral. 
       I4 = 0.D0
       IF(EXTPEAKING) THEN
c          print*,'externew2: PROB_E,PROB_EP',PROB_E,PROB_EP
c          print*,'externew2: PROBTOT_E',PROBTOT_E    !kp: PROBTOT_E culprit for NANs
         I4 = ( (PROBTOT_E-PROB_E)*I2/PROB_E + 
     >          (PROBTOT_EP-PROB_EP)*I3/PROB_EP )/2.D0
       ELSE
         EPMAX = E0/(1.D0 + (2.D0*E0/mtarg)*SIN2) - DELTA2 ! csk
         W2CUT = MP*MP + (1.D0/XCUT - 1.D0)*
     >       (4.D0*E0*E0*SIN2/(1.D0 + 2.D0*E0*SIN2/(MP*XCUT)))
csk    if (W2CUT .le. 0.0D0) W2CUT = 0.001 ! sek apr12
         
         IF(.NOT.TAIL_ON) EPMAX = (MP*MP + 2.D0*MP*E0 - 
     >          W2CUT)/(2.D0*MP + 4.D0*E0*SIN2)
         EPMIN = EPFINAL + DELTA
         
         if (EPMAX .gt. EPFINAL) then
           BT = BZ*TA
           VARMIN = ((EPMIN - EPFINAL)/EPMIN)**BT
           VARMAX = ((EPMAX - EPFINAL)/EPMAX)**BT
           IF(EPMAX.GT.EPMIN) CALL QUADMO3(I4,EXTERNI4NEW,VARMIN,
     >                       VARMAX,EPSS,NLVL, mtarg)
         endif
       ENDIF


       ANSWER = I1 + I2 + I3 + I4
c       print*,'externnew2: I1,I2,I3,I4,SIGNOM,I1+I2+I3+I4:',
c     >      I1,I2,I3,I4,SIGNOM,I1+I2+I3+I4 !kp: I4 culprit for NANs
c      WRITE(6,'(6(1PE13.5))') I1,I2,I3,I4,SIGNOM,I1+I2+I3+I4

       RETURN
       END
*******************************************************************************

       SUBROUTINE EXTERNI2NEW(VAR,DUM,ANSWER, mtarg)
*******************************************************************************
* 4/94, LMS.
*******************************************************************************
       IMPLICIT NONE
       INTEGER NOK, NBAD
       REAL*8 ANSWER, VAR, DUM, E0, EPFINAL, THETA, E,
csk     >        MP/0.939D0/, EP, IB, TB, TA, ELOSS_PROB, Z, !SEK apr12
     >        MP/0.938272D0/, EP, IB, TB, TA, ELOSS_PROB, Z,
     >        SIGMAR, SIGMARNEW, FZ, BT, SIGNOM, XCUT,
     >        DELTA, W2, BZ, BREMS, BT_INV, Y, NORMCOR
     
     	real*8 mtarg ! csk
     	
       COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT
       COMMON /KINLIST2/ DELTA, W2, TB, TA
       COMMON /B_ELOSS/ FZ, BZ
       LOGICAL FL_ION/.FALSE./

       BT = BZ*TA
       BT_INV = 1.D0/BT
       Y = VAR**BT_INV
       EP = EPFINAL/(1.D0 - Y)
       SIGNOM = SIGMAR(E0,EP,THETA)
       IB = BREMS(Y,Z,FZ,BZ)
! NORMCOR is just the analytic integration of the Bremsstrahlung probability
! function BREMSSTRAHLUNG (of which BREMS is a part) over all allowed E
! Needed to normalize the probability function.
       NORMCOR = (1.D0 - BT/(1.D0 + BT) + 
     >            3.D0*BT/(4.D0*(2.D0 + BT)) + 
     >            FZ*(1.D0 - BT/(1.D0 + BT)))

       ANSWER = IB*SIGNOM/NORMCOR
c      write(6,*) ep,signom, ANSWER

       RETURN
       END
*******************************************************************************

       SUBROUTINE EXTERNI3NEW(VAR,DUM,ANSWER, mtarg)
*******************************************************************************
* VAR = Y**BT
* 4/94, LMS.
*******************************************************************************
       IMPLICIT NONE
       INCLUDE 'instruct.inc'
       INTEGER NOK, NBAD
       REAL*8 ANSWER, VAR, DUM, E0, EPFINAL, THETA, E,
csk     >        MP/0.939D0/, EP, IB, TB, TA, ELOSS_PROB, Z,
     >        MP/0.938272D0/, EP, IB, TB, TA, ELOSS_PROB, Z,
     >        SIGMAR, SIGMARNEW, FZ, BT, SIGNOM, D, PB,
     >        DEPOLARIZATION, XCUT, DELTA, W2, BT_INV, Y,
     >        BREMS, B, BZ, NORMCOR
     
     	real*8 mtarg ! csk
     	
       COMMON /KINLIST/E0, EPFINAL, THETA, E, Z,XCUT
       COMMON /KINLIST2/ DELTA, W2, TB, TA
       COMMON /B_ELOSS/ FZ, BZ
       LOGICAL FL_ION/.FALSE./

       BT = BZ*TB
       BT_INV = 1.D0/BT
       Y = VAR**BT_INV
       E = E0*(1.D0 - Y)
       SIGNOM = SIGMAR(E,EPFINAL,THETA)
       IB = BREMS(Y,Z,FZ,B)
! NORMCOR is just the analytic integration of IB function over all allowed E
! Needed to normalize the probability function.
       NORMCOR = (1.D0 - BT/(1.D0 + BT) + 
     >            3.D0*BT/(4.D0*(2.D0 + BT)) + 
     >            FZ*(1.D0 - BT/(1.D0 + BT)))

       D = DEPOLARIZATION(E0,E,Z)
       PB = 1.D0 - D
       ANSWER = IB*SIGNOM/NORMCOR
       IF(FL_POL.AND.DEPOL) ANSWER = ANSWER*PB

       RETURN
       END
*******************************************************************************

       SUBROUTINE EXTERNI4NEW(VAR,DUM,ANSWER, mtarg)
*******************************************************************************
* 12/97, LMS.
*
* 4/2000 frw changed common block
*            COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT, EP, Y
*            to
*            COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT
*            which is the standard definition throughout
*            since EXTERNI5NEW (below) also uses EP, it got a new common block
*            of its own -- Y is NOT used anywhere else and thus deleted
*******************************************************************************

       IMPLICIT NONE
       INCLUDE 'instruct.inc'

       INTEGER NOK, NBAD, NLVL
       REAL*8 ANSWER, VAR, DUM, E0, EPFINAL, THETA, E,
csk     >        MP/0.939D0/, W2, IB, TB, TA, Z, I5, DELTA,
     >        MP/0.938272D0/, W2, IB, TB, TA, Z, I5, DELTA,
     >        EPS/5.D-4/, H1/0.01D0/, HMIN/1.0D-30/, EMIN, EMAX,
     >        PI/3.141592654D0/, SIN2, ENERGY, EPRIME, 
     >        ALITTLE/1.0D-14/, XCUT, BT, BT_INV, Y, EP, 
     >        W2CUT, VARMIN, VARMAX, DELTA2, FZ, BZ, 
     >        NORMCOR
     
     	real*8 mtarg ! csk
     	
       COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT
       COMMON /KINLIST2/ DELTA, W2, TB, TA
       COMMON /KINLIST3/ EP
       COMMON /B_ELOSS/ FZ, BZ
       EXTERNAL RKQC_EX2, EXTERNI5NEW

       DELTA2 = 1.D-5
       BT = BZ*TA
       BT_INV = 1.D0/BT
       Y = VAR**BT_INV
       EP = EPFINAL/(1.D0 - Y)
! NORMCOR is just the analytic integration of IB function over all allowed E
! Needed to normalize the probability function.
       NORMCOR = (1.D0 - BT/(1.D0 + BT) + 
     >            3.D0*BT/(4.D0*(2.D0 + BT)) + 
     >            FZ*(1.D0 - BT/(1.D0 + BT)))


       SIN2 = (SIN(THETA*PI/360.D0))**2
       EMAX = E0 - DELTA
       EMIN = EP/(1.D0 - (2.D0*EP/mtarg)*SIN2) + DELTA2 ! csk
       W2CUT = MP*MP + (1.D0/XCUT - 1.D0)*
     >     (4.D0*EP*EP*SIN2)/
     >     (1.D0 - 2.D0*EP*SIN2/(MP*XCUT))
       IF(.NOT.TAIL_ON) EMIN = (W2CUT - MP*MP + 
     >       2.D0*MP*EP)/(2.D0*MP - 4.D0*EP*SIN2)

       I5 = 0.D0
       IF(EMAX.GT.EMIN) THEN
         BT = BZ*TB
         VARMIN = ((E0 - EMAX)/E0)**BT
         VARMAX = ((E0 - EMIN)/E0)**BT
         CALL QUADMO4(I5,EXTERNI5NEW,VARMIN,VARMAX,EPS,NLVL, mtarg)
       ENDIF


       ANSWER = I5/NORMCOR
       RETURN
       END
*******************************************************************************

       SUBROUTINE EXTERNI5NEW(VAR,DUM,ANSWER, mtarg)
*******************************************************************************
* 4/94, LMS.
*
* 4/2000 frw changed common block
*            COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT, EP, YEP
*            to
*            COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT
*            which is the standard definition throughout
*            since EXTERNI4NEW (above) also uses EP, it got a new common block
*            of its own -- YEP is NOT used anywhere else and thus deleted
*******************************************************************************
       IMPLICIT NONE
       INCLUDE 'instruct.inc'
       INTEGER NOK, NBAD
       REAL*8 ANSWER, VAR, DUM, E0, EPFINAL, THETA, E,
csk     >        MP/0.939D0/, EP, IB_E, IB_EP, TB, TA, ELOSS_PROB,
     >        MP/0.938272D0/, EP, IB_E, IB_EP, TB, TA, ELOSS_PROB,
     >        Z, SIGMAR, SIGMARNEW, DELTA, FZ, BT, SIGNOM, W2, Q2,
     >        SIN2, JACOBIAN, BB, CC, PI/3.141592654D0/, ROOT, 
     >        DEPOLARIZATION, D, PB, XCUT, BZ, BT_INV, Y,
     >        BREMS, B, IB, NORMCOR, YEP
     
     	real*8 mtarg ! csk

       COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT
       COMMON /KINLIST2/ DELTA, W2, TB, TA
       COMMON /KINLIST3/ EP
       COMMON /B_ELOSS/ FZ, BZ
       LOGICAL FL_ION/.FALSE./
     	
       BT = BZ*TB
       BT_INV = 1.D0/BT
       Y = VAR**BT_INV
       E = E0*(1.D0 - Y)
       SIGNOM = SIGMAR(E,EP,THETA)
! NORMCOR is just the analytic integration of IB function over all allowed E
! Needed to normalize the probability function.
       NORMCOR = (1.D0 - BT/(1.D0 + BT) + 
     >            3.D0*BT/(4.D0*(2.D0 + BT)) + 
     >            FZ*(1.D0 - BT/(1.D0 + BT)))


       IB_E = BREMS(Y,Z,FZ,BZ)
       IB_EP= BREMS(YEP,Z,FZ,BZ)
       D = DEPOLARIZATION(E0,E,Z)
       PB = 1.D0 - D
       ANSWER = IB_E*IB_EP*SIGNOM/NORMCOR
       IF(FL_POL.AND.DEPOL) ANSWER = ANSWER*PB

       RETURN
       END
*******************************************************************************

ckp     9/16/12: Following two routines simply copied & pasted from extern.f to avoid using the whole extern.f
****************************************************************************************************
c     SUBROUTINE EXTERNI2(EPVAR,DUM,ANSWER)
      SUBROUTINE EXTERNI2(EPVAR,DUM,ANSWER,mtarg) ! csk
*******************************************************************************
* 4/94, LMS.
*******************************************************************************
       IMPLICIT NONE
       INTEGER NOK, NBAD
       REAL*8 ANSWER, EPVAR, DUM, E0, EPFINAL, THETA, E
       REAL*8 MP/0.938272D0/
       REAL*8 EP, IB, TB, TA, ELOSS_PROB, Z
       REAL*8 SIGMAR, SIGMARNEW, FZ, BT, SIGNOM, XCUT
       REAL*8 DELTA, W2

       real*8 mtarg             ! csk

       COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT
       COMMON /KINLIST2/ DELTA, W2, TB, TA
       COMMON /B_ELOSS/ FZ, BT
       LOGICAL FL_ION/.FALSE./

       EP = EPVAR
       SIGNOM = SIGMAR(E0,EP,THETA)
       IB = ELOSS_PROB(EP,EPFINAL,TA,Z,FL_ION)
       ANSWER = IB*SIGNOM
c       write(6,*) e0,ep,theta,signom,ib

       RETURN
       END
*******************************************************************************

       SUBROUTINE BINTE(ENVAR,DUM,ANSWER)
       REAL*8 ENVAR, DUM, ANSWER, E0, EPFINAL, THETA, E
       REAL*8 TB, TA, Z, IB, XCUT, DELTA, W2, ELOSS_PROB

       logical FL_ION

       COMMON /KINLIST/E0, EPFINAL, THETA, E, Z, XCUT
       COMMON /KINLIST2/ DELTA, W2, TB, TA

       IB = ELOSS_PROB(E0,ENVAR,TB,Z,FL_ION)
       ANSWER = IB
       RETURN
       END

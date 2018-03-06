       REAL*8 FUNCTION ELOSS_PROB(E0,E,T,Z,FL_ION)
******************************************************************************
* This subroutine returns the probability of finding an electron in the
* energy interval between E and E+DE at a depth of t (in radiation lengths)
* where E0 is the energy of the electron at t=0.
*
* For large E0 and/or large (E0-E) the loss due to ionization is completely
* negligible. It is included, here, as an option for future use (beyond
* E143), but is currently not implemented.
* 
* This function corresponds to I(E0,E,T) as discussed by Tsai (SLAC-PUB-848)
* where I(E0,E,T) = IB(E0,E,T) + II(E0,E,T). 
* IB is defined by Eq. B.43 and B.3. 
* There are two options for bremsstralung spectra, phi:
* The first is given by Eq. B.6 (complete screening). When this is used with
* our definition of IB, a small normalization correction (on the order 
* of 1.005) has been added for completeness. 
*
* The second is taken from Tsai, Rev. mod. phys., 46, 815 (1974) 
* (see, Eq. 3.83 and 3.84). This is an improved complete screening 
* formula. When used with our definition of IB, a much larger normalization 
* correction (0.98-0.99) is needed. Also, it is not clear whether this
* function, when used with IB, still satisfies the diffusion eq. given
* by Tsai, B.8. This should be the subjecy of further study. The
* diffusion equation describes the change in probability (with T)
* for a fixed energy bin due to the radiating in/out probabilities for
* that bin.
* 
* II is not taken directly from Tsai. It is defined to be the Landau
* probability distribution (L. Landau, J. Phys. USSR 8, 201 (1944)),
* which is valid for fast electrons (for more info, see A. Katramatou,
* SLAC NPAS-TN-86-8), and uses the approximate universal function,
* phi(lambda) as given by Tabata and Ito, NIM 158, 521 (1979).
*
* FL_ION: Flag which turns on/off (T/F) ionization inclusion. When turned on,
* I = IB + II with E0 replaced by E0-DELTA_0 everywhere (DELTA_0 is the
* most probable energy loss due to ionization). When turned off, I = IB.
*
* 2/94, LMS.
******************************************************************************
       IMPLICIT NONE

       REAL*8 E0, E, T, Z, DELTA_0, IB, II, I
       LOGICAL FL_ION

       IF(FL_ION) THEN
! Not yet implemented
c        CALL IONIZATION(E0,E0,T,Z,II)
       ELSE
         II = 0.D0
       ENDIF

       CALL BREMSSTRAHLUNG(E0,E,T,Z,IB)
        
       I = II + IB
       ELOSS_PROB = I

       RETURN
       END
******************************************************************************


       SUBROUTINE BREMSSTRAHLUNG(E0,E1,T,Z,IB)
**************************************************************************
* Returns bremsstrahlung probability, IB, for radiative corrections.
* Formulas come from Mo and Tsai, Rev. Mod. Phys. 41, 205 (1969) and
* Tsai, SLAC-PUB-848 and Rev. Mod. Phys. 46, 815 (1974).
*
* 2/94, LMS.
**************************************************************************
       IMPLICIT NONE

       REAL*8 E0, E1, T, Z, IB, B , FZ, 
     >        RHO, K, GAMFAC, GAMMA, BREMS, Y, NORMCOR, BT
       COMMON /B_ELOSS/ FZ, B

       K = E0 - E1
       Y = K/E0
       RHO = BREMS(Y,Z,FZ,B)/K    
       BT = B*T

C GAMFAC factor cancels out when NORMCOR is applied so don't calculate.
C      GAMFAC = 1.D0/GAMMA(1.D0 + BT)
       GAMFAC = 1.D0

       IB = (Y**BT)*RHO*T*B*GAMFAC
      
! NORMCOR is just the analytic integration of IB function over all allowed E.
       NORMCOR = GAMFAC*(1.D0 - BT/(1.D0 + BT) + 
     >            3.D0*BT/(4.D0*(2.D0 + BT)) + 
     >            FZ*(1.D0 - BT/(1.D0 + BT)))


       IB = IB/NORMCOR
       RETURN
       END
******************************************************************************


      FUNCTION BREMS (Y,Z,FZ,B)
********************************************************************************
*
*  Returns value of Bremstrung probability multiplied by k (photon energy
*     normalized to 1 at Y=0
*
*  Formula from Tsai; Rev of Mod Physics 46,(1974)815,  equation 3.83
*      Complete screening
*
*  Y = k/E : fraction of beam energy carried by photon
*  Z =     : atomic number of target
*
*      6/4/87 Steve Rock
*
* Note that the straggling parameter B is defined to be the limit as y goes 
* to zero of (4/3)*K*RHO in subroutine BREMSSTRAHLUNG where
* BREMS = K*RHO
*
* 2/92, LMS. Fixed error in Coulomb correction term.
* 4/94, FZ passed as an argument.
* 7/97, LMS. Fixed inconsistency problem between B and BREMS definitions.
********************************************************************************
      IMPLICIT NONE
      REAL*8 Z, Y, LRAD, LRADP, FC, BREMS, FZ, B,
     >       ALPHA/7.2973515D-3/
C     REAL*8 LOG1, LOG2, ETA

! TSAI page 846 Table B.2  and eq. 3.67, 3.68
      IF(Z.EQ.1.0)THEN
       LRAD = 5.31D0
       LRADP= 6.144D0
      ELSEIF(Z.EQ.2.0) THEN
       LRAD = 4.79D0
       LRADP=5.62D0
      ELSEIF(Z.EQ.3.0) THEN
       LRAD = 4.74D0
       LRADP= 5.805D0
      ELSEIF(Z.EQ.4.0) THEN
       LRAD = 4.71D0
       LRADP= 5.924D0
      ELSE
       LRAD = DLOG(184.15D0) - DLOG(Z)/3.D0
       LRADP= DLOG(1194.D0) -2.D0*DLOG(Z)/3.D0
      ENDIF

! Get coulomb correction from Tsai eq. 3.3
      FC = 1.202D0*(Z*ALPHA)**2 - 1.0369D0*(Z*ALPHA)**4 +
     >     1.008D0*(Z*ALPHA)**6/(1.D0+ (Z*ALPHA)**2)

! Tsai eq.(3.83) normalized to be 1 at Y=0
      FZ = (Z+1.D0)/(12.D0*(Z *(LRAD -FC) + LRADP))
      BREMS = (1.D0 -Y +.75D0 *Y**2) + (1.D0-Y)*FZ 
      B = 4.D0/3.D0*(1.D0 + FZ)

! Old value for B:
C      LOG1 = LOG(184.15D0/Z**(1.D0/3.D0))
C      LOG2 = LOG(1194.D0/Z**(2.D0/3.D0))
C      ETA = LOG2/LOG1
C      B = 4.D0/3.D0*(1.D0 + (Z + 1.D0)/(9.D0*(Z + ETA)*LOG1))

C     FZ = 0.D0
C     BREMS = (1.D0 -Y +.75D0 *Y**2)

      RETURN
      END
******************************************************************************

       FUNCTION GAMMA(X)
*****************************************************************************
* This subroutine calculates the GAMMA function of non-integer number, X, 
* using Arfken, Eq. 10.56.
*
* 11/92, LMS.
*****************************************************************************
       IMPLICIT NONE
       INTEGER N
       REAL*8 X, Y, ANSWER, XTRAFACT, GAMMA

       REAL*8 B(8)/-0.577191652D0,  0.988205891D0, -0.897056937D0,
     >              0.918206857D0, -0.756704078D0,  0.482199394D0,
     >             -0.193527818D0,  0.035868343D0/

       Y = X
       XTRAFACT = 1.D0
       DO WHILE(Y.GT.1.D0)
         XTRAFACT = XTRAFACT*Y
         Y = Y - 1.D0
       ENDDO
       DO WHILE(Y.LT.0.D0)
         XTRAFACT = XTRAFACT/(Y + 1.D0)
         Y = Y + 1.D0
       ENDDO

       ANSWER = 1.D0
       DO N = 1,8
         ANSWER = ANSWER + B(N)*Y**N
       ENDDO
         
       ANSWER = ANSWER*XTRAFACT
! ANSWER = X! and GAMMA(X) = (X-1)!
       GAMMA = ANSWER/X

       RETURN
       END
******************************************************************************


       FUNCTION DEPOLARIZATION(E0,E,Z)
******************************************************************************
* Subroutine returns the depolarization of a longitudinally polarized electron
* due to bremsstrahlung emmission for an electron with initial energy E0 and 
* final energy E.
* Formulae are from H. Olsen and L. C. Maximon, Pys. Rev. 114, 887 (1959).
* Calculations assume complete screening.
* This calculation reproduces curve shown in Fig. 5 very well.
*
* 8/96, LMS.
******************************************************************************
       IMPLICIT NONE
       INCLUDE 'radcon.inc'
       INTEGER I
       REAL*8 E0, E, Z, DEPOLARIZATION, LOG1, PSI1, PSI2, 
     >        K, XI1Z, FZ, A, SUM

C From Eq. 6.20
       SUM = 0.D0
       A = Z*ALPHA
       DO I = 1,20
         SUM = SUM + 1.D0/(I*(I*I + A*A))
       ENDDO
       FZ = A*A*SUM
 
C From Eqs. 9.11 and 9.13:
       LOG1 = LOG(183.D0/Z**(1.D0/3.D0))
       PSI1 = 4.D0*(LOG1 - FZ)
       PSI2 = PSI1 - 2.D0/3.D0
       K = E0-E
       XI1Z = 1.D0
       DEPOLARIZATION = K*K*(PSI1 - XI1Z*XI1Z*(PSI1 - 2.D0*PSI2/3.D0))/
     >         ((E0*E0 + E*E)*PSI1 - 2.D0*E0*E*PSI2/3.D0)
 
       RETURN
       END

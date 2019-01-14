*******************************************************************************
*     This file contains the up-to date models (as of 2009)  of
*     polarized and unpolarized structure functions
*     to be used in conjunction with "models.f" 
*     15-June-2009, SEK.
*******************************************************************************


********************************************************************************
      SUBROUTINE A1new(X,Q2,TARG,SIGNP,A1)     
********************************************************************************
*     Subroutine returns A1 based on fits to world data.
*     SIGNP is the ratio of neutron to proton corss sections at the desired
*     kinematics. It is supplied by calling routine.
*     
*     The input parameter IA1 is intended to be used for functional testing only.
*     IA1 = 1 yields nominal fits/models
*     IA1 = anything else is for testing alternate fits/models.
*     
*     7/97, LMS.
*     11/97, updated with FRW's most recent A1P global fit.
*     6/98  added IF structure to select external A1 fit when IA1=0
*     2000-2001, SEK: added to more possible external A1 fit functions: sek_fit
*     and sek_nce ("new century edition"), including EG1a data
********************************************************************************

      IMPLICIT NONE
      INCLUDE 'instruct.inc'

      LOGICAL INITIALIZE, Q2DEP/.TRUE./
      INTEGER J                             
      REAL*8 X, Q2, Wsq, AFIT(3), EXPON, CHT, SIGNP, SIGPN, 
     >     DIL_F2, A1, delA1, A1N, A1P, CHTERR/0.15/
      CHARACTER*1 TARG, TARGLAST

      A1 = 0.0D0
csk
      if (IA1.ge.4 .and. IA1.le.6 .and. X.gt.0.D0) then
         Wsq = 0.939D0*0.939D0 + Q2*(1.D0/X - 1.D0)
         call sek_me(Wsq, Q2, TARG, SIGNP, A1, delA1) ! 2006-7 refit including error
         if (IA1.eq.5) A1 = A1 + delA1
         if (IA1.eq.6) A1 = A1 - delA1
csk    New "Millenium Edition" based on new fit to ELSA, E142-155, EMC/SMC/Compass, HERMES,
csk        Hall A and EG1a/b
csk    and new definition of xi. Options 5 and 6 to evaluate uncertainty due to A1 model
csk    based on MINUIT error matrix. SEK 1-Nov-2006 updated 1-Oct-2007

      else
         write(6,*) ' '
         write(6,*) '***** ERROR:  invalid value for parameter IA1 !'
         write(6,*) ' '

      endif                     ! A1 model selection via IA1

      RETURN
      END
********************************************************************************


********************************************************************************
      SUBROUTINE ELASnew(IFFMOD,IPAULI,Q2,E0,TARG,ITAIL,M,W1,W2,G1,G2)
*******************************************************************************
*     This subroutine returns the structure fnctions W1, W2, G1, and G2,
*     evaluated under elastic conditions, i.e., in terms of nucleon form
*     factors. Formulae are from Eqs. 17 and 18, Kuchto and Shumeiko,
*     NP B219, (1983), 412. NO FERMI SMEARING DONE.
*     
*     2/93, LMS.
*     3/1/94, Updated. LMS.
******************************************************************************
      IMPLICIT NONE
      CHARACTER*3 TARG
      INTEGER IFFMOD
      INTEGER ITAIL, IPAULI
      REAL*8 Q2, E0, W1, W2, G1, G2, GEP, GMP, GEN, GMN, 
     >     TAU, M, GM, GE, Z, N, PS1, PS2

      Z = 0.D0
      N = 0.D0
      IF(TARG.EQ.'NH3') THEN
         Z = 1.D0
      ELSEIF(TARG.EQ.'ND3'.OR.TARG.EQ.'LID') THEN
         Z = 1.D0
         N = 1.D0
      ELSEIF(TARG.EQ.'HE3') THEN
         Z = 2.D0
         N = 1.D0
      ELSEIF(TARG.EQ.'NEU') THEN
         Z = 0.D0
         N = 1.D0
      ENDIF

      TAU = Q2/(4.D0*M*M)
      PS1 = 1.D0
      PS2 = 1.D0
      
      IF(TARG.EQ.'ND3') CALL PAULI_SUPPnew(IPAULI,2,Q2,E0,PS1,PS2)
      IF(TARG.EQ.'LID') CALL PAULI_SUPPnew(IPAULI,2,Q2,E0,PS1,PS2)
      IF(TARG.EQ.'HE3') CALL PAULI_SUPPnew(IPAULI,3,Q2,E0,PS1,PS2)


      IF(ITAIL.EQ.2) THEN       ! elastic proton or quasielastic.
C     CALL FFMODEL(IFFMOD,Q2,GEP,GMP,GEN,GMN)
         IFFMOD = 22
         CALL NewFORM(IFFMOD,Q2,GEP,GEN,GMP,GMN)


         W1 = 2.D0*PS1*M*TAU*(GMP*GMP*Z + GMN*GMN*N)
         W2 = 2.D0*PS2*M*(GEP*GEP*Z + GEN*GEN*N + 
     >        TAU*(GMP*GMP*Z + GMN*GMN*N))/(1.D0 + TAU)

         IFFMOD = 21
         CALL NewFORM(IFFMOD,Q2,GEP,GEN,GMP,GMN)

!     Only count polarized nucleons for G1 and G2.
         IF(TARG.EQ.'HE3') THEN
            Z = -0.028D0
            N = 0.86
ckag       these are the polarization numbers from the nuclear physics of 3He
         ENDIF
ckag         GMP = GMP*DSQRT(Z*PS2)
ckag         GEP = GEP*DSQRT(Z*PS2)
ckag         GMN = GMN*DSQRT(N*PS2)
ckag         GEN = GEN*DSQRT(N*PS2)
         IF(Z.LT.0) THEN
            GMP = GMP*DSQRT((-Z)*PS2)
            GEP = GEP*DSQRT((-Z)*PS2)
         ELSE
            GMP = GMP*DSQRT(Z*PS2)
            GEP = GEP*DSQRT(Z*PS2)
         ENDIF
         IF(N.LT.0) THEN
            GMN = GMN*DSQRT((-N)*PS2)
            GEN = GEN*DSQRT((-N)*PS2)
         ELSE
            GMN = GMN*DSQRT(N*PS2)
            GEN = GEN*DSQRT(N*PS2)
         ENDIF

         G1 = GMP*(GEP + TAU*GMP)/(M*(1.D0 + TAU))
     >        + GMN*(GEN + TAU*GMN)/(M*(1.D0 + TAU))
         G2 = GMP*(GEP - GMP)/(2.D0*M*(1.D0 + TAU))
     >        + GMN*(GEN - GMN)/(2.D0*M*(1.D0 + TAU))
      ENDIF


      IF(ITAIL.EQ.1) THEN       ! elastic nuclear (not proton).
         GE = 0.D0
         GM = 0.D0
         IF(TARG.EQ.'HE3') THEN
            CALL FFHE3new(Q2,GE,GM)
         ELSEIF (TARG.EQ.'ND3'.OR.TARG.EQ.'LID') THEN
            CALL FFDnew(Q2,GE,GM)
         ENDIF
         W1 = 2.D0*M*TAU*GM*GM
         W2 = 2.D0*M*(GE*GE + TAU*GM*GM)/(1.D0 + TAU)
         G1 = GM*(GE + TAU*GM)/(M*(1.D0 + TAU))
         G2 = GM*(GE - GM)/(2.D0*M*(1.D0 + TAU))
      ENDIF
      RETURN
      END
      
********************************************************************************
      SUBROUTINE FFDnew(Q2,GE,GM)
***************************************************************************
*     This subroutine returns elastic charge and magnetic form factors for
*     deuterium. 
*     Errors are included for future model dependent studies.
*     2/94, LMS.
***************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      LOGICAL ERROR/.FALSE./
      REAL*8 Q2, GE, GM, AQ, BQ, Q, S1, S2, S3, S4, S5, 
     >     BQERR, AQERR, TAU
      REAL*8 BQA/0.0046D0/, BQAERR/0.0006D0/, BQB/6.8D0/,
     >     BQBERR/0.24D0/, BQC/9.44D-09/, BQCERR/1.28D-08/,
     >     BQD/5.46D0/, BQE/1.6588D0/
      REAL*8 AQA/24.4819D0/, AQAERR/0.1913D0/, AQB/-75.0541D0/, 
     >     AQBERR/0.0425D0/, AQC/162.5866D0/, AQCERR/0.0437D0/,
     >     AQD/3.1238D0/, AQDERR/0.5446D0/, AQE/1.3093D0/,
     >     AQEERR/0.7254D0/

      AQ = 1.D0/(AQA*Q2**4 + AQB*Q2**3 + AQC*Q2**2 +
     >     AQD*Q2 + AQE)**2
      IF(ERROR) THEN
         S1 =  2.D0/(AQA*Q2**4 + AQB*Q2**3 + AQC*Q2**2 +
     >        AQD*Q2 + AQE)**3
         S2 = (AQAERR*Q2**4)**2 + (AQBERR*Q2**3)**2 +
     >        (AQCERR*Q2**2)**2 + (AQDERR*Q2)**2 +
     >        AQEERR**2
         AQERR = DSQRT(S1*S1*S2)
      ENDIF

      Q = DSQRT(Q2)
      BQ = BQA*EXP((-BQB)*Q2) + BQC*EXP((-BQD)*(Q - BQE)**2)
      IF(ERROR) THEN
         S1 = BQAERR*EXP((-BQB)*Q2)
         S2 = BQA*Q2*BQBERR*EXP((-BQB)*Q2)
         S3 = BQCERR*EXP((-BQD)*(Q - BQE)**2)
         BQERR = DSQRT(S1*S1 + S2*S2 + S3*S3)
      ENDIF
      TAU = Q2/(4.D0*MD*MD)

!     Note that A(Q2) = (GC(Q2))**2 + (8/9)*TAU**2*(GQ(Q2))**2 +
!     (2/3)*TAU*(1+TAU)(GM(Q2))**2 and 
!     B(Q2) = (4/3)*TAU*(1+TAU)**2*(GM(Q2))**2 where
!     GC is the charge form factor, GQ is the quadrupole form factor and
!     GM is the magnetic form factor. Here, it is assumed that GE and GM
!     follow the same functional form as given for elastic nucleon
!     scattering.
      if (BQ .gt. 0.0) then
         GM = DSQRT(BQ/(2.D0*TAU))
      else
         GM = 0.0D0
      endif
      if ((AQ*(1.D0+ TAU)).gt.(TAU*GM*GM)) then
         GE = DSQRT( AQ*(1.D0+ TAU) - TAU*GM*GM)
      else
         GE = 0.0D0
      endif

      RETURN
      END
********************************************************************************
      SUBROUTINE FFHE3new(Q2,GE,GM)
***************************************************************************
*     This subroutine returns elastic charge and magnetic form factors for
*     HE3. Low Q2 parameterization is from McCarthy, et al PRC 15, 1396 (1977).
*     High Q2 parameterization for the charge form factor is from Arnold, 
*     et al, PRL 40, 1429 (1978). The high Q2 parameterization is for a 
*     measured combination of the charge and magnetic form factors. Here, it
*     is assumed that the small magnetic form factor can be obtained from the
*     low Q2 parameterization evaluated at high Q2.
*     Errors are included for future model dependent studies.
*     
*     2/94, LMS.
***************************************************************************
      IMPLICIT NONE
      LOGICAL ERROR/.FALSE./
      REAL*8 Q2, Q2FM, GE, GM, FC, FM, FCERR, FMERR, S1, S2,
     >     S3, S4, S5, S6, Q, TAU, MU/-3.2D0/, M/2.80833D0/, AQ2,
     >     AQ2ERR, FCHIGH, FCHIGHERR, FRAC, Z,
     >     HC2/0.0389379D0/     ! (GeV fm)**2
      REAL*8 AFC/0.675D0/, AFCERR/0.008D0/, BFC/0.366D0/, 
     >     BFCERR/0.025D0/, CFC/0.836D0/, CFCERR/0.032D0/,
     >     DFC/-0.00678D0/, DFCERR/0.00083D0/, PFC/0.9D0/,
     >     PFCERR/0.16D0/, Q0FC/3.98D0/, Q0FCERR/0.09D0/
      REAL*8 AFM/0.654D0/, AFMERR/0.024D0/, BFM/0.456D0/,
     >     BFMERR/0.029D0/, CFM/0.821D0/, CFMERR/0.053D0/
      REAL*8 AA/0.034D0/, AAERR/0.004D0/, BA/2.72D0/, BAERR/0.09D0/


      TAU = Q2/(4.D0*M*M)
      Q2FM = Q2/HC2
      Q = DSQRT(Q2FM)
      IF(Q2.LT.0.8D0) THEN 
         FC = DABS(EXP((-AFC)*AFC*Q2FM) - 
     >        BFC*BFC*Q2FM*EXP((-CFC)*CFC*Q2FM)
     >        + DFC*EXP(-(((Q - Q0FC)/PFC)**2)))
         IF(ERROR) THEN
            S1 = 2.D0*AFC*Q2FM*AFCERR*EXP((-AFC)*AFC*Q2FM)
            S2 = 2.D0*BFC*Q2FM*BFCERR*EXP((-CFC)*CFC*Q2FM)
            S3 = 2.D0*CFC*BFC*BFC*Q2FM*Q2FM*CFCERR*EXP((-CFC)*CFC*Q2FM)
            S4 = DFCERR*EXP(-(((Q - Q0FC)/PFC)**2))
            S5 = 2.D0*DFC*(Q - Q0FC)/(PFC*PFC)*Q0FCERR*
     >           EXP(-(((Q - Q0FC)/PFC)**2))
            S6 = S5*(Q - Q0FC)*PFCERR/(PFC*Q0FCERR)
            FCERR = DSQRT(S1*S1 + S2*S2 + S3*S3 + S4*S4 + S5*S5 + S6*S6)
         ENDIF
      ENDIF

      FM = DABS(EXP((-AFM)*AFM*Q2FM) - BFM*BFM*Q2FM*EXP((-CFM)*CFM*Q2FM))
      IF(ERROR) THEN
         S1 = 2.D0*AFM*Q2FM*AFMERR*EXP((-AFM)*AFM*Q2FM)
         S2 = 2.D0*BFM*Q2FM*BFMERR*EXP((-CFM)*CFM*Q2FM)
         S3 = 2.D0*CFM*BFM*BFM*Q2FM*Q2FM*CFMERR*EXP((-CFM)*CFM*Q2FM)
         FMERR = DSQRT(S1*S1 + S2*S2 + S3*S3)
      ENDIF

      IF(Q2.GT.0.7D0) THEN
         AQ2 = AA*EXP((-BA)*Q2)
         IF(ERROR) THEN
            S1 = AAERR*EXP((-BA)*Q2)
            S2 = AQ2*Q2*BAERR
            AQ2ERR = DSQRT(S1*S1 + S2*S2)
         ENDIF
         FCHIGH = DSQRT(DABS(AQ2*AQ2*(1.D0 + TAU) - FM*FM*MU*MU*TAU))
         IF(ERROR) THEN
            S1 = AQ2*AQ2ERR*(1.D0 + TAU)/FCHIGH
            S2 = FM*FMERR*MU*MU*TAU/FCHIGH
            FCHIGHERR = DSQRT(S1*S1 + S2*S2)
         ENDIF
         IF(Q2.GE.0.8D0) THEN
            FC = FCHIGH
            FCERR = FCHIGHERR
         ELSE                   ! Require continuity over overlap region. 
            FRAC = (Q2 - 0.7D0)/0.1D0
            FC = FRAC*FCHIGH + (1.D0 - FRAC)*FC
            IF(ERROR) THEN
               S1 = FRAC*FCHIGHERR
               S2 = (1.D0 - FRAC)*FCERR
               FCERR = DSQRT(S1*S1 + S2*S2)
            ENDIF
         ENDIF
      ENDIF
      IF(ERROR.AND.Q2.LT.15.0) THEN
         FC = FC+FCERR
         FM = FM+FMERR
      ENDIF

!     Absorb Z from Mott cross section here.
      Z = 2.D0
      GE =  Z*DABS(FC)
      GM =  Z*MU*DABS(FM)
      RETURN
      END
*******************************************************************************

*******************************************************************************
      SUBROUTINE F1F2new(X,Q2,TARG,F1,F2,SIGNP,SUPPRESS)
*******************************************************************************
*     This subroutine returns F1 and F2, the inelastic structure functions for 
*     either proton or deuteron, or 3He targets (TARG = 'P','D', 'N', OR '3')
*     SIGNP is only defined for TARG = 'N'
*     
*     12/94, LMS.
*     
*     10/99, SEK. Implement ALL low-Q2 patches in this routine (except Resonances)
*     
*     12/11, SEK. Update with newest "cludge" for F2d, Rd from Eric Christy. Also,
*     increase number of options beyond 13: 14,15,16 will AVOID the
*     label c11   ad-hoc R numbers contained in "resaR_fit_final.dat" and 17-19 will
*     pick the above-mentioned new version of f1f209 (ending with dr)
*     3/14 SEK: Update with newest SF model by Christy/Kalantarians        
******************************************************************************
csk                Meaning of "SFChoice"
csk    10: 2007 version of F1n/F1p/F1A by Peter Bosted/Eric Christie (c) 2007, HERMES
csk    11: 2009 version of F1n/F1p/F1A by Peter Bosted/Eric Christie (c) 2009, HERMES
csk    12: Same version as 11, but with errors added to F2 (and proportionally F1)
csk    13: Same version as 11, but with errors subtracted from R (F2 unchanged)
csk    14: Same version as 11, but without tabulated R values in RR substituted
csk    15-16: Analog to 12-13, but without tabulated R values in RR substituted
csk    17: Newest kludge for Rd, F2d from Eric Christy December 2011. Should be only for D
csk    18-19: Analog to 15-16 for new kludge
csk    20: New 2014 fit by Christie/Kalantarians. Interpolate at very low Q2.
csk           21-22: Analog errors to 12-13

      IMPLICIT NONE

      INCLUDE 'radcon.inc'
      INCLUDE 'instruct.inc'

      REAL*8 X, Q2, F1, F2, W1R, W2R, W1I, W2I, NU, ST, SY, SLOPE, 
     >     DSLOPE, R, DR, WSQ, FRAC, SIN2, F2D, WD1, WD2, ERR1, 
     >     ERR2, ERR3, W1ERR, W2ERR, ERROR1, ERROR2, SIGNP, W, FMAX
      
      REAL*8 Xeff, Q2eff, pival, W1n, W2n, W1Rnres, W2Rnres, W1Dnres, 
     >     W2Dnres, SQRTWSQ, F1res, F2res, Rtemp, dRtemp
      REAL*8 A,B,X2,X3, Q2old, Q2new, Rold, Rnew, FracOld, FracNew, F1new, F2new
      REAL*8 Z_targ, A_targ
      REAL*8 sigtp,siglp, sigt,sigl, f2pm, f1pm, f2nm, f1nm, rn, xvaln(100)
      
      real x4,qsq4,f2smc,err_lo,err_hi
      real*8 SIGTOT, ESIGTOT,EF2
      integer iT, iT2, iflag    ! For f2allm; =1 if Q2>0
      
      integer I, J, IQ2eg1, IWeg1
      real*8 Rin(2), ResR(40,203,2), FRACWHI, FRACWLO, FraclgQ2hi, FraclgQ2lo
      

      CHARACTER*1 TARG
      LOGICAL GOOD, DORES, DOINEL, SUPPRESS, BlendOldNew, Filled /.FALSE./
      data xvaln / 
     &     0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01, !! AV18 w/ CJ Off-Shell ResNR v0
     &     0.14333E+01,0.13223E+00,0.21779E+00,0.10504E+00,0.19322E+00,
     &     0.16957E+00,0.38000E+00,0.77800E+01,0.46008E+01,0.24837E+01,
     &     0.18847E+01,0.30456E+00,0.52820E+02,0.28644E-08,0.28301E+01,
     &     0.58282E-01,0.44231E+03,0.38004E+00,0.39618E+01,0.47920E-01,
     &     0.29969E+03,0.28080E+00,0.24375E+01,0.00000E+00,0.10000E+01,
     &     0.10000E+01,0.20000E+01,0.51557E-01,0.22273E+04,0.46523E+00,
     &     0.86906E+01,0.20978E+03,0.20023E+00,0.17272E+01,0.14890E+00,
     &     -.14301E+03,0.11214E+01,0.23328E+01,-.10144E+00,-.93649E-02,
     &     0.19456E+00,0.19650E+01,0.38000E+00,0.15634E+01,0.14462E+00,
     &     0.10020E+01,0.99810E+00,0.99771E+00,0.10000E+01,0.10145E+01,
     &     0.10039E+01,0.10017E+01,0.10000E+01,0.10000E+01,0.10000E+01,
     &     0.10000E+01,0.10000E+01,0.70981E+02,0.22717E-01,0.71572E+01,
     &     0.00000E+00,0.00000E+00,0.10000E+01,0.20000E+01,0.00000E+00,
     &     0.17507E+04,0.51302E+00,0.74121E+02,0.00000E+00,0.83186E+02,
     &     0.32211E-02,0.75949E+01,0.00000E+00,0.00000E+00,0.10000E+01,
     &     0.10000E+01,0.00000E+00,0.40157E+01,0.16316E-02,0.84797E+00,
     &     0.00000E+00,0.21412E+03,0.55385E+01,0.18668E+03,0.10019E-03,
     &     0.10034E-01,0.31531E+01,0.00000E+00,0.00000E+00,0.18030E+03,
     &     0.43045E+02,0.39896E+00,0.00000E+00,0.36953E+00,0.00000E+00 /

      pival = 3.1415927
      F1 = 0.D0
      F2 = 0.D0
      SIGNP = 0.D0
      if(X .le. 0.0D0) return
      NU = Q2/(2.D0*MN*X)       ! NUCLEON mass
      WSQ = MN2 + 2.D0*MN*NU - Q2 
      
c11   if(x.gt.1.0D0.or.Q2.le.0.0D0.or.Q2.gt.1.0D4.or.x.le.0.0D0)then
      if(WSQ.le.0.0D0 .or. Q2.le.0.0D0 .or. Q2.gt.2.0D4)then
c     write(6,*) ' Wrong kinematics: ', x,Q2,WSQ,NU
         return
      endif
csk       
csk        Fit by Peter Bosted / Eric Christie 2009

      if (.not. Filled) then
         Filled = .TRUE.
         do IQ2eg1 = 1, 40
            do IWeg1 = 1, 203
               do J = 1, 2
                  ResR(IQ2eg1,IWeg1,J) = 0.0D0
               enddo
            enddo
         enddo
         OPEN (UNIT=7, FILE='resaR_fit_final.dat', STATUS='OLD')
         do I = 1, 3840
            read(7,*) IQ2eg1, IWeg1, (Rin(J), J=1,2)
            if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
            do J = 1,2
               ResR(IQ2eg1,IWeg1,J) = Rin(J)
            enddo
         endif
      enddo
      CLOSE(UNIT=7)
      endif
c11   NOTE: This prefilled array of R values is a kludge stemming from a time where the coded R didn't 
c11   work properly for all kinematics. It's not clear that it is still needed. For compatibility reasons
c11   we keep it for now but replace it for SFChoice = 14,15,16 and 17.


      DOINEL = .FALSE.
      DORES = .FALSE. 
      BlendOldNew = .FALSE.
      iT2 = 0
      iflag = 1

      if(SFChoice .lt. 10 .or. SFChoice .gt. 22) then ! 10 for legacy code F1F207, 20-22=2014
         write(6,*) ' Wrong SFChoice: ',AsymChoice, SFChoice, TARG
         return
      endif
      
      if(TARG .eq. 'P') then
         Z_targ = 1.0D0
         A_targ = 1.0D0
         iT = 1
      elseif(TARG .eq. 'D') then
         Z_targ = 1.0D0
         A_targ = 2.0D0
         iT = 2
      elseif(TARG .eq. 'N') then
         Z_targ = 0.0D0
         A_targ = 1.0D0
         iT = 3
      elseif(TARG .eq. '3') then
         Z_targ = 2.0D0
         A_targ = 3.0D0
         iT = 2
         iT2 = 1
      else
         write(6,*) ' Wrong Target ', TARG
         return
      endif
      if ((dabs(A_targ-1.0D0).gt.1.0D-5) .AND. SFChoice .gt. 19) then
         write(6,*) ' Latest models only for p and n ', TARG
         return
      endif
      

      IF(Q2.LT.12.D0.AND.WSQ.LT.10.0D0.AND..NOT.
     >     SUPPRESS) DORES = .TRUE. ! Note increased limits
      IF(.NOT.DORES.OR.WSQ.GT.6.0D0 .or. Q2.gt.8.0D0) DOINEL = .TRUE.
      
      if(DORES) then
         W = DSQRT(WSQ)
         if(WSQ .le. 1.155 .and. A_targ .lt. 1.5) return
c11   
         if(WSQ .le. 4.0D0) then
            Q2old = 0.08D0*(WSQ - 2.0D0)
            Q2new = 0.08D0*WSQ
         else
            Q2old = 0.16
            Q2new = 0.32
         endif
         
         if((SFChoice .gt. 19).and. (Q2.gt. Q2old)) then
            BlendOldNew = (Q2. le. Q2new)
            if (Z_targ .ge. 0.99d0) then
               call rescsp(WSQ,Q2,sigtp,siglp) 
               f2pm = sigtp+siglp
c     f2pm = scale*f2pm
               f2pm = f2pm*2.d0*X*(WSQ-MN2)/8.d0/pival/pival/alpha/0.3894d3
               f2pm = f2pm/(1.d0+4.d0*MN2*X*X/Q2) 

c     f1pm = scale*f1pm
               f1pm = sigtp*(WSQ-MN2)/8.d0/pival/pival/alpha/0.3894d3
               
               Rnew = 0.0d0
               if(sigtp .gt. 0.0) Rnew = siglp/sigtp
               F1new = f1pm
               F2new = f2pm
            else
               call rescsn(WSQ,Q2,xvaln,sigt,sigl)
               f2nm = sigt+sigl
               rn = sigl/sigt

c     f2nm = scale*f2nm
               f2nm = f2nm*2.d0*X*(WSQ-MN2)/8.d0/pival/pival/alpha/0.3894d3
               f2nm = f2nm/(1.d0+4.d0*MN2*X*X/Q2) 

c     f1nm = scale*f1nm
               f1nm = sigt*(WSQ-MN2)/8.d0/pival/pival/alpha/0.3894d3
               
               Rnew = rn
               F1new = f1nm
               F2new = f2nm
            endif
            if (BlendOldNew) then
               call F1F2IN09(Z_targ, A_targ, Q2, WSQ, F1res, F2res, R)
               FracOld = (Q2new-Q2)/(Q2new-Q2old)
               FracNew = (Q2-Q2old)/(Q2new-Q2old)
               F1 = F1res*FracOld + F1new*FracNew
               F2 = F2res*FracOld + F2new*FracNew
               R = R*FracOld + Rnew*FracNew
            else
               F1 = F1new
               F2 = F2new
               R = Rnew
            endif
         elseif(SFChoice .gt. 16) then
            if (W.gt.2.414159265)then
               call F1F2IN09(Z_targ, A_targ, Q2, WSQ, F1, F2, R)
            else
               call F1F2IN09dr(Z_targ, A_targ, Q2, WSQ, F1, F2, R)
               if(W.gt.2.1)then
                  call F1F2IN09(Z_targ, A_targ, Q2, WSQ, F1res, F2res, Rtemp)
                  F1 = F1res + 0.5*(F1-F1res)*(1.0+cos(10.0*(W-2.1)))
                  F2 = F2res + 0.5*(F2-F2res)*(1.0+cos(10.0*(W-2.1)))
                  R = Rtemp + 0.5*(R-Rtemp)*(1.0+cos(10.0*(W-2.1)))
               endif
            endif
         elseif (SFChoice .eq. 10) then
C     jixie            call F1F2IN07(Z_targ, A_targ, Q2, WSQ, F1, F2, R)
         else
            call F1F2IN09(Z_targ, A_targ, Q2, WSQ, F1, F2, R)
         endif

         Rtemp = R
         dRtemp = 0.1D0*dabs(R) + 0.02D0
         IWeg1 = W*1.0D2-0.5D0
         dR = dRtemp
         if (IWeg1 .lt. 107) goto 444
         if (((SFChoice .lt. 14).or.BlendOldNew).and.(IWeg1. le. 202)) then
         R = 0.0
         dR = 0.0
         FRACWHI = W*1.0D2 - 0.5D0 - IWeg1
         FRACWLO = 1.D0 - FRACWHI
         FraclgQ2hi =(DLOG(Q2/0.0084427D0)/DLOG(10.0D0))*13.0D0
         if(FraclgQ2hi.ge.1.0D0)then
            IQ2eg1 = FraclgQ2hi
            FraclgQ2hi = FraclgQ2hi - IQ2eg1
            FraclgQ2lo = 1.0D0 - FraclgQ2hi
            if (IQ2eg1 .gt. 39) then
               CALL R1998(X,Q2,R,DR,GOOD)
               goto 444
            endif
         else
            R = (FRACWLO*ResR(1,IWeg1,2) + FRACWHI*ResR(1,IWeg1+1,2))*Q2/0.01008D0
            dR = (FRACWLO*ResR(1,IWeg1,1) + FRACWHI*ResR(1,IWeg1+1,1))*Q2/0.01008D0
            dR = dabs(R-dR)+0.1D0*dabs(R) + 0.02D0
            goto 444
         endif
         R = FRACWLO*FraclgQ2lo*ResR(IQ2eg1,IWeg1,2) +
     >        FRACWHI*FraclgQ2lo*ResR(IQ2eg1,IWeg1+1,2) +
     >        FRACWLO*FraclgQ2hi*ResR(IQ2eg1+1,IWeg1,2) +
     >        FRACWHI*FraclgQ2hi*ResR(IQ2eg1+1,IWeg1+1,2)
         if(BlendOldNew) then
            R = FracOld*R + FracNew*Rnew
         endif
         dR = FRACWLO*FraclgQ2lo*ResR(IQ2eg1,IWeg1,1) +
     >        FRACWHI*FraclgQ2lo*ResR(IQ2eg1,IWeg1+1,1) +
     >        FRACWLO*FraclgQ2hi*ResR(IQ2eg1+1,IWeg1,1) +
     >        FRACWHI*FraclgQ2hi*ResR(IQ2eg1+1,IWeg1+1,1)
         dR = 0.5*dabs(R-dR)+0.1D0*dabs(R) + 0.02D0
      endif                     !(of existing kludge for SFChoice = 10-13)
 444  continue
      if (R .lt. 0.0D0) R = 0.0D0
      if ((R .lt. 0.02).and.(R .lt. Q2)) then ! Give R a minimum value to avoid A2Soffer = 0
         if (Q2 .lt. 0.02D0) then
            R = Q2
         else
            R = 0.02D0
         endif
      endif
      if ((SFChoice .eq. 13).or.(SFChoice .eq. 16).or.(SFChoice.eq.19).or.
     >     (SFChoice.eq.22)) then
         if (R .lt. dR) then
            R = R + dR          !
         else
            R = R - dR
         endif
      endif
      F1 = (4.D0*Mn*Mn*X*X/Q2+1.D0)*F2/2.0D0/X/(1+R)

      if(A_targ .gt. 1.5) then
         if((SFChoice .gt. 10).and.(SFChoice .lt. 17)) then
            call F1F2QE09(Z_targ, A_targ, Q2, WSQ, F1res, F2res)
         elseif (SFChoice .eq. 10) then
C     jixie                         call F1F2QE07(Z_targ, A_targ, Q2, WSQ, F1res, F2res)
         else
            call F1F2QE09dr(Z_targ, A_targ, Q2, WSQ, F1res, F2res)
         endif
c     if(F1res .lt. 0) write(6,*) Z_targ,A_targ,Q2,WSQ,F1res,F2res
         if (F1res .gt. 0.0) F1 = F1+F1res
         if (F2res .gt. 0.0) F2 = F2+F2res
      endif
      if ((SFChoice.eq.12).or.(SFChoice.eq.15).or.(SFChoice.eq. 18).or.(SFChoice.eq. 21))  then
      F1 = F1*1.03
      F2 = F2*1.03
      endif
      F1res = F1
      F2res = F2
      endif                     ! DORES
c     write(6,*) ' End of DORES: ', F1, F2
      
      if (DOINEL) then
         call SIGMATOT_PARAM(X,Q2,WSQ,iflag,SIGTOT,ESIGTOT,F2,EF2)
         if (WSQ .le. 1.4D0) then
            F2 = F2 * (WSQ**0.25 - 1.155**0.25)/(1.40**0.25 - 1.155**0.25)
            if (F2 .lt. 0.0D0) F2 = 0.0D0
         endif
         if ((SFChoice.eq.12).or.(SFChoice.eq.15).or.(SFChoice.eq. 18).or.(SFChoice.eq. 21))  then
         F2 = F2+EF2
      endif
      if(TARG .ne. 'P') then
         Xeff = X
         Q2eff = Q2
         if (Q2eff .lt. 0.4D0) then
            if (Q2eff .lt. 0.0001) RETURN
            Q2eff = 0.4D0
            Xeff = X*(1.D0 + 0.2715D0/Q2)/1.67875D0
            if (Xeff .ge. 1.0D0) RETURN
         endif
         X2 = Xeff*Xeff
         X3 = X2*Xeff
         A = 0.979 -1.692*X +2.797*X2 -4.313*X3 +3.075*X3*X
         B = (-.171)*X2 + .244*X3 ! replaced 10/22/97 by correct value on x3
         SIGNP = A *(1+X2/Q2eff)*(Q2eff/20.)**B
         if(SIGNP .le. 0.0) SIGNP = 0.0
         if(SIGNP .ge. 1.0) SIGNP = 1.0
         if(TARG .eq. 'N') F2 = F2*SIGNP
         if(TARG .eq. 'D') F2 = F2*(1.0+SIGNP)
         if(TARG .eq. '3') F2 = F2*(2.0+SIGNP)
      endif                     ! TARG .ne. 'P'
      CALL R1998(X,Q2,R,DR,GOOD)
      if ((SFChoice .eq. 13).or.(SFChoice .eq. 16).or.(SFChoice.eq.19).or.(SFChoice.eq.22)) then
      if (R .lt. dR) then
         R = R + dR             !
      else
         R = R - dR
      endif
      endif
      F1 = (4.D0*Mn*Mn*X*X/Q2+1.D0)*F2/2.0D0/X/(1+R)
      if (DORES) then
         FRAC = 0.0
         if (Q2 .gt. 8.0D0) FRAC = (Q2 - 8.0D0)/4.0D0
         if (WSQ .gt. 6.0D0) FRAC = FRAC + (1.0D0 - FRAC)*(WSQ - 6.0D0)/(10.0D0 - 6.0D0)
         FRAC = FRAC*1.570796327 ! Pi/2
         F1 = F1*(dsin(FRAC))**2 + F1res*(dcos(FRAC))**2
         F2 = F2*(dsin(FRAC))**2 + F2res*(dcos(FRAC))**2
      endif
      endif                     ! DOINEL
      return
      END

***************************************************************************
      
      SUBROUTINE G1G2DIS(X,Q2,TH,TARG,G1,G2,A1,A2,G1F1,IX)
***************************************************************************
*     This subroutine evaluates g1(x,Q2), g2(x,Q2), A1(x,Q2), A2(x,Q2) for
*     various model assumptions. The g2ww formula is given by
*     g2ww = -g1(x,q2) + integral over g1(x',q2)/x' from x to 1 limits of 
*     integration.
*     
*     Fits to A1 are used rather than fits to g1/F1 to make sure that the 
*     positivity constraint on A1 always holds.
*     
*     IPOL = 1: Use A1 fit and g2=g2ww
*     2. Use A1 fit and g2=g2ww + g2tw3 fit
*     3. Use A1 fit and g2=g2ww + g2tw3 fit/model
*     4: Use A1 fit and g2=0
*     5: Use A1 fit and A2 = 0
*     6: Use A1 fit and A_transverse = 0
*     7: Use Soffer Limit for A2 (SEK)
*     
*     11/94, LMS.
*     7/97, LMS. Updated to call A1MOD to evaluate fits to A1.
*     9/00, SEK. Updated to use A2 as input (starting point) for g2WW
*     3/04, SEK. Included new "AsymChoice" to choose different models
***************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'

      INTEGER J, NSTEPS, JLAST, IFIRST, NS


      REAL*8 X, Q2, G1F1, G2, SUM, G1, DELX, I(1000),Q2LAST, 
     >     A1, A2, F1, F2, R, DR, IX, XX, XMAX, XLAST,
     >     XPLO, XP, XPHI, GAMMA2, M/0.939D0/, ST, SY, !changed M
     >     SLOPE, DSLOPE, DXP, F2P, F2D, LIMIT, RGAMMA2,
     >     SIGNP, RADCON/1.745329252D-2/,
     >     INVF1, DELLOGX, LOGX, LOBIN, HIBIN, EXPW,
     >     GAM2, ERR, Q2eff, Xeff, W

      REAL*8 TH, SIN2, TAN2, ROOTQ2, EPSILON, ETA, ZETA, 
     >     E, EP, NU, G2BAR, tau, delov1ptau
      CHARACTER*1 TARG
      LOGICAL GOOD, SUPPRESS
      
      SUPPRESS = .FALSE.        ! Suppresses resonances in favor of scaling
!     IPOL = 1

      W = dsqrt(M**2 + Q2*(1.D0/x - 1.D0)) ! needed for IPOL=7

      NS = 20
      IF(.NOT.INTPEAKING.OR..NOT.EXTPEAKING) NS = 10
      NSTEPS = MAX(2,INT((1.D0 - X)*NS))
      if ((Q2 .lt. 10.0) .and. (A2 .ne. 0.0)) then
         XMAX = Q2/(4.0D0 - M*M + Q2)
      else
         XMAX = 1.0D0
      endif
      DELLOGX = (LOG(XMAX) - LOG(X))/DFLOAT(NSTEPS)
csk       EXPW = EXP(DELLOGX/2.D0)
csk       EXPW = EXPW - 1.D0/EXPW

      JLAST = NSTEPS+1
      I(JLAST) = 0.D0           !CSK Linda's first jump was too big
      if ((A2 .ne. 0.0) .and. (XMAX .lt. 1.0D0)) then
         CALL F1F2new(XMAX,Q2,TARG,F1,F2,SIGNP,SUPPRESS)
         tau = Q2/4.0D0/M**2/XMAX**2
         I(JLAST) = sqrt(tau)*F1*A2
      endif
      G2BAR = 0.0
      GAM2 = 4.D0*M*M/Q2
      
!     Following loop solves for g1 and g2ww from A1 (differential equation).
!     See E143 tech. note 82 by L. Stuart for specifics.
!     Integration is done on a logarithmic scale (dx varies logarithmically)
!     because it is faster and accuracy is not compromised.
      IF(IPOL.LE.3) THEN
         IFIRST = NSTEPS
         DO J=IFIRST,1,-1
            XX = X*EXP( (DFLOAT(J) - 0.5D0)*DELLOGX )
csk           DELX = XX*EXPW

            GAMMA2 = GAM2*XX*XX
csk           RGAMMA2 = DSQRT(GAMMA2)
            tau = 1/GAMMA2
            CALL F1F2new(XX,Q2,TARG,F1,F2,SIGNP,SUPPRESS)
            CALL A1new(XX,Q2,TARG,SIGNP,A1)

            IF(IPOL.EQ.2.OR.IPOL.EQ.3) THEN
               CALL G2TW3new(XX,Q2,TARG,IPOL-1,G2BAR)
            ENDIF
            delov1ptau = DELLOGX/(1.0 + tau)
            I(J) = (I(JLAST)*(1.0 + 0.5*delov1ptau)
     >           + (G2BAR + tau*A1*F1)*delov1ptau)/
     >           (1.0 - 0.5*delov1ptau)
            JLAST = J
         ENDDO
      ENDIF

      GAMMA2 = GAM2*X*X
      RGAMMA2 = DSQRT(GAMMA2)
      CALL F1F2new(X,Q2,TARG,F1,F2,SIGNP,SUPPRESS)
      CALL A1new(X,Q2,TARG,SIGNP,A1)
      Xeff = X
      Q2eff = Q2
      CALL R1998(Xeff,Q2eff,R,DR,GOOD) ! (used for DR only)

      if (x*F1 .gt. 0.0D0) then
         R = F2/F1/2.0/x*(1+GAMMA2)-1.0
      else
         R = 0.0D0
      endif
      

      INVF1 = 0.D0
      IF(F1.GT.0.D0) INVF1 = 1.D0/F1

      IF(IPOL.LE.3) THEN
         IX = I(1)
         IF(IPOL.EQ.2.OR.IPOL.EQ.3) THEN
            CALL G2TW3new(X,Q2,TARG,IPOL-1,G2BAR)
            I(1) = I(1) + G2BAR
         ENDIF
         G1 = (A1*F1 + GAMMA2*I(1))/(1.D0 + GAMMA2)
         G2 = -G1 + I(1)
         G1F1 = G1*INVF1
         A2 = RGAMMA2*INVF1*(G1 + G2)
      ELSEIF(IPOL.EQ.4) THEN
         G2 = 0.D0
         G1 = F1*A1
         A2 = RGAMMA2*INVF1*(G1 + G2)
         G1F1 = A1
         IX = 0.D0
      ELSEIF(IPOL.EQ.5) THEN
         A2 = 0.D0
         G1 = F1*A1/(1.D0 + GAMMA2)
         G2 = -G1
         G1F1 = G1*INVF1
         IX = 0.D0
      ELSEIF(IPOL.EQ.6) THEN
         SIN2 = (DSIN(TH*25.0/2.D0))**2
         TAN2 = SIN2/(1.D0 - SIN2)
         ROOTQ2 = DSQRT(Q2)
         NU = Q2/(2.D0*M*X)
         E = 0.5D0*(NU + DSQRT(NU*NU + Q2/SIN2))
         EP = E - NU
         EPSILON = 1.D0/(1.D0 + 2.D0*TAN2*(1.D0 + NU*NU/Q2))
         ETA = EPSILON*ROOTQ2/(E - EP*EPSILON)
         ZETA = ETA*(1.D0 + EPSILON)/(2.D0*EPSILON)

         G1 = F1*A1*(1.D0 + RGAMMA2*ZETA)/(1.D0 + GAMMA2)
         G2 = F1*A1*(ZETA/RGAMMA2 - 1.D0)/(1.D0 + GAMMA2)
         G1F1 = G1*INVF1
         A2 = RGAMMA2*INVF1*(G1 + G2)
         IX = 0.D0
      ELSEIF (IPOL .eq. 7) THEN
         A2 = DSQRT(0.5*(A1 + 1.0D0)*R) ! Soffer limit
         G1F1 = 1.0D0/(1.0D0 + GAMMA2)*(A1 + RGAMMA2*A2)
         G1 = G1F1*F1
         G2 = F1/(1.0D0 + GAMMA2)*(A2/RGAMMA2 - A1)
         IX = G1 + G2
      ELSE
         A1 = 0.D0
         A2 = 0.D0
         G1 = 0.D0
         G2 = 0.D0
         G1F1 = 0.D0
         IX = 0.D0
      ENDIF
!     At low Q2 <0.5 A2 can get too big. Constrain.
      LIMIT = DSQRT(R+DR)       
c     LIMIT = DSQRT(R)
      IF(DABS(A2).GT.LIMIT) THEN
         IF(A2.GT.LIMIT)  A2 = LIMIT
         IF(A2.LT.-LIMIT) A2 = -LIMIT
         G1 = F1*(A1 + RGAMMA2*A2)/(1.D0 + GAMMA2)
         G2 = F1*(A2/RGAMMA2 - A1)/(1.D0 + GAMMA2)
         G1F1 = G1*INVF1
         IX = G2 + G1
      ENDIF
      RETURN
      END      

********************************************************************************

      SUBROUTINE G1G2new(X,Q2,TH,TARG,G1,G2,A1,A2)
*******************************************************************************
*     This subroutine returns model nucleon spin structure functions and
*     asymmetries in the resonance region. It combines both the resonant
*     and nonresonant background components. If it is called with kinematics
*     outside the resonance region, it will return the deep inelastic results
*     calculated in the subroutine g1g2.for.
*     
*     11/94, LMS. 
*     Changed 10/29/95 SEK to allow 2nd resonance region and delta to have
*     arbitrary A1
*     6/18/97 SEK, Changed to allow AO parametrization to be used.
*     10/97, LMS, Implemented IPOLRES for old and new resonance model:
*     IPOLRES = 1: Improved model using A0 parameterization for resonance 
*     asymmetries (W2<2.0 GeV**2)
*     2: Old model using fits to E143 data to determine resonance 
*     asymmetries.
*     3: Extrapolate polarized DIS into resonance region.
*     5/98 FRW changed transition from RES to DIS to work better at large Q2
*     4/99 FRW changed transition from RES to DIS for 4.0 < W2 < 4.3
*     due to limit of validity range of SEK model (IPOLRES=1), W2 < 4.0
*     
*     ... Many changes by SEK, most recent 4-Jun-01
*     ... Completely new segment to implement updated A1p models SEK 24-Nov-2008
******************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'kpaVarChanges.inc' !Added by Jixie, to access changeA1

      REAL*8 X, Q2, G1, G2, A1, A2, W2, NU, TH, XMAX, A1n, A2n, A1p, A2p

      REAL*8 MP/0.939D0/, F1, F2, F1NR, F1R, F2R, F3R, SIGTOTERR, !changed MP
     >     SIG_NRES(2), SIG_RES1(2), SIG_RES2(2), SIG_RES3(2),
     *     SIGd_NRES(2), SIGd_RES1(2), SIGd_RES2(2), SIGd_RES3(2),
     *     SIGTOTd, SIGTOTn, SIGRESd, SIGRESn, SIGRESp,
     >     SIGTOT, R, R1, R2, R3, SNPB/0.40D0/, SNPR/1.D0/, 
     >     DUM1, GAMMA2, G1F1, RGAMMA2, TERM, DIL_F2, DUMMY, F1p, F1n

      INTEGER IWdmt, IQ2dmt, IQ2eg1, IWeg1
      REAL*8 W, DELW, FRACWLO, FRACWHI, FRACQ2LO, FRACQ2HI, FWLO, FWHI
      
      REAL*8 ResA1p(40,203,3), asymin(3), A1pRESnew, FraclgQ2hi, FraclgQ2lo,
     *     ResA1n(40,203,3), ResA2p(40,203,3),ResA2n(40,203,3),
     *     A2pRESnew, A1nRESnew, A2nRESnew, A2RESP, A2RESN

      INTEGER IW, IQ2, indexA1, indexA2
      INTEGER I, J

      CHARACTER*1 TARG, TARGP /'P'/, TARGN /'N'/
      LOGICAL FILLED /.FALSE./,goroper /.TRUE./,SUPPRESS /.FALSE./, LoQ2 /.FALSE./

      IPOLRES = 1               ! I dont know why I have to say this here explicitely

      A1 = 0.0
      A2 = 0.0
      G1 = 0.0
      G2 = 0.0
      if(AsymChoice .lt. 11) return
      if(AsymChoice .lt. 14) then
         indexA1 = 14 - AsymChoice
         indexA2 = 3
      else if(AsymChoice .le. 15) then
         indexA1 = 3
         indexA2 = 16 - AsymChoice
      else
         return
      endif
      
      IF(.NOT.FILLED.AND.IPOLRES.EQ.1) THEN
         FILLED = .TRUE.
         OPEN (UNIT=7, FILE='resa1p_fit_final.dat', STATUS='OLD')
         OPEN (UNIT=8, FILE='resa1n_fit_final.dat', STATUS='OLD')
         OPEN (UNIT=9, FILE='resa2p_fit_final.dat', STATUS='OLD')
         OPEN (UNIT=10, FILE='resa2n_fit_final.dat', STATUS='OLD')
         do I = 1, 3840
            read(7,*) IQ2eg1, IWeg1, (asymin(J), J=1,3)
            if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
            do J = 1,3
               ResA1p(IQ2eg1,IWeg1,J) = asymin(J)
            enddo
         endif
         read(8,*) IQ2eg1, IWeg1, (asymin(J), J=1,3)
         if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
         do J = 1,3
            ResA1n(IQ2eg1,IWeg1,J) = asymin(J)
         enddo
      endif
      read(9,*) IQ2eg1, IWeg1, (asymin(J), J=1,3)
      if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
      do J = 1,3
         ResA2p(IQ2eg1,IWeg1,J) = asymin(J)
      enddo
      endif
      read(10,*) IQ2eg1, IWeg1, (asymin(J), J=1,3)
      if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
      do J = 1,3
         ResA2n(IQ2eg1,IWeg1,J) = asymin(J)
      enddo
      endif
      enddo
      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      CLOSE(UNIT=10)
      do I = 1,40
         do J = 1,3
            ResA1p(I,107,J) = 1.0D0
            ResA1n(I,107,J) = 1.0D0
            ResA2p(I,107,J) = 0.0D0
            ResA2n(I,107,J) = 0.0D0
         enddo
      enddo
      ENDIF                     ! Filled
      
c     IF(IPOLRES.EQ.3) then
c     CALL G1G2DIS(X,Q2,TH,TARG,G1,G2,A1,A2,G1F1,DUM1)       
c     goto 555
c     ENDIF
      
      W2 = MP*MP + Q2*(1.D0/X - 1.D0)

      W = DSQRT(W2)
      IWeg1 = W*1.0D2-0.5D0
      if (IWeg1 .lt. 107) then
         A1 = 1.0
         return
      endif
      FRACWHI = W*1.0D2 - 0.5D0 - IWeg1
      FRACWLO = 1.D0 - FRACWHI
      FraclgQ2hi =(DLOG(Q2/0.0100786D0)/DLOG(10.0D0))*13.0D0
      LoQ2 = (FraclgQ2hi.lt.0.0D0)
      if(loQ2)then
         IQ2eg1 = 1
         FraclgQ2hi = 0.0D0
         FraclgQ2lo = 1.0D0
      else
         IQ2eg1 = FraclgQ2hi+1
         FraclgQ2hi = FraclgQ2hi - (IQ2eg1-1.0D0)
         FraclgQ2lo = 1.0D0 - FraclgQ2hi
      endif
      if(IQ2eg1.gt.40) then
         CALL G1G2DIS(X,Q2,TH,TARG,G1,G2,A1,A2,G1F1,DUM1)
         return
      endif
      
!     Check if we are in the resonance region or not. W2 cut inspired by data and AO model.
      if(IWeg1.gt.199)then
         IW = 199
         FWLO = 0.5
         FWHI = 0.5
         if(IQ2eg1.gt. 39)then
c     CALL G1G2DIS(X,Q2,TH,TARGP,G1,G2,A1,A2,G1F1,DUM1)
            A2 = 0.0
            A2RESP = FWLO*FraclgQ2lo*ResA2p(IQ2eg1,IW,indexA2) +
     >           FWHI*FraclgQ2lo*ResA2p(IQ2eg1,IW+1,indexA2) +
     >           FraclgQ2hi*A2
c     CALL G1G2DIS(X,Q2,TH,TARGN,G1,G2,A1,A2,G1F1,DUM1)
            A2RESN = FWLO*FraclgQ2lo*ResA2n(IQ2eg1,IW,indexA2) +
     >           FWHI*FraclgQ2lo*ResA2n(IQ2eg1,IW+1,indexA2) +
     >           FraclgQ2hi*A2
         else
            A2RESP = FWLO*FraclgQ2lo*ResA2p(IQ2eg1,IW,indexA2) +
     >           FWHI*FraclgQ2lo*ResA2p(IQ2eg1,IW+1,indexA2) +
     >           FWLO*FraclgQ2hi*ResA2p(IQ2eg1+1,IW,indexA2) +
     >           FWHI*FraclgQ2hi*ResA2p(IQ2eg1+1,IW+1,indexA2)
            A2RESN = FWLO*FraclgQ2lo*ResA2n(IQ2eg1,IW,indexA2) +
     >           FWHI*FraclgQ2lo*ResA2n(IQ2eg1,IW+1,indexA2) +
     >           FWLO*FraclgQ2hi*ResA2n(IQ2eg1+1,IW,indexA2) +
     >           FWHI*FraclgQ2hi*ResA2n(IQ2eg1+1,IW+1,indexA2)
            if(LoQ2)then
               A2RESP = A2RESP*DSQRT(Q2/0.0100786D0)
               A2RESN = A2RESN*DSQRT(Q2/0.0100786D0)
            endif
         endif
         
         XMAX = Q2/(4.0D0 - MP*MP + Q2) ! Limit of W=2
         CALL F1F2new(XMAX,Q2,TARGP,F1p,F2,DUM1,SUPPRESS)
         CALL F1F2new(XMAX,Q2,TARGN,F1n,F2,DUM1,SUPPRESS)
         if (TARG .eq. 'P') then
            A2 = A2RESP
         else if (TARG .eq. 'N') then
            A2 = A2RESN
         else if (TARG .eq. '3') then
            A2 = 
     *           (0.87D0*A2RESn*F1n - 2.D0*0.027D0*A2RESp*F1p)/(2.0D0*F1p + F1n)
         else
            A2 = 0.925*(F1p*A2RESp + F1n*A2RESn)/(F1p+F1n)
         endif
         A2p = A2RESP
         A2n = A2RESN
csk        This will be passed as A2_WW(xRR) to G1G2
         CALL G1G2DIS(X,Q2,TH,TARGP,G1,G2,A1p,A2p,G1F1,DUM1)
         CALL G1G2DIS(X,Q2,TH,TARGN,G1,G2,A1n,A2n,G1F1,DUM1)
      endif                     ! IWeg1 gt 199
      
      CALL F1F2new(X,Q2,TARGP,F1p,F2,DUM1,SUPPRESS)
      CALL F1F2new(X,Q2,TARGN,F1n,F2,DUM1,SUPPRESS)
      
      IF(IWeg1.gt.202) then        
         if (TARG .eq. 'P') then
            A2 = A2p
            A1 = A1p
         else if (TARG .eq. 'N') then
            A2 = A2n
            A1 = A1n
         else if (TARG .eq. '3') then
            A2 = (0.87D0*A2n*F1n - 2.D0*0.027D0*A2p*F1p)/(2.0D0*F1p + F1n)
            A1 = (0.87D0*A1n*F1n - 2.D0*0.027D0*A1p*F1p)/(2.0D0*F1p + F1n)
         else
            A2 = 0.925*(F1p*A2p + F1n*A2n)/(F1p+F1n)
            A1 = 0.925*(F1p*A1p + F1n*A1n)/(F1p+F1n)
         endif
         NU = Q2/(2.D0*MP*X)
         GAMMA2 = Q2/(NU*NU)
         RGAMMA2 = DSQRT(GAMMA2)
         TERM = 1.D0 + GAMMA2
         CALL F1F2new(X,Q2,TARG,F1,F2,DUM1,SUPPRESS)
         G1 = F1*(A1 + RGAMMA2*A2)/TERM
         G2 = F1*(A2/RGAMMA2 - A1)/TERM
         return
      endif
      
      A2 = 0.0
      if(IQ2eg1.gt. 39)then
         CALL G1G2DIS(X,Q2,TH,TARGP,G1,G2,A1,A2,G1F1,DUM1)
         A1pRESnew = FRACWLO*FraclgQ2lo*ResA1p(IQ2eg1,IWeg1,indexA1) +
     >        FRACWHI*FraclgQ2lo*ResA1p(IQ2eg1,IWeg1+1,indexA1) +
     >        FraclgQ2hi*A1
         A2pRESnew = FRACWLO*FraclgQ2lo*ResA2p(IQ2eg1,IWeg1,indexA2) +
     >        FRACWHI*FraclgQ2lo*ResA2p(IQ2eg1,IWeg1+1,indexA2) +
     >        FraclgQ2hi*A2
         A2 = 0.0
         CALL G1G2DIS(X,Q2,TH,TARGN,G1,G2,A1,A2,G1F1,DUM1)
         A1nRESnew = FRACWLO*FraclgQ2lo*ResA1n(IQ2eg1,IWeg1,indexA1) +
     >        FRACWHI*FraclgQ2lo*ResA1n(IQ2eg1,IWeg1+1,indexA1) +
     >        FraclgQ2hi*A1
         A2nRESnew = FRACWLO*FraclgQ2lo*ResA2n(IQ2eg1,IWeg1,indexA2) +
     >        FRACWHI*FraclgQ2lo*ResA2n(IQ2eg1,IWeg1+1,indexA2) +
     >        FraclgQ2hi*A2
      else
         A1pRESnew = FRACWLO*FraclgQ2lo*ResA1p(IQ2eg1,IWeg1,indexA1) +
     >        FRACWHI*FraclgQ2lo*ResA1p(IQ2eg1,IWeg1+1,indexA1) +
     >        FRACWLO*FraclgQ2hi*ResA1p(IQ2eg1+1,IWeg1,indexA1) +
     >        FRACWHI*FraclgQ2hi*ResA1p(IQ2eg1+1,IWeg1+1,indexA1)
         A2pRESnew = FRACWLO*FraclgQ2lo*ResA2p(IQ2eg1,IWeg1,indexA2) +
     >        FRACWHI*FraclgQ2lo*ResA2p(IQ2eg1,IWeg1+1,indexA2) +
     >        FRACWLO*FraclgQ2hi*ResA2p(IQ2eg1+1,IWeg1,indexA2) +
     >        FRACWHI*FraclgQ2hi*ResA2p(IQ2eg1+1,IWeg1+1,indexA2)
         A1nRESnew = FRACWLO*FraclgQ2lo*ResA1n(IQ2eg1,IWeg1,indexA1) +
     >        FRACWHI*FraclgQ2lo*ResA1n(IQ2eg1,IWeg1+1,indexA1) +
     >        FRACWLO*FraclgQ2hi*ResA1n(IQ2eg1+1,IWeg1,indexA1) +
     >        FRACWHI*FraclgQ2hi*ResA1n(IQ2eg1+1,IWeg1+1,indexA1)
         A2nRESnew = FRACWLO*FraclgQ2lo*ResA2n(IQ2eg1,IWeg1,indexA2) +
     >        FRACWHI*FraclgQ2lo*ResA2n(IQ2eg1,IWeg1+1,indexA2) +
     >        FRACWLO*FraclgQ2hi*ResA2n(IQ2eg1+1,IWeg1,indexA2) +
     >        FRACWHI*FraclgQ2hi*ResA2n(IQ2eg1+1,IWeg1+1,indexA2)
      endif
      if(LoQ2)then
         A2pRESnew = A2pRESnew*DSQRT(Q2/0.0100786D0)
         A2pRESnew = A2pRESnew*DSQRT(Q2/0.0100786D0)
      endif
      
      if (IWeg1.gt.199) then
         FWHI = (W - 2.005)/0.030
         if (FWHI .gt. 1.0) FWHI = 1.0
         if (FWHI .lt. 0.0) FWHI = 0.0
         FWLO = 1.0 - FWHI
         A1p = FWLO*A1pRESnew + FWHI*A1p
         A1n = FWLO*A1nRESnew + FWHI*A1n
         A2p = FWLO*A2pRESnew + FWHI*A2p
         A2n = FWLO*A2nRESnew + FWHI*A2n
      else
         A1p = A1pRESnew
         A1n = A1nRESnew
         A2p = A2pRESnew
         A2n = A2nRESnew
      endif
      
      if (TARG .eq. 'P') then
         A2 = A2p
         A1 = A1p
      else if (TARG .eq. 'N') then
         A2 = A2n
         A1 = A1n
      else if (TARG .eq. '3') then
         A2 = (0.87D0*A2n*dum1 - 2.D0*0.027D0*A2p)/(2.0D0 + dum1)
         A1 = (0.87D0*A1n*dum1 - 2.D0*0.027D0*A1p)/(2.0D0 + dum1)
      else
         A2 = 0.925*(F1p*A2p + F1n*A2n)/(F1p+F1n)
         A1 = 0.925*(F1p*A1p + F1n*A1n)/(F1p+F1n)
      endif
      
c     =========================== Jixie: see kpaVarChanges.inc
      A1 = A1 + changeA1
c     =========================== Jixie

      NU = Q2/(2.D0*MP*X)
      GAMMA2 = Q2/(NU*NU)
      RGAMMA2 = DSQRT(GAMMA2)
      TERM = 1.D0 + GAMMA2
      CALL F1F2new(X,Q2,TARG,F1,F2,DUM1,SUPPRESS)
      G1 = F1*(A1 + RGAMMA2*A2)/TERM
      G2 = F1*(A2/RGAMMA2 - A1)/TERM
      
 555  continue

      RETURN
      END      



********************************************************************************
      SUBROUTINE G2TW3new(X,Q2,TARG,IP,MOD)
*******************************************************************************
*     This program returns  G2 TWIST-3  model from a fit to the data.
*     11/94, LMS.
******************************************************************************
      IMPLICIT NONE

      INTEGER IP
      CHARACTER*1 TARG
      REAL*8 X, Q2, MOD

      IF(IP.LE.2) THEN
ckag            fits to e155x data
         IF(TARG.EQ.'P') MOD = 0.137D0*LOG(4.582D0*X)*(1.D0-X)**3
c     IF(TARG.EQ.'P') MOD = 0.14D0*LOG(4.7D0*X)*(1.D0-X)**3
         IF(TARG.EQ.'D') MOD = 0.0943D0*LOG(6.827D0*X)*(1.D0-X)**3
c     IF(TARG.EQ.'D') MOD = 0.06D0*LOG(6.8D0*X)*(1.D0-X)**3
C     IF(TARG.EQ.'P') MOD = 0.36D0*LOG(3.4D0*X)*(1.D0-X)**3
C     IF(TARG.EQ.'D') MOD = 1.26D0*LOG(3.37D0*X)*(1.D0-X)**3
         IF(TARG.EQ.'N'.OR.TARG.EQ.'3')
     >        MOD = 0.0943D0*LOG(6.827D0*X)*(1.D0-X)**3 -
     >        0.137D0*LOG(4.582D0*X)*(1.D0-X)**3
c     >       MOD = -1.27747D0*LOG(3.25584D0*X)*(1.D0-X)**3
c     >       MOD = (-0.43885D0)*(1.D0-X)**3 * (1.D0 - 0.333209D0/X)
         if(ip.eq.2) MOD = MOD*DSQRT(3.D0/Q2)
      ELSEIF(IP.EQ.3) THEN
         IF(TARG.EQ.'P') MOD = (1.D0-X)**3*X**1.037*
     >        (-0.682 + 3.510*X - 4.089*X**2)
         IF(TARG.EQ.'D') MOD = (1.D0-X)**3*X**2.779*
     >        (-21.726 + 88.233*X - 82.085*X**2)
      ENDIF
      RETURN
      END
********************************************************************************

********************************************************************************
      SUBROUTINE INELnew(E,EP,TH,TARG,W1,W2,G1,G2)
*******************************************************************************
*     This subroutine returns W1 and W2, the inelastic structure functions and
*     G1 and G2, the inelastic spin structure functions for either proton,
*     deuteron, neutron, or 3He (TARG = 'P', 'D', 'N', or '3')
*     
*     2/93, LMS.
*     6/93, LMS. Get G1 using E155 formulae
*     6/93, LMS. Added deuterium.
******************************************************************************
      IMPLICIT NONE

      INCLUDE 'radcon.inc'
      INCLUDE 'instruct.inc'

      REAL*8 E,EP,TH, X, Q2, F1, F2, W1, W2, NU, G1, G2, 
     >     G1X, G2X, WSQ, DUM1, DUM2, A1, A2, SIN2, SIGNP


      CHARACTER*1 TARG
      LOGICAL GOOD, SUPPRESS

      W1 = 0.D0
      W2 = 0.D0
      G1 = 0.D0
      G2 = 0.D0

      SIN2 = (DSIN(TH*RADCON/2.D0))**2
      Q2 = 4.D0*E*EP*SIN2
      NU = E - EP
      X = Q2/(2.D0*MN*NU)
      WSQ = MN2 + 2.D0*MN*NU - Q2

      IF(WSQ.LT.1.15D0.AND.TARG.EQ.'P') RETURN       
c     IF(WSQ.LT.1.15D0.AND.TARG.NE.'D') RETURN 
      
c     IF(WSQ.GT.4.0) RETURN  

      IF(FL_POL.AND..NOT.FL_UNPOL) GOTO 35

      SUPPRESS = .FALSE.
      IF(NORES) SUPPRESS = .TRUE.
      CALL F1F2new(X,Q2,TARG,F1,F2,SIGNP,SUPPRESS)
      W1 = F1/MN
      W2 = F2/NU


 35   CONTINUE
      IF(FL_POL) THEN

         CALL G1G2new(X,Q2,TH,TARG,G1X,G2X,A1,A2)

         G1 = G1X/(MN*MN*NU)
         G2 = G2X/(MN*NU*NU)

      ENDIF


      RETURN
      END      
********************************************************************************

********************************************************************************
      SUBROUTINE NewFORM(IG,QQG,GEP,GEN,GMP,GMN)
********************************************************************************
*     CALCULATE NUCLEON FORM FACTORS
*     
csk        CHANGED 3-Oct-2007 following Arrington's paper
*     IG = 21 - Up-to-date fits for TRUE form factors (-> asymmetries)
*     22 - Fudged fit for "nominal" form factors (-> cross sections)
*     QQG = INPUT Q SQUARED (GEV**2)
********************************************************************************
      IMPLICIT NONE
      
      INTEGER IG, INPUT, IOUT
      REAL*8 QQ, QQG, QG, TAU, GEP, GEN, GMP, GMN, GT, T1, T2,
     >     TOP, BOT, F1S, F1V, F2S, F2V, GD, RS, RV, F1E,
     >     F2E, F1M, F2M, F1, F2, F3, GES, GMS, GEV, GMV,
     >     F1RHO, F2RHO, F1P, F2P, QQP, C1, C2, C3, C4, 
     >     F2VK, F2SK, F1N, F2N, QQG1, QQG2, ALPH, RHO, STUFF
      REAL*8 GAM/0.25D0/, BR/0.672D0/, BW/1.102D0/, BF/0.112D0/, 
     >     AF/-0.052D0/, RMN2/0.88172D0/, RMW2/.6146D0/, 
     >     RMF2/1.0384D0/, RMR2/0.5852D0/, GAMR/.112D0/, 
     >     PI/3.14159D0/, RMPI/.139D0/,RMPI2/.019321D0 /
      
      REAL*8 RMUP/2.792782D0/, RMUN/-1.913148D0 /

!     CONVERT TO FM**-2
      QQ  = QQG/(.197328D0)**2
      TAU = QQG/(4.D0*RMN2)

csk        2007 parametrizations: Proton Arrington et al., Gmn CLAS and Mainz, Madey Gen

      GEP = 0.0D0
      GEN = 0.0D0
      GMP = 0.0D0
      GMN = 0.0D0

      if (IG .eq. 21) then      ! "True" proton form factors
         GEP = (1.0D0+TAU*(3.439+TAU*(-1.602+TAU*0.068)))/
     *        (1.0D0+TAU*(15.055+TAU*(48.061+TAU*(99.304+TAU*(0.012+TAU*8.65)))))
         GMP = RMUP*(1.0D0+TAU*(-1.465+TAU*(1.26+TAU*0.262)))/
     *        (1.0D0+TAU*(9.627+TAU*TAU*TAU*(11.179+TAU*13.245)))
      else if (IG .eq. 22) then ! "Fudged" proton form factors to get cross section right
         GEP = (1.0D0+TAU*(-1.651+TAU*(1.287+TAU*(-0.185))))/
     *        (1.0D0+TAU*(9.531+TAU*(0.591+TAU*TAU*TAU*4.994)))
         GMP = RMUP*(1.0D0+TAU*(-2.151+TAU*(4.261+TAU*0.159)))/
     *        (1.0D0+TAU*(8.647+TAU*(0.001+TAU*(5.245+TAU*(82.817+TAU*14.191)))))
      else
         write(6,*) ' Wrong Farm Factor Model'
         return
      endif
      GMN = RMUN/(1.0D0+3.26*QQG/(1-0.272*QQG/(1+0.0123*QQG/(1-2.52*QQG/(1+2.55*QQG)))))
      GD = 1.D0/(1.D0+QQG/.71D0)**2
      GEN = 0.888D0*(-RMUN)*TAU*GD/(1.0D0+3.21*TAU)        

 900  RETURN

      END
********************************************************************************

********************************************************************************

      SUBROUTINE PAULI_SUPPnew(IPAULI,IA,QSQ,E0,PS1,PS2)
!-----------------------------------------------------------------------
!     Gets Pauli Suppression factor for quasi-elastic scattering from
!     Several different Models.
!     The VanOrden Model is very Slow.
!     Modified from S. Rock's version. 8/96, LMS.
!----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INTEGER IPAULI, IA
      REAL*8 QSQ, E0, TH, PS1, PS2, FDEL, Q, TAU, FF, GE, GM, 
     >     W1, W2, SQF, PS_R4

!     no suppression
      PS1 = 1.D0
      PS2 = 1.D0

!     Stein
      IF(IPAULI.EQ.1) THEN
         FF = 0.D0                                                           
         IF (IA.EQ.2.AND.QSQ.LE.8.) THEN
            SQF  = DSQRT(QSQ/.197328D0**2)                                         
            FF= 1.58D0/SQF*(DATAN(SQF/0.93D0)   
     >           - 2.D0*DATAN(SQF/3.19D0)+DATAN(SQF/5.45D0)) 
            IF (FF.LE.0.) FF = 0.D0    
         ENDIF                                     
         IF(IA.EQ.3) THEN
            IF(IA.EQ.3) THEN
               CALL FFHE3new(QSQ,GE,GM)
               TAU = QSQ/(4.D0*MHE3**2)
               W1 = TAU*GM**2                                              
               W2 = (GE**2+W1)/(1.D0+TAU) 
               FF = DSQRT(W2)/2.D0 ! iZ of 3He
            ENDIF
         ENDIF
         PS2 = (1.D0-FF**2)     !old model from Stein

!     Tsai RMP 46,816(74) eq.B54
      ELSEIF(IPAULI.EQ.2) THEN
         TAU = QSQ/(4.D0*MN**2)
         Q = DSQRT(QSQ*TAU+QSQ)
         IF(Q.LT.2.D0*PFERMI) THEN
            PS2=3.D0*Q*(1.D0-0.08333D0*(Q/PFERMI)**2)/(4.D0*PFERMI)
            PS1 = PS2
         ENDIF

!     DeForest and Walecka, Adv. in Phys. 15, 1 (1966).
      ELSEIF(IPAULI.EQ.3) THEN
         Q = DSQRT(QSQ)/PFERMI
         IF(Q.LT.2.D0) THEN
            PS2=0.75D0*Q - Q*Q*Q/16.D0
            PS1 = PS2
         ENDIF

!     Van Orden
      ELSEIF(IPAULI.EQ.4) THEN
         CALL  Q_E_VANORDENnew(QSQ,E0,PS_R4)
         PS2= PS_R4
         PS1 = PS2
      ENDIF
      
      RETURN
      END
********************************************************************************
*******************************************************************************

      SUBROUTINE Q_E_VANORDENnew(Q2_ELAS_GEV,E_GEV,SUPPRESSION)
!---------------------------------------------------------------------------
C     This program compute the quasi-elastic cross section
C     based on the Van Orden calculation using the fermi gas model
!     input energy now in GeV
!     It returns the Suppression factor for the quasi-elastic peak.
!-----------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'radcon.inc'

      INTEGER IOPT, I0, I00, I, IOPT1
      REAL*8 E_GEV, THETA, Z, A, KF_GEV, SUPPRESSION
      REAL*8 SUM, AP, FMT, FACTOR, N, MT, TH, SINSQ,
     >     FMOTT, WMOTT, TANSQ, EP_ELAS, EF, RLR, RTR,
     >     SRL, RATIO, W, WBAR, QV, Q2BAR, F2, XETA,
     >     ZETA, T1, T2, E1, E, SL, ST, RS, RSSIG, XSECT,
     >     CROSS, QT, C, Q2_ELAS, E2, W0, W1, Q2, EINC, KF,
     >     KF3, Z1, A1, QV2, WSQ, MT4PI, MP2, Q2_ELAS_GEV
      REAL*8 MUP2/7.784D0/, MUN2/3.648D0/
      REAL*8 EBAR/1.D0/, WDEL/2.D0/ !$$ are these OK?
      REAL*8 AMASS/931.5D0/, MP

!-----------------------------------------------------------------------------------------
!     In this notation, W=nu
!---------------------------------------------------------------------------------------
      EINC = 1000.D0*E_GEV
      Q2_ELAS = 1.D6 * Q2_ELAS_GEV
      MP = MN*1000.D0

      SUM = 0.D0

      IF(IOPT.EQ.1)THEN
 1       EP_ELAS = EINC - Q2_ELAS/(2.D0*MP)
         IF(EP_ELAS.GT.0) THEN
            SINSQ = Q2_ELAS/(4.D0*EINC*EP_ELAS)
            TH = 2.D0*DASIN(DSQRT(SINSQ))
         ELSE
            EINC = EINC + 5.D0
            GO TO 1
         ENDIF

         FMOTT = ALPHA*ALPHA*DCOS(TH/2.D0)**2/4.D0/SINSQ/SINSQ*HC2*1.D-24
         WMOTT = FMOTT/EINC**2
         TANSQ = DTAN(TH/2.D0)**2
         QT = DSQRT(Q2_ELAS**2/(4.D0*MP2) + Q2_ELAS)
         IF(QT.GT. 2.D0*KF) THEN
            SUPPRESSION = 1.D0
            RETURN
         ENDIF
         W0=MAX(EINC-(EP_ELAS+2.D0*KF),2.D0)
         W1=MAX(EINC-(EP_ELAS-2.D0*KF),0.D0)
      ENDIF
      
      RLR = 0.D0
      RTR = 0.D0
      SRL = 0.D0
      RATIO = 1.D0

      W = W0
      I00 = INT(W0/WDEL)
      I0 = INT(W1/WDEL)
      DO 17 I=I00,I0
         WBAR = W-EBAR
         IF(WBAR.LE.0.) GO TO 15
         IF(IOPT.EQ.1) Q2 = 4.D0*EINC*(EINC-W)*SINSQ
         WSQ = W**2
         QV2 = Q2 + WSQ
         QV = DSQRT(QV2)
         Q2BAR = QV2 - WBAR**2
         E1 = DSQRT(QV2*MP2/Q2BAR+QV2/4.D0) - WBAR/2.D0
         IF(EF.LT.E1)GO TO 18   ! do not calculate 
         RATIO = Q2/(Q2+WSQ)
         F2 = 1.D0/(1.D0 + Q2/855.D0**2)**4
         XETA = Q2/(4.D0*MP2)
         ZETA = 1.D0 + XETA
!     T1=F2*Q2/2.*(((1.+2.79*XETA)/ZETA+(1.79/ZETA))**2+N/Z*3.65)
!     T1 = 2Mp**2 * DIPOLE *Tau*(MuP2 + N/Z * MuN2) = Tau*(Gmp**2 + Gmn**2 )
!     = 2Mp*Tau*(F1+(Mu-1)F2)**2 where
!     F1= DIPOLE*(1+TAU*Mu)/(1+Tau)    F2= DIPOLE/(1+Tau)
!     T2=2.*MP22*(((1.+2.79*XETA)/ZETA)**2+XETA*((1.79/ZETA)**2+
!     >  N/Z*1.91**2))*F2  !$$$*** I think the neutron term should be divided by ZETA
!     T2=2Mp**2 *DIPOLE* (Gep +Tau*Gmp)/(1+Tau)  + neutron
!     =2Mp**2 * (F1**2 +Tau*(MuP-1)F2**2)
         
!     Below is Steve's Redoing
         T1 = F2*Q2/2.D0*(MUP2 +N/Z*MUN2)
         T2 = 2.D0*MP2*F2*
     >        ( (1.D0+MUP2*XETA)/ZETA +N/Z* (0.D0 + MUN2*XETA)/ZETA)
         E2 = EF-WBAR
         E = E1
         IF(E2.GT.E1) E = E2

         RLR = (.75D0*Z/(KF3*QV))*(T2*((EF**3 - E**3)/3.D0
     >        + W*(EF**2 - E**2)/2.D0 + WSQ*(EF - E)/4.D0)/MP2
     >        - QV2*T1*(EF - E)/Q2)

         RTR = (.75D0*Z/(KF3*QV))*(2.D0*T1*(EF - E) 
     >        + T2*(Q2BAR*(EF**3 - E**3)/(3.D0*QV2) + Q2BAR*WBAR*(EF**2 - 
     >        E**2)/(2.D0*QV2) - (Q2BAR**2/(4.D0*QV2)
     >        + MP2)*(EF - E))/MP2)
 15      CONTINUE
         SL = RLR*MT4PI
         ST = RTR*MT4PI

         RS = RATIO*RATIO*SL+(0.5D0*RATIO+TANSQ)*ST
         RSSIG = WMOTT*RS
         XSECT = FACTOR*WMOTT*RS
         IF(SL.EQ.0.AND.ST.EQ.0.)GO TO 18
C     SRL=SRL+RLR
         IF(IOPT.EQ.1)THEN
            SUM = SUM + XSECT* WDEL
         ENDIF
 18      W = W + WDEL
 17   CONTINUE

      F2 = 1.D0/(1.D0 + Q2_ELAS/855.D0**2)**4
      XETA = Q2_ELAS/(4.D0*MP2)
      ZETA = 1.D0 + XETA
      
      CROSS = WMOTT*F2*1.D33 *
     >     (Z*((1.D0 + MUP2*XETA)/ZETA + 2.D0*TANSQ*MUP2*XETA) +
     >     N*((0.D0 + MUN2*XETA)/ZETA + 2.D0*TANSQ*MUN2*XETA))
      SUPPRESSION = SUM/CROSS
      RETURN


      ENTRY  Q_E_VANORDEN_INITnew(Z1,A1,KF_GEV,IOPT1)  
      Z = Z1
      A = A1
      KF= 1000.D0*KF_GEV
      KF3 = KF**3
      IOPT = IOPT1
      AP = ALPHA/PI
      FMT = A*AMASS
      FACTOR = 4.D0*PI/FMT*1.D33
      N = A - Z
      MP2 = (MN*1000.D0)**2
      MT = 931.5D0*(Z + N)
      EF = DSQRT(KF**2 + MP2)
      MT4PI = MT/4.D0/PI
      RETURN
      END

********************************************************************************
*     New routines needed for 2014 F2n and F2p
********************************************************************************

      SUBROUTINE RESCSP(w2,q2,sigtp,siglp)
      IMPLICIT none

      real*8 w2,q2,sigtp,siglp
      real*8 xvalp(100),xval1(50),xvalL(50)
      Integer i
      real*8 sigtdis,sigLdis

      data xvalp / 

     &     0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     &     0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     &     0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
     &     0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01,
     &     0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     &     -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     &     0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     &     0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     &     -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     &     0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     &     0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     &     0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     &     0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     &     0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     &     0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     &     0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     &     0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     &     0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     &     0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     &     0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /

      do i=1,50
         xval1(i) = xvalp(i)
         xvalL(i) = xvalp(50+i) 
         if(i.LE.12) xvalL(i) = xval1(i) !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
      
      call resmodp(1,w2,q2,xval1,sigTp)
      call resmodp(2,w2,q2,xvalL,sigLp)
      call disp(w2,q2,sigtdis,sigLdis)
      if(w2.GT.9.0) then 
c     write(6,*) "resmod: ",w2,q2,sigTp,sigLp
         call disp(w2,q2,sigtdis,sigLdis)
c     write(6,*) "dismod: ",w2,q2,sigtdis,sigLdis
         sigTp = sigtdis
         sigLp = sigLdis
      endif

      return
      end

      SUBROUTINE DISP(w2,q2,sigt,sigl)
      IMPLICIT none

      real*8 w2,q2,x,sigt,sigl,f1,f2,fL,r,dr
      Integer i
      real*8 mp2,pi,pi2,alpha,t1,t2
      logical goodfit

      mp2 = 0.938272
      mp2 = mp2*mp2
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.03599

      x = q2/(w2+q2-mp2)
      call f2allm(x,q2,f2)
      call r1990(x,q2,r,dr,goodfit)
      f1 = f2/2./x/(r+1.0)*(1.0+4.0*mp2*x*x/q2)
      fL = 2.*x*r*f1    

      sigt = 0.3894e3*f1*pi2*alpha*8.0/abs(w2-mp2)
      sigL = r*sigt

      return
      end

c     allm97, NMC published measured points Q2>0.75 GeV2
c     for values Q<1 use data of E665!
c     parameterization of F2 , according to
c     H.Abramowicz and A.Levy, hep-ph/9712415
c     
c     3*10-6 < x  < 0.85, W2>3GeV2
c     0.   < Q2 < 5000 GeV2, dof=0.97
c     
      
      SUBROUTINE f2allm(x,q2,f2a)

      IMPLICIT NONE
      
      REAL*8 x,q2,M22,f2a
      REAL*8 SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
      REAL*8 S11,A11,B11,M12,S21,A21,B21,M02,LAM2,Q02
      REAL*8 S12,S13,A12,A13,B12,B13,S22,S23,A22,A23
      REAL*8 B22,B23,w2,w,z
      REAL*8 ALFA,XMP2
C     POMERON
      PARAMETER (
     ,     S11   =   0.28067, S12   =   0.22291, S13   =   2.1979,
     ,     A11   =  -0.0808 , A12   =  -0.44812, A13   =   1.1709,
     ,     B11   =   0.60243**2, B12   =   1.3754**2, B13   =   1.8439,
     ,     M12   =  49.457 )
      
C     REGGEON
      PARAMETER (
     ,     S21   =   0.80107, S22   =   0.97307, S23   =   3.4942,
     ,     A21   =   0.58400, A22   =   0.37888, A23   =   2.6063,
     ,     B21   =   0.10711**2, B22   =   1.9386**2, B23   =   0.49338,
     ,     M22   =   0.15052 )
C     
      PARAMETER ( M02=0.31985, LAM2=0.065270, Q02=0.46017 +LAM2 )         
      PARAMETER ( ALFA=112.2, XMP2=0.8802)
      
C     
      W2=q2*(1./x -1.)+xmp2
      W=sqrt(w2)
      
C     
      IF(Q2.EQ.0.) THEN
         S=0.
         Z=1.

C     
C     POMERON
C     
         XP=1./(1.+(W2-XMP2)/(Q2+M12))
         AP=A11
         BP=B11
         SP=S11
         F2P=SP*XP**AP
C     
C     REGGEON
C     
         XR=1./(1.+(W2-XMP2)/(Q2+M22))
         AR=A21
         BR=B21
         SR=S21
         F2R=SR*XR**AR
C     
      ELSE
         S=LOG(LOG((Q2+Q02)/LAM2)/LOG(Q02/LAM2))
         Z=1.-X   
C     
C     POMERON
C     
         XP=1./(1.+(W2-XMP2)/(Q2+M12))
         AP=A11+(A11-A12)*(1./(1.+S**A13)-1.)
         BP=B11+B12*S**B13
         SP=S11+(S11-S12)*(1./(1.+S**S13)-1.)
         F2P=SP*XP**AP*Z**BP
C     
C     REGGEON
C     
         XR=1./(1.+(W2-XMP2)/(Q2+M22))
         AR=A21+A22*S**A23
         BR=B21+B22*S**B23
         SR=S21+S22*S**S23
         F2R=SR*XR**AR*Z**BR
         
C     
      ENDIF                                     
      
      
      f2a = q2/(q2+m02)*(F2P+F2R)
      
      
      RETURN
      END                                  

**********************************************************************

      SUBROUTINE RESCSN(w2,q2,xvaln,sigtn,sigln)
      IMPLICIT none

      real*8 w2,q2,sigtn,sigln
      real*8 xvaln(100),xval1(50),xvalL(50)
      integer i

      do i=1,50
         xval1(i) = xvaln(i)
         xvalL(i) = xvaln(50+i) 
         if(i.LE.12) xvalL(i) = xval1(i) !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
      
      call resmodn(1,w2,q2,xval1,sigtn)
      call resmodn(2,w2,q2,xvalL,sigLn)

      return
      end
      
********************************************************************


      SUBROUTINE RESMODP(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
         q20 = 0.05
      else
         q20 = 0.125
      endif 
      

CCCC  single pion branching ratios  CCCC

      br(1,1) = 1.0             !!!  P33(1232)       
      br(2,1) = 0.45            !!!  S11(1535)   
      br(3,1) = 0.65            !!!  D13(1520)
      br(4,1) = 0.65            !!!  F15(1680)
      br(5,1) = 0.4             !!!  S11(1650)
      br(6,1) = 0.65            !!!  P11(1440) roper 
      br(7,1) = 0.50            !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0             !!!  P33(1232)
      br(2,3) = 0.45            !!!  S11(1535) 
      br(3,3) = 0.0             !!!  D13(1520)
      br(4,3) = 0.0             !!!  F15(1680)
      br(5,3) = 0.1             !!!  S11(1650)
      br(6,3) = 0.0             !!!  P11(1440) roper   
      br(7,3) = 0.0             !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
         br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC  Meson angular momentum   CCCC



      ang(1) = 1.               !!!  P33(1232)
      ang(2) = 0.               !!!  S11(1535)
      ang(3) = 2.               !!!  D13(1520)
      ang(4) = 3.               !!!  F15(1680)
      ang(5) = 0.               !!!  S15(1650)
      ang(6) = 1.               !!!  P11(1440) roper   
      ang(7) = 3.               !!!  F37(1950)

      do i=1,7                  !!!  resonance damping parameter  !!!
         x0(i) = 0.215
c     x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
         br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
      
      
      dip = 1./(1.+q2/0.71)**2. !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC   Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6                  !!!  Read in resonance masses     !!!
         num = num + 1
         mass(i) = xval(i)
      enddo
      do i=1,6                  !!!  Read in resonance widths     !!!
         num = num + 1
         intwidth(i) = xval(num)
         width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then          !!!  Put in 4th resonance region  !!!
         mass(7) = xval(43)
         intwidth(7) = xval(44)
         width(7) = intwidth(7)
      else
         mass(7) = xval(47)
         intwidth(7) = xval(48)
         width(7) = intwidth(7) 
      endif

      do i=1,7
         kr(i) = (mass(i)**2-mp2)/2./mp
         kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
         epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
         ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
         epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
         ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
         eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
         petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

         pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &        *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c     !!!  1-pion decay mode


         pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &        *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &        **(ang(i)+2)      !!!  2-pion decay mode

         pwid(i,2) = W/mass(i)*pwid(i,2)


         pwid(i,3) = 0.         !!!  eta decay mode


         if(i.EQ.2.OR.i.EQ.5) then
            pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &           *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c     !!!  eta decay only for S11's 
         endif 



         pgam(i) = (kcm/kcmr(i))**2*
     &        (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

         pgam(i) = intwidth(i)*pgam(i)

         width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC   End resonance kinematics and Widths calculations   CCC


CCC   Begin resonance Q^2 dependence calculations   CCC
      

      do i=1,6
         do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
         enddo

         if(sf.EQ.1) then

            height(i) = rescoef(i,1)*
     &           (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &           (1.+q2/0.91)**rescoef(i,4)
            
         else

            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &           *exp(-1.*rescoef(i,3)*q2)

         endif

         
         height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then          !!!  4th resonance region  !!!

         height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
         height(7) = xval(49)/(1.+q2/0.91)**1. 

      endif
      height(7) = height(7)*height(7)



CCC   End resonance Q^2 dependence calculations   CCC

      
      do i=1,3                  !!!  Non-Res coefficients  !!!
         do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
         enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
         sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &        + (mass(i)*width(i))**2.)
         sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
         if(sf.eq.1) sigrsv(i) = sigr(i)
         sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC   Finish resonances / start non-res background calculation   CCC

      
      sig_nr = 0.

      if(sf.EQ.1) then

         do i=1,2  

            h_nr(i) = nr_coef(i,1)/     
     &           (q2+nr_coef(i,2))**
     &           (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
            sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
         enddo

         sig_nr = sig_nr*xpr(1)
         sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

         do i=1,1
            sig_nr = sig_nr + nr_coef(i,1)*
     &           (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &           /(1.-xb)/(q2+q20)
     &           *(q2/(q2+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
         enddo

      endif


      sig = sig_res + sig_nr
      
      if(w2.LT.1.159) sig = 0.0

 1000 format(8f12.5)

      RETURN 
      END 

*****************************************************************************

      SUBROUTINE RESMODN(sf,w2,q2,xval,sig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,xth(4),sig,xval(50),mass(7),width(7)
      REAL*8 height(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3)
      REAL*8 sig_res,sig_4L,sigtemp,slope,t,xpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.939565

c     mp = 0.938272

      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)
      if(w.LT.(mp+mpi)) wdif(1) = 0.0
      if(w.LT.(mp+mpi+2.*mpi)) wdif(2) = 0.0

c     write(6,*) "here"

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
         q20 = 0.05
      else
         q20 = 0.125
      endif 
      

CCCC  single pion branching ratios  CCCC

      br(1,1) = 1.0             !!!  P33(1232)       
      br(2,1) = 0.45            !!!  S11(1535)   
      br(3,1) = 0.65            !!!  D13(1520)
      br(4,1) = 0.65            !!!  F15(1680)
      br(5,1) = 0.4             !!!  S11(1650)
      br(6,1) = 0.65            !!!  P11(1440) roper 
      br(7,1) = 0.50            !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0             !!!  P33(1232)
      br(2,3) = 0.45            !!!  S11(1535) 
      br(3,3) = 0.0             !!!  D13(1520)
      br(4,3) = 0.0             !!!  F15(1680)
      br(5,3) = 0.1             !!!  S11(1650)
      br(6,3) = 0.0             !!!  P11(1440) roper   
      br(7,3) = 0.0             !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
         br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC  Meson angular momentum   CCCC



      ang(1) = 1.               !!!  P33(1232)
      ang(2) = 0.               !!!  S11(1535)
      ang(3) = 2.               !!!  D13(1520)
      ang(4) = 3.               !!!  F15(1680)
      ang(5) = 0.               !!!  S15(1650)
      ang(6) = 1.               !!!  P11(1440) roper   
      ang(7) = 3.               !!!  F37(1950)

      do i=1,7                  !!!  resonance damping parameter  !!!
         x0(i) = 0.215
c     x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
         br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
      
      
      dip = 1./(1.+q2/0.71)**2. !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/0.71)**1.

      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+2.*mpi)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LT.(mp+mpi)) xpr(1) = 0.0
      if(w.LT.(mp+mpi+2.*mpi)) xpr(2) = 0.0

      t = log(log((q2+m0)/0.330**2)/log(m0/0.330**2))

CCC   Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6                  !!!  Read in resonance masses     !!!
         num = num + 1
         mass(i) = xval(i)
      enddo
      do i=1,6                  !!!  Read in resonance widths     !!!
         num = num + 1
         intwidth(i) = xval(num)
         width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then          !!!  Put in 4th resonance region  !!!
         mass(7) = xval(43)
         intwidth(7) = xval(44)
         width(7) = intwidth(7)
      else
         mass(7) = xval(47)
         intwidth(7) = xval(48)
         width(7) = intwidth(7) 
      endif

      do i=1,7
         kr(i) = (mass(i)**2-mp2)/2./mp
         kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
         epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
         ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
         epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
         ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
         eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
         petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

         pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &        *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c     !!!  1-pion decay mode


         pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &        *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &        **(ang(i)+2)      !!!  2-pion decay mode

         pwid(i,2) = W/mass(i)*pwid(i,2)


         pwid(i,3) = 0.         !!!  eta decay mode


         if(i.EQ.2.OR.i.EQ.5) then
            pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &           *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c     !!!  eta decay only for S11's 
         endif 



         pgam(i) = (kcm/kcmr(i))**2*
     &        (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

         pgam(i) = intwidth(i)*pgam(i)

         width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC   End resonance kinematics and Widths calculations   CCC


CCC   Begin resonance Q^2 dependence calculations   CCC
      

      do i=1,6
         do j=1,4
            num = num + 1
            rescoef(i,j)=xval(num)
         enddo

         if(sf.EQ.1) then

            height(i) = rescoef(i,1)*
     &           (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))/
     &           (1.+q2/0.91)**rescoef(i,4)
            
         else

            height(i) = rescoef(i,1)*q2/(1.+rescoef(i,2)*q2)
     &           *exp(-1.*rescoef(i,3)*q2)

         endif

         
         height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then          !!!  4th resonance region  !!!

         height(7) = xval(45)*q2/(1.+xval(46)*q2)*exp(-1.*xval(47)*q2)

      else
         height(7) = xval(49)/(1.+q2/0.91)**1. 

      endif
      height(7) = height(7)*height(7)



CCC   End resonance Q^2 dependence calculations   CCC

      
      do i=1,3                  !!!  Non-Res coefficients  !!!
         do j=1,4
            num = num + 1
            nr_coef(i,j)=xval(num)
         enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7
         sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &        + (mass(i)*width(i))**2.)
         sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
         if(sf.eq.1) sigrsv(i) = sigr(i)
         sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w


CCC   Finish resonances / start non-res background calculation   CCC

      
      sig_nr = 0.

      if(sf.EQ.1) then

         do i=1,2  

            h_nr(i) = nr_coef(i,1)/     
     &           (q2+nr_coef(i,2))**
     &           (nr_coef(i,3)+nr_coef(i,4)*q2+xval(44+i)*q2**2)
            sig_nr = sig_nr +h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
         enddo

         sig_nr = sig_nr*xpr(1)
         sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

         do i=1,1
            sig_nr = sig_nr + nr_coef(i,1)*
     &           (1.-xb)**(nr_coef(i,2)+0.0*t )
c     &      (1.-xpr(1))**(nr_coef(i,2)+0.0*t )
c     &        *xpr(1)**(xval(41)+xval(42)*t)
     &           *xb**(xval(41)+xval(42)*t)
     &           *(1./(1.+q2/nr_coef(i,3)))**nr_coef(i,4)
         enddo

      endif


      sig_res = abs(sig_res)
      sig_nr = abs(sig_nr)

      sig = sig_res + sig_nr
      
      if(w.LT.(mp+mpi)) sig = 0.0

 1000 format(8f12.5)

      RETURN 
      END 
      
***********************************************************************************

!     File R1990.FORTRN.                                                       
!     Reference:  L.W.Whitlow, SLAC-Report-357,                                      
!     Ph.D. Thesis, Stanford University,                                 
!     March 1990.                                                        
!     For details see file HELP.DOCUMENT.                                            
      
!     Program contains 135 lines of Fortran code, of 72 characters each, with        
!     one subroutine.                                                                
      
      
      SUBROUTINE R1990(X,Q2,R,DR,GOODFIT)                                       
      
!     Model for R, based on a fit to world R measurements. Fit performed by         
!     program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details           
!     see Reference.                                                                
!     
!     Three models are used, each model has three free parameters.  The             
!     functional forms of the models are phenomenological and somewhat              
!     contrived.  Each model fits the data very well, and the average of            
!     the fits is returned.  The standard deviation of the fit values is            
!     used to estimate the systematic uncertainty due to model dependence.          
!     
!     Statistical uncertainties due to fluctuations in measured values have         
!     have been studied extensively.  A parametrization of the statistical          
!     uncertainty of R1990 is presented in FUNCTION DR1990.                         
!     
!     The three model forms are given by:                                           
!     
!     R_A = A(1)/LOG(Q2/.04)*FAC + A(2)/[Q2„4+A(3)„4]„.25 ;                     
!     R_B = B(1)/LOG(Q2/.04)*FAC + B(2)/Q2 + B(3)/(Q2**2+.3**2) ;               
!     R_C = C(1)/LOG(Q2/.04)*FAC + C(2)/[(Q2-Q2thr)„2+C(3)„2]„.5 ,              
!     ...where Q2thr = 5(1-X)„5 ;                     
!     where FAC = 1+12[Q2/(1+Q2)][.125„2/(.125„2+x„2)] gives the          
!     x-dependence of the logarithmic part in order to match Rqcd         
!     at high Q2.                                                         
!     
!     Each model fits very well.  As each model has its own strong points           
!     and drawbacks, R1990 returns the average of the models.  The                  
!     chisquare for each fit (124 degrees of freedom) are:                          
!     R_A: 110,    R_B: 110,    R_C: 114,    R1990(=avg): 108                   
!     
!     This subroutine returns reasonable values for R for all x and for all         
!     Q2 greater than or equal to .3 GeV.                                           
!     
!     The uncertainty in R originates in three sources:                             
!     
!     D1 = uncertainty in R due to statistical fluctuations of the data         
!     and is parameterized in FUNCTION DR1990, for details see             
!     Reference.                                                           
!     
!     D2 = uncertainty in R due to possible model dependence, approxi-          
!     mated by the variance between the models.                            
!     
!     D3 = uncertainty in R due to possible epsilon dependent errors            
!     in the radiative corrections, taken to be +/- .025.  See             
!     theses (mine or Dasu's) for details.                                 
!     
!     and the total error is returned by the program:                               
!     
!     DR = is the total uncertainty in R, DR = sqrt(D1„2+D2„2+D3„2).            
!     DR is my best estimate of how well we have measured R.  At           
!     high Q2, where R is small, DR is typically larger than R.  If        
!     you have faith in QCD, then, since R1990 = Rqcd at high Q2,          
!     you might wish to assume DR = 0 at very high Q2.                     
!     
!     NOTE:    In many applications, for example the extraction of F2 from          
!     measured cross section, you do not want the full error in R          
!     given by DR.  Rather, you will want to use only the D1 and D2        
!     contributions, and the D3 contribution from radiative                
!     corrections propogates complexely into F2.  For more informa-        
!     tion, see the documentation to dFRC in HELP.DOCUMENT, or             
!     for explicite detail, see Reference.                                 
!     
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      IMPLICIT NONE                                                             
      REAL*8 FAC,RLOG,Q2THR,R_A,R_B,R_C,R, D1,D2,D3,DR,DR1990,X,Q2                
      REAL*8 A(3), B(3), C(3)                                    
      LOGICAL GOODFIT                                                           
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      DATA A/ .06723, .46714, 1.89794 /                                  
      DATA B/ .06347, .57468, -.35342 /                                   
      DATA C/ .05992, .50885, 2.10807 /
      
      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(X**2+.125**2))                       
      RLOG  = FAC/LOG(Q2/.04)   !   <--- we use natural logarithms only!           
      Q2thr = 5.*(1.-X)**5                                                      
      
      R_A   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))                        
      R_B   = B(1)*RLOG + B(2)/Q2 + B(3)/(Q2**2+.3**2)                          
      R_C   = C(1)*RLOG + C(2)/SQRT((Q2-Q2thr)**2+C(3)**2)                      
      R     = (R_A+R_B+R_C)/3.                                                  
      
      D1    = DR1990(X,Q2)                                                      
      D2    = SQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.)                       
      D3    = .023*(1.+.5*R)                                                    
      IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3                               
      DR    = SQRT(D1**2+D2**2+D3**2)                                           
      
      GOODFIT = .TRUE.                                                          
      IF (Q2.LT..3) GOODFIT = .FALSE.                                           
      RETURN                                                                    
      END                                                                       
      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      FUNCTION DR1990(X,Q2)                                                     
      
!     Parameterizes the uncertainty in R1990 due to the statistical                 
!     fluctuations in the data.  Values reflect an average of the R-values          
!     about a neighborhood of the specific (x,Q2) value.  That neighborhood         
!     is of size [+/-.05] in x, and [+/-33%] in Q2.  For details, see               
!     Reference.                                                                    
!     
!     This subroutine is accurate over all (x,Q2), not only the SLAC deep           
!     inelastic range.  Where there is no data, for example in the resonance        
!     region, it returns a realistic uncertainty, extrapolated from the deep        
!     inelastic region (suitably enlarged).  We similarly estimate the              
!     uncertainty at very large Q2 by extrapolating from the highest Q2             
!     measurments.  For extremely large Q2, R is expected to fall to zero,          
!     so the uncertainty in R should not continue to grow.  For this reason         
!     DR1990 uses the value at 64 GeV for all larger Q2.                            
!     
!     XHIGH accounts for the rapidly diminishing statistical accuracy for           
!     x>.8, and does not contribute for smaller x.                                  
      
      
      IMPLICIT NONE                                                             
      REAL*8 U(10,10),DR1990,QMAX,Q,S,A,XLOW,XHIGH,X,Q2                           
      
      
      QMAX = 64.                                                                
      
      Q = MIN(Q2,QMAX)                                                          
      S = .006+.03*X**2                                                         
      A = MAX(.05,8.33*X-.66)                                                   
      
      XLOW  = .020+ABS(S*LOG(Q/A))                                              
      XHIGH = .1*MAX(.1,X)**20/(.86**20+MAX(.1,X)**20)                          
      
      DR1990 = SQRT(XLOW**2+XHIGH**2)                                           
      RETURN                                                                    
      END                                                                       
!     
!     End of file R1990.FORTRN.  135 Fortran lines.                                  


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
c     sk
      if (IA1.ge.4 .and. IA1.le.6 .and. X.gt.0.D0) then
         Wsq = 0.939D0*0.939D0 + Q2*(1.D0/X - 1.D0)
         call sek_me(Wsq, Q2, TARG, SIGNP, A1, delA1) ! 2006-7 refit including error
         if (IA1.eq.5) A1 = A1 + delA1
         if (IA1.eq.6) A1 = A1 - delA1
c     sk    New "Millenium Edition" based on new fit to ELSA, E142-155, EMC/SMC/Compass, HERMES,
c     sk	Hall A and EG1a/b
c     sk    and new definition of xi. Options 5 and 6 to evaluate uncertainty due to A1 model
c     sk    based on MINUIT error matrix. SEK 1-Nov-2006 updated 1-Oct-2007

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
c     kag       these are the polarization numbers from the nuclear physics of 3He
         ENDIF
c     kag         GMP = GMP*DSQRT(Z*PS2)
c     kag         GEP = GEP*DSQRT(Z*PS2)
c     kag         GMN = GMN*DSQRT(N*PS2)
c     kag         GEN = GEN*DSQRT(N*PS2)
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
******************************************************************************
      IMPLICIT NONE

      INCLUDE 'radcon.inc'
      INCLUDE 'instruct.inc'

      REAL*8 X, Q2, F1, F2, W1R, W2R, W1I, W2I, NU, ST, SY, SLOPE, 
     >     DSLOPE, R, DR, WSQ, FRAC, SIN2, F2D, WD1, WD2, ERR1, 
     >     ERR2, ERR3, W1ERR, W2ERR, ERROR1, ERROR2, SIGNP, W, FMAX
      
      REAL*8 Xeff, Q2eff, pival, W1n, W2n, W1Rnres, W2Rnres, W1Dnres, 
     >     W2Dnres, SQRTWSQ, F1res, F2res, Rtemp, dRtemp
      REAL*8 A,B,X2,X3
      REAL*8 Z_targ, A_targ
      
      real x4,qsq4,f2smc,err_lo,err_hi
      real*8 SIGTOT, ESIGTOT,EF2
      integer iT, iT2, iflag    ! For f2allm; =1 if Q2>0
      
      integer I, J, IQ2eg1, IWeg1
      real*8 Rin(2), ResR(40,203,2), FRACWHI, FRACWLO, FraclgQ2hi, FraclgQ2lo
      

      CHARACTER*1 TARG
      LOGICAL GOOD, DORES, DOINEL, SUPPRESS, Filled /.FALSE./

      pival = 3.1415927
      F1 = 0.D0
      F2 = 0.D0
      SIGNP = 0.D0
      if(X .le. 0.0D0) return
      NU = Q2/(2.D0*MN*X)       ! NUCLEON mass
      WSQ = MN2 + 2.D0*MN*NU - Q2 
      
c11   if(x.gt.1.0D0.or.Q2.le.0.0D0.or.Q2.gt.1.0D4.or.x.le.0.0D0)then
      if(WSQ.le.0.0D0 .or. Q2.le.0.0D0 .or. Q2.gt.1.0D4)then
c     write(6,*) ' Wrong kinematics: ', x,Q2,WSQ,NU
         return
      endif
c     sk       
c     sk	Latest fit by Peter Bosted / Eric Christie 2009

      if (.not. Filled) then
         Filled = .TRUE.
         do IQ2eg1 = 1, 40
	    do IWeg1 = 1, 203
               do J = 1, 2
                  ResR(IQ2eg1,IWeg1,J) = 0.0D0
               enddo
	    enddo
         enddo
         OPEN (UNIT=77, FILE='resaR_fit_final.dat', STATUS='OLD')
         do I = 1, 3840
            read(77,*) IQ2eg1, IWeg1, (Rin(J), J=1,2)
            if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
               do J = 1,2
                  ResR(IQ2eg1,IWeg1,J) = Rin(J)
               enddo
            endif
         enddo
         CLOSE(UNIT=77)
      endif
c11   NOTE: This prefilled array of R values is a kludge stemming from a time where the coded R didn't 
c11   work properly for all kinematics. It's not clear that it is still needed. For compatibility reasons
c11   we keep it for now but replace it for SFChoice = 14,15,16 and 17.


      DOINEL = .FALSE.
      DORES = .FALSE. 
      iT2 = 0
      iflag = 1

      if(SFChoice .lt. 11 .or. SFChoice .gt. 19) then
c     write(6,*) ' Wrong Choice: ',AsymChoice, SFChoice, TARG
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
c     write(6,*) ' Wrong Target ', TARG
         return
      endif

      IF(Q2.LT.10.D0.AND.WSQ.LT.9.0D0.AND..NOT.
     >     SUPPRESS) DORES = .TRUE.
      IF(.NOT.DORES.OR.WSQ.GT.6.0D0 .or. Q2.gt.4.0D0) DOINEL = .TRUE.
      
      if(DORES) then
         if(WSQ .le. 1.155 .and. A_targ .lt. 1.5) return
c11   
         if(SFChoice .lt. 17) then
	    call F1F2IN09(Z_targ, A_targ, Q2, WSQ, F1, F2, R)
         else
	    call F1F2IN09dr(Z_targ, A_targ, Q2, WSQ, F1, F2, R)
         endif
c     write(6,*) ' Called F1F2IN09'
         Rtemp = R
         dRtemp = 0.1D0*dabs(R) + 0.02D0
         W = DSQRT(WSQ)
         IWeg1 = W*1.0D2-0.5D0
         dR = dRtemp
         if (IWeg1 .lt. 107) goto 444
         if ((SFChoice .lt. 14).and.(IWeg1. le. 202)) then
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
     >           FRACWHI*FraclgQ2lo*ResR(IQ2eg1,IWeg1+1,2) +
     >           FRACWLO*FraclgQ2hi*ResR(IQ2eg1+1,IWeg1,2) +
     >           FRACWHI*FraclgQ2hi*ResR(IQ2eg1+1,IWeg1+1,2)
            dR = FRACWLO*FraclgQ2lo*ResR(IQ2eg1,IWeg1,1) +
     >           FRACWHI*FraclgQ2lo*ResR(IQ2eg1,IWeg1+1,1) +
     >           FRACWLO*FraclgQ2hi*ResR(IQ2eg1+1,IWeg1,1) +
     >           FRACWHI*FraclgQ2hi*ResR(IQ2eg1+1,IWeg1+1,1)
	    dR = 0.5*dabs(R-dR)+0.1D0*dabs(R) + 0.02D0
         endif                  !(of existing kludge for SFChoice = 11-13)
 444     continue
         if (R .lt. 0.0D0) R = 0.0D0
         if ((R .lt. 0.02).and.(R .lt. Q2)) then ! Give R a minimum value to avoid A2Soffer = 0
	    if (Q2 .lt. 0.02D0) then
               R = Q2
	    else
               R = 0.02D0
	    endif
         endif
         if ((SFChoice .eq. 13) .or.(SFChoice .eq. 16) .or. (SFChoice.eq.19)) then
	    if (R .lt. dR) then
               R = R + dR       !
	    else
               R = R - dR
	    endif
         endif
         F1 = (4.D0*Mn*Mn*X*X/Q2+1.D0)*F2/2.0D0/X/(1+R)

         if(A_targ .gt. 1.5) then
	    if(SFChoice .lt. 17) then
               call F1F2QE09(Z_targ, A_targ, Q2, WSQ, F1res, F2res)
	    else
               call F1F2QE09dr(Z_targ, A_targ, Q2, WSQ, F1res, F2res)
	    endif
c     if(F1res .lt. 0) write(6,*) Z_targ,A_targ,Q2,WSQ,F1res,F2res
	    if (F1res .gt. 0.0) F1 = F1+F1res
	    if (F2res .gt. 0.0) F2 = F2+F2res
         endif
         if ((SFChoice.eq.12) .or. (SFChoice.eq.15) .or. (SFChoice.eq. 18))  then
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
         if ((SFChoice.eq.12) .or. (SFChoice.eq.15) .or. (SFChoice.eq. 18))  then
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
         endif                  ! TARG .ne. 'P'
         CALL R1998(X,Q2,R,DR,GOOD)
         if ((SFChoice .eq. 13) .or.(SFChoice .eq. 16) .or. (SFChoice.eq.19)) then
	    if (R .lt. dR) then
               R = R + dR       !
	    else
               R = R - dR
	    endif
         endif
         F1 = (4.D0*Mn*Mn*X*X/Q2+1.D0)*F2/2.0D0/X/(1+R)
         if (DORES) then
	    FRAC = 0.0
	    if (Q2 .gt. 4.0D0) FRAC = (Q2 - 4.0D0)/6.0D0
	    if (WSQ .gt. 6.0D0) FRAC = FRAC + (1.0D0 - FRAC)*(WSQ - 6.0D0)/(9.0D0 - 6.0D0)
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
c     sk       EXPW = EXP(DELLOGX/2.D0)
c     sk       EXPW = EXPW - 1.D0/EXPW

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
c     sk           DELX = XX*EXPW

            GAMMA2 = GAM2*XX*XX
c     sk           RGAMMA2 = DSQRT(GAMMA2)
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
****************************************************************************** 7/30/13
*     
*     kpa: 7/30/13: added kpaVarChanges.inc to control the change of A1 in G1G2new(..) of newSFsSEK.f by 
*     a known amount such as 0.1. This same file is used in set_things_up.f to control 
*     the values of AsymChoice & SFchoice (which used be obtained from rcslacpol.file
*     in previous versions with each taking 11 for standard simulations), but our
*     final g1 extraction & systematic error estimations, we want to repeat the 
*     simulations by varying A1 to A1+0.1, AsymChoice from 11 to 12 & 15 and SFchoice 
*     from 11 to 12 & 13 (meaning 6 separate simulations - only one thing changed, not 
*     one change on top of another)
*
****************************************************************************** 7/30/13
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'kpaVarChanges.inc' !kpa: 7/30/13  (to control change in A1)
c      REAL*8 changeA1 

      REAL*8 X, Q2, G1, G2, A1, A2, W2, NU, TH, XMAX, A1n, A2n, A1p,A2p

      REAL*8 MP/0.939D0/, F1, F2, F1NR, F1R, F2R, F3R, SIGTOTERR, !changed MP
     >     SIG_NRES(2), SIG_RES1(2), SIG_RES2(2), SIG_RES3(2),
     *     SIGd_NRES(2), SIGd_RES1(2), SIGd_RES2(2), SIGd_RES3(2),
     *     SIGTOTd, SIGTOTn, SIGRESd, SIGRESn, SIGRESp,
     >     SIGTOT, R, R1, R2, R3, SNPB/0.40D0/, SNPR/1.D0/, 
     >     DUM1, GAMMA2, G1F1, RGAMMA2, TERM, DIL_F2, DUMMY, F1p, F1n

      INTEGER IWdmt, IQ2dmt, IQ2eg1, IWeg1
      REAL*8 W, DELW, FRACWLO, FRACWHI, FRACQ2LO, FRACQ2HI, FWLO, FWHI
      
      REAL*8 ResA1p(40,203,3),asymin(3),A1pRESnew,FraclgQ2hi,FraclgQ2lo,
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
         OPEN (UNIT=77, FILE='resa1p_fit_final.dat', STATUS='OLD')
         OPEN (UNIT=78, FILE='resa1n_fit_final.dat', STATUS='OLD')
         OPEN (UNIT=79, FILE='resa2p_fit_final.dat', STATUS='OLD')
         OPEN (UNIT=80, FILE='resa2n_fit_final.dat', STATUS='OLD')
         do I = 1, 3840
	    read(77,*) IQ2eg1, IWeg1, (asymin(J), J=1,3)
	    if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
               do J = 1,3
                  ResA1p(IQ2eg1,IWeg1,J) = asymin(J)
               enddo
	    endif
	    read(78,*) IQ2eg1, IWeg1, (asymin(J), J=1,3)
	    if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
               do J = 1,3
                  ResA1n(IQ2eg1,IWeg1,J) = asymin(J)
               enddo
	    endif
	    read(79,*) IQ2eg1, IWeg1, (asymin(J), J=1,3)
	    if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
               do J = 1,3
                  ResA2p(IQ2eg1,IWeg1,J) = asymin(J)
               enddo
	    endif
	    read(80,*) IQ2eg1, IWeg1, (asymin(J), J=1,3)
	    if((1.le.IQ2eg1).and.(IQ2eg1.le.40).and.(108.le.IWeg1).and.(IWeg1.le.203)) then
               do J = 1,3
                  ResA2n(IQ2eg1,IWeg1,J) = asymin(J)
               enddo
	    endif
         enddo
         CLOSE(UNIT=77)
         CLOSE(UNIT=78)
         CLOSE(UNIT=79)
         CLOSE(UNIT=80)
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
         
         XMAX = Q2/(4.0D0 - MP*MP + Q2)	! Limit of W=2
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
c     sk	This will be passed as A2_WW(xRR) to G1G2
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

c     =========================== kpa: 7/30/13  (see kpaVarChanges.inc)
      A1 = A1 + changeA1
c     =========================== kpa: 7/30/13
      
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
c     kag            fits to e155x data
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
*   June 2017 SEK: many changes, including treatment of D and elastic p
******************************************************************************
      IMPLICIT NONE

      INCLUDE 'radcon.inc'
      INCLUDE 'instruct.inc'

      REAL*8 E,EP,TH, X, Q2, F1, F2, W1, W2, NU, G1, G2, 
     >     G1X, G2X, WSQ, DUM1, DUM2, A1, A2, SIN2, SIGNP
           INTEGER YoniIndex    ! CSK


           CHARACTER*1 TARG
           LOGICAL GOOD, SUPPRESS

           W1 = 0.D0
           W2 = 0.D0
           G1 = 0.D0
           G2 = 0.D0

           SIN2 = (DSIN(TH*RADCON/2.D0))**2
           Q2 = 4.D0*E*EP*SIN2
           if(Q2 .lt. 0.0001) return
           NU = E - EP
           X = Q2/(2.D0*MN*NU)
           if (X .lt. 0.0001) return
           if((TARG.ne.'D').AND.(TARG.ne.'H')) then ! Implement pseudo proton el. SF
              if(X.gt.0.99) return
           else
              if(X.gt. 1.98) return
           endif
           WSQ = MN2 + 2.D0*MN*NU - Q2
           if(WSQ.le.0.0D0) return


            SUPPRESS = .FALSE.
            IF(NORES) SUPPRESS = .TRUE.

           IF(WSQ.LT.1.15D0.AND.TARG.EQ.'P') THEN
C   SEK2017: In this case we treat the proton with a fake "SF"
            CALL ProFake(X,Q2,F1,F2,SIGNP,G1X,G2X,A1,A2,suppress)
            RETURN
           ENDIF
c     IF(WSQ.LT.1.15D0.AND.TARG.NE.'D') RETURN 
           
c     IF(WSQ.GT.4.0) RETURN  
            IF(FL_POL.AND..NOT.FL_UNPOL) GOTO 35

C     SK
           if(TARG.ne.'D')then
              CALL F1F2new(X,Q2,TARG,F1,F2,SIGNP,SUPPRESS)
           else
              YoniIndex=2
              call DSFs(X,Q2,YoniIndex,F1,F2,SIGNP,G1X,G2X,A1,A2,suppress)
           endif
C     SK	
           W1 = F1/MN
           W2 = F2/NU


 35        CONTINUE
           IF(FL_POL) THEN

C     SK
              if(TARG.ne.'D')then       
                 CALL G1G2new(X,Q2,TH,TARG,G1X,G2X,A1,A2)
              else
                 YoniIndex=2
                 call DSFs(X,Q2,YoniIndex,F1,F2,SIGNP,G1X,G2X,A1,A2,suppress)
              endif
C     SK
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
c     sk	CHANGED 3-Oct-2007 following Arrington's paper
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

c     sk	2007 parametrizations: Proton Arrington et al., Gmn CLAS and Mainz, Madey Gen

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

 900          RETURN

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

        SUBROUTINE ProFake(X, Q2, F1, F2, R, G1, G2, A1, A2, suppress)
*******************************************************************************
*
* SEK 2017: Subroutine to fake elastic proton through equivalent SFs
* 
******************************************************************************
       IMPLICIT NONE

        include 'instruct.inc'

        REAL*8 X, Q2, W, G1n, G1p, G2n, G2p, F1n, F1p, F2n, F2p
        real*8  F1, F2, R, G1, G2, A1, A2, Nu, QoverNu, GEP,GEN,GMP,GMN
        real*8 deltaX
csk17
        integer IG
        CHARACTER*1 TARG, TARGP /'P'/, TARGN /'N'/
        logical suppress

        IG = 22
        deltaX = 0.01*0.939/Q2
        if(deltaX .lt. 0.01) deltaX = 0.01
        if(deltaX .gt. 0.5) deltaX = 0.5
csk17   Keep within W resolution of about +/- 5 MeV
        call NewFORM(IG,Q2,GEP,GEN,GMP,GMN)
        if((x.GT.(1.0 - deltaX)).AND.(X.LT.(1.0 + deltaX)))THEN
          F1 = 0.5D0*GMP**2/(2.0*deltaX)
        else
          F1 = 0.0
          F2 = 0.0
          R = 0.0
          G1 = 0.0
          G2 = 0.0
          A1 = 1.0
          A2 = 0.0
          RETURN
        endif
	    QoverNu = 4.0D0*0.93827D0*0.93827D0/Q2
        R = QoverNu*GEP**2/GMP**2
        F2 = 2.0D0*F1*(R + 1.0D0)/(QoverNu + 1.0D0)
        A1 = 1.0D0

        IG = 21
        call NewFORM(IG,Q2,GEP,GEN,GMP,GMN)
        A2 = dsqrt(QoverNu)*GEP/GMP

        G1 = F1*(A1+DSQRT(QoverNu)*A2)/(1.0D0 + QoverNu)
        G2 = F1*(A2/DSQRT(QoverNu)-A1)/(1.0D0 + QoverNu)

        RETURN
        END

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
      I00 = W0/WDEL
      I0 = W1/WDEL
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

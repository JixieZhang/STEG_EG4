      FUNCTION SIGMAR(ERC,EPRC,THETA)
*******************************************************************************
*     
*******************************************************************************
C     Jixie: 20190108, Remove FakeProton
C     Further kludged version June 2017 SEK to allow for all W radiated proton el
*******************************************************************************
C     
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      INCLUDE 'prinplot.inc'
      
      LOGICAL OK, EXTERNALv     ! If false only internal done.
      LOGICAL SAVEv
      INTEGER NPTS, ITAIL, NEXT
      REAL*8 SIG, DELINF, SIGMAR, THETA, ERC, EPRC
      CHARACTER*60 TITLE
      CHARACTER*8 TARGETS(3)/' PROTON ','NEUTRON ','DEUTERON'/
      COMMON /INFO/EXTERNALv, NEXT

      E = ERC
      EP = EPRC

      THETAD = THETA
!     print*,'sigmarNEW_SEK.f: THETA=',THETA   !kp: 6/1/12

!     Can only do one at a time for external (to save on cpu), so if both
!     EXTERNALv, FL_POL, and FL_UNPOL are all true then do internal only (default).
      IF(EXTERNALv.AND.FL_POL.AND.FL_UNPOL) EXTERNALv = .FALSE.
      NPTS = NEXT

      PP_SIGA(NPTS) = 0.D0
      PP_SIGP(NPTS) = 0.D0
      PP_SIGELTAILA(NPTS) = 0.D0
      PP_SIGELTAILP(NPTS) = 0.D0
      PP_SIGQTAILA(NPTS) = 0.D0
      PP_SIGQTAILP(NPTS) = 0.D0
      PP_SIGINTAILA(NPTS) = 0.D0
      PP_SIGINTAILP(NPTS) = 0.D0

c     print*,'sigmar: EPRC=',EPRC

!     If polarized target nucleus is other than proton, calculate the
!     elastic radiative tail. For proton, this tail is calculated in
!     SIGRN after GETREADY(ITAIL=2) call is made. SIGRN calculates 
!     the (quasi)elastic tail from target of mass M, with NO 
!     smearing corrections at the moment for any target.
      IF(TARG.NE.'NH3') THEN
         ITAIL = 1
         CALL GETREADY(ITAIL)
         IF(TAIL_ON) CALL SIGRN(ITAIL,SIG,PP_SIGELTAILA(NPTS),
     >        PP_SIGELTAILP(NPTS))
      ENDIF
      ITAIL = 2
      CALL GETREADY(ITAIL)
      PP_X(NPTS) = XS
      PP_Y(NPTS) = YS
      OK = .TRUE.


CSK       IF(W2.LT.MPITHR) OK = .FALSE.
CSKold       IF((W2.LT.MPITHR).and.(TARG.NE.'ND3')) OK = .FALSE.
CSK


      if(TARG .eq. 'ND3') then 
         OK = ((W2 .gt. 0.0D0).and.(W2 .gt. (0.882D0-0.499D0*Q2))) ! 
CSK September 2012 - needs to be checked for Q^2 > 1.76
CSK First condition is strictly speaking too limiting, but for now I want to prevent mischief	   
         if(OK) MPITHR = 0.882D0-0.499D0*Q2
         if(MPITHR .le. 0.0D0) MPITHR = 0.001D0
      else
         IF(W2.LT.MPITHR) OK = .FALSE.
      endif


      IF(OK) CALL SIG0(W2,1.D0,1.D0,SIG,PP_SIGA(NPTS),PP_SIGP(NPTS))


      SIG0KINA = PP_SIGA(NPTS)
      SIG0KINP = PP_SIGP(NPTS)
      SIG0KIN = SIG0KINA + SIG0KINP

      CALL DELTAVR(PP_DVR(NPTS),DELINF)


      IF(TARG.EQ.'NH3'.AND.TAIL_ON) CALL SIGRN(ITAIL,SIG, 
     >     PP_SIGELTAILA(NPTS),PP_SIGELTAILP(NPTS))
     
CSK       IF(TARG.NE.'NH3'.AND.TAIL_ON) CALL SIGRN(ITAIL,SIG,
CSK     >      PP_SIGQTAILA(NPTS),PP_SIGQTAILP(NPTS))
CSK
      IF((TARG.NE.'NH3'.and.TARG.NE.'ND3').AND.TAIL_ON) 
     >     CALL SIGRN(ITAIL,SIG,PP_SIGQTAILA(NPTS),PP_SIGQTAILP(NPTS))
CSK

!     Skip special case
      SAVEv = FL_POL
      IF(POLTYPE.EQ.'TRAN'.AND.IPOL.EQ.5.AND.FL_POL) THEN
         FL_POL = .FALSE.
         IF(W2.GT.4.3D0) PP_SIGP(NPTS) = 0.D0
      ENDIF
      IF(INTPEAKING) THEN
         IF(OK) CALL SIGRF_FAST(SIG,PP_SIGINTAILA(NPTS),
     >        PP_SIGINTAILP(NPTS))
      ELSE
         IF(OK) CALL SIGRF(SIG,PP_SIGINTAILA(NPTS),
     >        PP_SIGINTAILP(NPTS))
      ENDIF
      FL_POL = SAVEv

!     Convert cross sections to nb/sr/GeV
      PP_SIGA(NPTS)   = PP_SIGA(NPTS)*HC2*EP*1000.D0
      PP_SIGP(NPTS)   = PP_SIGP(NPTS)*HC2*EP*1000.D0
      PP_SIGELTAILA(NPTS) = PP_SIGELTAILA(NPTS)*HC2*EP*1000.D0
      PP_SIGELTAILP(NPTS) = PP_SIGELTAILP(NPTS)*HC2*EP*1000.D0
      PP_SIGQTAILA(NPTS) = PP_SIGQTAILA(NPTS)*HC2*EP*1000.D0
      PP_SIGQTAILP(NPTS) = PP_SIGQTAILP(NPTS)*HC2*EP*1000.D0
      PP_SIGINTAILA(NPTS) = PP_SIGINTAILA(NPTS)*HC2*EP*1000.D0
      PP_SIGINTAILP(NPTS) = PP_SIGINTAILP(NPTS)*HC2*EP*1000.D0

c     print*,'sigmar.f L106: NPTS, pp_siga, HC2, EP', !kp: 6/3/12
c     >      NPTS, pp_siga(NPTS), HC2, EP

!     If MULTISOFT then exponentiate. 
!     See Kuchto and Shumeiko, NPB219, 412 (1983). Eq. 39.
      IF(MULTISOFT) then
         PP_SIGRADP(NPTS) = EXP(DELINF)*
     >        (PP_SIGP(NPTS)*(1.D0 + PP_DVR(NPTS) - DELINF) + 
     >        PP_SIGINTAILP(NPTS)) + PP_SIGELTAILP(NPTS) +
     >        PP_SIGQTAILP(NPTS)
         PP_SIGRADA(NPTS) = EXP(DELINF)*
     >        (PP_SIGA(NPTS)*(1.D0 + PP_DVR(NPTS) - DELINF) + 
     >        PP_SIGINTAILA(NPTS)) + PP_SIGELTAILA(NPTS) +
     >        PP_SIGQTAILA(NPTS)
      ELSE
         PP_SIGRADP(NPTS) = PP_SIGP(NPTS)*(1.D0 + PP_DVR(NPTS)) + 
     >        PP_SIGINTAILP(NPTS) + PP_SIGELTAILP(NPTS) +
     >        PP_SIGQTAILP(NPTS)
         PP_SIGRADA(NPTS) = PP_SIGA(NPTS)*(1.D0 + PP_DVR(NPTS)) + 
     >        PP_SIGINTAILA(NPTS) + PP_SIGELTAILA(NPTS) +
     >        PP_SIGQTAILA(NPTS)
      ENDIF

      IF(PP_SIGP(NPTS).NE.0.D0) PP_DELTP_I(NPTS) = 
     >     (PP_SIGRADP(NPTS)/PP_SIGP(NPTS) - 1.D0)*100.D0
      IF(PP_SIGA(NPTS).NE.0.D0) PP_DELTA_I(NPTS) = 
     >     (PP_SIGRADA(NPTS)/PP_SIGA(NPTS) - 1.D0)*100.D0

 40   IF(PP_SIGA(NPTS).NE.0.D0) PP_A0(NPTS)  = 
     >     PP_SIGP(NPTS)/PP_SIGA(NPTS)

c     PP_ARAD_I(NPTS) = PP_A0(NPTS)*(1.D0 + 
c     >    PP_DELTP_I(NPTS)/100.D0)/(1.D0 + PP_DELTA_I(NPTS)/100.D0)
      PP_ARAD_I(NPTS) = PP_SIGRADP(NPTS)/PP_SIGRADA(NPTS)

      IF(PP_ARAD_I(NPTS).NE.0.D0) PP_RAT1_I(NPTS) = 
     >     PP_A0(NPTS)/PP_ARAD_I(NPTS)
      IF(PP_SIGA(NPTS).NE.0.D0) PP_RAT2(NPTS) = 
     >     (PP_SIGELTAILA(NPTS) + PP_SIGQTAILA(NPTS))/
     >     PP_SIGA(NPTS)*100.D0
      IF(PP_SIGP(NPTS).NE.0.D0) PP_RAT3(NPTS) = 
     >     (PP_SIGELTAILP(NPTS) + PP_SIGQTAILP(NPTS))/
     >     PP_SIGP(NPTS)*100.D0
      IF(.NOT.EXTERNALv.AND.PP_SIGA(NPTS).NE.0.D0) THEN
         PP_ADIFF_I(NPTS) = 100.D0*(PP_SIGP(NPTS)*(PP_SIGINTAILA(NPTS) + 
     >        PP_SIGELTAILA(NPTS) + PP_SIGQTAILA(NPTS)) - 
     >        PP_SIGA(NPTS)*(PP_SIGINTAILP(NPTS) + 
     >        PP_SIGELTAILP(NPTS) + PP_SIGQTAILP(NPTS)))/
     >        (PP_SIGA(NPTS)*PP_SIGA(NPTS)*(1.D0 +PP_DVR(NPTS)) + 
     >        PP_SIGA(NPTS)*(PP_SIGINTAILA(NPTS) +
     >        PP_SIGELTAILA(NPTS) + PP_SIGQTAILA(NPTS))) 
      ENDIF

      IF(FL_POL) SIGMAR = PP_SIGRADP(NPTS)
      IF(FL_UNPOL) SIGMAR = PP_SIGRADA(NPTS)

      RETURN
      END
********************************************************************************

      SUBROUTINE DELTAVR(DVR,DELINF)
******************************************************************************
*     This subroutine calculates Eq. 33 of Kuchto, NP B219, 412 (1983).
*     
*     2/93, LMS.
*     6/96, DER. Added hadron vac. pol.
*****************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'

      REAL*8 DVR, DELINF, SPENCE, CC1, CC2, CC3, 
     >     CC4, CC5, DELVACH, HH1, HH2, HH3, MLEP2
      LOGICAL DONE/.FALSE./

      IF(.NOT.DONE) THEN
         DONE = .TRUE.
         CC1 = 13.D0/6.D0
         CC2 = 2.D0/3.D0
         CC3 = 48.D0/9.D0
         CC4 = 8.D0*MTAU2/3.D0
         CC5 = PI2/6.D0
      ENDIF

!     Following needed for correct implementation of Eq. 33 for both incident 
!     muons and electrons.
      IF(LEPTON.EQ.'MUON') THEN
         MLEP2 = ME2
      ELSE
         MLEP2 = MMU2
      ENDIF

!     h-vac constants come from Burkhardt, H.; TASSO NOTE No.192
      IF(Y .LT. 1.D0) THEN
         HH1 = -1.345D-9
         HH2 = -2.302D-3
         HH3 =  4.091
      ELSEIF(Y .LT. 64.D0) THEN
         HH1 = -1.512D-3
         HH2 = -2.822D-3
         HH3 =  1.218
      ELSE
         HH1 = -1.1344D-3
         HH2 = -3.0680D-3
         HH3 =  9.9992D-1
      ENDIF

      DELVACH = ((-2.D0)/AOP) * (HH1 + HH2*LOG(1.D0 + HH3*Y))

      DELINF = TR*LOG((NUMAX*NUMAX)/(SP*XP))
      DVR    = DELINF + AOP*(CC1*LITLM + CC2*LOG(Y/MLEP2) - CC3 
     >     + CC2*(Y + 2.D0*MTAU2)*LMTAU + 
     >     CC4/Y*(1.D0 - 2.D0*MTAU2*LMTAU) + DELVACH - 
     >     0.5D0*(LOG(SP/XP))**2
     >     + SPENCE((S*X - M2*Y)/(SP*XP)) - CC5)
      
      RETURN
      END
********************************************************************************
      FUNCTION DOT(V1,V2)
***************************************************************************
*     This function is the dot product between 4-vector V1 and V2.
*     
*     2/93, LMS.
**************************************************************************
      IMPLICIT NONE
      INTEGER I
      REAL*8 DOT, V1(4), V2(4), SIGN(4)/1.D0, 1.D0, 1.D0, -1.D0/

      DOT = 0.D0
      DO I = 1,4
         DOT = DOT + SIGN(I)*V1(I)*V2(I)
      ENDDO
      RETURN
      END
********************************************************************************
      SUBROUTINE GETREADY(ITAIL)
*******************************************************************************
*     This subroutine initializes variables, constants, etc..
*     
*     1/93, LMS
*******************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      INTEGER ITAIL
      REAL*8 DOT, LAMBDAMT, RTLAMBDAMT

      IUNPOL = 12
      IF(LEPTON.EQ.'MUON') THEN
         ML = MMU
      ELSE                      ! Assume incident electron
         ML = ME
      ENDIF
!     For deuterium and helium, the efFective pion threshold will be less than
!     that for free proton because of Fermi smearing effects.
      IF(TARG.NE.'NH3') THEN
         MPITHR = 1.1D0
      ELSE
         MPITHR = (MPI + MN)**2
      ENDIF

      IF(ITAIL.EQ.1) THEN
         IF(TARG.EQ.'ND3') M = MD
         IF(TARG.EQ.'LID') M = MD
         IF(TARG.EQ.'HE3') M = MHE3
      ELSEIF(ITAIL.EQ.2) THEN
         M      = MN
      ELSE
         RETURN
      ENDIF
      AOP    = ALPHA/PI         
      EU2    = EU*EU
      ED2    = ED*ED
      ML2    = ML*ML
      M2     = M*M
      MMU2   = MMU*MMU
      MTAU2  = MTAU*MTAU
      PI2    = PI*PI

      THETAR   = THETAD*RADCON

!     Define four-vectors.
c     -kag
      IF(POLTYPE.EQ.'LONG'.OR.POLTYPE.EQ.'long') THEN
         SCAT_PHIR = 90.D0*RADCON
         TS_BETA = 90.D0*RADCON
         TS_ALPHA = 0.D0
      ELSEIF(POLTYPE.EQ.'TRAN'.OR.POLTYPE.EQ.'tran') THEN
         SCAT_PHIR = 90.D0*RADCON
         TS_BETA = 90.D0*RADCON
         TS_ALPHA = 90.D0*RADCON
      ENDIF
c     -kag
      K1(1)    = 0.D0
      K1(2)    = 0.D0
      K1(3)    = E
      K1(4)    = E

c     -kag       K2(1)    = 0.D0
c     -kag       K2(2)    = EP*SIN(THETAR)
c     -kag       K2(3)    = EP*COS(THETAR)
c     -kag       K2(4)    = EP
      K2(1)    = EP*SIN(THETAR)*COS(SCAT_PHIR)
      K2(2)    = EP*SIN(THETAR)*SIN(SCAT_PHIR)
      K2(3)    = EP*COS(THETAR)
      K2(4)    = EP

      P1(1)    = 0.D0
      P1(2)    = 0.D0
      P1(3)    = 0.D0
      P1(4)    = M

      XI(1)    = 0.D0
      XI(2)    = 0.D0
      XI(3)    = PL*E/ML
      XI(4)    = PL*E/ML

      ETA(1)   = PN*SIN(TS_ALPHA)*COS(TS_BETA)
      ETA(2)   = PN*SIN(TS_ALPHA)*SIN(TS_BETA)
      ETA(3)   = PN*COS(TS_ALPHA)
      ETA(4)   = 0.D0
c     -kag       ETA(1)   = 0.D0
c     -kag       ETA(4)   = 0.D0
c     -kag       IF(POLTYPE.EQ.'LONG'.OR.POLTYPE.EQ.'long') THEN
c     -kag         ETA(2) = 0.D0
c     -kag         ETA(3) = PN
c     -kag       ELSEIF(POLTYPE.EQ.'TRAN'.OR.POLTYPE.EQ.'tran') THEN
c     -kag         ETA(2) = PN
c     -kag         ETA(3) = 0.D0
c     -kag       ENDIF

      Q(1)     = K1(1) - K2(1)
      Q(2)     = K1(2) - K2(2)
      Q(3)     = K1(3) - K2(3)
      Q(4)     = K1(4) - K2(4)

      QDOTXI   = DOT(Q,XI)
      QDOTETA  = DOT(Q,ETA)
      XIDOTETA = DOT(XI,ETA)
      P1DOTXI  = DOT(P1,XI)
      P1DOTQ   = DOT(P1,Q)
      K1DOTETA = DOT(K1,ETA)
      K2DOTETA = DOT(K2,ETA)
      K2DOTXI  = DOT(K2,XI)
      S        = -(2.D0*DOT(P1,K1))
      X        = -(2.D0*DOT(P1,K2))
      Y        = DOT(Q,Q)


      Q2        = 4.D0*E*EP*(SIN(THETAR/2.D0))**2
      W2        = M2 + S - X - Y
      XS        = Q2/(S - X)
c     print*,'GETREADY(): XS(BjrknX),Q2,S,X,E,M=',XS,Q2,S,X,E,M  !kp:4/20/12
      YS        = (E - EP)/E
      SX        = S - X
      SSP       = S + X
      LAMBDAS   = S*S - 4.D0*ML2*M2
      RTLAMBDAS = SQRT(LAMBDAS)
      LAMBDAX   = X*X + 4.D0*ML2*M2
      RTLAMBDAX = SQRT(LAMBDAX)
      LAMBDAM   = Y*Y + 4.D0*ML2*Y
      RTLAMBDAM = SQRT(LAMBDAM)
      NN        = 2.D0*ALPHA*ALPHA/RTLAMBDAS
      LM        = LOG((RTLAMBDAM + Y)/(RTLAMBDAM - Y))/RTLAMBDAM

      LAMBDAMT  = Y*Y + 4.D0*MTAU2*Y
      RTLAMBDAMT= SQRT(LAMBDAMT)
      LMTAU     = LOG((RTLAMBDAMT + Y)/(RTLAMBDAMT - Y))/RTLAMBDAMT
      LITLM     = LOG(Y/ML2)
      TR        = AOP*(LITLM - 1.D0)

      YM        = Y + 2.D0*ML2
      J0        = 2.D0*(YM*LM - 1.D0)

      SP        = X + Y
      XP        = S - Y
      NUMAX     = W2 - MPITHR
      LAMBDASP  = SP*SP - 4.D0*ML2*W2
      LAMBDAXP  = XP*XP - 4.D0*ML2*W2
      RTLAMBDASP= SQRT(LAMBDASP)
      RTLAMBDAXP= SQRT(LAMBDAXP)
      LSP       = LOG((SP + RTLAMBDASP)/(SP - RTLAMBDASP))/RTLAMBDASP
      LXP       = LOG((XP + RTLAMBDAXP)/(XP - RTLAMBDAXP))/RTLAMBDAXP

      AZ        = SX*SX + 4.D0*M2*Y
      DELTA     = 2.D0*(Y*(S*X - M2*Y) - ML2*AZ)
      C1        = X*SX - 2.D0*M2*Y
      C2        = 2.D0*ML2*SX - X*Y
      C1HAT     = (-S)*SX - 2.D0*M2*Y
      C2HAT     = 2.D0*ML2*SX +S*Y

      RETURN
      END
********************************************************************************
      SUBROUTINE RR(RP,RA,R,ITAIL)
*****************************************************************************
*     This subroutine computes the  the polarized and unpolarized contributions
*     to the function R given in Eq. 12 of Kuchto and Shumeiko, NPB219 (1983), 412.
*     
*     Assume that VARDEF has been called defined before this routine.
*     
*     2/93, LMS.
*     Added more exact formulas for R1 and R2 as given by A.6-A.9. Previous R1 
*     and R2 are the approximate formulas (A.11-A.14) and may only be valid 
*     for the longitudinal scattering case. If the logical OLDR1R2 is TRUE then
*     the approximate values are used, otherwise the more exact formulas 
*     are used. 2/94, LMS.
*****************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      INCLUDE 'radvar.inc'

      INTEGER II
      INTEGER ITAIL             ! 1 elastic, 2 quasielastic, 3 inelastic.

      LOGICAL OLDR1R2/.FALSE./
      CHARACTER*1 TGT
      REAL*8 RP, RA, R, XX, QQ, WW1, WW2, GG1,SIN2, E1, E2, 
     >     GG2, TWOMNU, S1, S2, R1, R2, R1TILDE, R2TILDE,
     >     SUM1, SUM2, KDOTETA, NU, TERM, R1_0, R2_0

      RA = 0.D0
      RP = 0.D0
!     Get structure functions.
      TWOMNU = MX2 - M2 + T
      XX = T/TWOMNU

      NU = TWOMNU/(2.D0*M)
      SIN2 = (SIN(THETAR/2.D0))**2
      TERM = T/(4.D0*SIN2)
      E1 = (NU + SQRT(NU*NU + 4.D0*TERM))/2.D0
      E2 = E1 - NU
      IF(TARG.EQ.'NH3') TGT = 'P'
      IF(TARG.EQ.'HE3') TGT = '3'
      IF(TARG.EQ.'ND3') TGT = 'D'
      IF(TARG.EQ.'LID') TGT = 'D'
      IF(TARG.EQ.'NEU') TGT = 'N'
      IF(ITAIL.EQ.1.OR.ITAIL.EQ.2) CALL ELASnew(IFFMOD,IPAULI,T,E1,
     >     TARG,ITAIL,M,WW1,WW2,GG1,GG2)
      IF(ITAIL.EQ.3) CALL INELnew(E1,E2,THETAD,TGT,WW1,WW2,GG1,GG2)

      IF(FL_POL) THEN
         IF(OLDR1R2) THEN
            SUM1 = 0.D0
            SUM2 = 0.D0
            DO II = 1, 2
               IF(II.EQ.1) KDOTETA = K1DOTETA
               IF(II.EQ.2) KDOTETA = K2DOTETA
               SUM1 = SUM1 + KDOTETA*((-T)*(2.D0*FIR + TD*ID) +
     >              2.D0*ML2*(Y*AAHAT(2,II) + T*AA(2,II)) -
     >              (T + Y)*(AA(1,II) + AAHAT(1,II)))
               SUM2 = SUM2 + KDOTETA*(X*AAHAT(1,II) - S*AA(1,II) +
     >              2.D0*ML2*S*(AA(2,II) + AAHAT(2,II)) -
     >              YM*SSP*AAD(II))
            ENDDO

!     Following are from Eq. A.11-A.14
            R1TILDE = M*PL*(SUM1 + QDOTETA*(Y*(I(1) + IHAT(1)) +
     >           2.D0*ML2*TD*IHAT(2)))
            R2TILDE = M*PL*T*((S*QDOTETA - TT*K1DOTETA)*(2.D0*FIR + 
     >           I(1)) - IHAT(1)*(X*QDOTETA - TT*K2DOTETA) -
     >           U*YM*QDOTETA*ID + SUM2)

            R1 = R1TILDE + 2.D0*ML2*M*PL*(2.D0*FIR*QDOTETA + 
     >           K1DOTETA*(TPOS*AD(1) + (TNEG - 2.D0*YM)*AAD(1)
     >           + T*(AHAT(2,1) - AAHAT(2,1)) +
     >           4.D0*ML2*(AA(2,1) +AAHAT(2,1)) -
     >           2.D0*ML2*(A(2,1) + AHAT(2,1))) +
     >           K2DOTETA*(TPOS*(AD(4) + AAD(1)) -
     >           2.D0*YM*AAD(2) + T*(AHAT(2,4) - AAHAT(2,1)) +
     >           2.D0*ML2*(AA(2,2) + AAHAT(2,2) -
     >           AA(2,1) - AAHAT(2,1) - A(2,4) - AHAT(2,4))))
            R2 = R2TILDE - 2.D0*ML2*M*PL*T*(K1DOTETA*(X*(AHAT(2,1) - 
     >           AAHAT(2,1)) - S*AD(1) + (S - TT)*AAD(1)) +
     >           K2DOTETA*(X*(AHAT(2,4) + AAHAT(2,1)) +
     >           TT*AAHAT(2,1) - S*(AAD(1) + AD(4))))
         ELSE
            SUM1 = 0.D0
            SUM2 = 0.D0
            DO II = 1, 2
               IF(II.EQ.1) KDOTETA = K1DOTETA
               IF(II.EQ.2) KDOTETA = K2DOTETA
               SUM1 = SUM1 + KDOTETA*(2.D0*ML2*(AA(2,II) + 
     >              AAHAT(2,II)) + YD*AAD(II))
               SUM2 = SUM2 + KDOTETA*(K2DOTXI*SX*AAD(II) +
     >              2.D0*P1DOTXI*(ML2*(AA(2,II) +
     >              AAHAT(2,II)) - YM*AAD(II)))
            ENDDO

!     Following are from Eq. A.6-A.9
            R1_0 = (-2.D0)*ML*M*(2.D0*FIR*(K2DOTXI*QDOTETA + 
     >           T*XIDOTETA)+K2DOTXI*(QDOTETA*TD*ID + SUM1))
            R2_0 = (-2.D0)*ML*M*T*(FIR*(2.D0*P1DOTXI*QDOTETA + 
     >           TT*XIDOTETA) - K2DOTXI*QDOTETA*ID*U + SUM2)
            R1 = R1_0 + 2.D0*ML*M*(
     >           K2DOTXI*(K1DOTETA*(TNEG*AAD(2)
     >           + TPOS*AD(4) + TM*(AHAT(2,4) - AAHAT(2,2))
     >           - 2.D0*ML2*(A(2,4) - AA(2,2))) + K2DOTETA*
     >           (TPOS*(AAD(2) + AD(2)) + T*(AHAT(2,2)
     >           - AAHAT(2,2)) - 2.D0*ML2*(AA(2,2) +AAHAT(2,2)
     >           + A(2,2) + AHAT(2,2)))) 
     >           + P1DOTXI*(K1DOTETA*(TNEG*AAD(3)
     >           + TPOS*AD(5) + TM*(AHAT(2,5) - AAHAT(2,3))
     >           - 2.D0*ML2*(A(2,5) - AA(2,3))) + K2DOTETA*
     >           (TPOS*(AAD(3) + AD(6)) + T*(AHAT(2,6)
     >           - AAHAT(2,3)) - 2.D0*ML2*(AA(2,3) +AAHAT(2,3)
     >           + A(2,6) + AHAT(2,6)))) 
     >           + XIDOTETA*(TPOS*A0D + TM*A0HAT(2) - 2.D0*ML2*A0(2)))
            R2 = R2_0 - 2.D0*ML*M*T*(
     >           K2DOTXI*(K1DOTETA*((S - TT)*AAD(2) -S*AD(4) +
     >           X*(AHAT(2,4) -AAHAT(2,2))) + K2DOTETA*
     >           ((X + TT)*AAHAT(2,2) + X*AHAT(2,2) - S*
     >           (AAD(2) + AD(2))))
     >           + P1DOTXI*(K1DOTETA*((S - TT)*AAD(3) -S*AD(5) +
     >           X*(AHAT(2,5) -AAHAT(2,3))) + K2DOTETA*
     >           ((X + TT)*AAHAT(2,3) + X*AHAT(2,6) - S*
     >           (AAD(3) + AD(6))))
     >           + XIDOTETA*(X*A0HAT(2) - S*A0D))

         ENDIF
         RP = M*GG1*R1 + GG2*R2/M
      ENDIF

      IF(FL_UNPOL) THEN
         S1 = I(0) + (TM*YM + 0.5D0*TD*TD)*ID - 
     >        ML2*TM*(I(2) + IHAT(2))
         S2 = ID*(T*(LAMBDAS + X*X - SX*TT) - 
     >        M2*(T*T + Y*Y) + 2.D0*ML2*(2.D0*S*X + 
     >        SX*TT - TT*TT)) - X*TT*I(1) - S*TT*IHAT(1) 
     >        + 2.D0*ML2*((M2*T - S*(S - TT))*I(2) +
     >        (M2*T - X*(X + TT))*IHAT(2)) - 2.D0*M2*I(0)

         RA = 2.D0*M*WW1*S1 + WW2*S2/(2.D0*M)
C     write(6,*) s1,s2,t,td,id

      ENDIF
      R = RP + RA

      RETURN
      END
********************************************************************************

      SUBROUTINE SIG0(MX2,XX1,XX2,SIG,SIGA,SIGP)
******************************************************************************
*     This subroutine returns the deep inelastic Born cross section, including 
*     spin-dependent terms for the kinematics described by variables in 
*     RADVAR.INC Formula is from Kuchto, NP B219, 412 (1983).
*     
*     1/93, LMS.
*     2/93, LMS. Added MX2, XX1, and XX2 arguments for evaluating the cross
*     section at various kinematics. See Eq. 27. For MX2 = W2, XX1
*     and XX2 are equal to one.
*****************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'

      INTEGER I, IQ, IW
      CHARACTER*1 TGT
      REAL*8 SIG, SIGA, SIGP, WW1, WW2, G1, G2, TEST
      REAL*8 MX2, XX1, XX2, STEMP, XTEMP, YTEMP, XSTEMP,
     >     NTEMP, K1P(4), K2P(4), QQ(4), QDXI, QDETA,W2_K, 
     >     P1DQ, DOT, QQQ, XIP(4), XIDETA, P1DXI, E1, E2
      EXTERNAL DOT

      DO I = 1,4
         K1P(I) = K1(I)*XX1
         K2P(I) = K2(I)*XX2
         QQ(I)  = K1P(I) - K2P(I)
         XIP(I) = K1P(I)*PL/ML
      ENDDO
      QDXI = DOT(QQ,XIP)
      QDETA= DOT(QQ,ETA)
      P1DQ = DOT(P1,QQ)
      XIDETA = DOT(XIP,ETA)
      P1DXI = DOT(P1,XIP)

      STEMP = S*XX1
      XTEMP = X*XX2
      YTEMP = M2 + STEMP - XTEMP - MX2
      XSTEMP = YTEMP/(STEMP - XTEMP)
      NTEMP = 2.D0*ALPHA*ALPHA/SQRT(STEMP*STEMP - 4.D0*ML2*M2)
      W2_K = M2 + STEMP - XTEMP - YTEMP

      IF(TARG.EQ.'NH3') TGT = 'P'
      IF(TARG.EQ.'HE3') TGT = '3'
      IF(TARG.EQ.'ND3') TGT = 'D'
      IF(TARG.EQ.'LID') TGT = 'D'
      IF(TARG.EQ.'NEU') TGT = 'N'
      CALL INELnew(E*XX1,EP*XX2,THETAD,TGT,WW1,WW2,G1,G2)

      IF(FL_UNPOL) SIGA = NTEMP/(YTEMP*YTEMP)*(2.D0*M*WW1*
     >     (YTEMP - 2.D0*ML2) + WW2*(STEMP*XTEMP - M2*YTEMP)/M)
      IF(FL_POL) SIGP = NTEMP/(YTEMP*YTEMP)*(4.D0*ML*
     >     (M2*G1*(QDXI*QDETA - YTEMP*XIDETA) - 
     >     YTEMP*G2*(P1DXI*QDETA - P1DQ*XIDETA)))

      SIG = SIGA + SIGP

      RETURN
      END
********************************************************************************
      SUBROUTINE SIGRF(SIG,SIGRFA,SIGRFP)
**************************************************************************
*     This subroutine calculates the polarized and unpolarized contributions
*     to the finite lepton bremsstrahlung cross section defined by Eq. 14 in 
*     Kuchto and Shumeiko, NP B219 (1983), 412.                      
*     
*     2/93, LMS.
**************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      INCLUDE 'radvar.inc'

      INTEGER NOK, NBAD, NLVL
      REAL*8 SIG, SIGRFA, SIGRFP, MX2MIN, MX2MAX,
     >     ALITTLE/1.D-8/
      EXTERNAL RKQC, SRFA, SRFP

      MX2MIN = MPITHR
      MX2MAX = W2 - ALITTLE

      SIGRFA = 0.D0
      SIGRFP = 0.D0
      IF(MX2MIN.LT.MX2MAX) THEN
c     IF(FL_UNPOL) CALL ODEINT(SIGRFA,1,MX2MIN,MX2MAX,EPS,H1,
c     >                            HMIN,NOK,NBAD,SRFA,RKQC)
c     IF(FL_POL) CALL ODEINT(SIGRFP,1,MX2MIN,MX2MAX,EPS,H1,
c     >                          HMIN,NOK,NBAD,SRFP,RKQC)
c     Replace with faster integration routine
         IF(FL_UNPOL) CALL QUADMO1(SIGRFA,SRFA,MX2MIN,MX2MAX,
     >        EPS,NLVL)
         IF(FL_POL) CALL QUADMO1(SIGRFP,SRFP,MX2MIN,MX2MAX,
     >        EPS,NLVL)
      ENDIF
      SIG = SIGRFA + SIGRFP

      RETURN
      END


      SUBROUTINE SRFA(MX2VAR,DUM,RESULT)
**************************************************************************
*     This subroutine calculates the unpolarized contribution to the inner
*     integral of the finite lepton bremsstrahlung cross section defined by 
*     Eq. 14 in Kuchto and Shumeiko, NP B219 (1983), 412.                      
*     
*     2/93, LMS.
**************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'

      INTEGER NOK, NBAD, NLVL
      REAL*8 MX2VAR, DUM, RESULT, TMIN, TMAX, MX2LAST/0.D0/, 
     >     LASTRESULT/0.D0/
      EXTERNAL RKQC_INT, RFA

      MX2 = MX2VAR
      IF(MX2.NE.0.D0.AND.MX2.EQ.MX2LAST) THEN
         RESULT = LASTRESULT
         RETURN
      ENDIF
      MX2LAST = MX2
      TMIN = ((W2 - MX2)*(SX - SQRT(AZ)) + 2.D0*MX2*Y)/
     >     (2.D0*W2)
      TMAX = ((W2 - MX2)*(SX + SQRT(AZ)) + 2.D0*MX2*Y)/
     >     (2.D0*W2)

      RESULT = 0.D0
C     IF(TMIN.LT.TMAX) CALL ODEINT_INT(RESULT,1,TMIN,TMAX,
C     >         EPS,H1,HMIN,NOK,NBAD,RFA,RKQC_INT)
c     Replace with faster integration routine
      IF(TMIN.LT.TMAX) CALL QUADMO2(RESULT,RFA,TMIN,TMAX,
     >     EPS,NLVL)
      LASTRESULT = RESULT
      RETURN
      END

      SUBROUTINE SRFP(MX2VAR,DUM,RESULT)
**************************************************************************
*     This subroutine calculates the polarized contribution to the inner
*     integral of the finite lepton bremsstrahlung cross section defined by 
*     Eq. 14 in Kuchto and Shumeiko, NP B219 (1983), 412.                      
*     
*     2/93, LMS.
**************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'

      INTEGER NOK, NBAD, NLVL
      REAL*8 MX2VAR, DUM, RESULT, TMIN, TMAX, MX2LAST/0.D0/, 
     >     LASTRESULT/0.D0/
      EXTERNAL RKQC_INT, RFP

      MX2 = MX2VAR
      IF(MX2.NE.0.D0.AND.MX2.EQ.MX2LAST) THEN
         RESULT = LASTRESULT
         RETURN
      ENDIF
      MX2LAST = MX2
      TMIN = ((W2 - MX2)*(SX - SQRT(AZ)) + 2.D0*MX2*Y)/
     >     (2.D0*W2)
      TMAX = ((W2 - MX2)*(SX + SQRT(AZ)) + 2.D0*MX2*Y)/
     >     (2.D0*W2)

      RESULT = 0.D0
C     IF(TMIN.LT.TMAX) CALL ODEINT_INT(RESULT,1,TMIN,TMAX,
C     >        EPS,H1,HMIN,NOK,NBAD,RFP,RKQC_INT)
c     Replace with faster integration routine
      IF(TMIN.LT.TMAX) CALL QUADMO2(RESULT,RFP,TMIN,TMAX,
     >     EPS,NLVL)


      LASTRESULT = RESULT
      RETURN
      END


      SUBROUTINE RFA(TVAR,DUM,RESULT)
***************************************************************************
*     This subroutine computes the unpolarized part of integrand in Eq. 14 of 
*     Kuchto and Shumeiko, NP B219 (1983), 412. Note that all the kinematic 
*     variables have been absorbed.
*     
*     2/93, LMS.
****************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'

      INTEGER ITAIL/3/          ! 1 elastic, 2 quasielastic, 3 inelastic.
      REAL*8 TVAR,DUM, RESULT, RP, RA, R, TLAST/0.D0/, 
     >     LASTRESULT/0.D0/, W2CUT

      T = TVAR
      IF(T.NE.0.D0.AND.T.EQ.TLAST) THEN
         RESULT = LASTRESULT
         RETURN
      ENDIF
      TLAST = T
      CALL VARDEF
      CALL RR(RP,RA,R,ITAIL)
      RESULT = RA*NN*AOP/(T*T) - AOP*FIR*SIG0KINA

      W2CUT = MN*MN + (1.D0/XCUT - 1.D0)*T
      IF(.NOT.TAIL_ON.AND.MX2.LT.W2CUT) RESULT =
     >     (-AOP)*FIR*SIG0KINA

      LASTRESULT = RESULT
      RETURN
      END

      SUBROUTINE RFP(TVAR,DUM,RESULT)
***************************************************************************
*     This subroutine computes the polarized part of integrand in Eq. 14 of 
*     Kuchto and Shumeiko, NP B219 (1983), 412. Note that all the kinematic 
*     variables have been absorbed.
*     
*     2/93, LMS.
****************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'

      INTEGER ITAIL/3/          ! 1 elastic, 2 quasielastic, 3 inelastic.
      REAL*8 TVAR,DUM, RESULT, RP, RA, R, TLAST/0.D0/, 
     >     LASTRESULT/0.D0/, W2CUT

      T = TVAR
      IF(T.NE.0.D0.AND.T.EQ.TLAST) THEN
         RESULT = LASTRESULT
         RETURN
      ENDIF
      TLAST = T
      CALL VARDEF
      CALL RR(RP,RA,R,ITAIL)
      RESULT = RP*NN*AOP/(T*T) - AOP*FIR*SIG0KINP

      W2CUT = MN*MN + (1.D0/XCUT - 1.D0)*T
      IF(.NOT.TAIL_ON.AND.MX2.LT.W2CUT) RESULT =
     >     (-AOP)*FIR*SIG0KINP

      LASTRESULT = RESULT
      RETURN
      END
********************************************************************************
      SUBROUTINE SIGRF_FAST(SIG,SIGRFA,SIGRFP)
***************************************************************************
*     This subroutine calculates the polarized and unpolarized contributions
*     to the finite lepton bremsstrahlung cross section defined by Eq. 27 in 
*     Kuchto and Shumeiko, NP B219 (1983), 412. This formula assumes the
*     peaking approximation.                     
*     
*     2/93, LMS.
**************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'
      include 'forDebug1.inc'   !8/29/12 

      INTEGER II, NOK, NBAD, NREGION, NLVL
      PARAMETER (NREGION=21)
      REAL*8 SIG, SIGRFA, SIGRFP, MX2MIN, MX2MAX, XSAVE(NREGION),
     >     ALITTLE/1.D-8/, XMIN, XMAX, DELX, SIGSUMP, XX1,
     >     SIGSUMA, SIGP(NREGION), SIGA(NREGION), W2LO, 
     >     W2HI, XLO, XHI, SIN2
      EXTERNAL RFPEAKA, RFPEAKP, RKQC

      MX2MIN = MPITHR
      MX2MAX = W2 - ALITTLE

      SIGRFA = 0.D0
      SIGRFP = 0.D0
      IF(MX2MIN.LT.MX2MAX) THEN
c     IF(FL_UNPOL) CALL ODEINT(SIGRFA,1,MX2MIN,MX2MAX,EPS,H1,
c     >                            HMIN,NOK,NBAD,RFPEAKA,RKQC)
c     IF(FL_POL) CALL ODEINT(SIGRFP,1,MX2MIN,MX2MAX,EPS,H1,
c     >                          HMIN,NOK,NBAD,RFPEAKP,RKQC)
c     Replace with faster integration routine
         IF(FL_UNPOL) CALL QUADMO1(SIGRFA,RFPEAKA,MX2MIN,MX2MAX,
     >        EPS,NLVL)
         IF(FL_POL) CALL QUADMO1(SIGRFP,RFPEAKP,MX2MIN,MX2MAX,
     >        EPS,NLVL)
      ENDIF
      SIG = SIGRFA + SIGRFP       


      RETURN
      END


      SUBROUTINE RFPEAKA(MX2VAR,DUM,RESULT)
**************************************************************************
*     This subroutine calculates the unpolarized contribution to the inner
*     integral of the finite lepton bremsstrahlung cross section defined by 
*     Eq. 27 in Kuchto and Shumeiko, NP B219 (1983), 412.                      
*     
*     2/93, LMS.
**************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'

      REAL*8 MX2VAR, DUM, RESULT, MX2, NU, NU1, NU2, XX1, 
     >     XX2, T1, T2, SIG1, SIG2, SIGA1, SIGA2, 
     >     SIGP1, SIGP2, MX2LAST/0.D0/, LASTRESULT/0.D0/,
     >     W2CUT1, W2CUT2, QS1, QS2

      MX2 = MX2VAR
      IF(MX2.NE.0.D0.AND.MX2.EQ.MX2LAST)THEN
         RESULT = LASTRESULT
         RETURN                 !Keep previous result
      ENDIF
      MX2LAST = MX2
      NU  = W2 - MX2
      NU1 = NU/(S - Y)
      NU2 = NU/(X + Y)
      XX1 = 1.D0 - NU1
      XX2 = 1.D0 + NU2
      T1  = AOP*(0.5D0*(1.D0 + XX1*XX1)*LITLM - XX1)
      T2  = AOP*(0.5D0*(1.D0 + 1.D0/(XX2*XX2))*LITLM - 1.D0/XX2)
      CALL SIG0(MX2,XX1,1.D0,SIG1,SIGA1,SIGP1)
      CALL SIG0(MX2,1.D0,XX2,SIG2,SIGA2,SIGP2)
      RESULT = (T1*SIGA1 + XX2*T2*SIGA2 -
     >     2.D0*TR*SIG0KINA)/NU

      IF(.NOT.TAIL_ON) THEN
         QS1 = Q2*E/(EP + Q2/(2.D0*MN*XCUT))
         QS2 = Q2*EP/(E - Q2/(2.D0*MN*XCUT))
         W2CUT1 = MN*MN + (1.D0/XCUT - 1.D0)*QS1
         W2CUT2 = MN*MN + (1.D0/XCUT - 1.D0)*QS2
         RESULT = (-2.D0)*TR*SIG0KINA/NU
         IF(MX2.GT.W2CUT1) RESULT=RESULT + XX2*T2*SIGA2/NU
         IF(MX2.GT.W2CUT2) RESULT=RESULT + T1*SIGA1/NU
      ENDIF

      LASTRESULT = RESULT

      RETURN
      END


      SUBROUTINE RFPEAKP(MX2VAR,DUM,RESULT)
**************************************************************************
*     This subroutine calculates the polarized contribution to the inner
*     integral of the finite lepton bremsstrahlung cross section defined by 
*     Eq. 27 in Kuchto and Shumeiko, NP B219 (1983), 412.                      
*     
*     2/93, LMS.
**************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'

      REAL*8 MX2VAR, DUM, RESULT, MX2, NU, NU1, NU2, XX1, 
     >     XX2, T2, T1P, SIG1, SIG2, SIGA1, SIGA2, 
     >     SIGP1, SIGP2, MX2LAST/0.D0/, LASTRESULT/0.D0/,
     >     W2CUT1, W2CUT2, QS1, QS2

      MX2 = MX2VAR
      IF(MX2.NE.0.D0.AND.MX2.EQ.MX2LAST)THEN
         RESULT = LASTRESULT
         RETURN                 !Keep previous result
      ENDIF
      MX2LAST = MX2
      NU  = W2 - MX2
      NU1 = NU/(S - Y)
      NU2 = NU/(X + Y)
      XX1 = 1.D0 - NU1
      XX2 = 1.D0 + NU2
      T2  = AOP*(0.5D0*(1.D0 + 1.D0/(XX2*XX2))*LITLM - 1.D0/XX2)
      T1P = AOP*(0.5D0*(1.D0 + XX1)*LITLM - 1.D0 - NU1*XX1)
      CALL SIG0(MX2,XX1,1.D0,SIG1,SIGA1,SIGP1)
      CALL SIG0(MX2,1.D0,XX2,SIG2,SIGA2,SIGP2)
      RESULT = (T1P*SIGP1 + XX2*T2*SIGP2 -
     >     2.D0*TR*SIG0KINP)/NU

      IF(.NOT.TAIL_ON) THEN
         QS1 = Q2*E/(EP + Q2/(2.D0*MN*XCUT))
         QS2 = Q2*EP/(E - Q2/(2.D0*MN*XCUT))
         W2CUT1 = MN*MN + (1.D0/XCUT - 1.D0)*QS1
         W2CUT2 = MN*MN + (1.D0/XCUT - 1.D0)*QS2

         RESULT = (-2.D0)*TR*SIG0KINP/NU
         IF(MX2.GT.W2CUT1) RESULT=RESULT + XX2*T2*SIGP2
         IF(MX2.GT.W2CUT2) RESULT=RESULT + T1P*SIGP1
      ENDIF

      LASTRESULT = RESULT

      RETURN
      END
********************************************************************************
      SUBROUTINE SIGRN(ITAIL,SIG,SIGRNA,SIGRNP)
**************************************************************************
*     This subroutine calculates the polarized and unpolarized contributions
*     to the elastic tail cross section defined by Eq. 19 in Kuchto and
*     Shumeiko, NP B219 (1983), 412. The multisoft photon correction is given
*     by Eq. 41.                     
*     
*     2/93, LMS.
**************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'

      INTEGER NOK, NBAD, ITAIL, ITL, NLVL
      REAL*8 SIG, SIGRNA, SIGRNP, TMIN, TMAX, SOFTCORR
      COMMON /TAIL/ITL
      EXTERNAL RKQC, RNA, RNP

      ITL = ITAIL
      MX2 = M2
      TMIN = ((W2 - MX2)*(SX - SQRT(AZ)) + 2.D0*MX2*Y)/
     >     (2.D0*W2)
      TMAX = ((W2 - MX2)*(SX + SQRT(AZ)) + 2.D0*MX2*Y)/
     >     (2.D0*W2)

      SIGRNA = 0.D0
      SIGRNP = 0.D0
      IF(TMIN.LT.TMAX) THEN
c     IF(FL_UNPOL) CALL ODEINT(SIGRNA,1,TMIN,TMAX,EPS,H1,
c     >                            HMIN,NOK,NBAD,RNA,RKQC)
c     IF(FL_POL) CALL ODEINT(SIGRNP,1,TMIN,TMAX,EPS,H1,
c     >                          HMIN,NOK,NBAD,RNP,RKQC)
c     Replace with faster integration routine
         IF(FL_UNPOL) CALL QUADMO1(SIGRNA,RNA,TMIN,TMAX,EPS,NLVL)
         IF(FL_POL) CALL QUADMO1(SIGRNP,RNP,TMIN,TMAX,EPS,NLVL)
      ENDIF

      IF(MULTISOFT) THEN
         SOFTCORR = ((W2 - MX2)**2/(S*XP))**TR
         SIGRNA = SIGRNA*SOFTCORR
         SIGRNP = SIGRNP*SOFTCORR
      ENDIF

      SIG = SIGRNA + SIGRNP
      RETURN
      END


      SUBROUTINE RNA(TVAR,DUM,RESULT)
***************************************************************************
*     This subroutine computes the unpolarized part of R for elastic
*     scattering, where R = R(t) in Eq. 19 in Kuchto and Shumeiko, NP B219 
*     (1983), 412. Note that all the kinematic variables have been absorbed.
*     
*     2/93, LMS.
****************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'

      INTEGER ITAIL
      REAL*8 TVAR,DUM, RESULT, RP, RA, R, TLAST/0.D0/, 
     >     LASTRESULT/0.D0/
      COMMON /TAIL/ITAIL

      T = TVAR
      IF(T.NE.0.D0.AND.T.EQ.TLAST) THEN
         RESULT = LASTRESULT
         RETURN
      ENDIF
      TLAST = T

      CALL VARDEF
      CALL RR(RP,RA,R,ITAIL)
      RESULT = RA*NN*AOP/(T*T)

      LASTRESULT = RESULT
      RETURN
      END

      SUBROUTINE RNP(TVAR,DUM,RESULT)
***************************************************************************
*     This subroutine computes the polarized part of R for elastic
*     scattering, where R = R(t) in Eq. 19 in Kuchto and Shumeiko, NP B219 
*     (1983), 412. Note that all the kinematic variables have been absorbed.
*     
*     2/93, LMS.
****************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'

      INTEGER ITAIL
      REAL*8 TVAR,DUM, RESULT, RP, RA, R, TLAST/0.D0/, 
     >     LASTRESULT/0.D0/
      COMMON /TAIL/ITAIL

      T = TVAR
      IF(T.NE.0.D0.AND.T.EQ.TLAST) THEN
         RESULT = LASTRESULT
         RETURN
      ENDIF
      TLAST = T

      CALL VARDEF
      CALL RR(RP,RA,R,ITAIL)
      RESULT = RP*NN*AOP/(T*T)

      LASTRESULT = RESULT
      RETURN
      END
********************************************************************************
      FUNCTION SPENCE(X)                                                        
*******************************************************************************
*     Spence function. This function is evaluated using the "recipe" given by 
*     TSAI, SLAC-PUB-848 (1971), pp. 14-15.
*     
*     1/93, LMS.
******************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      
      INTEGER I
      REAL*8 X, SPENCE, Y, SUM, TERM           
      
      IF(X.EQ.1.D0) THEN
         SPENCE = PI2/6.D0
      ELSEIF(X.EQ.-1.D0) THEN
         SPENCE = (-PI2)/12.D0
      ELSE
         IF(ABS(X).LT.1.D0) Y = X
         IF(ABS(X).GT.1.D0) Y = 1.D0/X
         TERM = Y
         SUM = Y
         I = 2
         DO WHILE(ABS(TERM).GT.ABS(Y*1.0D-9))
            TERM = TERM*Y*FLOAT((I-1)*(I-1))/FLOAT(I*I)
            SUM = SUM + TERM 
            I = I + 1
         ENDDO
         IF(ABS(X).LT.1.D0) THEN
            SPENCE = SUM
         ELSEIF(X.GT.1.D0) THEN
            SPENCE = (-0.5D0)*(LOG(ABS(X))**2) + PI2/3.D0 - SUM
         ELSEIF(X.LT.-1.D0) THEN
            SPENCE = (-0.5D0)*(LOG(ABS(X))**2) - PI2/6.D0 - SUM
         ENDIF
      ENDIF
      RETURN
      END
********************************************************************************
      SUBROUTINE VARDEF
*******************************************************************************
*     This subroutine defines variables needed for integration purposes.
*     The final goal here is to evaluate integral given by Eq. 14 in
*     Kuchto and Shumeiko, NPB219 (1983), 412. Most of these equations are 
*     found on page 433.
*     
*     2/93, LMS
*******************************************************************************
      IMPLICIT NONE
      INCLUDE 'instruct.inc'
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'radvar.inc'

      INTEGER N, II
      REAL*8 TERM1, TERM2, TERM3, TA0(2), TA0HAT(2)

!     Assume T and MX2, the integration variables have been previously defined.
      TT       = T + MX2 - M2
      TD       = T - Y
      U        = SX - TT
      BZ       = S*Y*U + X*(SX*T - Y*TT) - 2.D0*M2*Y*TD
      CZ       = (X*T - Y*(S - TT))**2 + 
     >     4.D0*ML2*(U*(SX*T - Y*TT) - M2*TD*TD)
      BZHAT    = (-X)*Y*U - S*(SX*T - Y*TT) - 2.D0*M2*Y*TD
      CZHAT    = ((-S)*T + Y*(X + TT))**2 + 
     >     4.D0*ML2*(U*(SX*T - Y*TT) - M2*TD*TD)
      
!     Changed I(-1) = 0.D0, and IHAT(-1) = 0.D0 definitions as per instructions
!     from Nikolai Shumeiko, 5/27/94. LMS.
      I(0)     = 1.D0/SQRT(AZ)
      I(1)     = 1.D0/SQRT(CZ)
      I(2)     = I(1)*BZ/CZ 
      I(-1)    = I(0)*BZ/AZ
      IHAT(-1) = (-I(0))*BZHAT/AZ
      IHAT(0)  = I(0)
      IHAT(1)  = 1.D0/SQRT(CZHAT)
      IHAT(2)  = (-IHAT(1))*BZHAT/CZHAT

      TERM1    = U*C2 + TD*LAMBDAX
      TERM2    = U*C2HAT + TD*(LAMBDAS + C1HAT)
      TERM3    = U*LAMBDAM + TD*C2

      DO N = 0,2
         AA(N,1) = (C1*I(N-1) - TERM1*I(N))/DELTA
         AA(N,2) = (C1HAT*I(N-1) + TERM2*I(N))/DELTA
         AA(N,3) = (Y*SSP*I(N-1) - TERM3*I(N))/DELTA

         AAHAT(N,1) = (C1*(IHAT(N-1) - TD*IHAT(N)) - 
     >        TERM1*IHAT(N))/DELTA
         AAHAT(N,2) = (C1HAT*(IHAT(N-1) - TD*IHAT(N)) + 
     >        TERM2*IHAT(N))/DELTA
         AAHAT(N,3) = (Y*SSP*(IHAT(N-1) - TD*IHAT(N)) - 
     >        TERM3*IHAT(N))/DELTA
      ENDDO

      IF(FL_POL) THEN
         DO N = 1,2
            TA0(N)    = AA(N-1,1) + AA(N-1,2) + TD*AA(N,1) + U*AA(N,3)
            TA0HAT(N) = AAHAT(N-1,1) + AAHAT(N-1,2) - TD*AAHAT(N,2) + 
     >           U*AAHAT(N,3)
            A0(N) = TA0(N)/2.D0
            A0HAT(N) = TA0HAT(N)/2.D0

            A(N,1)    = (C1*AA(N-1,1) - TERM1*AA(N,1) - 
     >           LAMBDAX*TA0(N))/DELTA
            A(N,2)    = (C1HAT*AA(N-1,2) + TERM2*AA(N,2) - 
     >           LAMBDAS*TA0(N))/DELTA
            A(N,4)    = (C1HAT*AA(N-1,1) + TERM2*AA(N,1) +
     >           (LAMBDAS + C1HAT)*TA0(N))/DELTA
            A(N,5)    = (Y*SSP*AA(N-1,1) - TERM3*AA(N,1) - 
     >           C2*TA0(N))/DELTA
            A(N,6)    = (Y*SSP*AA(N-1,2) - TERM3*AA(N,2) +
     >           C2HAT*TA0(N))/DELTA

            AHAT(N,1) = (C1*(AAHAT(N-1,1) - TD*AAHAT(N,1)) - 
     >           TERM1*AAHAT(N,1) - LAMBDAX*TA0HAT(N))/DELTA
            AHAT(N,2) = (C1HAT*(AAHAT(N-1,2) - TD*AAHAT(N,2)) + 
     >           TERM2*AAHAT(N,2) - LAMBDAS*TA0HAT(N))/DELTA
            AHAT(N,4) = (C1HAT*(AAHAT(N-1,1) - TD*AAHAT(N,1)) + 
     >           TERM2*AAHAT(N,1) + 
     >           (LAMBDAS + C1HAT)*TA0HAT(N))/DELTA
            AHAT(N,5) = (Y*SSP*(AAHAT(N-1,1) - TD*AAHAT(N,1)) - 
     >           TERM3*AAHAT(N,1) - C2*TA0HAT(N))/DELTA
            AHAT(N,6) = (Y*SSP*(AAHAT(N-1,2) - TD*AAHAT(N,2)) - 
     >           TERM3*AAHAT(N,2) + C2HAT*TA0HAT(N))/DELTA

         ENDDO       
      ELSE
         DO N = 1,2
            DO II = 1,6
               A(N,II) = 0.D0
               AHAT(N,II) = 0.D0
            ENDDO
         ENDDO
      ENDIF

      TPOS   = T + 2.D0*YM
      TNEG   = T - 2.D0*YM

      ID     = (I(1) - IHAT(1))/TD
      FIR    = YM*ID - ML2*(I(2) + IHAT(2))
      TM     = T - 2.D0*ML2
      YD     = TM - YM
      A0D = (A0(1) - A0HAT(1))/TD
      DO II = 1,6
         AD(II) = (A(1,II) - AHAT(1,II))/TD 
         IF(II.LE.3) AAD(II) = (AA(1,II) - AAHAT(1,II))/TD 
      ENDDO

      RETURN
      END
********************************************************************************

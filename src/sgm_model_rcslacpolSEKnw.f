      subroutine sgm_model_rcslacpol(NPTS,TMOD,EE,EEP, !E = Ebeam, EP = E' = E_prime and THETAD = theta in degrees
     >     THETADD,XXCUT,SCAT_PHIRR,sig)


      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      INCLUDE 'prinplot.inc'
      
      include 'binning.inc'     !kp: from cs_init.f
c     ================================================= 2/27/13 =========
c     
c     I suspected the value of 7.0 given to Z before calling INTEG4(TMOD,NPTS,Q2OUT, EB_INDEX)
c     I thought may be I should give it a value of 1.0 for deuteron rather than 7.0 for N in NH3.
c     
c     
c     ===================================================================
      real Ebeam, sig
      real sgm_model,degrad,ppi
      data degrad/0.0174532925e0/,ppi/3.14159265e0/ ! ppi used: "Symbol 'pi' at (1) already has basic type of REAL"

      INTEGER EB_INDEX          !Jixie: EB_INDEX will be signed according to input beam erngy and target type
c     REAL*8 CSRAD              !kp: radiated cross-section to be returned by this function.
      REAL*8 EE,EEP,THETADD,XXCUT,SCAT_PHIRR !kp: added because of "Error: COMMON attribute conflicts with DUMMY attribute in 'e' at (1)", see radvarfix.inc
      INTEGER DBG_MY_KINE
      LOGICAL EXTERNALv
      INTEGER IP, NPTS, TMOD,NEXT
      REAL*8 ENERGY, Z, EINT, SANSWER, TA, Q2OUT, ERC, 
     >     EPRC, THRC, ANSWER, SIGMAR, SIGMARNEW, MULT, 
     >     XCUTRC, W2BIN, TB
c     CHARACTER*8 TARGETS(5)/'PROTON  ','3HELIUM ',
c     >      'ND3     ','LiD     ','NEUTRON '/
      real*8 Ess,thh,Epp
      real*8 nu,qq2,xx1,maa,mpp
      parameter (maa=12*0.93827) !kp: carbon mass in GeV
      parameter (mpp=0.93827)
      
      real*8 mtarg              ! csk
      
      COMMON /KINLIST/ERC,EPRC,THRC,EINT,Z,XCUTRC
c     COMMON /INFO/EXTERNALv,NEXT


      
      THETAD = THETADD          !kp: 4/3/12


      sig = 0.0
      Ebeam=EE

      Ess=dble(Ebeam)
      thh=dble(THETADD*degrad)
      Epp=dble(EEP)             !kp: electron momentum
c     Kinematic tests
c     sgm_model=0.e0
      nu=Ebeam-Epp              !kp: energy transfer
c     kp3/2/13      if(nu.lt.0.e0) return     !The originial sgm_model.f file had this line disabled.
      if(nu.lt.0.e0) return

      qq2=2.e0*Ebeam*Epp*(1.e0-cos(thh)) !kp: Q^2 = - 4-mom-transfer-squared = 4*Eb*Ep*sin^2(th/2)
c     x=q2/(2.e0*maa*nu)                     !kp: Lorentz invariant Bjorken scaling variable = Q^2/2Pq = Q^2/2M*nu
c     xx1=qq2/(2.e0*mpp*nu)
c     kp3/2/13      xx1=qq2/(2.e0*mn*nu)      !kp: average nucleon mass = 0.938272 (defined in radcon.inc)
c     if(xx1.gt.1.e0) return     !The originial sgm_model.f file had this line disabled.
c     kp3/2/13      if(xx1.gt.1.9e0) return     !kp: 4/19/12 for deut/ND3 target!The originial sgm_model.f file had this line disabled.
      xx1=qq2/(2.e0*mn*nu)
      if(xx1.gt.1.9e0) return



      
      
c     print*,'start calculation'
      SCAT_PHIR = SCAT_PHIRR
      
      ENERGY = Ebeam
      ERC = Ebeam
      EPRC = EEP
      XCUTRC = XXCUT
      THRC = THETADD

      
c     IP = 1                   !Don't see any use of it


c     print*, 'NPTS: ThDeg',NPTS,THRC !kp: 

      IF(ZTEST) Z = 1.D0
      IF(EXPER.EQ.'E80 '.OR.EXPER.EQ.'E130') THEN
         IF(.NOT.ZTEST) Z = 7.D0 !Z of Nitrogen in NH3 or ND3?
         TB = 0.04294D0
         TA = 0.04414D0
         CALL INTSIMPLE(TB,TA,NPTS,Q2OUT)
c     Other If blocks for other experiments removed     !kp: 3/27/12

c     ELSEIF(EXPER.EQ.'EG1B') THEN    !kp: 9/15/12
c     IF(.NOT.ZTEST) Z = 7.D0    !kp: 9/15/12
c     CALL INTEG1B(TMOD,NPTS,Q2OUT)  !kp: Disabled on 9/15/12 to avoid compiling unused files rleg1b.f & integ1b.f 

c     ELSEIF(EXPER.EQ.'EG4') THEN
      ELSEIF(EXPER.EQ.'EG4A') THEN !A added for convenience of not having to modify a lot of things
         IF(.NOT.ZTEST) Z = 7.D0
c     IF(.NOT.ZTEST) Z = 1.D0 !kp: 2/27/13
         CALL INTEG4(TMOD,NPTS,Q2OUT, EB_INDEX)
      ENDIF 





c     print*,'finished integration.'

!     kp: Above integ routines return values of PP_SIGRADA(NPTS), PP_SIGRADANOTAIL(NPTS), PP_SIGRADP(NPTS) (and possibly other things too, haven't checked that yet) which are defined in prinplot.inc as:       PP_SIGRADA(PP_NMAX),     ! Unpolarized radiated cross section.

      PP_ARAD_E(NPTS) = PP_SIGRADP(NPTS)/PP_SIGRADA(NPTS)
      PP_RAT1_E(NPTS) = PP_A0(NPTS)/PP_ARAD_E(NPTS)
      IF(PP_SIGP(NPTS).NE.0.D0) THEN
         PP_DELTP_E(NPTS) = (PP_SIGRADP(NPTS)/PP_SIGP(NPTS)-1.0)
     >        *100.0
      ENDIF
      PP_DELTA_E(NPTS) = (PP_SIGRADA(NPTS)/PP_SIGA(NPTS)-1.0)*100.0
      PP_ADIFF_E(NPTS) = (PP_A0(NPTS) - PP_ARAD_E(NPTS))*100.D0

      PP_FRC_E(NPTS) = PP_SIGRADANOTAIL(NPTS)/PP_SIGRADA(NPTS) 
      PP_ARC_E(NPTS) = (PP_A0(NPTS) - 
     >     PP_ARAD_E(NPTS)/PP_FRC_E(NPTS))*100.D0

      IF(POLTYPE.EQ.'TRAN'.AND.IPOL.EQ.5) PP_DELTP_I(NPTS) = 0.D0

      
c     call PRINT_SOME_OUTPUTS(NPTS,Q2OUT)    





c     Obtain 2-differential cross section
c     sig=sgm_model(Ebeam,theta_el,p_el) !c      write(31,33) i_p_el,i_th_h,i_ph_h,i_p_h,sig
      sig=PP_SIGRADA(NPTS)      !8/3/12 For unpolarized radiated cross-section
c     sig=PP_SIGA(NPTS)      !8/3/12 For unpolarized Born cross-section
      print*,'x, sig = ',xx1,sig



      return
      END
********************************************kp: Program Ends here **********






      SUBROUTINE INTSIMPLE(TB,TA,NPTS,Q2OUT)
*******************************************************************************
*     This routine is for evaluating external+internal RC's for a single TB and TA
*     
*     7/97, LMS.
*******************************************************************************
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      INCLUDE 'prinplot.inc'
      
      LOGICAL EXTFLAG
      INTEGER NPTS, IS, NEXT
      REAL*8 TB, TA, SIGRU, SIGRUNT, SIGRP, Q2OUT,SANSWER, 
     >     SIGMAR, MULT

      COMMON /INFO/EXTFLAG,NEXT

      NEXT = NPTS
      FL_POL = .TRUE.
      FL_UNPOL = .TRUE.
      EXTFLAG = .FALSE.         !turn off momentarily.
      SANSWER = SIGMAR(E,EP,THETAD)
      Q2OUT = Q2

      EXTFLAG = .TRUE.
      NEXT = 500

      CALL RADCROSS(TB,TA,SIGRU,SIGRUNT,SIGRP)

      PP_SIGRADA(NPTS) = SIGRU
      PP_SIGRADANOTAIL(NPTS) = SIGRUNT
      PP_SIGRADP(NPTS) = SIGRP

      RETURN
      END
*******************************************************************************


      SUBROUTINE RADCROSS(TB,TA,SIGRU,SIGRUNT,SIGRP)
*******************************************************************************
*     This routine is for evaluating fully radiated cross sections for some 
*     TB and TA
*     
*     7/97, LMS.
*******************************************************************************
      IMPLICIT NONE
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      
      REAL*8 TB, TA, TB1, TA1, SIGRU, SIGRP, SIGRUNT, MULT
      
      real*8 MP/0.938272D0/, mtarg ! csk

      MULT = 1.D0
      IF(TARGNEG.OR.TARGPOS) THEN
         IF(TARGNEG) MULT = 0.95D0
         IF(TARGPOS) MULT = 1.05D0
      ENDIF
c     kag	temporary fix for no external
c     MULT = 0
      TB1 = TB*MULT
      TA1 = TA*MULT
c     sk
      mtarg = MP
      if (TARG .eq. 'ND3') mtarg = 1.876D0
c     sk	   
      
c     print*,'RADCROSS(): TB= ',TB

!     Tails are turned off for calculation of dilution rc which is the ratio of
!     rates from DIS to DIS+tails. Tail is defined by events with x > XCUT which
!     is read from the kinematic input file. 
      IF(TB.GT.0.D0) THEN       ! Full external calculation
         FL_POL = .FALSE.
         FL_UNPOL = .TRUE.
c     kp:         IF(IEXTERNAL.EQ.1) CALL EXTERN2(TB1,TA1,SIGRU)  !Commented out again: 9/8/12
         IF(IEXTERNAL.EQ.2) CALL EXTERN2NEW(TB1,TA1,SIGRU, mtarg)
         IF(XCUT.LT.1.0) THEN
            TAIL_ON = .FALSE.
c     kp:           IF(IEXTERNAL.EQ.1) CALL EXTERN2(TB1,TA1,SIGRUNT)
            IF(IEXTERNAL.EQ.2) CALL EXTERN2NEW(TB1,TA1,SIGRUNT,mtarg)
            TAIL_ON = .TRUE.
         ENDIF
         FL_POL = .TRUE.
         FL_UNPOL = .FALSE.
c     kp:         IF(IEXTERNAL.EQ.1) CALL EXTERN2(TB1,TA1,SIGRP)
         IF(IEXTERNAL.EQ.2) CALL EXTERN2NEW(TB1,TA1,SIGRP,mtarg)
      ELSE                      ! External after scattering only
         FL_POL = .FALSE.
         FL_UNPOL = .TRUE.
c     kp:         IF(IEXTERNAL.EQ.1) CALL EXTERN1(TA1,SIGRU)   !Commented out again: 9/8/12
c     kp:         IF(IEXTERNAL.EQ.2) CALL EXTERN1NEW(TA1,SIGRU)   !kp: I think, SEK forgot to add 'mtarg'? 9/14/12
         IF(IEXTERNAL.EQ.2) CALL EXTERN1NEW(TA1,SIGRU,mtarg) !kp: 9/14/12
         IF(XCUT.LT.1.0) THEN
            TAIL_ON = .FALSE.
c     kp:           IF(IEXTERNAL.EQ.1) CALL EXTERN1(TA1,SIGRUNT)
            IF(IEXTERNAL.EQ.2) CALL EXTERN1NEW(TA1,SIGRUNT,mtarg)
            TAIL_ON = .TRUE.
         ENDIF
         FL_POL = .TRUE.
         FL_UNPOL = .FALSE.
c     kp:         IF(IEXTERNAL.EQ.1) CALL EXTERN1(TA1,SIGRP)
         IF(IEXTERNAL.EQ.2) CALL EXTERN1NEW(TA1,SIGRP,mtarg)
      ENDIF

      RETURN
      END
*******************************************************************************


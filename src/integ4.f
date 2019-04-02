      SUBROUTINE INTEG4(IMODE,NPTS,Q2OUT,EB_INDEX)
*****************************************************************************
*     This program was made of a copy of integ1b.f       Jan, 12, 2012, K. Adhikari
*     I simply replaced eg1b lines with the ones corresponding to eg4 ones. 
*     
*     Remember to change the value of Eb (1 or 2 for 1.339 or 2.999 GeV ND3 runs 
*     respectively). The EG4 routine RLEG4 needs that unlike RLEG1B
*     
*     Also, beware of the target label mismatch - for ND3, the eg4 label is 2, for this code it is 1.
*     For that reason, I am using another index I_EG4 which will have 2 for ND3 and so on
*******************************************************************************
*     
*     
*     This program integrates over the CLAS EG4 target for the external 
*     corrections. It uses the subroutine 
*     
*     RLEG4(RADA,RADB,X,Y,Z,THA,PHIA,I_EG4, EB_INDEX) instead of RLEG1B(RADA,RADB,X,Y,Z,THA,PHIA,I)
*     
*     RADA & RADB are target radiation lengths for 
*     EG1 target before (RADB) and after (RADA) scattering.
*     
*     (THA,PHIA) are the horizontal and vertical angles in radians 
*     for scattered electrons. PHIA is zero at positive x axis and 
*     goes in clockwise direction. THETA is 0 at downstream Z and 
*     goes in counterclockwise dirextion.
*     
*     (X,Y,Z) defines the scattering point in cm.  The EG1 Proton
*     and deuteron targets are centered at -55.0 cm in Z.  
*     
*     Z is defined as downstream, X is to the left and Y straight up
*     
*     I = 1(ND3 target)
*     2(NH3 target)
*     3(carbon target)
*     4(empty target)
*     
*     IMODE defines the integration space
*     =1 Use center of target to represent integration space.
*     2 Integration over z of the target only.
*     3 Integrate over raster area at z=target center.
*     4 Full integration over target: VERY CPU consuming.
*     
*     Based on T.J. Liu's radlength.f subroutine written for the 
*     E143 target model from May94. Modified by  Al Tobias Jun97
*     RADA was rewritten to fit the HallB target system.
*     April 2000 RHF
**************************************************************

      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      INCLUDE 'prinplot.inc'
      
      INCLUDE 'kpaVarChanges.inc' !kpa: 8/6/13  (to control/change in RA by 20% for systematics)

      INTEGER EB_INDEX, I_EG4   !kp: EB_INDEX = 1 for 1.339 GeV, 2 for 1.998 GeV; I_EG4=2 for ND3

      INTEGER ITARG             !jixie, add to call xiaochao's subroutine
      REAL*8  PACKF             !jixie, add to call xiaochao's subroutine
      
      LOGICAL EXTFLAG, SAVEv, frwOK
      INTEGER NPTS, IS, NEXT, IMODE, I, J, K, L,
     >     PCOIL, DCOIL, IR, NPT, IPT, NPTTOT, NPTAREA,
     >     LSAVE, NZSTEP/10/, NRSTEP/7/, IDEN/8/
      REAL*8 TA, TLEN, Zv, SIGRU, SIGRP, Q2OUT, THD, !kp: THD added
     >     SANSWER, SIGRUTOT, SIGRPTOT, SIGMAR, MULT, 
     >     THET, RAD, DR, RMAX/0.7D0/, ANGLE,
     >     SIGRUSAVE(50), SIGRPSAVE(50), SIGRUNTSAVE(50),
     >     TBSAVE(50), TASAVE(50), SIGRUNT, SIGRUNTTOT,
     >     XTARG, YTARG, ELLCUT,
     >     chamberl, zcenter
c     -kag added

      REAL*8 XSHIFT/0.D0/, YSHIFT/0.D0/
      REAL*8 TB/0.026D0/, TA2/0.047D0/, TA5/0.040D0/ !E155

      COMMON /INFO/EXTFLAG,NEXT

      SIGRUTOT = 0.D0
      SIGRUNTTOT = 0.D0
      SIGRPTOT = 0.D0
      NEXT = NPTS

      FL_POL = .TRUE.
      FL_UNPOL = .TRUE.
      EXTFLAG = .FALSE.         !turn off momentarily.

      write(6,'(A,5F12.5)')'integ4.f L83: E,EP,THETAD,Q2: ',
     > E,EP,THETAD,Q2  
      
      SANSWER = SIGMAR(E,EP,THETAD)
      Q2OUT = Q2
      EXTFLAG = .TRUE.
      NEXT = 500
c     NEXT = 300000            !kp: 9/17/12 (got the following crash report)
c     Subscript out of range on file line 31, procedure sigmarNEW. Attempt to access the 300000-th element of variable pp_siga. Abort
c     we have PARAMETER (PP_NMAX = 505*505) in prinplot.inc, so giving the 504*504 should always work (as long as we are upto 500*500 calcs)
      NEXT = 504*504            !kp: 9/17/12

      chamberl = 1.0D0
c     zcenter = -55.0D0  !EG1b value
      zcenter = -100.930D0
c     kag	test change
c     kag       zcenter = zcenter - 0.5D0
c     chamber length in z and  target center in z (all in cm)
c     RMAX is the raster radius in cm (chamber radius is 0.73 cm)

      THET = THETAD*RADCON      !kp: this line missing, added on 1/13/12 (Thanks to Dr. Kuhn)
      THD = THETAD
c     print*,'integ4.f: thetad=',thetad

      I=2                       !for EG1, 1 is for ND3, 2 is for NH3
      I_EG4=1                   !added by Jixie. For EG4,  I_EG4 =1 is for NH3
      IF(TARG.EQ.'ND3') I = 1   
      IF(TARG.EQ.'ND3') I_EG4 = 2
      
!     kp: Replace calls to RLEG1B(..) for TA and TB with a new fn that uses this (http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Testfiles/Fortran/WdMake/RadLCppRootFortran/DropBox_RadLCppRootFortran/ChkCppFns/averageRADAgraphsEG4_20M_pol2pol3fits.gif) parameterization for TA (as a functionof theta) and a constant value for TB (assuming the scattering from the target center). The plot is for RADA (averaged over the target volume) vs theta. RADB = 0.01081 (http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Testfiles/Fortran/WdMake/RadLCppRootFortran/radLengthVsTheta.gif)
      IMODE=1
      IF(IMODE.EQ.1) THEN
c     kp: THET doesn't seem to have any value assigned. I guess due to some typo, so I am adding THET = THETAD*RADCON above 1/13/12
c     CALL RLEG4(TA,TB,0.D0,0.D0,zcenter, THET,0.D0,I_EG4, EB_INDEX)
         TB=0.01081
         if(THD.le.4.0) THEN
            TA=0.0249+(1.21d-4)*THD+(8.14d-6)*(THD**2.0) !kp: 4/7/12
         elseif((THD.gt.4.0).and.(THD.le.12.0)) THEN
            TA=0.0236+0.000637*THD-(2.09d-5)*(THD**2.0) !kp: 4/7/12
         elseif((THD.gt.12.0).and.(THD.le.52.0)) THEN
            TA=0.0268+(1.43d-4)*THD
     >           -(3.81d-6)*(THD**2.0)+(8.19d-8)*(THD**3.0) !kp: 4/7/12
         endif
c     TB=0.0D0              !kp: 5/25/12 (for Debug, turn it off during a normal run)
c     TA=0.0D0              !kp: 5/25/12 (for Debug, turn it off during a normal run)
c     TA=0.000001D0         !8888888888           kp: 1/7/13  (for Debug, to turn on all material in GSIM)
c     TB=4.0*TB             !kp:6/7/12 (for Debug)
c     TA=4.0*TA             !kp:6/7/12 (for Debug)

c     TA = 1.2*TA           !kp:8/6/13 (for Systematics)
         TA = TA*factorRA       !kp:8/6/13 kpaVarChanges.inc (for Systematics) (alternative to above line)
C     
C     Added by jixie:
C     The above TA,TB are used for ND3 target. For NH3, call xiaochao's subroutine
C     Here I always use ITARG=1 (1.0cm top NH3) for the inelastic for energies 3.0, 2.3, 2.0 and 1.3 GeV. 
C     For 1.1 GeV only the bottom cell was used (ITARG=11).
C     For 2.3 GeV or 3.0 GeV, and long target, only the top cell(1.0cm) was used (ITARG=1).
C     Jx: 20190214, add short target (0.5cm) case
C     Jx: 20190326, add bottom long target (1cm) cases for 2.3, 2.0 and 1.3 GeV
         IF(I_EG4 .EQ. 1) THEN
            ITARG = defaultTarget
            IF (EB_INDEX .EQ. 1) THEN
               ITARG = 11
            ELSEIF (EB_INDEX .GE. 4 .AND. defaultTarget .NE. 5 ) THEN   
               ITARG = 1
            ENDIF
            CALL rleg4_simp(TA,TB,PACKF,zcenter,THET,ITARG,EB_INDEX)
         ENDIF  
         write(*,'("TB,TA,PF=",3F10.5)') TB,TA,PACKF

         CALL RADCROSS(TB,TA,SIGRU,SIGRUNT,SIGRP)
         
         PP_SIGRADA(NPTS) = SIGRU
         PP_SIGRADANOTAIL(NPTS) = SIGRUNT
         PP_SIGRADP(NPTS) = SIGRP
      ELSEIF(IMODE.EQ.2) THEN
         TLEN = chamberl/FLOAT(NZSTEP)
         DO IS = 1,NZSTEP
            Zv = zcenter - chamberl/2.0D0 + TLEN/2.D0 + FLOAT(IS-1)*TLEN
            CALL RLEG4(TA,TB,0.D0,0.D0,Zv, THET,0.D0,I_EG4, EB_INDEX)
            
            CALL RADCROSS(TB,TA,SIGRU,SIGRUNT,SIGRP)

            SIGRUTOT = SIGRUTOT + SIGRU
            SIGRUNTTOT = SIGRUNTTOT + SIGRUNT
            SIGRPTOT = SIGRPTOT + SIGRP
         ENDDO
         PP_SIGRADA(NPTS) = SIGRUTOT/FLOAT(NZSTEP)
         PP_SIGRADANOTAIL(NPTS) = SIGRUNTTOT/FLOAT(NZSTEP)
         PP_SIGRADP(NPTS) = SIGRPTOT/FLOAT(NZSTEP)
      ELSEIF(IMODE.GE.3) THEN
         IF(IMODE.EQ.3) NZSTEP = 1
         TLEN = chamberl/FLOAT(NZSTEP)
         NPTAREA = 0.D0
         DO IS = 1,NZSTEP
            Zv = zcenter - chamberl/2.0D0 + TLEN/2.D0 + FLOAT(IS-1)*TLEN
            LSAVE = 0
            DR = RMAX/FLOAT(NRSTEP - 1)
            NPTTOT = 1 + IDEN*((NRSTEP-1)*NRSTEP)/2
            DO IR = 1,NRSTEP
               RAD = FLOAT(IR-1)*DR
               NPT = 1
               IF(IR.GT.1) NPT = IDEN*(IR-1)
               DO IPT = 1,NPT
                  ANGLE = 2.D0*PI*FLOAT(IPT-1)/FLOAT(NPT)
                  XTARG = RAD*COS(ANGLE)
                  YTARG = RAD*SIN(ANGLE)
                  CALL RLEG4(TA,TB,xtarg,ytarg,Zv,THET,0.D0,I_EG4,
     >             EB_INDEX)
                  NPTAREA = NPTAREA + 1
                  DO L = 1, LSAVE
                     IF(ABS(TBSAVE(L)-TB).LT.0.00001.AND.
     >                    ABS(TASAVE(L)-TA).LT.0.00001) THEN
                        SIGRU = SIGRUSAVE(L)
                        SIGRUNT = SIGRUNTSAVE(L)
                        SIGRP = SIGRPSAVE(L)
                        SAVEv = .FALSE.
                        GO TO 50
                     ENDIF
                  ENDDO

                  SAVEv = .TRUE.

                  CALL RADCROSS(TB,TA,SIGRU,SIGRUNT,SIGRP)

 50               SIGRUTOT = SIGRUTOT + SIGRU
                  SIGRUNTTOT = SIGRUNTTOT + SIGRUNT
                  SIGRPTOT = SIGRPTOT + SIGRP
                  IF(SAVEv) THEN
                     IF(LSAVE.LT.50) LSAVE = LSAVE + 1
                     TBSAVE(LSAVE) = TB
                     TASAVE(LSAVE) = TA
                     SIGRUSAVE(LSAVE) = SIGRU
                     SIGRUNTSAVE(LSAVE) = SIGRUNT
                     SIGRPSAVE(LSAVE) = SIGRP
                  ENDIF
 55            ENDDO
            ENDDO
            WRITE(6,'('' LSAVE= '',i3)')LSAVE
         ENDDO
         PP_SIGRADA(NPTS) = SIGRUTOT/FLOAT(NPTAREA)
         PP_SIGRADANOTAIL(NPTS) = SIGRUNTTOT/FLOAT(NPTAREA)
         PP_SIGRADP(NPTS) = SIGRPTOT/FLOAT(NPTAREA)

      ENDIF
      RETURN
      END
*******************************************************************************

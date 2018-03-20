      SUBROUTINE rleg4_simp(RADA,RADB,PF,Zv,THETA,ITARG,IEb)
**********************
C     X. Zheng, Feb 2018, simplified version to calculate rad length for EG4 
**********************
*     
*     Outputs:
*     
*     RADA & RADB are target radiation lengths for 
*     EG4 target before (RADB) and after (RADA) scattering.
*     
*     PF is the packing factor used in the calculation
*     
*     Inputs:
*     
*     simplified version:
*     input THETA is the scattering angle in radians
*     input Zv is vertex Z, not used at all. Always assume the center of target to be Z
*     
*     The EG4 Proton
*     and deuteron targets are centered at -100.93 cm in Zv.  
*     
*     Input
*     ITARG = target ID
*     1 (NH3 long top target,1cm) ! density 0.917 g/cm3
*     2 (ND3 target, 1cm)         ! density 1.056 g/cm3
*     3 (empty cut with helium, 1cm) ! helium density 0.145 g/cm3 at 1K
*     4 (long carbon, 1cm)          ! carbon density 2.166 g/cm3
*     5 (NH3 short target, 0.5cm)
*     6 (short carbon, 0.5cm)
*     7 (long carbon no helium, 1cm)
*     8 (empty cup without helium, 1cm)
*     9 (short carbon without helium, 0.5cm) ! carbon density 2.166 g/cm3
*     10 (not used)
*     11 (NH3 bottom long target, 1cm)
*     
*     Input
*     IEb is energy ID, and when combined with ITARG, is used to assign the packing factor values
*     = 1 (Beam energy 1.1 GeV)
*     = 2 (Beam energy 1.3 GeV)
*     = 3 (2.0 GeV)
*     = 4 (2.3 GeV)
*     = 5 (3.0 GeV)
*     
*     Other tags:
*     IPF is packing factor ID
*     = 1 using Sarah's values
C     https://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_28,_2011
*     = 2 using XZ's exclusive-channel analysis values
C     see page 66 of excl channel analysis note
C     
C     DEBUG is flag for printing debugging info
*     
**************************************************************
      
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      REAL*8 PF_INI(2,5,11) ! iflag, iEb, itarg
      integer IPF/1/            ! flag for selecting packing factor values.
!     IPF=1 is using Sarah's, IPF=2 is using XZ's exclusive results
      
      logical DEBUG/.true./

      DEBUG=.false.

      if ((ITARG.ge.7).and.(ITARG.le.10)) then
         print *,'this routine works only for targets with helium'
         goto 200
      endif
      
      PI=acos(-1.)
      r2d = 180./PI
      if (DEBUG) print *,'RLEG4_simp check inputs: ',
     >     theta*r2d,Zv,ITARG,IEb

C     initialize packing factors
      do l=1,5
         do j=1,11
            do k=1,2
               PF_INI(k,l,j)=0.
            enddo
         enddo
      enddo
C     IPF=1 using Sarah's values https://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_28,_2011
      PF_INI(1,5,1)=0.782
      PF_INI(1,4,1)=0.682
      PF_INI(1,4,5)=0.720
      PF_INI(1,3,1)=0.716
      PF_INI(1,3,11)=PF_INI(1,3,1) ! not available from Sarah's
      PF_INI(1,2,1)=0.717
      PF_INI(1,2,5)=0.602
      PF_INI(1,2,11)=0.657
      PF_INI(1,1,11)=0.625 

C     Preliminary PACKING FRACTION OF ND3 (from S. Phillips) http://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_28%2C_2011
      PF_INI(1,2,2)=0.624
      PF_INI(1,3,2)=0.764
      
C     IPF=2 using XZ's exclusive values (see analysis note)
      PF_INI(2,5,1)=0.65
      PF_INI(2,4,1)=0.65
      PF_INI(2,4,5)=0.30
      PF_INI(2,3,1)=0.65
      PF_INI(2,3,11)=0.65
      PF_INI(2,2,1)=0.70
      PF_INI(2,2,5)=0.35
      PF_INI(2,2,11)=0.70
      PF_INI(2,1,11)=0.75
      
C     for ND3 no excl results available so keep using Sarah's
C     Preliminary PACKING FRACTION OF ND3 (from S. Phillips) http://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_28%2C_2011
      PF_INI(2,2,2)=0.624
      PF_INI(2,3,2)=0.764

      PF=PF_INI(IPF,IEb,ITARG)
      IF (DEBUG) print *,'packing factor=',PF

C     DENSITIES OF VARIOUS TARGET MATERIALS (rho g/cm3)	

      RHO_LHE=.145d+0           !DENSITY OF LHE at 1K
      RHO_N15H3=.917d+0         !DENSITY OF NH3 (15NH3)
      RHO_N15D3= 1.056d+0       !DENSITY OF ND3 (15ND3)
      RHO_CARBON=2.265d+0       !DENSITY OF CARBON        
      RHO_AL=2.7d+0             !DENSITY OF ALUMINUM
      RHO_KAPTON=1.42d+0        !DENSITY OF KAPTON
      RHO_PCTFE=2.2d+0          !DENSITY OF PCTFE (insert teflon rho)
      RHO_AIR=.00129d+0         !DENSITY OF AIR 
      RHO_MYLAR= 1.39d+0        !DENSITY OF MYLAR (ignore aluminum)
      RHO_CEREX= 1.58d+0        !DENSITY OF CEREX (calc by mylar ratio)
      RHO_CU = 8.96d+0          !DENSITY OF COPPER
      RHO_FEP =2.2d+0           !DENSITY OF FEP (insert teflon rho) 
      RHO_NMR =6.391d+0         !DENSITY OF NMR = 0.38*FEPD+0.62*CUD
      
C     https://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_8,_2010
c$$$  Sarah,
c$$$  Here is a message from Oscar Rondon regarding the density of solid ammonia.
c$$$  Bottom line:
c$$$  14NH3 = 0.867 g/cc
c$$$  14ND3 = 1.007 g/cc
c$$$  15NH3 = 0.917 g/cc
c$$$  15ND3 = 1.057 g/cc
      
c     http://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_14%2C_2011#Target_Composition 
      RHO_N14H3D=0.867d+0       !Density of 14NH3
      RHO_N14D3D=1.007d+0       !Density of 14ND3
      RHO_N15H3=0.917d+0
      RHO_N15D3=1.057d+0
      RHO_NH3=0.98*RHO_N15H3+0.02*RHO_N14H3 !Composition not 100% pure 15NH3 rather the isotop percentages were 98% 15NH3 & 2% 14NH3
      RHO_ND3=0.97*RHO_N15D3+0.02*RHO_N14D3+0.01*RHO_N15H3 !Composition not 100% pure 15ND3 rather the isotop percentages were 97% 15ND3, 2% 14ND3 & 1% 15NH3

      if (DEBUG) print *,'check density:',RHO_NH3,RHO_ND3

C     XZ: using values from the Hall C table:
C     https://userweb.jlab.org/~jones/rss/rsstgt.htm
      RHO_NH3=0.917d+0
      RHO_ND3=1.057d+0
      if (DEBUG) print *,'check density:',RHO_NH3,RHO_ND3

C     RADIATION CONSTANTS OF VARIOUS TARGET MATERIALS (Xo g/cm2)

      X0_LHE=94.32d+0           !X0 OF LHE, PDG
      X0_CARBON=42.70d+0        !X0 OF CARBON , PDG
      X0_AL=24.01d+0            !X0 OF ALUMINUM, PDG
      X0_KAPTON=40.56d+0        ! X0 OF KAPTON, PDG
      X0_PCTFE=34.84d+0         ! X0 OF PCTFE (teflon), PDG
      X0_AIR=36.62d+0           !X0 OF AIR, PDG
      X0_MYLAR=39.95d+0         !X0 OF MYLAR, PDG
      X0_CEREX= 45.4d+0         !X0 FOR CEREX
      X0_CU= 12.86d+0           !X0 FOR CU, PDG
      X0_FEP=37.84d+0           !X0 FOR FEP (teflon)
      X0_NMR=0.38*X0_FEP+0.62*X0_CU !X0 FOR NMR

C     using values from the Hall C table:
C     https://userweb.jlab.org/~jones/rss/rsstgt.htm
      X0_NH3=43.255d+0          !X0 OF NH3 
      X0_ND3=50.500d+0          !X0 OF ND3

C     XZ: calculate X0 using weighted average <- note I can't find 15N and using 14N is probably incorrect
c$$$  X0_N = 37.99d+0         ! X0 of 14N_2
c$$$  X0_H = 63.04d+0         ! X0 of H_2
c$$$  X0_D = 125.97d+0        ! X0 of D_2
c$$$  X0_NH3 = (14./17)*X0_N+(3./17)*X0_H  !X0 OF NH3, gives 40.85
c$$$  X0_ND3 = (14./20)*X0_N+(6./20)*X0_D  !X0 OF ND3, gives 48.06

      if (DEBUG) print *,'check X0:',X0_NH3,X0_ND3

C     DIMENSIONS OF EG4 TARGET CHAMBERS/Cells

      if ((ITARG.eq.5).or.(ITARG.eq.6).or.(ITARG.eq.9)) then
         CHAMBERL= 0.5d+0       !INNER LENGTH OF TARGET CHAMBER IN CM 
      else                      ! long targets
         CHAMBERL= 1.0d+0       !INNER LENGTH OF TARGET CHAMBER IN CM 
      endif
      CHAMBERR= 0.725d+0        !INNER RADIUS OF TARGET CHAMBER IN CM
c     Table B-B: Inner Diameter (ID) is 14.5 and Outer one (OD) is 15.0mm so, the inner radius in cm is 1.45/2=0.725
      WALLT= 0.025d+0           !WALL THICKNESS OF CHAMBER (AVERAGE)
c     I think its equal to (OD-ID)/2 = 0.05/2 = 0.025

C     XZ now assign target thickness and density
      
c     carbon: (Dr. Ripani's measurements: thin disk (Target# 6) = 1.08±0.02 mm & thick disk (Target# 4) = 2.16±0.05 mm)
      
      NMRT=0.00838d+0           !EFFECTIVE THICKNESS OF NMR COILS=83.8 microns

      CARBONL=0.
      
      if ((ITARG.eq.4).or.(ITARG.eq.7)) then
         CARBONL=0.216d+0       !THICKNESS OF CARBON TARGET in cm
         DTARG1=CARBONL/2.
         DTARG2=CARBONL/2.
         RHO_TARG=RHO_CARBON
         X0_TARG=X0_CARBON
      elseif ((ITARG.eq.6).or.(ITARG.eq.9)) then
         CARBONL=0.108d+0       ! thin carbon target in cm
         DTARG1=CARBONL/2.
         DTARG2=CARBONL/2.
         RHO_TARG=RHO_CARBON
         X0_TARG=X0_CARBON
      elseif ((ITARG.eq.1).or.(ITARG.eq.5).or.(ITARG.eq.11)) then
         DTARG1=CHAMBERL*PF/2.
         DTARG2=CHAMBERL*PF/2.
         RHO_TARG=RHO_NH3
         X0_TARG=X0_NH3
      elseif (ITARG.eq.2) then
         DTARG1=CHAMBERL*PF/2.
         DTARG2=CHAMBERL*PF/2.
         RHO_TARG=RHO_ND3
         X0_TARG=X0_ND3
      else                      ! nothing in the target!
         RHO_TARG=0.
         X0_TARG=1000000.
      endif
      
      BANJOL = 2.13             ! total thickness of banjo in cm, see http://www.jlab.org/~adhikari/myHome/secure/Corrections/Mcor/CheckMomCorEg4/MT_highPnTheta/RootFiles/BanjoWindows.html

      BANJOL = 2.19             ! should be 2.19 cm, see dwg.No.66840-E-04410 (item IV in table A). WE can use either values. The uncertainty in this number is probably more than 0.05cm anyways!
      
      if ((ITARG.eq.7).or.(ITARG.eq.8).or.(ITARG.eq.9)) then ! these are empty targets
         DHE1=0.
         DHE2=0.
      else
         DHE1=(BANJOL)/2.-DTARG1
         DHE2=(BANJOL)/2.-DTARG2
      endif
      
      if (DEBUG) then
         print *,'target length in cm',CHAMBERL
         print *,'carbon thickness in cm',CARBONL
         print *,'target effective thickness in cm',DTARG1,DTARG2
         print *,'helium total thickness in cm',DHE1,DHE2
      endif

C     now start calculating thicknesses for Before and After
      RADB=0.
      RADA=0.
      
C     Z-dimensions/Locations (in cm) of EG4 beam line materials 
C     and their thicknesses in cm
      ZVVENTR=-166.1d+0         ! upstream vacuum vessel entrance window, 50microns Al
      DVVENTR=50./1.d+4
      RADB=RADB+DVVENTR*RHO_AL/X0_AL
      if (DEBUG)
     >     write(*,'("RADB1=",2F12.6)')RADB,DVVENTR*RHO_AL/X0_AL
C     
C     These two shields are noted in the drawing as "not in beamline"
C     DWG # 66840-E-04410 (eg4a Run Winter 2006)
C     
      Z77KSLD1=-127.1d+0        ! 77-K upstream heat-shield, 14 microns Al
      D77KSLD1=14./1.d+4
c$$$  RADB=RADB+D77KSLD1*RHO_AL/X0_AL
c$$$  if (DEBUG) 
c$$$  >       write(*,'("RADB2=",2F12.6)')RADB,D77KSLD1*RHO_AL/X0_AL
      
      Z4KSLD1=-121.0d+0         ! 4-Kelvin Shield (Upstream), 14 microns Al
      D4KSLD1=14./1.d+4
c$$$  RADB=RADB+D4KSLD1*RHO_AL/X0_AL
c$$$  if (DEBUG) 
c$$$  >       write(*,'("RADB3=",2F12.6)')RADB,D4KSLD1*RHO_AL/X0_AL
      
      ZTARGCNTR=-100.93d+0      ! Zv of Target (chamber) center (Mid point of Bajo-windows) in cm
      
      ZBANJOENTR=ZTARGCNTR-BANJOL/2. ! Bajo Entrance Window , 71 microns Al
      DBANJOENTR=71./1.d+4
      RADB=RADB+DBANJOENTR*RHO_AL/X0_AL
      if (DEBUG) 
     >     write(*,'("RADB4=",2F12.6)')RADB,DBANJOENTR*RHO_AL/X0_AL
      
      ZTARGENTR=ZTARGCNTR-CHAMBERL/2 ! Z location (cm) of Target Chamber Entrance Window, 25 microns kapton
      DTARGENTR=25./1.d+4
      RADB=RADB+DTARGENTR*RHO_KAPTON/X0_KAPTON
      if (DEBUG) write(*,'("RADB5=",2F12.6)')
     >     RADB,DTARGENTR*RHO_KAPTON/X0_KAPTON

      RADB=RADB+DHE1*RHO_LHE/X0_LHE
      if (DEBUG) write(*,'("RADB6=",2F12.6)')RADB,DHE1*RHO_LHE/X0_LHE
      RADB=RADB+DTARG1*RHO_TARG/X0_TARG
      if (DEBUG) write(*,'("RADB7=",5F12.6)')
     >     RADB,DTARG1*RHO_TARG/X0_TARG,DTARG1,RHO_TARG,X0_TARG

      if (DEBUG) write(*,'("\n")')
      
      RADA=RADA+DTARG2*RHO_TARG/X0_TARG
      if (DEBUG) 
     >     write(*,'("RADA1=",2F12.6)')RADA,DTARG2*RHO_TARG/X0_TARG

      RADA=RADA+DHE2*RHO_LHE/X0_LHE
      if (DEBUG) 
     >     write(*,'("RADA2=",2F12.6)')RADA,DHE2*RHO_LHE/X0_LHE

      ZTARGEXIT=ZTARGCNTR+CHAMBERL/2 ! Z location (cm) of Target Chamber Exit Window, 25 microns kapton
      DTARGEXIT=25./1.d+4
      RADA=RADA+DTARGEXIT*RHO_KAPTON/X0_KAPTON
      if (DEBUG) 
     >     write(*,'("RADA3=",2F12.6)')
     >     RADA,DTARGEXIT*RHO_KAPTON/X0_KAPTON
      
      ZBANJOEXIT=-99.86d+0      ! Banjo Exit Window, 71 microns Al
      ZBANJOEXIT=ZBANJOENTR+BANJOL
      DBANJOEXIT=71./1.d+4
      RADA=RADA+DBANJOEXIT*RHO_AL/X0_AL
      if (DEBUG) 
     >     write(*,'("RADA4=",2F12.6)')RADA,DBANJOEXIT*RHO_AL/X0_AL
      
      if (DEBUG) print *,'Banjo:',ZBANJOENTR,ZBANJOEXIT,BANJOL 
      if (DEBUG) print *,'Targt:',ZTARGENTR,ZTARGEXIT,CHAMBERL 

      ZFOIL=-94.1d+0            ! insulating foil, aluminized mylar, 6 microns
      DFOIL=6./1.d+4
      RADA=RADA+DFOIL*RHO_MYLAR/X0_MYLAR
      if (DEBUG) 
     >     write(*,'("RADA5=",2F12.6)')RADA,DFOIL*RHO_MYLAR/X0_MYLAR
      
      Z77KSLD2=-92.6d+0         ! heat-shield foil, 25 microns Al
      D77KSLD2=25./1.d+4
      RADA=RADA+D77KSLD2*RHO_AL/X0_AL
      if (DEBUG) 
     >     write(*,'("RADA6=",2F12.6)')RADA,D77KSLD2*RHO_AL/X0_AL
      
C     now superinsulation:

c$$$  https://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_8,_2010
c$$$  
c$$$  Each layer of superinsulation consists of one piece of aluminized mylar
c$$$  (~6 um thick) and two pieces of Cerex.  Cerex is a lightweight 
c$$$  fabric spun from nylon and comes in various weights. 
c$$$  The thickness of aluminum on the mylar is incredibly thin (~250
c$$$  angstroms).   Dave Kashy (I think) made the following measurements 
c$$$  of the *areal* density of super-insulation:
c$$$  one piece Al mylar: 0.88 mg/cm2
c$$$  one piece Cerex:  1.0 mg/cm2


      
      ZSUPERINS1=-93.5d+0       ! 10-layered super-insulation, each layer is made of one sheet of aluminized mylar, 0.88mg/cm^2, and 2 sheets of cerex, 1.0 mg/cm^2/sheet 
      DSUPERINS1a=10*0.88/1000./RHO_MYLAR
      DSUPERINS1b=10*2*1.0/1000./RHO_CEREX
      RSUPERINSIR1 = 1.25d+0    !RADIUS OF HOLE CUT IN ALL SHEETS EXCEPT ONE (M. Zarecky Table-A, inside diameter 2.5cm)
      if (DEBUG) print *,'min angle for super ins 1:',
     >     atan(RSUPERINSIR1/(ZSUPERINS1-ZTARGCNTR))*180./PI, ! this gives 9.55 deg
     >     RSUPERINSIR1,(ZSUPERINS1-ZTARGCNTR)
      IF (THETA.gt.atan(RSUPERINSIR1/(ZSUPERINS1-ZTARGCNTR))) then
         RADA=RADA+DSUPERINS1a*RHO_MYLAR/X0_MYLAR
         if (DEBUG) write(*,'("RADA7a=",5F12.6)')
     >        RADA,DSUPERINS1a*RHO_MYLAR/X0_MYLAR,
     >        DSUPERINS1a,RHO_MYLAR,X0_MYLAR
         RADA=RADA+DSUPERINS1b*RHO_CEREX/X0_CEREX
         if (DEBUG) write(*,'("RADA7b=",5F12.6)')
     >        RADA,DSUPERINS1b*RHO_CEREX/X0_CEREX,
     >        DSUPERINS1b,RHO_CEREX,X0_CEREX
      ENDIF
      
      ZSUPERINS2=-90.5d+0       ! 20-layered super-insulation, each layer is made of one sheet of aluminized mylar, 0.88mg/cm^2, and 2 sheets of cerex, 1.0 mg/cm^2/sheet 
      DSUPERINS2a=20*0.88/1000./RHO_MYLAR
      DSUPERINS2b=20*2*1.0/1000./RHO_CEREX
      RSUPERINSIR2 = 1.25d+0    !RADIUS OF HOLE CUT IN ALL SHEETS EXCEPT ONE (M. Zarecky Table-A, inside diameter 2.5cm)
      if (DEBUG) print *,'min angle for super ins 2:',
     >     atan(RSUPERINSIR2/(ZSUPERINS2-ZTARGCNTR))*180./PI, ! this gives 6.83 deg
     >     RSUPERINSIR2,(ZSUPERINS2-ZTARGCNTR)
      IF (THETA.gt.atan(RSUPERINSIR2/(ZSUPERINS2-ZTARGCNTR))) then
         RADA=RADA+DSUPERINS2a*RHO_MYLAR/X0_MYLAR
         if (DEBUG) write(*,'("RADA8a=",5F12.6)')
     >        RADA,DSUPERINS2a*RHO_MYLAR/X0_MYLAR,
     >        DSUPERINS2a,RHO_MYLAR,X0_MYLAR
         RADA=RADA+DSUPERINS2b*RHO_CEREX/X0_CEREX
         if (DEBUG) write(*,'("RADA8b=",5F12.6)')
     >        RADA,DSUPERINS2b*RHO_CEREX/X0_CEREX,
     >        DSUPERINS2b,RHO_CEREX,X0_CEREX
      ENDIF

      ZVVEXIT1=-78.3d+0         ! downstream vacuum vessel beam-line exit window, 50 microns Al, for theta<5 deg
      DVVEXIT1=50./1.d+4
      
      ZVVEXIT2=-78.3d+0         ! for theta>5deg, there is the 280 microns Al vacuum vessel exit window axial
      DVVEXIT2=280./1.d+4
      IF (THETA.le.(5.*PI/180.)) then
         RADA=RADA+DVVEXIT1*RHO_AL/X0_AL
         if (DEBUG) 
     >        write(*,'("RADA9=",2F12.6)')RADA,DVVEXIT1*RHO_AL/X0_AL
      ENDIF

      RADA=RADA/cos(THETA)
      if (DEBUG) write(*,
     >     '("RADA10=",F12.6," cos(th)=",F12.6)')RADA,cos(THETA)

      IF (THETA.gt.(5.*PI/180.)) then ! calculate thickness of the 280 micron Al window
         xC=108.83
         R1=86.54
         R2=R1+DVVEXIT2

         tan2=tan(THETA)*tan(THETA)
         x1=(xC-sqrt(xC*xC-(1+tan2)*(xC*xC-R1*R1)))/(1+tan2)
         x2=(xC-sqrt(xC*xC-(1+tan2)*(xC*xC-R2*R2)))/(1+tan2)
         y1=x1*tan(THETA)
         y2=x2*tan(THETA)
         RADA=RADA+sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*RHO_AL/X0_AL
         if (DEBUG) 
     >        write(*,'("RADA11=",2F12.6)')RADA,
     >        sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))*RHO_AL/X0_AL
      ENDIF
C     last one: air gap to drift chamber, 78.3 cm
C     comment this out because this is in GEANT
C     
C     RADA=RADA+78.3*RHO_AIR/X0_AIR
C     if (DEBUG) 
C     >     write(*,'("RADA12=",2F12.6)')RADA,78.3*RHO_AIR/X0_AIR
      
      ZVVEXIT3=0d+0             ! downstream beam pipe vacuum window, 71 microns Al, way downstream and not considered below
      DVVEXIT3=71./1.d+4
      
      if (DEBUG) write(*,'("RADB,A=",2F10.5)')RADB,RADA
 200  END
      

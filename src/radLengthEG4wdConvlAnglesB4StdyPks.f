      SUBROUTINE rleg4(RADA,RADB,X,Y,Zv,THA,PHIA,I,Eb)
**********************
*     
*     3/2/13: It seems that the Z variable here (for vertex Z) may conflict with
*     Z-variable (for atomic # of targets) from common-block(s) in various 
*     RCSLACPOL files. So, I changed the Z here to Zv.
*     I will have to do the same in integ4.f too.
*     
*************************************************************
*     RADA & RADB are target radiation lengths for 
*     EG4 target before (RADB) and after (RADA) scattering.
*     
*     Considers only the case of ND3 polarized target (1 cm long)
*     
*     (THA,PHIA) are the conventional spherical polar and azimuthal angles
*     in radians, calculated as THA = acos(cz) and PHIA = atan2(cy,cz)
*     
*     (X,Y,Zv) defines the scattering point in cm.  The EG4 Proton
*     and deuteron targets are centered at -100.93 cm in Zv.  
*     
*     Zv is defined as downstream, X is to the left and Y straight up
*     
*     I = 2(ND3 target)
*     
*     Eb = 1 (Beam energy 1.339 GeV)
*     = 2 (Beam energy 1.989 GeV)
*     
*     
*     Written by K. Adhikari, based on the following code from K. Griffioen 
*     http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/SimStuffs/RadCorStuff/rleg1KGriffioenRadCor.f
*     
*     
*     
*     Some of the drawings and related derivations that I made to understand
*     the geometries and the effective thickness calculations can be found
*     at the following links:
*     OVC window parameters: http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/SimStuffs/RadCorStuff/OVC_exit.svg
*     OVC tilt and effective thickness: http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/SimStuffs/RadCorStuff/OVC_tilt.svg
*     Calculation of RADXY: http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/SimStuffs/RadCorStuff/radialDistInFoil.svg
*     To find out more about EG4 beam line materials, please look at the following web page that I made for this purpose:
*     http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/SimStuffs/RadCorStuff/rleg1KGriffioenRadCor.html
*     
*     If there are two variable names that are mostly identical except for the last letter 1 or 2, such as FOURKSLDR1, & FOURKSLDR2
*     then 1 indicates upstream and 2 indicates downstream.
*     
*     
*     kpa: 2/17/12  Added d+0 at the end of every assignment expression for real*8 variables to ensure they are double precision 
*     without default round-off that caused me initial troubles. (http://www-classes.usc.edu/engr/ce/108/text/fbk01.htm)
**************************************************************
      
      IMPLICIT NONE
      INTEGER I, Eb, N
      REAL*8 RADA, RADB, Zv, X, Y, THA, PHIA 
      REAL*8 TARGETENTRZ, TARGETEXITZ, TARGETCNTRZ, CARBONEXITZ 
      REAL*8 CHAMBERL, CHAMBERR, WALLT, CARBONL

      REAL*8 BEAMENTERRL, HEATSLDBRL, BANJOENTRL
      REAL*8 TARGETENTRL, HELIUMENTRL, TARGETEXITRL, OVCEXITRL
      REAL*8 HELIUMEXITRL, BANJOEXITRL, FOURKSLDRL1, SEV7KSLDRL1
      REAL*8 FOURKSLDRL2, SEV7KSLDRL2, SUPERINSRL1, SUPERINSRL2

      REAL*8 BANJOENTRZ, BANJOEXITZ, FOURKSLDZ1, FOURKSLDZ2
      REAL*8 SUPERINSZ1, SUPERINSZ2, SEV7KSLDZ1, SEV7KSLDZ2, OVCEXITZ
      REAL*8 BANJOEXITR, SUPERINSR1, SUPERINSR2, FOURKSLDR2
      REAL*8 SEV7KSLDR2, OVCEXITR, SUPERINSIR1, SUPERINSIR2
      
      REAL*8 LHED, NH3D, ND3D, CARBOND, ALD, KAPTOND, NMRD
      REAL*8 PCTFED, AIRD, MYLARD, HEGASD, POLYD, CEREXD

      REAL*8 LHEXO, NH3XO, ND3XO, CARBONXO, ALXO, CEREXXO 
      REAL*8 KAPTONXO, PCTFEXO, AIRXO, MYLARXO, POLYXO, NMRXO
      REAL*8 NH3D14, ND3D14, BANJOL, SI1LAYER, TRKFACTOR

      REAL*8 ROVC, YOVC, VOVC, THETAOVC, ZOVC, SIN2TH, COS2TH

      REAL*8 FND3, FNH3


      REAL*8 LT, RADXY, R, HYPOT, RTRI, ROUT, RWALL, RIN 
      REAL*8 LIN, LWALL, LOUT, POLTARGETRL, THETATRI
c     REAL*8 EMPTYTARGETRL, CARBONTARGETRL
      REAL*8 NMRT, RNMR, LNMR
      REAL*8 CUXO, CUD, FEPXO, FEPD

      DIMENSION RPART(6,2)      ! ARRAY WHICH HOLDS THE Zv and RADIUS  VALUES
! OF THE PARTICLE AT EACH OUTGOING TARGET LAYER
      REAL*8 RPART,RADXY_EG4
      real*8 r2d, ABSZH
      real*8 DZI, TERM1, TERM2C, TERM2W, TERM2N, RSPHI
      real*8 SINTH,COTTH,TCC,TWW,TNN,DIST_SIC !,DBGXYZ
      INTEGER PLUSMINUS         !Plus or Minus
      real*8 RADA_4DBG,RADA_4DBG1,RADA_4DBG2,RADA_4DBG3,RADA_4DBG4
      real*8 RADA_W1, RADA_W2, RADA_W3, RADA_W4, RADA_W5, RADA_W6
      r2d = 57.2957795130823229d+0
c     print *,'kp: ',tha*r2d,phia*r2d,x,y,z





c     kp: added on 2/15/12 (to avoid crash because of COT(THA) or 1.0/sin(THA) below
      IF(THA.GT.0.0d+0) THEN
         CONTINUE
      ELSE 
c     WRITE(*,*)  'THETA out of range.'
         GO TO 200
      ENDIF




C     DIMENSIONS OF EG4 TARGET CHAMBERS/Cells

      CHAMBERL= 1.0d+0          !INNER LENGTH OF TARGET CHAMBER IN CM  (This is for "long" targets; 0.5 cm for "short" targets)
      CHAMBERR= 0.725d+0        !INNER RADIUS OF TARGET CHAMBER IN CM
c     Table B-B: Inner Diameter (ID) is 14.5 and Outer one (OD) is 15.0mm so, the inner radius in cm is 1.45/2=0.725
      WALLT= 0.025d+0           !WALL THICKNESS OF CHAMBER (AVERAGE)
c     I think its equal to (OD-ID)/2 = 0.05/2 = 0.025
      CARBONL=0.227d+0          !THICKNESS OF CARBON TARGET
c     = 0.108 (short) & = 0.216 (long) (Only type 4 (long) used during ND3 run period
c     (Dr. Ripani's measurements: thin disk (Target# 6) = 1.08±0.02 mm & thick disk (Target# 4) = 2.16±0.05 mm)
      NMRT=0.00838d+0           !EFFECTIVE THICKNESS OF NMR COILS=83.8 microns



C     Z-dimensions/Locations (in cm) of EG4 beam line materials 
      
      SEV7KSLDZ1=-127.1d+0      ! 77-K upstream heat-shield
      FOURKSLDZ1=-121.0d+0      ! 4-Kelvin Shield (Upstream)
      BANJOENTRZ=-102.00d+0     ! Bajo Entrance Window
      TARGETCNTRZ=-100.93d+0    ! Zv of Target (chamber) center (Mid point of Bajo-windows) in cm
      TARGETENTRZ=TARGETCNTRZ-CHAMBERL/2 ! Z location (cm) of Target Chamber Entrance Window
      TARGETEXITZ=TARGETCNTRZ+CHAMBERL/2 ! Z location (cm) of Target Chamber Exit Window
      BANJOEXITZ=-99.86d+0      ! Bajo Exit Window
      FOURKSLDZ2=-93.5d+0       ! 4-Kelvin Shield (Downstream), same as that of 10-layered super-insulation
      SUPERINSZ1=-93.5d+0       ! 10-layered super-insulation (together with 4-K downstream shield)
      SEV7KSLDZ2=-92.6d+0       ! 77-K downstream heat-shield (in between two super-insulations)
      SUPERINSZ2=-90.5d+0       ! 10-layered super-insulation 
      OVCEXITZ=-78.3d+0         ! OVC (Outer Vacuum Can or Vacuum Vessel) Exit Window (Axial) (actually the center/bottom of the curved window)
      
      BANJOL = abs(BANJOENTRZ-BANJOEXITZ)

c     print *,TARGETENTRZ,TARGETEXITZ,BANJOL 

c     Outer radii of all layers (beyond the target exit and upto OVC exit)
c     Calculated by taking the Z distance from the targetface to the layer and multiplying it by tan(55)
c     This is done since acceptance of the magnet is 50 degrees (55 just gives a little room for alignment errors).

      BANJOEXITR = 2.24d+0      !radius of acceptance for banjo exit window = (BANJOEXITZ - TARGETENTRZ)*tan(55*pi/180)   
      FOURKSLDR2 = 11.33d+0     !radius of acceptance for 4-Kelvin Shield (Downstream) = (FOURKSLDZ2 - TARGETENTRZ)*tan(55*pi/180)
      SUPERINSR1 = 11.33d+0     !radius of acceptance for 10-layered super-insulation = (SUPERINSZ1 - TARGETENTRZ)*tan(55*pi/180)
      SEV7KSLDR2 = 12.61d+0     !radius of acceptance for 77-K downstream heat-shield = (SEV7KSLDZ2 - TARGETENTRZ)*tan(55*pi/180)
      SUPERINSR2 = 15.61d+0     !radius of acceptance for 20-layered super-insulation = (SUPERINSZ2 - TARGETENTRZ)*tan(55*pi/180)
      OVCEXITR = 33.03d+0       !radius of acceptance for OVC window = (OVCEXITZ - TARGETENTRZ)*tan(55*pi/180)

      SUPERINSIR1 = 1.25d+0     !RADIUS OF HOLE CUT IN ALL SHEETS EXCEPT ONE (M. Zarecky Table-A, inside diameter 2.5cm)
      SUPERINSIR2 = 1.25d+0     !RADIUS OF HOLE CUT IN ALL SHEETS EXCEPT ONE (M. Zarecky Table-A, inside diameter 2.5cm)


C     DENSITIES OF VARIOUS TARGET MATERIALS (rho g/cm3)	

      LHED=.145d+0              !DENSITY OF LHE
      NH3D=.917d+0              !DENSITY OF NH3 (15NH3)
      ND3D= 1.056d+0            !DENSITY OF ND3 (15ND3)
      CARBOND=2.265d+0          !DENSITY OF CARBON        
      ALD=2.7d+0                !DENSITY OF ALUMINUM
      KAPTOND=1.42d+0           !DENSITY OF KAPTON
      PCTFED=2.2d+0             !DENSITY OF PCTFE (insert teflon rho)
      AIRD=.00129d+0            !DENSITY OF AIR 
      MYLARD= 1.39d+0           !DENSITY OF MYLAR (ignore aluminum)
      CEREXD= 1.58d+0           !DENSITY OF CEREX (calc by mylar ratio)
      CUD = 8.96d+0             !DENSITY OF COPPER
      FEPD =2.2d+0              !DENSITY OF FEP (insert teflon rho) 
      NMRD =6.391d+0            !DENSITY OF NMR = 0.38*FEPD+0.62*CUD

c     http://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_14%2C_2011#Target_Composition 
      NH3D14=0.867d+0           !Density of 14NH3
      ND3D14=1.007d+0           !Density of 14ND3
      NH3D=0.98*NH3D+0.02*NH3D14 !Composition not 100% pure 15NH3 rather the isotop percentages were 98% 15NH3 & 2% 14NH3
      ND3D=0.97*ND3D+0.02*ND3D14+0.01*NH3D !Composition not 100% pure 15ND3 rather the isotop percentages were 97% 15ND3, 2% 14ND3 & 1% 15NH3


C     RADIATION CONSTANTS OF VARIOUS TARGET MATERIALS (Xo g/cm2)

      LHEXO=94.32d+0            !X0 OF LHE
      NH3XO=43.26d+0            !X0 OF NH3
      ND3XO=50.5002d+0          !X0 OF ND3 
      CARBONXO=42.70d+0         !X0 OF CARBON 
      ALXO=24.01d+0             !X0 OF ALUMINUM
      KAPTONXO=40.56d+0         !XO OF KAPTON
      PCTFEXO=34.84d+0          !XO OF PCTFE (teflon0
      AIRXO=36.66d+0            !XO OF AIR
      MYLARXO=39.95d+0          !XO OF MYLAR 
      CEREXXO= 45.4d+0          !XO FOR CEREX
      CUXO= 12.8616d+0          !XO FOR CU
      FEPXO=37.84d+0            !XO FOR FEP (teflon)
      NMRXO=0.38*FEPXO+0.62*CUXO !XO FOR NMR




C     Radiation Lengths (precisely - Radiation Exponents)
c     Known radiation exponents defined for windows etc. Calculated by taking X0/rho (rho = density,X0=rad. length in g/cm^2, L=X0/rho=rad. length in cm) 
c     to get radiation length=L (Wikipedia, PDG - definition (page 17), tables). Then take thickness of material divided by this L to get answer. 

c     BEAMENTERRL =0.0008         ! eg1: upstream vacuum vessel entrance window = 71 microns Alum (50 microns - table A-I.)
      BEAMENTERRL =0.00056d+0   ! eg4:  == 0.0008*50/71 or x/(X0*rho)=0.0050/(24.01/2.6999)
c     FOURKSLDRL1 =0.00022        ! eg1: INNER WINDOW ON 4KSHIELD = 20 MICRONS OF ALUM 
      FOURKSLDRL1 =0.000157d+0  ! eg4: == 0.0014/(24.01/2.6999)    Alum - 14 microns - table A-IIA.) 
c     HUNKSLDRL1  =0.00022        ! eg1: INNER PART OF 100KSHIELD = 20 MICRONS OF ALUM (100K shield replaced by 77K shield; 14 microns; Table A-II)
      SEV7KSLDRL1 =0.000157d+0  ! eg4: == 0.0014/(24.01/2.6999)    Alum - 14 microns - table A-II.) 
c     BANJOENTRL  =0.0008         ! eg1: banjo entrance window = 71 microns Alum (same - table A-III.)
      BANJOENTRL =0.000798d+0   ! eg4: == 0.0008 as before or x/(X0*rho)=0.0071/(24.01/2.6999) 
c     HELIUMENTRL =0.000384       ! eg1: includes (kp: only) 0.25 cm of He4 before target cell entrance window (He4 inside target chamber is considered differently)
c     HELIUMENTRL = 0.0008763   ! eg4: 0.0008763= 0.57/(94.32/0.145)  or         L_He=(2.14-1.0)/2=0.57cm  (2.14cm Banjo, 1.0cm or 0.5cm target cell)
c     HELIUMENTRL = 0.0012606   ! eg4: = 0.82/(94.32/0.145)          L_He=(2.14-0.5)/2=0.82cm  (2.14cm Banjo, 0.5cm target cell)
      HELIUMENTRL = ((BANJOL - 
     >     CHAMBERL)/2)/(LHEXO/LHED) ! eg4:  Same as HELIUMEXITRL assuming the Banjo & target center coincide.
c     TARGETENTRL =0.0000124      ! eg1: entrance window for target cell = 11  microns Alum  
c     = 0.0000882 ???      == (0.0025/(40.58/1.42)    (kapton: 25 microns, - M. Zarecky table B-C.)
      TARGETENTRL  = 0.0001155d+0 ! eg4: == (0.0033/(40.58/1.42)   (kapton: 33 microns C. Keith).)
c     TARGETEXITRL =0.000158      ! eg1: kapton exit window for target cell = 45 microns (Kapton:  25 microns, - table B-C. (33 microns C. Keith))
      TARGETENTRL  = 0.0001155d+0 ! eg4: == (0.0033/(40.58/1.42)   (kapton: 33 microns C. Keith).)
c     HELIUMEXITRL =0.000384      ! eg1: INCLUDES 0.25 CM OF HE4 AFTER TARGET CELL EXIT WINDOW 
c     Unlike HELIUMENTRL Shouldn't it depend on the scattering point and the angles? Yes, but, we use this constant with those factors later on.
      HELIUMEXITRL = ((BANJOL - 
     >     CHAMBERL)/2)/(LHEXO/LHED) ! eg4:  Same as HELIUMENTRL assuming the Banjo & target center coincide.
c     BANJOEXITRL =0.0008         ! eg1: BANJO EXIT WINDOW = 71 MICRONS OF ALUM (same - table A-VI.)
      BANJOEXITRL =0.000798d+0  ! eg4: == 0.0008 as before or x/(X0*rho)=0.0071/(24.01/2.6999)  
c     FOURKSLDRL2 =0.0012         ! eg1: INNER AND OUTER OVERLAP + 0.04mm ALUM/MYLAR TAPE  (?????)
c     FOURKSLDRL3 =0.00084        ! eg1: OUTER PART OF 4KSHIELD = 75 MICRONS OF ALUM (?????)
      FOURKSLDRL2 = 0.0000672d+0 ! eg4: there is a 6 microns Aluminium very close to superinsulation1 (so its is combined in SUPERINSRL1 for simplicity)
c     HUNKSLDRL2  =0.00136        ! eg1: INNER AND OUTER OVERLAP + 0.04mm OF MYLAR TAPE (?????)
c     HUNKSLDRL3  =0.001          ! eg1: OUTER PART OF 100KSHIELD = 3X30 MICRONS (?????)
      SEV7KSLDRL2 = 0.000281d+0 !eg4: = 25 microns, Aluminium (Table A-VIB or B-J = 0.0025/(24.01/2.6999)
      SI1LAYER = 0.0000881d+0   !eg4: = 0.00088/39.95 + 3*0.001/45.4 => radLength of one layer of super-insulation (= 1 layer of Alum/mylar and 3 ply of Cerex)
c     SUPERINSRL1 = 0.000001         ! eg1: INNER RADIUS IS JUST ONE SHEET OF ALUM/MYLAR   
      SUPERINSRL1 = 10*SI1LAYER + 0.001/45.4 ! eg4: (10 layers: Table A-VII) The last terms due to 1 PLY upstream 4K shield (Table B-K)
c     SUPERINSRL2 = 0.0000362        ! eg1: OUTER PART IS 13 SHEETS OF ALUM/MYLAR & CEREX INSULATION 
      SUPERINSRL2 = 20*SI1LAYER ! eg4: (20 layers: Table A-VIII)  =0.0000362*20/13 ?????
c     OVCEXITRL   = 0.0031          ! eg1: OVC EXIT WINDOW = 280 MICRONS OF ALUM (???) 
      OVCEXITRL = 0.0031d+0     ! eg4: ( same??? Table B-M?)  






C     PACKING FRACTIONS OF LID and NH3 TARGET MATERIALS (ratio < 1)

      IF (I .EQ. 2) THEN
         IF (Eb .EQ. 1) THEN
            FND3=0.624d+0       !Preliminary PACKING FRACTION OF ND3 (from S. Phillips) http://clasweb.jlab.org/rungroups/eg4/wiki/index.php/October_28%2C_2011
         ELSE IF (Eb .EQ. 2) THEN
            FND3=0.764d+0       !PACKING FRACTION OF ND3 ??????
         ELSE
C     ERROR MESSEGE FOR WRONG FLAG VALUE FOR Eb
            WRITE(*,*) 'ERROR:INVALID FLAG (Eb=1or2 only)' 
            GO TO 200 
         ENDIF
      ELSE 
C     ERROR MESSEGE FOR WRONG FLAG VALUE FOR I
         WRITE(*,*) 'ERROR:INVALID FLAG (I=2 only)' 
         GO TO 200 
      ENDIF

c     FNH3=.58         !PACKINGH3 ????





C     CREATE AND FILL ARRAY FOR SCATTERED PARTICLE LOCATION
      
      IF (I.EQ.2) THEN
         RPART(1,1) = TARGETEXITZ !Z dim of target exit window
      ENDIF
      RPART(2,1) = BANJOEXITZ   !Z dim of banjo exit window
      RPART(3,1) = SUPERINSZ1   !Z dim of 10-layered superinsulation (and also of dwonstream 4k shield which is one PLY of Alum/Mylar)
      RPART(4,1) = SEV7KSLDZ2   !Z dim of upstream 77K shield
      RPART(5,1) = SUPERINSZ2   !Z dim of 20-layered super insulation
      RPART(6,1) = OVCEXITZ     !Z dim of OVC exit window

      N=1
      DO WHILE (N.LT.7)
c     RPART(N,2) = RADXY(THA,PHIA,X,Y,Zv,RPART(N,1))
         RPART(N,2) = RADXY_EG4(THA,PHIA,X,Y,Zv,RPART(N,1))
c     print *,tha,phia,x,y,z
c     print *,RPART(N,1)
         N=N+1
      END DO
      
c     print *,rpart

C     INITIALIZE OUTPUTS

      RADA=0.0d+0
      RADB=0.0d+0
      
      RADA_4DBG=0.0d+0          !kp: 2/8/2012
      RADA_4DBG1=0.0d+0         !kp: 2/8/2012
      RADA_4DBG2=0.0d+0         !kp: 2/8/2012
      RADA_4DBG3=0.0d+0         !kp: 2/8/2012
      RADA_4DBG4=0.0d+0         !kp: 2/8/2012

      RADA_W1=0.0d+0            !kp: 2/13/2012, 
      RADA_W2=0.0d+0            !kp: 2/13/2012, 
      RADA_W3=0.0d+0            !kp: 2/13/2012, 
      RADA_W4=0.0d+0            !kp: 2/13/2012, 
      RADA_W5=0.0d+0            !kp: 2/13/2012, 
      RADA_W6=0.0d+0            !kp: 2/13/2012


C     TESTS TO SEE IF GIVEN SCATTERING POINT LIES WITHIN TARGET
C     NEED FIND OUT IF Zv=0 IS ACTUALLY AT BACK OF TARGET OR CENTER
C     OF NH3 TARGET (RIGHT NOW I AM ASSUMING CENTER IS Zv=-55)
      

      R=(X**2+Y**2)**0.5d+0
      
!kp: Following print lines just for debugging purpose
c     print*,'Radius R: ',R
c     print*,'z,zMax:',R,Zv,TARGETCNTRZ+CHAMBERL/2
c     print*,'zCtr,chL,chL/2: ',TARGETCNTRZ,
c     <       CHAMBERL,CHAMBERL/2


      IF (R .GT. CHAMBERR) THEN
         WRITE(*,*)  'ERROR:SCATTERING POINT OUTSIDE TARGET[1,2,4]'
         GO TO 200
      ENDIF
      IF (Zv .GT. (TARGETCNTRZ+CHAMBERL/2)) THEN
         WRITE(*,*)  'ERROR:SCATTERING POINT OUTSIDE TARGET[1,2,4]'
         GO TO 200
      ENDIF
      IF (Zv .LT. (TARGETCNTRZ-CHAMBERL/2)) THEN
         WRITE(*,*)  'ERROR:SCATTERING POINT OUTSIDE TARGET[1,2,4]'
         GO TO 200
      ENDIF

      



c     =========RADB contribution from different components============

C     RADB CONTRIBUTION FROM BEAM ENTRANCE WINDOW IN OVC,HEAT SHIELD WINDOWS (4k & 77k), BANJO ENTRANCE WINDOW,
C     LIQUID HELIUM (between Banjo and target entrance windows) AND THE TARGET ENTRANCE WINDOWS.

      RADB=BEAMENTERRL+FOURKSLDRL1
     >     +SEV7KSLDRL1+BANJOENTRL
     <     +HELIUMENTRL+TARGETENTRL
      
C     LT IS THE DISTANCE THE ELECTRON TRAVERSED IN THE TARGET BEFORE SCATTERING
      
      LT = ABS(Zv-TARGETENTRZ)

C     FOLLOWING RADB FROM ND3+HELIUM CONTRIBUTION
      IF (I .EQ. 2) THEN
         RADB=RADB+LT*FND3*ND3D/ND3XO
     >        +LT*(1-FND3)*LHED/LHEXO

C     ERROR MESSEGE FOR WRONG FLAG VALUE FOR I

      ELSE 
c     WRITE(*,*) 'ERROR:INVALID FLAG (I=1to4 only)' 
         WRITE(*,*) 'ERROR:INVALID FLAG (I=2 only)' 
         GO TO 200 
      ENDIF





C     ========== Determine the radiation length of target and helium in target per unit of target length ============

c     In fact, here we calculate the inverse (rho/X0) of radiation length (X0/rho in cm units), 
c     which when multiplied by the material thickness gives us the radiation exponent (= thickness/rad-length) 
c     which is a dimensionless quantity. Remember that RADB and RADA both are the total of all the radiation 
c     exponents before and after the scattering and therefore they are dimesionless (Remember that X must be a 
c     dimesionless number in a expression such as exp(X), sin(X), tan(X) etc).

      IF (I .EQ. 2) THEN
         POLTARGETRL=CHAMBERL*FND3*(ND3D/ND3XO) !Contribution only due to the ammonia beads
     >        + CHAMBERL*(1.0d+0-FND3)*(LHED/LHEXO) !Contribution due to the Helium inside the target chamber
      ENDIF






c     ============== Start Calculating the RADA with target exit window (includes target and helium inside the chamber) ===============
c     DZI=0.0                 !kp: Zc - Zv with Zc=Zv-coordinate of the track intersection (Ic) with the chamber cylinder
c     TERM1=0.0
      TERM2C=0.0d+0
      TERM2W=0.0d+0
      TERM2N=0.0d+0
      RSPHI=0.0d+0
      TCC=0.0d+0
      TWW=0.0d+0
      TNN=0.0d+0
c     PLUSMINUS=0.0         !kp: Must be given +1 or -1 value to it. Otherwise, it will fail it's purpose.
      COTTH=1.0/TAN(THA)
      SINTH=SIN(THA)

      
      TERM1=-( X*COS(PHIA)
     <     + Y*SIN(PHIA) )
c     RSPHI=X*Y*SIN(2*PHIA)   !kp: Old wrong line
      RSPHI=( X*SIN(PHIA) 
     <     -Y*COS(PHIA))**2
c     DZI=COTTH*(TERM1        !kp: Zc - Zv with Zc=Zv-coordinate of the track intersection (Ic) with the chamber cylinder
c     <       + (CHAMBERR**2 
c     <       - RSPHI)**0.5 )

c     IF(DZI.LT.0.0) THEN   !kp: PLUSMINUS is to select + or - sign for the squared-root above in DZI
c     PLUSMINUS=-1
c     ELSE
c     PLUSMINUS=+1
c     ENDIF
c     DIST_SIC(SINTH,TERM1,RSPHI,RR,PLUSMINUS)
      TCC=DIST_SIC(THA,TERM1
     >     ,RSPHI,CHAMBERR)     !kp: Track length from S to the intersection on the chamber cylinder
      TWW=DIST_SIC(THA,TERM1
     >     ,RSPHI,CHAMBERR+WALLT) !kp: Track length from S to the intersection on the wall cylinder
      TNN=DIST_SIC(THA,TERM1,RSPHI
     >     ,CHAMBERR+WALLT+NMRT) !kp: Track length from S to the intersection on the NMR cylinder
c     LIN=TCC                 !kp: Track length through chamber content
c     LWALL=TWW-TCC           !kp: Track length through chamber wall
c     LNMR=TNN-TWW            !kp: Track length through NMR layer
c     LOUT=HYPOT-TNN          !kp: Track length through Helium outside the chamber
      
c     ============= Above part added on 2/14/2012 ==========================================
C     R IS THE RADIUS FROM THE X,Y=0  POINT IN X-Y PLANE 
C     OF THE SCATTERED PARTICLE. 

C     TRKFACTOR below is the conversion factor that converts the Z-component of a given length into a corresponding
C     length along the particle track (by multiplying Z-length by TRKFACTOR we get the length along the path)
C     The alternative to the old factor 1.0//ABS((COS(THA)*COS(PHIA))) which seemed to use small polar angle (theta) approximation
c     TRKFACTOR = sqrt(1.0/(COS(THA)**2) + tan(PHIA)**2)     !This is the case when THA and PHIA are the horizontal and vertical angles
      TRKFACTOR = 1.0d+0/cos(THA)
c     TRKFACTOR = 1.0/ABS((COS(THA)*COS(PHIA)))                !Old value for the case when THA and PHIA are the horizontal and vertical angles

      R = (X**2 + Y**2)**0.5d+0 !Radial distance of beam from Z-axis
c     RTRI = RPART(1,2)-R                            !Old line
c     RTRI =  ((TARGETEXITZ-Zv)*tan(THA)
      HYPOT = (TARGETEXITZ-Zv)/cos(THA) !Hypotenuse of the triangle made by the track, beam path and the foil plane
c     THETATRI = ASIN(RTRI/HYPOT)                    !Old line
      THETATRI = THA
      ABSZH = ABS(TARGETEXITZ-Zv) !kp Horizontal Z-lenght of track inside the target


c     Comment out the following CALL & print lines while not debugging
c     CALL DBGXYZ(X,Y,Zv,THA,PHIA,TCC)
c     CALL DBGXYZ(X,Y,Zv,THA,PHIA,TWW)
c     CALL DBGXYZ(X,Y,Zv,THA,PHIA,TNN)
c     print*, 'radperL4CWN:',POLTARGETRL
c     <       ,PCTFED/PCTFEXO,NMRD/NMRXO
c     c        print*, 'Z, R_perp_TgtExit: ',RPART(1,1),RPART(1,2)
c     print*, 'HYPOT, R_perp_TgtExit: ',HYPOT,RPART(1,2)
c     c        print*, 'rC,tWall,tNMR:',CHAMBERR,WALLT,NMRT
c     print*, '3Rs:',CHAMBERR,CHAMBERR+WALLT
c     <       ,CHAMBERR+WALLT+NMRT
      
c     kp: Moved the following block down because it seemed that HYPOT didn't have a value (distorting the rada in 4th or ELSE block below
      LIN=TCC                   !kp: Track length through chamber content
      LWALL=TWW-TCC             !kp: Track length through chamber wall
      LNMR=TNN-TWW              !kp: Track length through NMR layer
      LOUT=HYPOT-TNN            !kp: Track length through Helium outside the chamber
      





c     IF (I.NE.3) THEN
      IF (RPART(1,2).LT.CHAMBERR) THEN !i.e., when the track passed through the exit window.
c     >    RADA = RADA   + POLTARGETRL*ABS(TARGETEXITZ-Zv)/ABS((COS(THA)*COS(PHIA)))   !Old line           
         RADA = RADA
     >        + POLTARGETRL*HYPOT
c     >          + POLTARGETRL*ABSZH*TRKFACTOR      !Contribution from the content inside the chamber
     >        + TARGETEXITRL*TRKFACTOR !contribution from the chamber wall
c     print *,trkfactor,'ra',rada            !kp
c     >          ,poltargetrl,abs(targetexitz-z)
         RADA_4DBG1=RADA
      ELSE IF (RPART(1,2).LT.(CHAMBERR
     >        + WALLT)) THEN    !i.e., when the track passed through the PCTFED or Kel-F wall & exited out the kapton cap).
         LWALL=HYPOT-TCC
         RADA = RADA + POLTARGETRL*LIN 
     >        + (PCTFED/PCTFEXO)*LWALL 
     >        + TARGETEXITRL*TRKFACTOR !Assumes that the exit window extends upto the outer wall surface
c     print *, '2: rada = ',rada !kp
         RADA_4DBG2=RADA
      ELSE IF (RPART(1,2).LT.(CHAMBERR  
     >        + WALLT + NMRT)) THEN !i.e., when the track passes through cell wall & NMR layer & exits out the side face. 
         LNMR=HYPOT-TWW
         RADA = RADA + POLTARGETRL*LIN 
     >        + (PCTFED/PCTFEXO)*LWALL 
     >        + (NMRD/NMRXO)*LNMR
c     print *, '3: rada = ',rada !kp
         RADA_4DBG3=RADA
      ELSE 
         RADA = RADA + POLTARGETRL*LIN 
     >        + (PCTFED/PCTFEXO)*LWALL
     >        + (NMRD/NMRXO)*LNMR 
     >        + (LHED/LHEXO)*LOUT
c     print *, LOUT, LMNR, LWALL
         RADA_4DBG4=RADA
      ENDIF
c     ENDIF





      RADA_4DBG = RADA
c     RADA = 0.0              !for Debug
c     RADA = RADA_4DBG4
      RADA_W1 =  RADA
c     print *, 'CLAS is ',rada






c     ============== Add Helium and Banjo window contribution to RADA ===============
c     IF (I .NE. 3) THEN                 Not-Carbon target 
      IF (RPART(2,2).GT.BANJOEXITR) THEN !if qv  is greater than radius/range of the exit
c     WRITE(*,*) 'PARTICLE OUTSDIE OF BANJO EXIT WINDOW RANGE'
         RADA=-1.0
         GO TO 200
      ELSE
c     R = (X**2+Y**2)**0.5      !This line already exists above at the start of RADA calculation
c     HYPOT = ((RPART(2,2)-R)**2+
c     >        (BANJOEXITZ-Zv)**2)**0.5 !Hypotenuse of the triangle made by the track, beam path and the Banjo-exit plane
         HYPOT=(BANJOEXITZ-Zv)*TRKFACTOR !kp: 2/14/12
c     LIN = ((RPART(1,2)-R)**2+
c     >        (TARGETEXITZ-Zv)**2)**0.5 !Hypotenuse of the triangle made by the track, beam path and the target-exit plane
         LIN=(TARGETEXITZ-Zv)*TRKFACTOR !kp: 2/14/12YPOT=(BANJOEXITZ-Zv)*TRKFACTOR1 !kp: 2/14/12
         LOUT = HYPOT - LIN     !Track length between the planes of target-exit and Banjo-exit (filled with 4He)
         RADA = RADA + BANJOEXITRL*TRKFACTOR !Contribution due to Banjo exit
     >        + LOUT*LHED/LHEXO !Contribution due to the 4He in between the target and Banjo exit windows.
c     Sum of the contributions from RADA above, & the lengths qr and Sq respectively . 
      ENDIF
c     ENDIF


c     RADA = 0.0
      RADA_W2=RADA-RADA_W1






c     ============== Add Super-insulation-1 (4K shield included as it has the same Z) contribution to RADA ===============

      IF (RPART(3,2) .GT. SUPERINSR1) THEN
c     WRITE(*,*) 'PARTICLE OUTSDIE OF SUPERINSULATION-1 RANGE'
         RADA=-2.0
         GO TO 200
      ELSE
         IF (RPART(3,2) .LT. SUPERINSIR1) THEN
c     RADA = RADA + SUPERINSIR1/ABS(COS(THA)*COS(PHIA))
	    RADA = RADA + 
     >           (FOURKSLDRL2+SI1LAYER)*TRKFACTOR !Just due to the 6 micron (Al) down 4k-shield and  1 layer of Sup-ins-1 that covers the beam-hole
         ELSE 
	    RADA = RADA + 
     >           (FOURKSLDRL2+SUPERINSRL1)*TRKFACTOR !Due to upstream 4k-shield & all 10 layers of Sup-ins-1
         ENDIF
      ENDIF


c     RADA = 0.0                !for Debug
      RADA_W3=RADA-RADA_W1-RADA_W2





c     ============== Add 77K shield  contribution to RADA ===============

      IF (RPART(4,2) .GT. SEV7KSLDR2) THEN
c     WRITE(*,*) 'PARTICLE OUTSDIE OF the upstream 77k Shield RANGE'
         RADA=-3.0
	 GO TO 200
      ELSE
         RADA = RADA + SEV7KSLDRL2*TRKFACTOR  
      ENDIF


c     RADA = 0.0                !for Debug
      RADA_W4=RADA-RADA_W1-RADA_W2-RADA_W3




c     ============== Add Super-insulation-2  contribution to RADA ===============

      IF (RPART(5,2) .GT. SUPERINSR2) THEN
c     WRITE(*,*) 'PARTICLE OUTSDIE OF SUPERINSULATION-2 RANGE'
         RADA=-4.0
         GO TO 200
      ELSE
         IF (RPART(5,2) .LT. SUPERINSIR2) THEN
c     RADA = RADA + SUPERINSIR1/ABS(COS(THA)*COS(PHIA))
	    RADA = RADA + SI1LAYER*TRKFACTOR !Just due to the 1 layer of Sup-ins-2 that covers the beam-hole
         ELSE 
	    RADA = RADA + SUPERINSRL2*TRKFACTOR !Due to all 20 layers of Sup-ins-2
         ENDIF
      ENDIF


c     RADA = 0.0                !for Debug
      RADA_W5=RADA-RADA_W1-RADA_W2-RADA_W3
     >     -RADA_W4





c     ============== Add OVC (==Outer Vacuum Can/Vessel) window contribution to RADA ===============

      IF (RPART(6,2).GT.OVCEXITR) THEN
c     WRITE(*,*) 'PARTICLE OUTSDIE OF OVC WINDOW RANGE'
         RADA=-5.0
         GO TO 200
      ELSE
         ROVC = 93.98           !Radius of curvature of OVC exit window (93.98 cm == 37 inches (Source Dave Kashy))
         ZOVC = ABS(OVCEXITZ - Zv) !Distance of OVC exit from the target center
         SIN2TH = sin(THA)**2
         COS2TH = cos(THA)**2
c     YOVC           !radial distance of the track hit point on the OVC from the Z-axis
c     VOVC below is Z-distance of the track hit point from the vertical plane (perpendicular to Z-axis) that stands at Z=OVCEXITZ i.e., 
c     VOVC is the ABS(Z_hit - OVCEXITZ) where Z_hit is the Z-coordinate of the point where the track hits the curved OVC exit window
         VOVC = (ROVC*COS2TH - ZOVC*SIN2TH) -
     >        sqrt((ROVC*COS2TH - ZOVC*SIN2TH)**2 
     >        - (ZOVC**2)*SIN2TH)
         THETAOVC = asin((ROVC-VOVC)/ROVC) ! angle of tilt of the curvature (tangent plane at the hit point) with the Z-axis
	 RADA = RADA + 
     >        (OVCEXITRL/sin(THETAOVC))*TRKFACTOR !Effective thickness along Z-axis at the hit point = OVCEXITRL/cos(THETAOVC)
      ENDIF


c     kp: For Debug
      RADA_W6=RADA-RADA_W1-RADA_W2-RADA_W3
     >     -RADA_W4-RADA_W5




C     ADD PATHLENGTH THROUGH REGION TWO OF CLAS
C     THIS WAS CALCULATED BY ASSUMING ONE PATHLENGTH
C     THROUGH CLAS = 2.5 METERS.  O.5 OF THAT IS THROUGH
C     THE DRIFT CHAMBERS WHICH ARE COMPOSED OF 90% ARGON AND
C     10% C02.  THE OTHER 2M IS ASSUMED TO BE AIR.

c     print *, 'CLAS is ',rada

      RADA = RADA + 0.01143d+0
c     RADA = RADA - RADA_W1
c     RADA = RADA_4DBG4
c     RADA =  RADA_W1           !kp: For Debug 2/13/2012
c     print *, 'CLAS is ',rada

 200  END




      REAL*8 FUNCTION DIST_SIC(THA,TERM1,RSPHI,RR) !kp Distance between S(X,Y,Zv) & Ic(track intersection with cylinder of RR)
      IMPLICIT NONE
      REAL*8 THA,TERM1,RSPHI,RR,PLUSMINUS !,DIST_SIC
      REAL*8 DZI,SINTH,COTTH,TSQRT
c     print *,'Cylinder R: ',RR
      SINTH=SIN(THA)
      COTTH=1.0/TAN(THA)
      TSQRT=(RR**2 - RSPHI)**0.5d+0
      DZI=COTTH*(TERM1+TSQRT)   !kp: Zc - Zv with Zc=Zv-coordinate of the track intersection (Ic) with the chamber cylinder
      IF(DZI.LT.0.0d+0) THEN    !kp: PLUSMINUS is to select + or - sign for the squared-root above in DZI
         PLUSMINUS=-1
      ELSE
         PLUSMINUS=+1
      ENDIF
      
c     DIST_SIC=ABS((TERM1+PLUSMINUS*TSQRT)/SINTH)
      DIST_SIC=(TERM1 + PLUSMINUS*TSQRT)/SINTH
c     print*,'DIST_SIC+',DIST_SIC
c     DIST_SIC=(TERM1 - PLUSMINUS*TSQRT)/SINTH
c     print*,'DIST_SIC-',DIST_SIC
      RETURN
      END FUNCTION

      

      SUBROUTINE DBGXYZ(X,Y,Zv,TH,PHI,TT) !kp This prints X,Y,Zv of the intersection point (for Debug purpose only) with a given cylinder
      IMPLICIT NONE
      REAL*8 TH,PHI,X,Y,Zv,TT,XI,YI,ZI
      XI=X+TT*SIN(TH)*COS(PHI)
      YI=Y+TT*SIN(TH)*SIN(PHI)
      ZI=Zv+TT*COS(TH)
      print *,'Zv, Xi,Yi,Zi,TT: ',Zv,XI,YI,ZI,TT
      END




      

C     Function Definition: REAL*8 FUNCTION RADXY_EG4(THA,PHIA,X,Y,Zv,ZPLANE)
C     RADXY IS DEFINED AS THE DISTANCE OF THE SCATTERED PARTICLE FROM THE XY ORIGIN IN THE LAYER SPECIFIED IN THE XY PLANE.

      REAL*8 FUNCTION RADXY_EG4(THA,PHIA,X,Y,Zv,ZPLANE) !kp The new version (replaced the old one (below, renamed) with this one 2/8/2012)
      IMPLICIT NONE
******************************This what was before ******
*     REAL*8 THA,PHIA,ZPLANE,Zv,X,Y,XPLANE,YPLANE
*     XPLANE=(((ZPLANE-Zv)*TAN(THA))+ X)       ! X coordinate of the track in the given plane of the foil?
*     YPLANE=(((ZPLANE-Zv)*TAN(PHIA))+ Y)      ! Y coordinate of the track in the plane of the foil?
*     RADXY=SQRT((XPLANE**2)+(YPLANE**2))     ! the radial distance of the track from the (X=0,Y=0) point in the plane/foil?
*********************************************************
c     If the beam line is at (X=0,Y=0), is RADXY the distance of the track from the beam line in the plane/foil?
c     http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/SimStuffs/RadCorStuff/radialDistInFoil.svg

      REAL*8 THA,PHIA,XPLANE,YPLANE,ZPLANE,Zv,X,Y,RADXY
c     REAL*8 RADXY_EG4
c     RADXY = ABS(ZPLANE-Zv)*TAN(THA) !Distance between points B and T
c     print *, 'radxy ',radxy,tha !kp: 2/8/2012
      XPLANE=X+(ZPLANE-Zv)*TAN(THA)*COS(PHIA)
      YPLANE=Y+(ZPLANE-Zv)*TAN(THA)*SIN(PHIA)
c     RADXY=SQRT((XPLANE**2)+(YPLANE**2))
      RADXY=((XPLANE**2)+(YPLANE**2))**0.5d+0
      RADXY_EG4=RADXY
      RETURN
      END FUNCTION



C     Function Definition: REAL*8 FUNCTION RADXY(THA,PHIA,X,Y,Zv,ZPLANE)
C     RADXY IS DEFINED AS THE DISTANCE OF THE SCATTERED PARTICLE FROM THE XY ORIGIN IN THE LAYER SPECIFIED IN THE XY PLANE.

      REAL*8 FUNCTION RADXY_OLD2(THA,PHIA,X,Y,Zv,ZPLANE) !kp The new version (replaced the old one (below, renamed) with this one 2/8/2012)
      IMPLICIT NONE
      REAL*8 THA,PHIA,ZPLANE,Zv,X,Y !,RADXY_OLD2
      RADXY_OLD2 = ABS(ZPLANE-Zv)*TAN(THA) !Distance between points B and T
c     print *, 'radxy ',radxy,tha !kp: 2/8/2012
      RETURN
      END FUNCTION









      REAL*8 FUNCTION RADXY_OLD1(THA,PHIA,X,Y,Zv,ZPLANE) !kp The old version (replaced this with above one 2/8/2012)
      IMPLICIT NONE
      REAL*8 THA,PHIA,ZPLANE,Zv,X,Y,PHIORIGIN, DPHI, PI
      REAL*8 DISTOB, DISTBT, RADXY !, RADXY_OLD1
      PI = acos(-1.0d+0)
      PHIORIGIN = atan2(Y,X)
c     if the phi of the scattering point is -ve, then add 2*pi to make it always fall between (0,2pi) rather than (pi,-pi)
      IF (PHIA.LT.0.0d+0) THEN   
         PHIA = PHIA + 2*PI
      ENDIF 
c     if PHIA is -ve, then add 2*pi to make it always fall between (0,2pi) rather than (pi,-pi)
      IF (PHIORIGIN.LT.0.0d+0) THEN   
         PHIORIGIN = PHIORIGIN + 2*PI
      ENDIF   
      
c     Above two angles are the angles made by following two radius vectors with the X-axis in the XY plane of the foil/layer specified
c     First is the one joining the origin (O) in this plane to the projected beam spot (let's call it point B) on the same plane.
c     Second is the one joining the origin (O) in this plane to the point (let's call it T) where the track meets the plane.
      DPHI = ABS(PHIA-PHIORIGIN) ! This is the angle between above two radius vectors.

c     Let's call the scattering point S(X,Y,Zv), then ST makes an angle of THA with SB (assuming SB is parallel to Z-axis)
      DISTOB = sqrt(X**2 + Y**2) !Distance between points O and B
      DISTBT = (ZPLANE-Zv)*TAN(THA) !Distance between points B and T
      RADXY = DISTOB*cos(DPHI) 
     >     + sqrt(DISTBT**2 - 
     >     (DISTOB*sin(DPHI))**2)
      RADXY_OLD1 = RADXY
c     print *, 'radxy ',radxy,dphi,tha !kp: 2/8/2012
c     print *, 's:', DISTBT**2 - (DISTOB*sin(DPHI))**2
      RETURN
      END FUNCTION



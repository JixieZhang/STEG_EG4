      subroutine SET_THINGS_UP(Ebeam,TMOD)
*******************************************************************************
*    This is a subroutine adapted from the rcslacpol program (some useless stuff 
*          may still be in the file that can be removed)
*    
*     This routine is supposed to initialize the values of the global/common variables
*             before the kinematic loop to calculate the cross-sections starts.
*
*    This routine will be called from STEG
*
*    Don't forget to update the value of the size (PP_NMAX) of the arrays (common 
*        variables) defined in prinplot.inc according to the # of bins in cs-map 
*        to be produced.
******************************************************************************* 7/30/13
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
       INCLUDE 'radcon.inc'
       INCLUDE 'radvarfix.inc'
       INCLUDE 'instruct.inc'
       INCLUDE 'prinplot.inc'

       include 'binning.inc'    !kp: from cs_init.f
       
       INCLUDE 'kpaVarChanges.inc' !kpa: 7/30/13  (to control change in A1)
c      REAL*8 changeA1 

c     =============kp: Variable declarations from cs_init.f =======================
      real Ebeam,cth_el_min,cth_el_max,
     &     p_el_min,p_el_max,cs_max,cs_int,sig
      real theta_el,p_el,t_current
      real v_p_el(0:n_p_el),v_th_el(0:n_th_el)
      real SIMPS,sgm_model,degrad,ppi,delta_phi
      integer i_th_el,i_p_el,ithe,ipe
      logical flag_acc
      data degrad/0.0174532925e0/,ppi/3.14159265e0/   ! ppi used: "Symbol 'pi' at (1) already has basic type of REAL"
c      common/inpt/ithe,ipe
c     ============ Variables needed to add the stuff to be done by read_map program of STEG =============
      real*8 csmx(0:n_th_el,0:n_p_el),fmx
      real*8 ctha,cthb,pa,pb,acc,acc_tot,norm
      integer STAT

c     ===== ===== My own Additions =====
       INTEGER EB_INDEX         !Jixie: EB_INDEX will be signed according to input beam erngy and target type
       REAL*8 CSRAD             !kp: radiated cross-section to be returned by this function.
       REAL*8 EE,EEP,THETADD,XXCUT,SCAT_PHIRR  !kp: added because of "Error: COMMON attribute conflicts with DUMMY attribute in 'e' at (1)", see radvarfix.inc
       INTEGER DBG_MY_KINE



c     ===========kp: Variables from RCSLACPOL =====
       LOGICAL EXTERNALv, HiPres
       INTEGER IP, NPTS, I, IPLOT, IQ, IW, IT, ITAIL, NEXT,
     >      LEN1, LEN2, LEN3, IS, TMOD
       REAL*8 ENERGY, Z, EINT, SANSWER, TA, Q2OUT, ERC, 
     >      EPRC, THRC, ANSWER, SIGMAR, SIGMARNEW, MULT, 
     >      XCUTRC, W2BIN, TB

       CHARACTER*40 InputFile, InFile, OutFile
       CHARACTER*40 FILENAME1, FILENAME2, FILENAME3
       CHARACTER*8 TARGETS(5)/'PROTON  ','3HELIUM ',
     >      'ND3     ','LiD     ','NEUTRON '/
       CHARACTER*6 EXPERIMENT
       CHARACTER*60 TITLE
       CHARACTER*9 DATENOW/'   040315'/
       CHARACTER*8 TIMENOW, CLOCK_, DATE
       CHARACTER*80 FileText(30)
       EXTERNAL DATE
cFRW   EXTERNAL CLOCK_
       COMMON /KINLIST/ERC,EPRC,THRC,EINT,Z,XCUTRC
       COMMON /INFO/EXTERNALv,NEXT


       common/inpt/ithe,ipe     !cs_init.f


       EXTERNALv = .TRUE.
c       EXTERNALv = .FALSE.          !kp: 9/7/12 executable steg4DirectXSplotsPol_ExtRadOff made with this line
       FL_POL = .TRUE.              !kp: If T do "polarized" analysis.  (instruct.inc)
       FL_UNPOL = .TRUE.            !kp: If T do "unpolarized" analysis.
       
cfrw   read command line flag for output format
c       HiPres = .true.          !     write(6,*) ' --- writing Brief, High Precision output ---'
       HiPres = .false.         ! use default if none specified







cfrw   get input file name from command line
       InFile = 'rcslacpol.file'
       OPEN (UNIT=27,FILE=InFile, STATUS='OLD', err=9911)
       print*,'Finished reading parameter file rcslacpol.file ...'



c     kp: What are the following things and where is EXTPEAKING declaration? ============================ 
       read(27,*) INTPEAKING     !kp: in 'instruct.inc' -  LOGICAL INTPEAKING:  ! Peaking internal (otherwise "exact")
       read(27,*) EXTPEAKING    !kp: 'grep' couldn't find it's declaration
       READ(27,'(1X,A3)') TARG            !kp: reading 1 column (1X) of character (ASCII) info? (example - NH3)
       READ(27,'(1X,D5.3)') PFERMI        !kp: reading 1 column (1X) of Double precision #?   (example - 0.000)
       READ(27,'(1X,A4)') POLTYPE         !kp: Look at rcslacpol.file 
       READ(27,'(1X,A6)') EXPERIMENT
       READ(27,'(1X,I2)') IFFMOD          !kp: reading 1 column (1X) of an integer in?
       READ(27,'(1X,I2)') IPOL
       READ(27,'(1X,I2)') IPOLRES
       READ(27,'(1X,I2)') IPAULI
       READ(27,'(1X,I2)') IA1         
       READ(27,'(1X,I2)') IEXTERNAL   !c        READ(27,'(1X,I2,//////////)') AsymChoice
       READ(27,'(1X,I2,/////////)') AsymChoice
       READ(27,'(1X,I2,//////)') SFchoice
       print*, AsymChoice, SFchoice
       READ(27,'(A40)') FILENAME1                !kp: kine.dat
       READ(27,'(A40)') FILENAME2                !kp: test.out
       READ(27,'(A40,/)') FILENAME3              !kp: test.top
c       print*, FILENAME1, FILENAME2, FILENAME3  !C       READ(27,*) DYNAMIC_MODEL                 !kp: !Model flags ? 
       READ(27,*) NMC
       READ(27,*) NORES
       READ(27,*) BODEK1
       READ(27,*) ERROR
       READ(27,*) TARGNEG
       READ(27,*) TARGPOS
       READ(27,*) DELTEST
       READ(27,*) ZTEST
       READ(27,*) MULTISOFT
c     kp: TS_ALPHA & TS_BETA defined in radvarfix.inc (look at this file, it has more things in it) as follows:
c     TS_ALPHA,                    ! target spin vector's polar angle (rad)
c     TS_BETA,                     ! target spin vector's azimuthal angle (rad)


******************************************************************************* 7/30/13
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
       AsymChoice=AsymChoiceKpa
       SFchoice=SFchoiceKpa
****************************************************************************** 7/30/13



       IF(POLTYPE.EQ.'ARBI') THEN !kp: In the file rcslacpol.file (also see radvarfix.inc), POLTYPE == LONG (meaning longitudinal polzn)
          READ(27,*) TS_ALPHA    !kp: polar angle of Target Spin?
          READ(27,*) TS_BETA
          TS_ALPHA = TS_ALPHA*RADCON !kp: h2model_thia.f:   DATA RADCON/0.017453/ = acos(0)/90 = deg2rad
          TS_BETA = TS_BETA*RADCON
          print*, TS_ALPHA, TS_BETA
       ENDIF



c Jixie: set EB_INDEX according to the target types and beam energies:
c For ND3: EB_INDEX = 1 for 1.339 GeV, 2 for 1.998 GeV;
c For NH3: EB_INDEX = 1 (1.1 GeV), 2 (1.3 GeV), 3 (2.0 GeV), 4 (2.3 GeV), 5 (3.0 GeV)

       EB_INDEX = 0   ! DEFAULT VALUE     
       IF(TARG .EQ. 'ND3') THEN
          IF(ABS(Ebeam-1.339D0) .LT. 0.1) THEN
             EB_INDEX = 1             
          ELSEIF(ABS(Ebeam-1.998D0) .LT. 0.1) THEN
             EB_INDEX = 2             
          ENDIF
       ELSE       
          IF(ABS(Ebeam-1.0539D0) .LT. 0.1) THEN
             EB_INDEX = 1             
          ELSEIF(ABS(Ebeam-1.3369D0) .LT. 0.1) THEN
             EB_INDEX = 2              
          ELSEIF(ABS(Ebeam-1.9911D0) .LT. 0.1) THEN
             EB_INDEX = 3              
          ELSEIF(ABS(Ebeam-2.2596D0) .LT. 0.1) THEN
             EB_INDEX = 4              
          ELSEIF(ABS(Ebeam-2.999D0) .LT. 0.1) THEN
             EB_INDEX = 5             
          ENDIF   
       ENDIF
c Jixie: make sure spit some error message if EB_INDEX is not signed properly    
       IF(EXPERIMENT(1:4) .EQ. 'EG4A' .AND. EB_INDEX .EQ. 0) THEN
          PRINT*,'Error: could not sign value to EB_INDEX!',
     >     ' unknown Ebeam (',Ebeam,')'
          CALL ABORT
       ENDIF



c     kp: From what I saw, it seems everything from each line (not just the values such as 'T') is stored in the arrays.
c     kp: It looks like the input file is re-read and the input info are stored in the array FileText(30) for later use.
       REWIND(27)                !kp: REWIND - positions the file associated with the specified unit to its initial point.
       DO I = 1,26
          READ(27,'(A80)') FileText(I)
      write(6,'(''>'',A80,''<'')') FileText(I)
       ENDDO


c     kp:  Following readout, just to make sure the input is read in correctly? 
      write(6,*) 'INTPEAKING', INTPEAKING        !kp: write(6,*) is for screen output (equivalent to 'cout' ..
      write(6,*) 'EXTPEAKING', EXTPEAKING
      write(6,*) 'TARG ', TARG
      write(6,*) 'PFERMI', PFERMI
      write(6,*) 'POLTYPE ', POLTYPE
      write(6,*) 'EXPERIMENT ', EXPERIMENT
      write(6,*) 'IFFMOD', IFFMOD
      write(6,*) 'IPOL', IPOL
      write(6,*) 'IPOLRES', IPOLRES
      write(6,*) 'IPAULI', IPAULI
      write(6,*) 'IA1', IA1
      write(6,*) 'IEXTERNAL', IEXTERNAL
      write(6,*) 'FILENAME1', FILENAME1
      write(6,*) 'FILENAME2', FILENAME2
      write(6,*) 'FILENAME3', FILENAME3
c      write(6,*) 'DYNAMIC_MODEL', DYNAMIC_MODEL
      write(6,*) 'NMC', NMC
      write(6,*) 'NORES', NORES
      write(6,*) 'BODEK1', BODEK1
      write(6,*) 'ERROR',  ERROR
      write(6,*) 'TARGNEG', TARGNEG
      write(6,*) 'TARGPOS', TARGPOS
      write(6,*) 'DELTEST', DELTEST
      write(6,*) 'ZTEST', ZTEST
      write(6,*) 'MULTISOFT', MULTISOFT
      write(6,*) 'Ebeam', Ebeam        !Jixie print this for debug
      write(6,*) 'EB_INDEX', EB_INDEX  !Jixie print this for debug







! Choose electron as incident particle
      LEPTON = 'ELEC'
      
!Integration parameters: DO NOT CHANGE
      EPS = 1.0D-4              !kp: radcon.inc:  common vars EPS, H1, HMIN => REAL*8 EPS, H1, HMIN, 
      H1 = 0.01D0
      HMIN = 1.0D-30

!     Default is tails on. Flag will be set appropriately in code when needed.
      TAIL_ON = .TRUE.          ! instruct.inc:! TAIL_ON: If .FALSE. turn off radiative tails.
!      TAIL_ON = .FALSE.          ! instruct.inc:! TAIL_ON: If .FALSE. turn off radiative tails.  5/25/12
      
! Following flag is true for normal running
      DEPOL = .TRUE.            !kp: instruct.inc:   ! DEPOL: If .TRUE. Apply electron bremsstrahkung depolarization correction.










      EXPER = EXPERIMENT(1:4)   ! instruct.inc:       COMMON /INSTRUCTc/ EXPER  (So, EXPER carries the first 4 chars of EXPERIMENT
C      print *,'EXPER =',EXPER
      IF(EXPER.EQ.'E155'.OR.EXPER.EQ.'E143'.OR.EXPER
     <     .EQ.'EG1A'.OR.EXPER.EQ.'EG1B'.OR.EXPER.EQ.'EG4A') THEN
         READ(EXPERIMENT(6:6),'(I1)') TMOD !kp: EXPERIMENT in rcslacpol.file has ' EG1B_1 ', so TMOD would be different from 1, may be 0.
         IF(TMOD.LT.1.OR.TMOD.GT.4) TMOD=1
         WRITE(6,'('' TMOD = '',I1)') TMOD !kp: TMOD=IMODE (grep) - one input to INTEG1(TMOD,NPTS,Q2OUT) etc
      ENDIF
!     ================
!    In SUBROUTINE INTEG1(IMODE,NPTS,Q2OUT) etc -   IMODE defines the integration space
*        =1 Use center of target to represent integration space.
*         2 Integration over z of the target only.
*         3 Integrate over raster area at z=target center.
*         4 Full integration over target: VERY CPU consuming.  
!     ================

      TMOD=1                    !kp: 4/4/12: I forced the mode 1 so as to avoid CPU-consuming run of the program









        




      IF(IPAULI.EQ.4) THEN      !kp:   instruct.inc:     >     IPAULI,         ! Pauli suppression model index. (see right above or in file newSFs.f
          IF(TARG.EQ.'ND3') 
     >         CALL Q_E_VANORDEN_INITnew(1.D0,2.D0,PFERMI,1) !kp: newSFs.f:  ENTRY  Q_E_VANORDEN_INITnew(Z1,A1,KF_GEV,IOPT1) (more below ENDIF)
          IF(TARG.EQ.'LID') 
     >         CALL Q_E_VANORDEN_INITnew(1.D0,2.D0,PFERMI,1)       
          IF(TARG.EQ.'HE3') 
     >         CALL Q_E_VANORDEN_INITnew(2.D0,3.D0,PFERMI,1) 
       ENDIF
!     newSFs.f:  SUBROUTINE ELASnew(IFFMOD,IPAULI,Q2,E0,TARG,ITAIL,M,W1,W2,G1,G2) -> returns the structure fnctions W1, W2, G1, and G2,
*              evaluated under elastic conditions, i.e., in terms of nucleon form factors.
!     newSFs.f:  SUBROUTINE PAULI_SUPPnew(IPAULI,IA,QSQ,E0,PS1,PS2) -> Gets Pauli Suppression factor for quasi-elastic scattering from
!              Several different Models. The VanOrden Model is very Slow.
!             IPAULI=1   for Stein model
!             IPAULI=2   for Tsai RMP 46,816(74) eq.B54
!             IPAULI=3   for DeForest and Walecka, Adv. in Phys. 15, 1 (1966).
!             IPAULI=4   for Van Orden model
!     





c       print*,'kp', '  bp', INTPEAKING, EXTPEAKING, TARG
c       print*, PFERMI,POLTYPE,' ',EXPERIMENT,IFFMOD,IPOL,IPOLRES
       print*, IPAULI,IA1,IEXTERNAL,AsymChoice,SFchoice
c       print*,'kp', '  bp'
c       return




       CLOSE(UNIT=27)

       goto 9123

c error handling block

9911   write(6,*) ' '
       write(6,*) ' = = = = =   E R R O R   = = = ='
       write(6,*) ' '
       write(6,*) ' could not open INPUT file'
       write(6,*) ' >', InFile, '<'
       write(6,*) ' '
       goto 9123



9123   continue

c       STOP
       return
       END

      program stegAllvKPA
      IMPLICIT NONE
      INCLUDE 'radcon.inc'
      INCLUDE 'radvarfix.inc'
      INCLUDE 'instruct.inc'
      INCLUDE 'prinplot.inc'
      
      include 'bcs.inc'
      include 'binning.inc'
      include 'csmap.inc'
      include 'regen.inc'

      INCLUDE 'kpaVarChanges.inc'   !Jixie: to acess defaultTarget
c     ======================================= 2/27/13============================
c     
c     First saved the original copy as ~/SingleProg2Stage/BackUp/stegAllmyOwnV1_27_13NoAvgXSinBinMultScB4_2_27_13.f & 
c     made the following modifications:
c     After seeing http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/SimStuffs/ClasSimComp/XsKPAvSEK/plotXSvW_StegSEKnw.gif
c     I wondered what could be limiting the kinematics the low & high Q2 bins.
c     I suspected XCUT value & so I am going to disable the line "XCUT = 1.0/(0.0658/(E*EP*sin2)+1.0)" below
c     so that would mean that XCUT will take the earlier assigned value of XCUT=A*0.9900
c     
c     
c     =============kp: Variable declarations from cs_init.f =======================
      real Ebeam,r_beam,x_beam,y_beam,sth_el_min,sth_el_max,p_el_min
     &     ,p_el_max,cs_maxm,cs_int,sig,theta_el_min,theta_el_max
     &     ,ddsth,ddp
      real theta_el,p_el,t_current,t_length,t_offset,degrad,raddeg
      real v_p_el(0:n_p_el),v_th_el(0:n_th_el),phi_el_min,phi_el_max
      real SIMPS,sgm_model,ppi,delta_phi
      integer i_th_el,i_p_el,ithe,ipe
      data degrad/0.0174532925e0/,ppi/3.14159265e0/ ! ppi used: "Symbol 'pi' at (1) already has basic type of REAL"
c     common/inpt/ithe,ipe
c     ============ Variables needed to add the stuff to be done by read_map program of STEG =============
c     real*8 csmx(0:n_th_el,0:n_p_el),fmx
      real*8 ctha,cthb,pa,pb,acc,acc_tot,norm
      integer STAT

c     ===== ===== My own Additions =====
      INTEGER EB_INDEX,seed     !Jixie: EB_INDEX will be signed according to input beam erngy and target type
      REAL*8 CSRAD              !kp: radiated cross-section to be returned by this function.
      REAL*8 EE,EEP,THETADD,XXCUT,SCAT_PHIRR !kp: added because of "Error: COMMON attribute conflicts with DUMMY attribute in 'e' at (1)", see radvarfix.inc
      INTEGER DBG_MY_KINE
c     By jixie: add this function to check infinity and NAN
      INTEGER IsInfinityOrNAN
      
c     By jixie: add this to check if file exsits
      logical FileExist

c     Jixie: PACKF, ITARG will be used to call xiaochao's subroutine to calculate TA, TB, PACKF      
      REAL*8  PACKF             
      INTEGER ITARG    ! 1/11 for top/bottom long target(1.0cm) and 5 for short target(0.5cm)                   
      
c     Jixie: the following will be used to call energy loss subroutine 
      REAL*8  p_el_d,theta_el_d
      INTEGER istick   ! 2 for long target(1.0cm) and 1 for short target(0.5cm)                   
      
c     integer arg_count         !kp 4/11/12 (crucial to put both cs_map and event generation function in the same prog)
c     character(len=255) cmd
      INTEGER csmapORevGen      !kp: 1 for csmap generation, 2 for event generation.

c     ===========kp: Variables from RCSLACPOL =====
      LOGICAL EXTERNALv, HiPres,flag_acc,flag_remap
      INTEGER IP, NPTS, I, IPLOT, IQ, IT, ITAIL, NEXT, !IW,
     >     LEN1, LEN2, LEN3, IS, TMOD, Nstory
      REAL*8 ENERGY, Z, EINT, SANSWER, TA, Q2OUT, ERC, 
     >     EPRC, THRC, ANSWER, SIGMAR, SIGMARNEW, MULT, 
     >     XCUTRC, W2BIN, TB
      character*1 answ
      character*100 filebos,csmapfile,accfile !kp: csmapfile,accfile 4/23/12
      CHARACTER*40 InputFile, InFile, OutFile
      CHARACTER*40 FILENAME1, FILENAME2, FILENAME3
      CHARACTER*8 TARGETS(5)/'PROTON  ','3HELIUM ',
     >     'ND3     ','LiD     ','NEUTRON '/
      CHARACTER*6 EXPERIMENT
      CHARACTER*60 TITLE
      CHARACTER*9 DATENOW/'   040315'/
      CHARACTER*8 TIMENOW, CLOCK_, DATE
      CHARACTER*80 FileText(30)
      EXTERNAL DATE

!Jixie: add p0_e,th0_e into the ntuple to keep original p_gev and theta_rad
      real p0_el, theta0_el
        
      real DNUP(9),mp,c,rran,d_s,r_d_s !,mpi,pi,me
      real phi_el,el_flag,sin2  !kp: 9/22/12 (added sin2 for XCUT as a function of kinematics (inside the kinematics loop)
      real vert_z,vert_x,vert_y,vert_r,vert_phi
      real vert_el_x,vert_el_y,vert_el_z
      real xmin(2),xmax(2),dxx(2),xx(2),fxmax,phi_sec
      real sgs, Qsq, www,Ep_elas,Ep_cut !kp:5/11/12: Qsq,www added to generate only inelastic events using overall map made earlier
      real THD,thSig,thCor,phCor,GausGenr !kp: 1/27/13: added for the purpose of applying mult. scat. corr.
      integer ipart,NW,IBID,MEVT,RUN
      integer NVAR,lrecl,maxpages,nparz,ntot,isec
      integer istat,icycle,ISTATUS1,NBANK
      integer IMCEV,IMCVX,IMCTK,IHEAD,NR,IPAWC
      integer znt,iphe,ireg_th_el,nfail,nbins(2),ireg(2)
      character*8 ncode(9)
      parameter(maxpages=1000,nparz=10000,ntot=100000)
      common/pawc/ipawc(maxpages*128)
      common/inpt/ithe,iphe,ipe !both in cs_init.f and steg.f of steg_Linux
      
c     !add p0_e,th0_e into the ntuple to keep original p_gev and theta_deg, in total 9 variables now
      data NVAR/9/,lrecl/1024/
      data ncode/'p_e','th_e','ph_e','id_e','x_e','y_e','z_e',
     &           'p0_e','th0_e'/
      
      data mp/0.938272029e0/,c/29.9792458e0/ !,pi/3.14159265359e0/,me/0.510998918E-03/     !kp: pi,me conflict from STEG & rcslacpol 4/15/12
c     FRW   EXTERNAL CLOCK_
      COMMON /KINLIST/ERC,EPRC,THRC,EINT,Z,XCUTRC
      COMMON /INFO/EXTERNALv,NEXT


!     kp:    The following block for command line argument stuff didn't seem to work with g77 (worked with gfortran, at least in Ubuntu)
!ccc  kp: http://gcc.gnu.org/onlinedocs/gfortran/GET_005fCOMMAND.html#GET_005fCOMMAND (2/7/2012) (These three lines just for fun)
c     call get_command(cmd)
ccc   print *, 'The command used was: ', cmd
!cc   kp: http://gcc.gnu.org/onlinedocs/gfortran/COMMAND_005fARGUMENT_005fCOUNT.html (2/7/2012)
c     arg_count = command_argument_count()
ccc   cc     print *, 'arg_count =', arg_count
c     if(arg_count.lt.1) then
c     print*, 'Enter values of csmapORevGen: 1 or 2'
c     print*, '1 for cs-map, 2 for event generation.'
c     return
c     endif 

ccc   call getarg(1,InputFile)
ccc   READ(InputFile,*) th      
c     call getarg(1,InputFile)
c     READ(InputFile,*) csmapORevGen






      degrad=ppi/180.e0
      raddeg=180.e0/ppi
      flag_acc=.false.
      flag_remap=.false.
      el_flag=1.e0

c     check if input file exist 
      FileExist=.false.
      inquire(file='stegAllvPREG.dat', exist=FileExist) 
      if(FileExist) then
        print *,'Reading arguments from file stegAllvPREG.dat' 
      endif
ccccc c ... Read Input from file 'stegAllvPREG.dat'
      if (FileExist) then
         open(1,file='stegAllvPREG.dat',status='old')     
         read(1,*) Ebeam
         print*,'Used beam energy [GeV]:', Ebeam
         read(1,*) t_current
         print*,'Used torus current [A]:', t_current
         read(1,*) r_beam
         print*,'Used beam spot radius [cm]:', r_beam
         read(1,*) x_beam, y_beam
         print*,'Used beam offset: bx and by [cm]:', x_beam, y_beam
         read(1,*)  t_length
         print*,'Used target length [cm]:', t_length
         read(1,*) t_offset
         print*,'Used target offset [cm]:', t_offset
         read(1,*) theta_el_min,theta_el_max
         print*,'Used theta_min and theta_max for detected electron [degree]:', theta_el_min, theta_el_max
         sth_el_min=sin(theta_el_min*degrad)
         sth_el_max=sin(theta_el_max*degrad)
         read(1,*) phi_el_min, phi_el_max
         print*,'Used phi_min and phi_max for detected electron [degree]:', phi_el_min, phi_el_max
         phi_el_min=phi_el_min*degrad
         phi_el_max=phi_el_max*degrad
         read(1,*) p_el_min, p_el_max
         print*,'Used p_min and p_max for detected electron [GeV]:',  p_el_min, p_el_max
         read(1,*) answ
         print*,'Used acceptance cut mode: Y/N  ', answ
         if(answ.eq.'Y') flag_acc=.true.
         read(1,*) answ
         print*,'Do you want to remap: Y/N  ' , answ !kp: not needed in cs_map generation
         if(answ.eq.'Y') flag_remap=.true.
         read(1,*) Nstory
         print*,'Used number of events:' , Nstory !kp: not needed in cs_map generation    
         read(1,'(A100)') csmapfile
         print*,'Used cs-map filename:' , csmapfile
c     read(*,'(A100)') siga_map
c     read(*,'(A100)') sigeltaila_map
c     read(*,'(A100)') sigqtaila_map
c     read(*,'(A100)') sigintaila_map
c     read(*,'(A100)') sigrada_map
c     read(*,'(A100)') sigradanotail_map
c     read(*,'(A100)') delta_i_map
c     read(*,'(A100)') delta_e_map
c     read(*,'(A100)') sigp_map
c     read(*,'(A100)') sigeltailp_map
c     read(*,'(A100)') sigqtailp_map
c     read(*,'(A100)') sigintailp_map
c     read(*,'(A100)') sigradp_map
c     read(*,'(A100)') deltp_i_map
c     read(*,'(A100)') deltp_e_map
         read(1,'(A100)') accfile
         print*,'Used acc-map filename:', accfile
         read(1,'(A100)') filebos
         print*,'Used output filename:', filebos !kp: not needed in cs_map generation
         read(1,*) csmapORevGen !kp: Used 1 or 2 from command line. 
         print*,'Used 1 OR 2 for csmap OR event Generation', csmapORevGen !kp:4/11/12: crucial to combine two functions in one program
         close(1)      
         print *,'Finish reading arguments from file stegAllvPREG.dat' 
      else
ccccc get Input from shell
         print*,'Enter beam energy [GeV]:'
         read*,Ebeam
         print*,'Enter torus current [A]:'
         read*,t_current
         print*,'Enter beam spot radius [cm]:'
         read*,r_beam
         print*,'Enter beam offset: bx and by [cm]:'
         read*,x_beam,y_beam
         print*,'Enter target length [cm]:'
         read*,t_length
         print*,'Enter target offset [cm]:'
         read*,t_offset
         print*,'Enter theta_min and theta_max for detected electron [degree]:'
         read*,theta_el_min,theta_el_max
         sth_el_min=sin(theta_el_min*degrad)
         sth_el_max=sin(theta_el_max*degrad)
         print*,'Enter phi_min and phi_max for detected electron [degree]:'
         read*,phi_el_min,phi_el_max
         phi_el_min=phi_el_min*degrad
         phi_el_max=phi_el_max*degrad
         print*,'Enter p_min and p_max for detected electron [GeV]:'
         read*,p_el_min,p_el_max
         print*,'Enter acceptance cut mode: Y/N' 
         read*,answ
         if(answ.eq.'Y') flag_acc=.true.
         print*,'Do you want to remap: Y/N' !kp: not needed in cs_map generation
         read*,answ
         if(answ.eq.'Y') flag_remap=.true.
         print*,'Enter number of events:' !kp: not needed in cs_map generation
         read*,Nstory
         print*,'Enter cs-map filename:'     
         read(*,'(A100)') csmapfile
c     read(*,'(A100)') siga_map
c     read(*,'(A100)') sigeltaila_map
c     read(*,'(A100)') sigqtaila_map
c     read(*,'(A100)') sigintaila_map
c     read(*,'(A100)') sigrada_map
c     read(*,'(A100)') sigradanotail_map
c     read(*,'(A100)') delta_i_map
c     read(*,'(A100)') delta_e_map
c     read(*,'(A100)') sigp_map
c     read(*,'(A100)') sigeltailp_map
c     read(*,'(A100)') sigqtailp_map
c     read(*,'(A100)') sigintailp_map
c     read(*,'(A100)') sigradp_map
c     read(*,'(A100)') deltp_i_map
c     read(*,'(A100)') deltp_e_map
         print*,'Enter acc-map filename:'     
         read(*,'(A100)') accfile
         print*,'Enter output filename:' !kp: not needed in cs_map generation
         read(*,'(A100)') filebos


         print*,'Enter 1 OR 2 for csmap OR event Generation.' !kp:4/11/12: crucial to combine two functions in one program
         read*,csmapORevGen     !kp: Enter 1 or 2 from command line. 
      endif


      write(6,'(A,5F12.5)'),' sth_el_min,sth_el_max,p_el_min,p_el_max,Ebeam :',
     +       sth_el_min,sth_el_max,p_el_min,p_el_max,Ebeam

c     =========kp: Initialize o/p vars and check if min & max of p & th are nonsensical or not ===== cs_init.f ==

      cs_maxm=0.0e+0
      cs_int=0.0e+0
      if(sth_el_min.ge.sth_el_max) return
      if(p_el_min.ge.p_el_max) return




c     kp: 5/15/2012:   As long as I understand (http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/v115/top.html), 
c         the following two lines are only to initialize the random number generator RANLUX(RRAN,1) by feeding 
c         it with some seed (recalculated within the local subroutine rran_init()) through the call of 
c         RLUXGO(4,idum,0,0) (idum is the recalculated seed, the number returned by GetPID() doesn't really 
c         go as the seed to the generator, but idum does depend on that number.)
      seed = GetPId()
      call rran_init(seed) 


c     Jixie: ============== For rcslacpol ======
c     Move SET_THINGS_UP() here because generating event mode also need to call
c     Xiaochao's routine to calculate target thickness: TA,TB and PF. 
c     Later on PF and target type are also used to calculate energy loss. 
      call SET_THINGS_UP(Ebeam,TMOD,EB_INDEX)

c     Jixie: ============== For rcslacpol ======


      if(csmapORevGen.eq.1) then
         goto 1111              !kp: proceed with Cross-section map generation
      elseif(csmapORevGen.eq.2) then
         goto 2222              !kp: use previously produced CSmap to generate events
      else
         print*,'Eneter either 1 or 2 for csMap or evnt. Gen.'
         goto 9123              !kp: end the program doing nothing
      endif


 1111 continue







      


c     kpa: ============== For rcslacpol ======
c     Jixie moved SET_THINGS_UP above but keep XCUT here, since XCUT is used only in creating xs map
c     Jixie changed XCUT to XCUT=A*0.99      
      XCUT=0.99D0               ! for NH3 
      IF (TARG.EQ.'ND3') THEN
         XCUT=2.D0*0.99D0       ! for ND3
      ENDIF
      
      SCAT_PHIR =0.0000
c     kpa: ============== For rcslacpol ======






c     kpa: ============== For rcslacpol ======
      NPTS = 0                  !kp: Looping/iteration variable?
      PL = -1.D0
      PN = 1.D0
c     kpa: ============== For rcslacpol ======
      print*,'# of Th & p bins: ', n_th_el,n_p_el






c     open(unit=84,file='cs_max_new.dat',status='unknown',IOSTAT=STAT)
c     open(unit=85,file='acc_fc.dat',status='unknown',IOSTAT=STAT)
      open(unit=84,file=csmapfile,status='unknown',IOSTAT=STAT)
      open(unit=85,file=accfile,status='unknown',IOSTAT=STAT)
      norm=0.e0
      acc_tot=0.e0

c     kp:  ================== ================== ============================     
      do i_th_el=1,n_th_el      !kp: Start of Do-loop in i_th_el
         ddsth=(sth_el_max-sth_el_min)/float(n_th_el)
         theta_el=asin(sth_el_min+(i_th_el-0.5D0)*ddsth) !kp:Added 0.5 to calculated the cross-section in the middle of the bin
         do i_p_el=1,n_p_el     !kp: Start of Do-loop in i_p_el
            ddp=(p_el_max-p_el_min)/float(n_p_el)
            p_el=p_el_min+(i_p_el-0.5D0)*ddp !kp:Added 0.5 to calculated the cross-section in the middle of the bin
            
            
            SCAT_PHIR = SCAT_PHIR*RADCON !radcon = radian converion factor = rad2deg or deg2rad?
            ENERGY = Ebeam
            ERC = Ebeam
            EPRC = p_el
            XCUTRC = XCUT
            THRC = theta_el/degrad
            SCAT_PHIR=0.000
            
            NPTS = NPTS + 1
            NEXT = NPTS
            IP = 1              !Don't see any use of it though


c     Following three lines added on 6/3/12
            E=Ebeam
            EP=p_el
            THETAD=THRC
c     kp     ================================ 9/22/12 ======
C     Jixie:  I disable these 3 lines. Sebastian suggeated to keep XCUT=A*0.99
C     sin2 = (1.0-cos(THRC*degrad))/2.0
C     XCUT = 1.0/(0.0658/(E*EP*sin2)+1.0)
C     XCUTRC = XCUT
c     kp     ===============================================

c     ========= kp: cs_init.f  ====== Integrate over given kinematic region  ===========
c     Obtain 2-differential cross section
c     sig=sgm_model(Ebeam,theta_el,p_el) !c      write(81,33) i_p_el,i_th_h,i_ph_h,i_p_h,sig
c     sig=PP_SIGRADA(NPTS)
            TMOD=1
            call sgm_model_rcslacpol(NPTS,TMOD,ENERGY, !E = Ebeam, EP = E' = E_prime and THETAD = theta in degrees
     >           EPRC,THRC,XCUTRC,SCAT_PHIR,sig,EB_INDEX)     
c     By Jixie: sometimes the output PP_SIGA(NPTS),PP_SIGRADA(NPTS),PP_SIGP(NPTS),PP_SIGRADP(NPTS) is infinity or NaN
c     I have to check them before print them into the output file -----20180205
            if( IsInfinityOrNAN(PP_SIGA(NPTS)) .gt. 0 .or. 
     >           IsInfinityOrNAN(PP_SIGRADA(NPTS)) .gt. 0  .or. 
     >           IsInfinityOrNAN(PP_SIGP(NPTS)) .gt. 0  .or.   
     >           IsInfinityOrNAN(PP_SIGRADP(NPTS)) .gt. 0  ) then
               PP_SIGA(NPTS) = 0.0D0 
               PP_SIGRADA(NPTS) = 0.0D0 
               PP_SIGP(NPTS) = 0.0D0 
               PP_SIGRADP(NPTS) = 0.0D0 
               sig = 0.D0
            endif
C     By Jixie: sometimes the radiated unpolarized XS is negative, if it happens, set them to zero   -----20180205
            if(  PP_SIGRADA(NPTS) .lt. 0.0D0 ) then              
               print*, "  Warning: radiated unpolarized XS is negative!!! sig=",sig
               PP_SIGA(NPTS) = 0.0D0 
               PP_SIGRADA(NPTS) = 0.0D0 
               PP_SIGP(NPTS) = 0.0D0 
               PP_SIGRADP(NPTS) = 0.0D0 
               sig = 0.D0
            endif         


c     By Jixie: 20180310   (cutting off events below W=0.95, or E'>EP_ELAS-EP_BIN_WIDTH)
c     regular E' = (M^2+2*M*E-W^2)/(4*E*sin^2(Theta/2)+2*M)
c     elastic scattering: E' = E/(1+E/M*(1-cosTh))
            Qsq=2.e0*Ebeam*p_el*(1.e0-cos(theta_el))       
            www=sqrt(mp**2 + 2.e0*mp*(Ebeam-p_el) - Qsq) !mp = 0.938272029e0 (defined above)
            Ep_elas = Ebeam/(1.e0+Ebeam/mp*(1.e0-cos(theta_el)))
            Ep_cut = Ep_elas - ddp
            if(www .lt. 0.95e0 .or. p_el .gt. Ep_cut) then              
               print*, "  Warning: under elastic peak, reset to zero !!! sig=",sig
               PP_SIGA(NPTS) = 0.0D0 
               PP_SIGRADA(NPTS) = 0.0D0 
               PP_SIGP(NPTS) = 0.0D0 
               PP_SIGRADP(NPTS) = 0.0D0 
               sig = 0.D0
            endif         


            write(6,'(A,I5,4F12.5)'),'steg.f: NPTS,th,p,PP_X(NPTS),sig=',
     >           NPTS,THRC,p_el,PP_X(NPTS),sig


            write(85,34) acc
            acc_tot=acc_tot+acc
            norm=norm+sig
c     write(84,34) sig
            write(84,44) EPRC,THRC,PP_SIGA(NPTS),PP_SIGRADA(NPTS), 
     >           PP_SIGP(NPTS),PP_SIGRADP(NPTS) 
     
c     write(20,12) PP_X(NPTS), Q2OUT, THETAD, PP_SIGA(NPTS), 
c    >  PP_SIGELTAILA(NPTS), PP_SIGQTAILA(NPTS), PP_SIGINTAILA(NPTS),  
c    >  PP_SIGRADA(NPTS), PP_SIGRADANOTAIL(NPTS)), PP_SIGP(NPTS), 
c    >  PP_SIGELTAILP(NPTS), PP_SIGQTAILP(NPTS), PP_SIGINTAILP(NPTS), 
c    >  PP_SIGRADP(NPTS)
c     12     format(F6.3,1x,F6.3,1x,F7.3,1x,11(E12.5,1x))

c     print*,'ip= ',i_p_el,' ith= ',i_th_el
         enddo                  !kp: End of Do-loop in i_p_el
      enddo                     !kp: End of Do-loop in i_th_el


      write(84,34) norm         !kp: One more # added at the o/p file so the total # of rows in the file is nTh*nP+1
      close(84)
      write(85,34) acc_tot
      close(85)
c     ============ from read_map.f =====================

      goto 9122                 !kp: This line added to skip the event generation part in order to generate the map.

      
















 2222 continue                  !Begin Event generation using the previously produced cs-map.
      print*, 'starting event generation!'
c     Importance sampling bins
      xmin(1)=sth_el_min
      xmax(1)=sth_el_max
      xmin(2)=p_el_min
      xmax(2)=p_el_max
      nbins(1)=n_th_el
      nbins(2)=n_p_el
      dxx(1)=abs(xmax(1)-xmin(1))/float(nbins(1))
      dxx(2)=abs(xmax(2)-xmin(2))/float(nbins(2))
      norm=0.0d+0
      acc_tot=0.0d+0

c     Init maps of cross section
      call read_map(norm,acc_tot,flag_acc)
      
      print*,'1.0+norm =',1.0+norm,' acc_tot =',acc_tot
      
c     Open Ntuple
      call hlimit(128*maxpages)
      call hropen(1,'ESCA','steg.hbook','N',lrecl,istat)
      call hbookn(1,' ',nvar,'//ESCA',nparz,ncode)
c     BOS bank parameters
      NW   = 11                 ! number of colums in event bank
      IBID = 12                 ! BOS output device number
      MEVT = 0                  ! event number
      RUN  = 0
c     print*, '1'
      CALL BOS(IWW,Nbcs)
      CALL BKFMT('HEAD','I')
      CALL BKFMT('MCEV','I')
      CALL BKFMT('MCTK','(6F,5I)')
      CALL BKFMT('MCVX','(4F,I)') ! MC vertex parameters
      CALL BLIST(IWW,'E=','HEAD')
c     print*, '2'
      CALL BLIST(IWW,'E+','MCEV')
      CALL BLIST(IWW,'E+','MCVX')
      CALL BLIST(IWW,'E+','MCTK')
      CALL FPARM('OPEN UNIT=12 FILE="' //filebos// '" WRITE RECL=32760'//
     >     ' ACTION=WRITE STATUS=NEW FORM=BINARY')
c     print*, '3'




      

c     This file will be used only temporarily - just to test the effects of multiple scattering
c     The outuput produced by this file will be read in & analyzed/plotted by another program
c     open(61,file='op2testMultScatEffects.dat',status='unknown') !kp: 2/5/13



c     Event Loop               !kp: regen is used by sum_uni_distr(..)
      regen=1
      nfail=0
      DO WHILE (MEVT.le.Nstory)
c     print*, '1'
c     Extract random kinematic point

         

c     print*,'dbg1: ireg(1) ireg(2)',ireg(1),ireg(2) 
c     Importance sampling on polar angle and momentum of electron
c     call sum_uni_distr(xmin,xmax,dx,nbins,cs_max,norm,x,ireg)
         call sum_uni_distr(xmin,xmax,dxx,nbins,fprob,norm,xx,ireg)
c     print*, '2'
         theta_el=asin(xx(1))
         p_el=xx(2)

         theta0_el = theta_el
         p0_el = p_el   




c     kp: 5/11/12 (added to generate exclusively inelastic events (cutting off events below W=1.0)
c     Disable/Comment out these three lines when we need full spectrum in W (not just inelastic events)
c     Qsq=2.e0*Ebeam*p_el*(1.e0-cos(theta_el))      !kp: Q^2 = - 4-mom-transfer-squared = 4*Eb*Ep*sin^2(th/2)
c     www=sqrt(mn**2 + 2*mn*(Ebeam-p_el) - Qsq)     !!kp: average nucleon mass = 0.939 (defined in radcon.inc)
c     if(www<1.0) cycle
c     kp: 5/11/12 (added to generate exclusively inelastic events (cutting off events below W=1.0)




         



c     Extract electron azimuthal angle within CLAS acceptance
         if(flag_acc) then
            call acceptance_el(Ebeam,p_el,(theta_el*raddeg),
     <           t_current,delta_phi)
            if(delta_phi.eq.0.e0) cycle
            isec=1+int(rran()*6.e0)
            phi_el=(float(isec-1)*60.e0+2.e0*(rran()
     <           -0.5e0)*delta_phi)*degrad
            if(phi_el.lt.0.0) phi_el=2.e0*pi+phi_el
         else
            phi_el=phi_el_min+rran()*(phi_el_max-phi_el_min)
         endif
c     goto  8012             !kp:4/15/12


c     8012 continue                  !kp: 4/7/12






c     All lines below with ckpppp were to have the executable which is similar to stegAll_upEb2 except for a slight diff. in 
c     event sampling.


c     ===== 1/27/13  This block added for the purpose of applying/adding multiple scattering effect ====== 
         THD=theta_el*raddeg    !kp: raddeg=1.0/radcon
c     =========The following block for TB & TA is a copy from integ4.f file ============
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
         
C     Added by jixie:
C     The above TA,TB are used for ND3 target. For NH3, call xiaochao's subroutine
C     Here I always use ITARG=1 (1.0cm top NH3) for the inelastic for energies 3.0, 2.3, 2.0 and 1.3 GeV. 
C     For 1.1 GeV only the bottom cell(1.0cm) was used (ITARG=11).
C     For 2.3 GeV or 3.0 GeV, and long target, only the top cell(1.0cm) was used (ITARG=1).
C     Jx: 20190214, add short target (0.5cm) case
C     Jx: 20190326, add bottom long target (1cm) cases for 2.3, 2.0 and 1.3 GeV
         IF(TARG .EQ. 'NH3') THEN
            ITARG = defaultTarget
            IF (EB_INDEX .EQ. 1) THEN
               ITARG = 11
            ELSEIF (EB_INDEX .GE. 4 .AND. defaultTarget .NE. 5 ) THEN   
               ITARG = 1
            ENDIF
            CALL rleg4_simp(TA,TB,PACKF,-100.93D0,theta_el,ITARG,
     >           EB_INDEX)
         ENDIF  
c     Eq. 27.14 from http://pdg.lbl.gov/2009/reviews/rpp2009-rev-passage-particles-matter.pdf for width of mult. scat. ang. distb. 
c     
c     th0 = thSig = sigma of Gaussian= (13.6 MeV/(beta*c*p))*z*sqrt(x/X0)*(1+0.038*ln(x/X0))        -----   (Eq. 27.14)
c     
c     beta = p/E=p/sqrt(p*p+m*m), & for electron it's nearly 1.0 (even for 0.2 GeV electron, the beta is about 0.9999968.., so I'll use 1.0)
c     Z or deuteron is 1 (just as for Hydrogen), c=1, or I will express cp in GeV (the usual unit that we use in data analysis)
         thSig=(13.6/(1000.0*p_el))*sqrt(TA)*(1.0+0.038*log(TA)) ! -----   (Eq. 27.14)  !beta=c=z=1, factor 1000 is to convert p to MeV
         thCor=GausGenr(0.0,thSig) !Generating a Gaussian random # with mean=0 & sigma=thSig
         theta_el = theta_el + thCor !As a correction, adding deflection due to the multiple scattering 
         thCor=GausGenr(0.0,thSig) ! new Gaussian random number for correction to phi_el
C     Jixie: turn off MSC for phi angle
C     phi_el = phi_el + thCor/sin(THD*radcon) !SEK: theta_new=theta_thrown+theta_x    &    phi_new=phi_thrown+theta_y/sin(theta_thrown)

c     !kp: remember a file with ID 61 is opened above outside the event loop
c     kpppp          write(61,41) p_el,(theta_el-thCor)*raddeg,theta_el*raddeg,
c     kpppp     <         (phi_el-thCor/sin(THD*radcon))*raddeg,phi_el*raddeg
c     =================1/27/13 ===========================================================================


C     Added by jixie: 20180406  
C     Calculate most probable ionization energy loss, original code from Mikhail Osipenko  
      IF (ITARG .EQ. 1 ) THEN
         istick = 2
      ELSE
         istick = 1
      ENDIF
      p_el_d = p_el
      theta_el_d = theta_el
      call eloss_ion_prob(p_el_d,theta_el_d,istick,PACKF)
      p_el = p_el_d  



c     print*,'MEVT=',MEVT
c     Accepted event
         MEVT=MEVT+1
c     Extracting vertex position
         vert_r=rran()*r_beam   ! Vertex radius
c     + accordiong e6a target radius=3.5-6 mm
c     + generating outside target volume (on a cilider surface
c     + but Gagik implemented 1cm radius
c     vert_r = 1.
         vert_phi=rran()*2.e0*pi ! Vertex phi
         vert_x=x_beam+vert_r*cos(vert_phi)
         vert_y=y_beam+vert_r*sin(vert_phi)
         vert_z=t_offset+(rran()-0.5e0)*t_length ! Target lenght + offset
c     No detached vertices in this reaction
         vert_el_x=vert_x
         vert_el_y=vert_y
         vert_el_z=vert_z
c     Booking Ntuple
         dnup(1) =p_el
         dnup(2) =theta_el*raddeg
         dnup(3) =phi_el*raddeg
         dnup(4) =el_flag
         dnup(5) =vert_el_x
         dnup(6) =vert_el_y
         dnup(7) =vert_el_z
         dnup(8) =p0_el
         dnup(9) =theta0_el*raddeg
c     Booking BOS bank
c     general banks
         iHEAD = NBANK('HEAD',0,8,1) ! BOS HEADER
         IWW(iHEAD+1)=2         ! Version Number
         IWW(iHEAD+2)=RUN       ! Run Number
         IWW(iHEAD+3)=MEVT      ! Event Number
         IWW(iHEAD+4)=0         ! Event Time
         IWW(iHEAD+5)=-1        ! Event Type
         IWW(iHEAD+6)=0         ! ROC: sync status is 0-OK, > 0bit pattern of offending ROC's
         IWW(iHEAD+7)=7         ! Event Classification from DAQ: 1-15 Physics Events
         IWW(iHEAD+8)=1         ! Level 1 Trigger Latch Word (16 bits)
         NR=MEVT
         iMCEV = NBANK('MCEV',0,2,1)
         IWW(iMCEV+1) = int(rran()*100000) ! first geant random number seed for event
         IWW(iMCEV+2) = int(rran()*100000) ! second seed
         iMCVX = NBANK('MCVX',0,5,1)
         rw(iMCVX+1) = vert_x   ! x of vertex
         rw(iMCVX+2) = vert_y   ! y
         rw(iMCVX+3) = vert_z   ! z
c     rw(iMCVX+4) = 0e0           ! secs of flight
         rw(iMCVX+4) = 1.e-9*(vert_z-t_offset)/c ! secs of flight
         IWW(iMCVX+5) = 0       ! vertex flag
c     Particle bank
         IMCTK = NBANK('MCTK',0,11,2)
c     Electron
         ipart = 11*(1-1)
         rw(imctk+ipart+1)=sin(theta_el)*cos(phi_el) ! x dir cosine at track origin
         rw(imctk+ipart+2)=sin(theta_el)*sin(phi_el) ! y dir cosine
         rw(imctk+ipart+3)=cos(theta_el) ! z dir cosine
         rw(imctk+ipart+4)=p_el ! momentum
         rw(imctk+ipart+5)=me   ! mass
         rw(imctk+ipart+6)=-1.e0 ! charge
         IWW(imctk+ipart+7)=11  ! track Particle Data Group id
         IWW(imctk+ipart+8)=0   ! track flag
         IWW(imctk+ipart+9)=1   ! beginning vertex number
         IWW(imctk+ipart+10)=0  ! ending vertex number
         IWW(imctk+ipart+11)=0  ! parent track
c     Increment event numer
         
c     Jixie: do not print this message event by event to speed up
C     nfail is not used any more            
!        print *,'event number: ',MEVT
!        print *,'number of trials: ',nfail
         nfail=0
         
         if((MEVT/1000*1000).eq.MEVT) print *,' event number ',MEVT
         CALL FWBOS(IWW,IBID,'E',iSTATUS1)
         CALL BDROP(IWW,'E')
         CALL BGARB(IWW)
         if(MEVT.eq.99) then
C     print *, 'jixie debug: call BOSTA'
            call BOSTA		! calling BOS statistic to check
C     print *, 'jixie debug: call BOSBK'
C     By Jixie: this line causes problem! Since it does not do anything to the bos banks, this line can be skipped       
C     call BOSBK(IWW)
         endif
C     print *, 'jixie debug: call HFN'
         CALL HFN(1,dnup)
      ENDDO                     ! End event loop

      
      close(61)                 !kp 2/5/13 Closing the file 'op2testMultScatEffects.dat'






      CALL FWBOS(IWW,IBID,'0',iSTATUS1)
      print *,' end write status ',iSTATUS1
      CALL FCLOS()
      write(*,*) 'Ok, RUN completed after ',Nstory,' events!'
      call hrout(1,icycle,' ')
      call hrend('ESCA')
c     Save new map
      if(flag_remap) call write_map(norm)


















 9122 continue                  !kp: This line added to skip the event generation part in order to generate the map.


c     ========= kp: cs_init.f  ============
 33   format(i5,i5,i5,i5,1pe11.4,1pe11.4) 
 34   format(1pe11.4)
 44   format(f18.9,f18.9,f22.9,f22.9,f22.9,f22.9) !kp:9/24/12
 35   format(1pe11.4$)          ! Missing comma in FORMAT statement at (^)
c     ^
c     35   format(1pe11.4)
c     ======== kp: read_map.f in read_map dir ====
c     10   format(i3)              !kp: Label 10 already defined at (1) when redefined at (2)
 40   format(f6.3,f8.5,f9.6)
 41   format(f10.2,f10.2,f10.2,f10.2,f10.2)
 20   format('Integrated cross section= ',1pe11.4)

      goto 9123


c     error handling block

 9911 write(6,*) ' '
      write(6,*) ' = = = = =   E R R O R   = = = ='
      write(6,*) ' '
      write(6,*) ' could not open INPUT file'
      write(6,*) ' >', InFile, '<'
      write(6,*) ' '
      goto 9123






 9123 continue

c     STOP
      return
      END




!     ===========================================================================    
c     Jixie: Using both user input 'num', process id 'pid' and system time
c     to form a positive interger 'idum' to initialize RANLUX()
c     https://sites.ifi.unicamp.br/mabernal/files/2015/02/RANLUX_James_CPC94.pdf
c     http://luscher.web.cern.ch/luscher/ranlux/
      subroutine rran_init(num)
      implicit none
      integer*4 num
      integer*4 idum,unixtime,pid,GetPId    
      character*28 ctime
      
      call getunixtime(unixtime)
      call getasciitime(unixtime,ctime)
    
      pid=GetPId()             
      
      write(6,*) unixtime,num,pid,
     &     float(unixtime-123736761)/float(unixtime+123736761),
     &     float(num-697899)/float(num+697899),
     &     float(pid-3835)/float(pid+3835)
      
      idum=int(float(num)*abs(
     &     float(unixtime-123736761)/float(unixtime+123736761)
     &     -float(num-697899)/float(num+697899)
     &     +float(pid-3835)/float(pid+3835)))      
      
      write(6,*) 'seed:',idum,' from start time ',ctime
      if(idum.lt.0 .or. idum.ne.idum) stop
      if(idum.eq.num) stop
      if(int(float(idum)/1000.0).lt.1) stop
      
      CALL RLUXGO(4,idum,0,0)  ! this line initialize RANLUX()
      return
      end

      FUNCTION RRAN()
      IMPLICIT NONE
      REAL RRAN
      CALL RANLUX(RRAN,1)
      RETURN
      END


c     Following function made & added on 1/27/13 for a Gaussian generator to make use for Multiple-scattering
c     Code source http://www.design.caltech.edu/erik/Misc/Gaussian.html
c     Returns a random # of gaussian distribution with mean & sigma
      function GausGenr(mean,sigma)
      IMPLICIT NONE
      real GausGenr,w,x1,x2,y1,y2,rran,mean,sigma
      w=1.0                     !kp: Without this the first # came out as NAN
      do while (w >= 1.0)
         x1 = 2.0 * rran() - 1.0
         x2 = 2.0 * rran() - 1.0
         w = x1 * x1 + x2 * x2
      end do
c     print*,w,log(w), -2.0*log(w), -2.0*log(w)/w, sqrt((-2.0*log(w))/w)
      w = sqrt( (-2.0 * log( w ) ) / w )
c     print*,w,
      y1 = x1 * w
      y2 = x2 * w
c     First the gaussian # of mean=0 & sigma=1
      GausGenr=y1               !kp: one could also use y2, I think 
c     Now, the gaussian # with given values for mean & sigma 
      GausGenr=sigma*y1+mean
      return
      end


C     Added by Jixie: check if it is inifinity(return 1) or NAN(return 2)
      integer function IsInfinityOrNAN(A)
      IMPLICIT NONE
      real*8, intent(in) :: A
      integer ret
      ret = 0
      if (A > HUGE(A)) then
         print*, "    Catch an INFINITY!!!"
         ret=1
      endif
      if(isnan(A)) then
         print*, "    Catch an NaN!!!"
         ret=2
      endif    
      IsInfinityOrNAN = ret    
      return 
      end function
      

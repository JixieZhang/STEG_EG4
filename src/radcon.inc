*******************************************************************************
* RADCON.INC : Include file of useful parameters and constants for 
* radiative corrections.
*
* 1/93, LMS.
******************************************************************************
       REAL*8 ME, ME2, MN, MN2, MMU, MTAU, PI, ALPHA, EU, ED, 
     >        MPI, MD, MHE3, RADCON, HC2

       PARAMETER (ME     = 0.000511D0,       ! Electron mass (GeV)
     >            ME2    = 2.61121D-07,      ! ME*ME
     >            MN     = 0.939D0,          ! Nucleon mass (GeV)
     >            MN2    = 0.881721D0,       ! MN*MN
     >            MMU    = 0.105658D0,       ! Muon mass (GeV)
     >            MTAU   = 1.7841D0,         ! Tauon mass (GeV)
     >            PI     = 3.141592654D0,   
     >            ALPHA  = 0.00729735D0,     ! Fine structure constant. 
     >            EU     = 0.66666667D0,     ! Up quark fractional charge. 
     >            ED     =-0.33333333D0,     ! Down quark fractional charge.
     >            MD     = 1.87561339D0,     ! Deuteron mass (GeV).
     >            MPI    = 0.1349739D0,      ! PI0 mass (GeV)
     >            MHE3   = 2.80841D0,        ! 3He nucleus mass (GeV).
     >            RADCON = 1.745329252D-2,   ! PI/180.0
     >            HC2    = 389.379D0)        ! (HBAR C)**2 = GeV**2*microbarns

       REAL*8 M,                         ! Target mass.
     >        ML,                        ! Incident lepton mass.
     >        AOP,                       ! ALPHA/PI.
     >        EU2,                       ! EU**2.
     >        ED2,                       ! ED**2.
     >        M2,                        ! M**2.
     >        ML2,                       ! ML*ML
     >        MMU2,                      ! MMU**2.
     >        MTAU2,                     ! MTAU**2.
     >        MPITHR,                    ! (MN +MPI)**2 = pion threshold.
     >        PI2                        ! PI**2.

! Integration variables
       REAL*8 EPS, H1, HMIN

!Fermi momentum
       REAL*8 PFERMI

       COMMON /RADCONST/ M, ML, AOP, EU2, ED2, M2, ML2,
     >                 MMU2, MTAU2, MPITHR, PI2, EPS, H1, 
     >                 HMIN, PFERMI


***************************************************************************
* This include file is for storing cariables used in printing and plotting
* results.
* 2/94, LMS.
***************************************************************************
       INTEGER PP_NMAX
c       PARAMETER (PP_NMAX = 500) 
c       PARAMETER (PP_NMAX = 501*501)   ! Added by KPA (3/28/2012)	
       PARAMETER (PP_NMAX = 505*505)   ! Added by KPA (9/17/2012)	

       REAL*8 PP_SIGA(PP_NMAX),        ! Unpolarized Born cross section.
     >        PP_SIGELTAILA(PP_NMAX),  ! Unpolarized elastic Rad. tail.
     >        PP_SIGQTAILA(PP_NMAX),   ! Unpolarized quasielastic rad. tail.
     >        PP_SIGINTAILA(PP_NMAX),  ! Unpolarized inelastic rad. tail.
     >        PP_SIGRADA(PP_NMAX),     ! Unpolarized radiated cross section.
     >        PP_SIGRADANOTAIL(PP_NMAX), ! Unpolarized radiated cross section 
                                         ! W/O background tails.
     >        PP_DELTA_I(PP_NMAX),     ! Unpolarized delta RC. Internal.
     >        PP_DELTA_E(PP_NMAX),     ! Unpolarized delta RC. External.
     >        PP_SIGP(PP_NMAX),        ! Polarized Born cross section.
     >        PP_SIGELTAILP(PP_NMAX),  ! Polarized elastic rad. tail.
     >        PP_SIGQTAILP(PP_NMAX),   ! Polarized quasielastic rad. tail.
     >        PP_SIGINTAILP(PP_NMAX),  ! Polarized inelastic tail.
     >        PP_SIGRADP(PP_NMAX),     ! Polarized radiated cross section.
     >        PP_DELTP_I(PP_NMAX),     ! Polarized delta RC. Internal.
     >        PP_DELTP_E(PP_NMAX),     ! Polarized delta RC. External.
     >        PP_X(PP_NMAX),           ! X scaling variable (Bjorken).
     >        PP_Y(PP_NMAX),           ! Y scaling variable ( (E-EP)/E ).
     >        PP_DVR(PP_NMAX),         ! DELTA_VR.
     >        PP_A0(PP_NMAX),          ! Born asymmetry (raw).
     >        PP_ARAD_I(PP_NMAX),      ! Radiated asymmetry. Internal.
     >        PP_ARAD_E(PP_NMAX),      ! Radiated asymmetry. External.
     >        PP_RAT1_I(PP_NMAX),      ! A0/A1. Internal.
     >        PP_RAT1_E(PP_NMAX),      ! A0/A1. External.
     >        PP_RAT2(PP_NMAX),        ! Unpolarized radiated/Born tail.
     >        PP_RAT3(PP_NMAX),        ! Polarized radiated/Born tail.
     >        PP_ADIFF_I(PP_NMAX),     ! (A0 - A1)*100 asymmetry difference.
     >        PP_ADIFF_E(PP_NMAX) ,    ! (A0 - A1)*100 asymmetry difference.
     >        PP_FRC_I(PP_NMAX),       ! Multiplicative RC dilution.
     >        PP_ARC_I(PP_NMAX),       ! Additive RC
     >        PP_FRC_E(PP_NMAX),       ! Multiplicative RC dilution.
     >        PP_ARC_E(PP_NMAX)        ! Additive RC

       COMMON /PRINPLOT/PP_SIGA, PP_SIGELTAILA, PP_SIGQTAILA,
     >                  PP_SIGINTAILA, PP_SIGRADA, PP_SIGRADANOTAIL,
     >                  PP_DELTA_I, PP_DELTA_E, PP_SIGP, 
     >                  PP_SIGELTAILP, PP_SIGQTAILP, PP_SIGINTAILP,
     >                  PP_SIGRADP, PP_DELTP_I, PP_DELTP_E, PP_X,
     >                  PP_Y, PP_DVR, PP_A0, PP_ARAD_I, PP_ARAD_E,
     >                  PP_RAT1_I, PP_RAT1_E, PP_RAT2, PP_RAT3, 
     >                  PP_ADIFF_I, PP_ADIFF_E, PP_FRC_I, PP_ARC_I,
     >                  PP_FRC_E, PP_ARC_E

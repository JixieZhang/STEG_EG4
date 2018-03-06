*******************************************************************************
* RADVARFIX.INC : Include file of useful variables for radiative corrections 
* Which remain fixed for a kinematic point.
*
* 1/93, LMS.
******************************************************************************
       CHARACTER*4 POLTYPE               ! Nucleon polarization type,
                                         ! ='LONG'  'TRAN' or 'ARBI'.
       CHARACTER*4 LEPTON                ! Lepton type,
                                         ! ='ELEC' or 'MUON'.
       CHARACTER*3 TARG                  ! TARG = 'NH3', 'HE3', or 'ND3'

       REAL*8 E,                         ! Incident electron energy.
     >        EP,                        ! Scattered electron energy.
     >        THETAD,                    ! Scattering angle (degrees).
     >        THETAR,                    ! Scattering angle (radians).
     >        K1(4),                     ! Incident lepton 4-momentum.
     >        K2(4),                     ! Scattered lepton 4-momentum.
     >        P1(4),                     ! Initial target 4-momentum.
     >        XI(4),                     ! lepton polarization vector.
     >        ETA(4),                    ! Target polarization vector.
     >        Q(4),                      ! K1 - K2.
     >        QDOTXI,                    ! Dot product.
     >        QDOTETA,                   ! Dot product.
     >        XIDOTETA,                  ! Dot product.
     >        P1DOTETA,                  ! Dot product.
     >        P1DOTXI,                   ! Dot product.
     >        P1DOTQ,                    ! Dot product.
     >        K1DOTETA,                  ! Dot product.
     >        K2DOTETA,                  ! Dot product.
     >        K2DOTXI,                   ! Dot product.
     >        Q2,                        ! 4-momentum transfer squared.
     >        XS,                        ! X scaling variable (Bjorken).
     >        YS,                        ! Y scaling variable.
     >        PL,                        ! Degree of longitudinal lepton pol.
     >        PN,                        ! Degree of nucleon polarization.
     >        S,                         ! -2*M*E
     >        X,                         ! -2*M*E' 
     >        Y,                         ! Q**2
     >        W2,                        ! Missing mass squared.
     >        SX,                        ! S - X.
     >        LAMBDAS,                   ! S*S - 4*ML2*M2
     >        RTLAMBDAS,                 ! SQRT(LAMBDAS)
     >        LAMBDAX,                   ! X*X - 4*ML2*M2
     >        RTLAMBDAX,                 ! SQRT(LAMBDAX)
     >        NN,                        ! 2*ALPHA**2/RTLAMBDAS
     >        LAMBDAM,                   ! Y*Y + 4*ML2*Y
     >        RTLAMBDAM,                 ! SQRT(LAMBDAM)
     >        YM,                        ! Y + 2*ML2
     >        LM, LMTAU, TR, LITLM,
     >        J0,                        ! 2*(YM*LM - 1)
     >        SP,                        ! X + Y.
     >        XP,                        ! S - Y.
     >        NUMAX,                     ! W2 - MPITHR
     >        LAMBDASP,                  ! SP*SP - 4*ML2*W2
     >        LAMBDAXP,                  ! XP*XP - 4*ML2*W2
     >        RTLAMBDASP,                ! SQRT(LAMBDASP)
     >        RTLAMBDAXP,                ! SQRT(LAMBDAXP)
     >        LSP, LXP,
     >        AZ,                        ! SX*SX + 4*M2*Y.
     >        SSP,                       ! S + X.
     >        DELTA,                     ! 2[Y(S*X - M2*Y) - ML2*AZ].
     >        C1, C2, C1HAT, C2HAT,
     >        SIG0KIN,                   ! Born cross section.
     >        SIG0KINA,                  ! Unpolarized part of SIG0KIN
     >        SIG0KINP,                  ! Polarized part of SIG0KIN
     >        XCUT,                      ! x Integration cuttoff for
     >                                   ! inelastic tail background calc.
     >        TS_ALPHA,			 ! target spin vector's polar angle (rad)
     >        TS_BETA,			 ! target spin vector's azimuthal angle (rad)
     >        SCAT_PHIR			 ! electron scattering azimuthal angle (rad)

       COMMON /RADVARFIX2/ POLTYPE, LEPTON, TARG
       COMMON /RADVARFIX3/ E, EP, THETAD, 
     >                    THETAR, K1, K2, P1, XI, ETA, Q, QDOTXI, 
     >                    QDOTETA, XIDOTETA, P1DOTETA, P1DOTXI, P1DOTQ, 
     >                    K1DOTETA, K2DOTETA, K2DOTXI, Q2, XS, YS, PL, 
     >                    PN, S, X, Y, W2, SX, LAMBDAS, RTLAMBDAS, 
     >                    LAMBDAX, RTLAMBDAX, NN, LAMBDAM, RTLAMBDAM, 
     >                    YM, LM, LMTAU, TR, LITLM, J0, SP, XP, NUMAX, 
     >                    LAMBDASP, LAMBDAXP, RTLAMBDASP, RTLAMBDAXP, 
     >                    LSP, LXP, AZ, SSP, DELTA, C1, C2, C1HAT, 
     >                    C2HAT, SIG0KIN, SIG0KINA, SIG0KINP, XCUT,
     >                    TS_ALPHA, TS_BETA, SCAT_PHIR


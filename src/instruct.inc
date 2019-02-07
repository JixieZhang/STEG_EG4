*******************************************************************************
* INSTRUCT.INC : Include file of instructional flags and indices.
*
* 1/93, LMS.
******************************************************************************
       LOGICAL INTPEAKING,               ! Peaking internal (otherwise "exact")
     >         FL_POL,                   ! If T do "polarized" analysis.
     >         FL_UNPOL                  ! If T do "unpolarized" analysis.

! If NMC = .TRUE., F2NMC is called, otherwise F2GLOB is called.
! If NORES is .TRUE., F2GLOB is extrapolated to the resonance region. 
! If BODEK1 is .TRUE., BODEK is called, otherwise H2MODEL is called.
! If NMC = .FALSE. and ERROR = .TRUE., W1 and W2 are increased by their errors
! TARGNEG: If .TRUE., TA and TB are decreased by 4%
! TARGPOS: If .TRUE., TA and TB are increased by 4%
! ZTEST: If .TRUE., Z changed to 1
! DELTEST: If .TRUE., change integration DELTA from 0.01 to 0.05
! EXTPEAKING: If .TRUE., The external corrections use peaking approx.
! TAIL_ON: If .FALSE. turn off radiative tails.
! DEPOL: If .TRUE. Apply electron bremsstrahkung depolarization correction.
! MULTISOFT: If .TRUE. correct for multi-soft photons via exponentiation.
!            Also, radiative tail definitions change slightly.


       LOGICAL NMC, NORES, BODEK1, ERROR, TARGNEG, TARGPOS, 
     >     DELTEST, ZTEST, EXTPEAKING, TAIL_ON, DEPOL, MULTISOFT
 
       INTEGER IUNPOL,                   ! F2 model index.
     >         IPOL,                     ! G1, G2 model index.
csk    IPOL = 1: standard version of A2 in the DIS region
csk    IPOL = 2: extra term g2tw3 in A2 in DIS
csk    IPOL = 4: g2 = 0 om DIS
csk    IPOL = 5: A2 = 0 in DIS
csk    IPOL = 7: A2 = Soffer bound in DIS
     >         IFFMOD,                   ! Form factor model index.
     >         IPAULI,                   ! Pauli suppression model index.
     >         IPOLRES,                  ! Polarized resonance region model.
     >         IA1,                      ! Neutron A1 model index.
csk     IA1 = 4: sek_me (2006/7 combined Doaa/sek fit to world data including EG1)
csk     IA1 = 5: sek_me; add correlated fit error to A1
csk     IA1 = 6: sek_me; subtract correlated fit error to A1
     >         IEXTERNAL                 ! External code testing index

csk       add new variables to have a simple way of chosing alternative models
      
       integer AsymChoice
csk    1 = standard: Use aparms.dat for A1, A2 in RR, IPOL = 1 -> g2WW in DIS
csk    2: standard except aparms0.dat for A1 (a1rr_ref.dat) in g2g1res
csk    3: standard except aparms0.dat for A2 (a2rr_ref.dat) in g2g1res
csk    4: standard except aparms0.dat for both A1 and A2 (can skip)
csk    5: standard except A2 = 0 for both RR and DIS (IPOL = 5)
csk    6: standard except A2 = a2rr_res.dat & A1 = a1rr_res.dat in the res. region
csk    7: standard except A2 contains extra term g2tw3 in DIS (IPOL = 2)
csk    8: standard except use new model for A2 in RR based on DMT
csk    9: standard except use Soffer limit for A2 throughout
csk    10: standard except use Soffer limit divided by W for A2 throughout
csk    11: New Standard Resonance Model 2008-09 (c) Guler/Kuhn
csk    12: Preliminary version v1 of above - 2008-09 (c) Guler/Kuhn 
csk    13: Old Standard 1 reimplemented in 2008-09 (c) Guler/Kuhn
csk NOTE: 11-13 work both with old and new version of strucfunc but give different results for A2
csk       14-15 only work with new version of strucfunc and newSF.f instead of models.f
csk    14: Alternative A2 resonance model: A2_MAID
csk    15: Old version of A2 resonance model 

       integer SFchoice
csk    1 = standard: Use R1998 throughout, F2NMC_NEW for DIS, H2MODEL for RR
csk    2: standard except use R1998 + DR in DIS and H2MODEL
csk    3: standard except use F2GLOB instead of F2NMC in DIS
csk    4: standard except use RRicco instead of R1998 in H2MODEL (RR)
csk    5: standard except use Bodek instead of H2MODEL in the RR
csk    6: standard except use H2HallC instead of H2MODEL for F1 and R in the RR
csk    7: standard except use R1998 down to Q2=0 (using S. Rock's version of R1998)
csk    8: standard except use new F1p and Flp fit by Eric Christy (Fall 2004)
csk    9: Like 8 except NMC-param. of F2n/p
csk    10: New version of F1n/F1p/F1d by Peter Bosted/Eric Christie (c) 2006
csk    11: Newest version of F1n/F1p/F1A by Peter Bosted/Eric Christie (c) 2007, HERMES
csk    12: Same version as 11, but with errors added to F2 (and proportionally F1)
csk    13: Same version as 11, but with errors subtracted from R (F2 unchanged)
csk    14: Same version as 11, but without tabulated R values in RR substituted
csk    15-16: Analog to 12-13, but without tabulated R values in RR substituted
csk    17: Newest kludge for Rd, F2d from Eric Christy December 2011. Should be only for D
csk    18-19: Analog to 15-16 for new kludge
csk    20: New version by Christie/Kalantarians 2014
csk    21-22 corresponding error estimates, analog to 12-13

      
      
! EXP is the eperiment name for target integration purposes: "E142", "E143",
! "E154", "E155", or "NONE" where "NONE" is for no target integration.
       CHARACTER*4 EXPER


cfrw  this common block was split first on Bernoulli (ODU) 4/8/99
      COMMON /INSTRUCTL/ INTPEAKING, FL_POL, FL_UNPOL,
     >      NMC, NORES, BODEK1, ERROR, TARGNEG,
     >      TARGPOS, DELTEST, ZTEST, TAIL_ON, EXTPEAKING, 
     >      DEPOL, MULTISOFT

       COMMON /INSTRUCTi/ IUNPOL, IPOL, IFFMOD, IPAULI, IPOLRES, 
     >      IA1, IEXTERNAL, AsymChoice, SFchoice

       COMMON /INSTRUCTc/ EXPER

! -*-Mode: Fortran; compile-command: "f77 -o r1990.o -c r1990.f"; -*-
      SUBROUTINE R1998(X,Q2,R,DR,GOODFIT)                               
                                                                        
!----------------------------------------------------------------       
! X      : Bjorken x                                                    
! Q2     : Q squared in (GeV/c)**2                                      
! R      :                                                              
! DR     : Absolute error on R                                          
! GOODFIT:  = .TRUE. if the X,Q2 are within the range of the fit.       
!-------------------------------------------------------------------    
! Model for R, based on a fit to world R measurements. Fit performed by 
! program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details   
! see Reference.                                                        
!                                                                       
! Three models are used, each model has three free parameters.  The     
! functional forms of the models are phenomenological and somewhat      
! contrived.  Each model fits the data very well, and the average of    
! the fits is returned.  The standard deviation of the fit values is    
! used to estimate the systematic uncertainty due to model dependence.  
!                                                                       
! Statistical uncertainties due to fluctuations in measured values have 
! have been studied extensively.  A parametrization of the statistical  
! uncertainty of R1990 is presented in FUNCTION DR1990.                 
!                                                                       
!                                                                       
! Each model fits very well.  As each model has its own strong points   
! and drawbacks, R1998 returns the average of the models.  The          
! chisquare for each fit (237 points with 6 parameters) are:  
!          ALL DATA  #PTS=237         |     X<  0.07 #PTS= 28
! FIT  #PARAM CHISQ  CHISQ/DF PROB(%) |  CHISQ  CHISQ/DF   PROB(%)
! R1998         217.4   0.94    73.1        28.7   1.06    37.5
! R1998a   6    219.4   0.95    69.8        28.8   1.07    37.1
! R1998b   6    217.7   0.94    72.6        27.8   1.03    42.2
! R1998c   6    221.9   0.96    65.5        30.9   1.15    27.4                

!                                                                       
! This subroutine returns reasonable values for R for all x and for all 
! Q2 greater than or equal to .3 GeV.                                   
!                                                                       
! The uncertainty in R originates in three sources:                     
!                                                                       
!     D1 = uncertainty in R due to statistical fluctuations of the data 
!          as reflected in the error of the fit. 
!          It is parameterized in FUNCTION DR1990, for details see     
!          Reference.                                                   
!                                                                       
!     D2 = uncertainty in R due to possible model dependence, approxi-  
!          mated by the variance between the models.                    
!                                                                       
!     D3 = uncertainty in R due to possible epsilon dependent errors    
!          in the radiative corrections, taken to be ± .025.  See     
!          theses (mine or Dasu's) for details.  This is copied from R1990                       
!                                                                       
! and the total error is returned by the program:                       
!                                                                       
!     DR = is the total uncertainty in R, DR = sqrt(D12+D22+D32). 670
!          DR is my best estimate of how well we have measured R.  At   
!          high Q2, where R is small, DR is typically larger than R.  If
!          you have faith in QCD, then, since R1990 = Rqcd at high Q2,  
!          you might wish to assume DR = 0 at very high Q2.             
!                                                                       
! NOTE:    In many applications, for example the extraction of F2 from  
!          measured cross section, you do not want the full error in R  
!          given by DR.  Rather, you will want to use only the D1 and D2
!          contributions, and the D3 contribution from radiative        
!          corrections propogates complexely into F2.  For more informa-
!          tion, see the documentation to dFRC in HELP.DOCUMENT, or     
!          for explicite detail, see Reference.                         
!                                                                       
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IMPLICIT NONE                                                     
      REAL*8 QP2,FAC,RLOG,Q2THR,R_A,R_B,R_C,R, D1,D2,D3,DR,DR1998,X,Q2
      REAL*8 Q2_SAVE

      REAL*8 Q2_LIMIT/.2/
      REAL*8 A(6) /4.8520D-02,  5.4704D-01,  2.0621D+00,
     >          -3.8036D-01,  5.0896D-01, -2.8548D-02/
      REAL*8 B(6) /4.8051D-02,  6.1130D-01, -3.5081D-01, 
     >          -4.6076D-01,  7.1697D-01, -3.1726D-02/
      REAL*8 C(6) /5.7654D-02,  4.6441D-01,  1.8288D+00,
     >           1.2371D+01, -4.3104D+01,  4.1741D+01/


      LOGICAL GOODFIT                                                   
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      GOODFIT = .TRUE.                                                  
      IF ((X.LT.0.02).OR.(Q2.LT.0.3)) GOODFIT = .FALSE.    

      Q2_SAVE = Q2             ! Bugs fixed in Low Q2 manipulation.
      IF(Q2.LT.Q2_LIMIT) THEN   ! added 11/15/01 to avoid low Q2 problems caused by RLOG
       Q2 = Q2_LIMIT
      ENDIF
                                 
      FAC = 1.D0+12.D0*Q2/(Q2+1.D0)*(.125D0**2 /(.125D0**2+X**2))
      RLOG  = FAC/DLOG(Q2/.04D0)!   <--- we use natural logarithms only!      


!Model A
      QP2 = (A(2)/(Q2**4+A(3)**4)**.25) * (1.D0+A(4)*X +A(5)*X**2) !1.07   .030
      R_A   = A(1)*RLOG +QP2*X**A(6)                                           
!Model B
      QP2 = (1.D0+B(4)*X+B(5)*X**2)*(B(2)/Q2 +B(3)/(Q2**2+.09D0)) !1.06   .042
      R_B   =  B(1)*RLOG+QP2*X**B(6)   
!Model C
      Q2thr =C(4)*(X) +c(5)*x**2 +c(6)*X**3  
      QP2 =  C(2)/DSQRT((Q2-Q2thr)**2+C(3)**2) 
      R_C   =  C(1)*RLOG+QP2    

      R     = (R_A+R_B+R_C)/3.D0                                          


      
      D1    = DR1998(X,Q2)                                              

      D2    = DSQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.D0)               
      D3    = .023D0*(1.D0+.5D0*R)                                            
              IF (Q2 .LT. 1.0D0 .OR. X.LT. 0.1D0) D3 = 1.5D0*D3                       


      DR    = DSQRT(D1**2+D2**2+D3**2)                                   
                                                                             
      Q2 = Q2_SAVE  ! restore Q2  added 11/2/01 to fix bug at low Q2 
      IF(Q2.LT.Q2_LIMIT) THEN   ! added 11/15/01 to avoid low Q2 problems caused by RLOG
       R = Q2/Q2_LIMIT * R
      ENDIF   
      RETURN                                                            
      END                                                               
                                                                        
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                                        
      REAL*8 FUNCTION DR1998(X,Q2)                                             
                                                                        
! Parameterizes the uncertainty in R1990 due to the statistical         
! fluctuations in the data.  Values reflect an average of the R-values  
! about a neighborhood of the specific (x,Q2) value.  That neighborhood 
! is of size [+/-.05] in x, and [+/-33%] in Q2.  For details, see       
! Reference.                                                            
!                                                                       
! This subroutine is accurate over all (x,Q2), not only the SLAC deep   
! inelastic range.  Where there is no data, for example in the resonance
! region, it returns a realistic uncertainty, extrapolated from the deep
! inelastic region (suitably enlarged).  We similarly estimate the      
! uncertainty at very large Q2 by extrapolating from the highest Q2     
! measurments.  For extremely large Q2, R is expected to fall to zero,  
! so the uncertainty in R should not continue to grow.  For this reason 
! DR1990 uses the value at 64 GeV for all larger Q2.                    
!                                                                       
! XHIGH accounts for the rapidly diminishing statistical accuracy for   
! x>.8, and does not contribute for smaller x.                          
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      REAL*8 Q2,X                   

      DR1998 = .0078D0 -.013D0*X +(.070D0 
     >     -.39D0*X+.70D0*X**2)/(1.7D0+Q2)
      RETURN                                                            
      END                                              
                 
!======================================================
      SUBROUTINE R1998abc_doc(unit)
      IMPLICIT NONE
      INTEGER UNIT
      
      WRITE(UNIT,'(/''R1998 COMMENTS''/
     >'' R forced to go linearly to 0 from value at Q2=0.2'')')
      RETURN
      END

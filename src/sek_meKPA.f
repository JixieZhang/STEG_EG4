* KPA 4/16/12: This file is the same as seek_me.f except for the fact that the input data files
*   (disa1n_fit_error.dat, disa1n_fit_param.dat,disa1p_fit_error.dat, disa1p_fit_param.dat) will be 
*   opened outside sek_me(), in order to make the program faster by making it open/read/close the 
*   files only once rather than multiple times.
*
*  I realized that the flag 'Firstcallsek' makes the code open and read the output files only once
*    unlike what I had suspected. So, I didn't modify the code.
*
      SUBROUTINE sek_me(Wsq, Q2, targ, signp, A1, delA1) !"Final Millenium Edition"
      
Cc    FRW  changed 10/8/97 from hard-coded constants to constants read from file
Cc    FRW  changed 11/22/97 to include target type '3' for He-3
Cc    FRW  changed 10/13/98 to handle NaNQ due to F1d=0
Cc    FRW  changed 11/10/98 to use mixed Q2 parameter
Cc    SEK  new century version for Sebastian's "xi"-based fit to A1p and A1n (DIS)
Cc    DNS  Millenium Edition for new fit by Doaa Teama, Nevzat Guler and SEK based on
Cc    DNS  ELSA photon, E143/155, EMC/SMC, HERMES and EG1a/b data. Includes error
Cc    DNS  estimate based on MINUIT error matrix. 30-Oct-2006 SEK Final version 1-Aug-2007 Doaa/SEK
      
      IMPLICIT NONE
      
      REAL*8 X, Q2, Wsq, A1, delA1, signp, XI, nu
      REAL*8 py/3.141592654/
      real*8 A1p, A1n, delA1p, delA1n
      real*8 F1p, F1n, F1d , A, B, C, D
      real*8 paramp(7),errparp(7,7), deriv(7)
      real*8 paramn(7),errparn(7,7)
      real*8 wd/0.05D0/
      real*8 M2/0.880354345984D0/
      real*8 scale
      real*8 phoQ2
      real*8 dummy1, dummy2
      integer i, j, k, n, m
      integer NUMP, NUMN  ! number of parameters for proton and neutron fits
      CHARACTER*1 targ    ! one of P, N, D, 3
      logical suppress, Firstcallsek /.TRUE./
      
      save Firstcallsek
      
Cc ======================================================================================
Cc ======================= Begin to read parameters =====================================
Cc ======================================================================================
Cc
Cc    NOTE:
Cc    Below we read the parameters from a file.
Cc

Cc    Proton (7 parameter fit disa1p 01S1, change NUMP if number of parameters changes)
Cc    ========
      NUMP = 7

Cc    Neutron (7 parameter fit disa1n 01S1, change NUMN if number of parameters changes) 
Cc    ========
c      NUMN = 7
c      NUMN = 6                  !kpa: 3/31/12
      NUMN = 7                  !kpa: 8/19/13  Otherwise "C = dsin(py* XI**paramn(7))" at 138th line below wont work due to paramn(7)

      if(Firstcallsek) then

        OPEN(UNIT=55,FILE='disa1p_fit_param.dat',STATUS='OLD')
        DO N=1,NUMP
          READ(55,*) paramp(N)
        ENDDO
        CLOSE(UNIT=55)     
        OPEN(UNIT=55,FILE='disa1p_fit_error.dat',STATUS='OLD')
        DO N=1,NUMP
          READ(55,*) (errparp(N,M),M=1,NUMP)
        ENDDO
        CLOSE(UNIT=55)

c        print*,'kpdbg: sek_me file open'

        OPEN(UNIT=55,FILE='disa1n_fit_param.dat',STATUS='OLD')
        DO N=1,NUMN
          READ(55,*) paramn(N)
        ENDDO
        CLOSE(UNIT=55)     
        OPEN(UNIT=55,FILE='disa1n_fit_error.dat',STATUS='OLD')
        DO N=1,NUMN
          READ(55,*) (errparn(N,M),M=1,NUMN)
        ENDDO
        CLOSE(UNIT=55)

        Firstcallsek = .FALSE.

      endif
	   

Cc ======================================================================================
Cc ======================== End to read parameters ======================================
Cc ======================================================================================

      scale = 1.D0 - 1.5D0*wd   ! F1F2 puts F1d, F2d as per nucleUS!
      A1 = 0.0D0
      delA1 = 1.0D0
      
Cc    sek introduce new variable which extrapolates smoothly to photon point
      phoQ2 = 0.27150
      nu = (Wsq - 0.93827D0*0.93827D0 + Q2)/1.87654D0
      if (nu .gt. 0.0D0 .and. Q2 .ge. 0.0) then
         XI = (Q2 + phoQ2)/0.93827D0/(nu + dsqrt(nu*nu+Q2))
         X = Q2/2.0D0/0.93827D0/nu
      else
         write (6,*) ' Bad nu, Wsq or Qsq!', nu, Wsq, Q2
         return
      endif
      
      if ( (targ.eq.'P').or.(targ.eq.'D').or.(targ.eq.'3') ) then
         A = (paramp(1)+paramp(2)*datan(paramp(3)*paramp(3)*Q2))
         B = (paramp(4) + paramp(5)*datan(paramp(6)*paramp(6)*Q2))
         C = dsin(py* XI**paramp(7))
         
         A1p =  XI**A * (1.0 + B* C)
         
         deriv(1) = dlog(XI)*A1p
         deriv(2) = deriv(1)*datan(paramp(3)*paramp(3)*Q2)
         deriv(3) = deriv(1)*2*paramp(2)*paramp(3)*Q2
     *        /(paramp(3)*paramp(3)*paramp(3)*paramp(3)*Q2*Q2 + 1)
         deriv(4) = XI**A * C
         deriv(5) = deriv(4)* datan(paramp(6)*paramp(6)*Q2)
         deriv(6) = deriv(4)*2*paramp(5)*paramp(6)*Q2
     *        /(paramp(6)*paramp(6)*paramp(6)*paramp(6)*Q2*Q2 + 1)
         deriv(7) = py * XI**paramp(7) * dcos(py*XI**paramp(7))* 
     *        XI**A *B*dlog(XI)
         
         delA1p = 0.0D0
         do j = 1, 7
            do k = 1, 7
               delA1p = delA1p + errparp(J,K)*deriv(J)*deriv(K)
            enddo
         enddo
      endif
      
      
      if ( (targ=='N').or.(targ=='D').or.(targ=='3') ) then
         A = (paramn(1)+paramn(2)*datan(paramn(3)*paramn(3)*Q2))
         B = (paramn(4) + paramn(5)*datan(paramn(6)*paramn(6)*Q2))
         C = dsin(py* XI**paramn(7))
c         C = dsin(py* XI**paramn(7))  !kp: 8/19/13: SEK found that (7) in paramn, which doesn't match with the NUMN = 6 value (corresponding to the # of parameters or columns in neutron-input-files such as disa1n_fit_error.dat, so this 7th parameter taking random value may explain why we have the unexpectedly steeply rising g1 (in the model) above W=2 such as shown in http://wwwold.jlab.org/Hall-B//secure/eg4/adhikari/Analysis/SimStuffs/Extracted_g1/g1_stdAndExtracted_StdD57_nStdD62C194S139w40BinsEb8.gif (black line is for the model).
         
         A1n =  XI**A * (1.0 + B* C)
         
         deriv(1) = dlog(XI)*A1n
         deriv(2) = deriv(1)*datan(paramn(3)*paramn(3)*Q2)
         deriv(3) = deriv(1)*2*paramn(2)*paramn(3)*Q2
     *        /(paramn(3)*paramn(3)*paramn(3)*paramn(3)*Q2*Q2 + 1)
         deriv(4) = XI**A * C
         deriv(5) = deriv(4)* datan(paramn(6)*paramn(6)*Q2)
         deriv(6) = deriv(4)*2*paramn(5)*paramn(6)*Q2
     *        /(paramn(6)*paramn(6)*paramn(6)*paramn(6)*Q2*Q2 + 1)
         deriv(7) = py * XI**paramn(7) * dcos(py*XI**paramn(7))* 
     *        XI**A *B*dlog(XI)
         
         delA1n = 0.0D0
         do j = 1, 7
            do k = 1, 7
               delA1n = delA1n + errparn(J,K)*deriv(J)*deriv(K)
            enddo
         enddo
      endif

      
      if (targ=='P') then
         
         A1 = A1p
         delA1 = dsqrt(delA1p)
         
         
      elseif (targ=='N') then
         A1 = A1n
         delA1 = dsqrt(delA1n)
         
      elseif (targ=='D') then
         suppress = .false.
         if ((M2 - Q2 + Q2/X) < 4.D0) suppress = .true. ! sek_me only called for DIS
         
         call F1F2new(X, Q2, 'P', F1p, dummy1, dummy2, suppress)
         call F1F2new(X, Q2, 'N', F1n, dummy1, dummy2, suppress)
         
         if ((F1p .eq. 0.D0).and.(F1n .eq. 0.D0)) then
            A1 = 0.D0
         else
            A1 = scale * (F1p*A1p + F1n*A1n) / (F1p+F1n)
            delA1 = scale/(F1p+F1n) *dsqrt(F1p*F1p*delA1p*delA1p
     >           +F1n*F1n*delA1n*delA1n)
         endif
         
      elseif (targ=='3') then
         suppress = .false.
         call F1F2new(X, Q2, '3', dummy1, dummy2, signp, suppress)
         A1 = (0.87D0*A1n*signp - 2.D0*0.027D0*A1p) / (signp + 2.D0)
	 delA1 = 0.87D0*delA1n*signp/(signp + 2.0D0)
      else
         
         write(6,'(''%%%%% illegal Target in SEK_ME: >'',
     >        A1, ''< !!!'')') targ
         
      endif
      
      return
      end
      

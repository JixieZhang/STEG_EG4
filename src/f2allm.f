*******************************************************************
*                                                                 *
*    NEW PARAMETERIZATION OF SIGMA_TOT AND F2                     *
*    BASED ON THE ALLM FUNCTIONAL FORM                            *
*                                                                 *
* Written by:  D. Gabbert and L. De Nardo                         *
* E-mail:      dgabbert@mail.desy.de, ldenardo@mail.desy.de       *
* Date:        April, 2007                                        *
* Ref: For publications please cite it as DESY report 07-107,
* hep-ph/0708.3196v1, DIS 2007.
*                                                                 *
*******************************************************************     
*                                                                 *
* USAGE                                                           *
* CALL SIGMATOT_PARAM(xbj,qsqr,wsqr,flag,sigtot,esigtot,f2,ef2)   *
*                                                                 *
* INPUT                                                           *
* xbj, qsqr, wsqr: kinematics                                     *
* flag             =1 (use only xbj and qsqr  as input)           * 
*                  =2 (use only qsqr and wsqr as input)           * 
*                                                                 *
* OUTPUT                                                          *
* sigtot, esigtot: sigma_tot and error on sigma_tot               *
* F2, eF2        : F2 and error on F2                             *
*                                                                 *
* COMMENTS                                                        *
* The flag=2 setting is necessary for Q2=0                        *
* F2 and eF2 will be meaningful only for qsqr>0,                  *
* for qsqr=0 values for F2 and eF2 will be set to -999.           *
*                                                                 *
* More details can be found in sigmatot.ps                        *
* and the DIS 2007 proceedings                                    *
*                                                                 *
*******************************************************************     

      SUBROUTINE SIGMATOT_PARAM(XBJ,QSQR,WSQR,FLAG,
     +                             SIGTOT,ESIGTOT,F2,EF2)


      IMPLICIT NONE

      DOUBLE PRECISION PAR(23)
      DOUBLE PRECISION EPAR(23)
      DOUBLE PRECISION COV(23,23)
      DOUBLE PRECISION DADP(23)
      SAVE COV,PAR,EPAR

*     VARIABLES

      DOUBLE PRECISION XBJ,QSQR,WSQR,SIGTOT,ESIGTOT,F2,EF2     
      INTEGER FLAG
      DOUBLE PRECISION A, DA

    

*     CONSTANTS     
      DOUBLE PRECISION PI,E,MP,ALPHA
      DOUBLE PRECISION HBARC, MPBARN
      PARAMETER (PI=3.14159)      
      PARAMETER (E=2.7182818)
      PARAMETER (MP=0.9382723)
      PARAMETER (ALPHA=0.00729735257)
*     HBARC = H BAR * C IN GEV*M
      PARAMETER (HBARC=1.97327E-16)
*     MPBARN= M^2 PER BARN
      PARAMETER (MPBARN=1.E-28)

*     COVARIANCE MATRIX ELEMENTS
      INTEGER I,J

*     QUANTITIES FOR THE FUNCTIONAL FORM
      DOUBLE PRECISION T, XP, XR
      DOUBLE PRECISION AP, BP, CP 
      DOUBLE PRECISION AR, BR, CR 
      DOUBLE PRECISION FP, FR
      DOUBLE PRECISION NU

*     VARIABLES FOR READING FILES
      INTEGER IOFLAG, IOSTAT
      INTEGER NEXTUN
      INTEGER LUN1,LUN2
      CHARACTER*120 LINE
      LOGICAL FIRSTCALL
      DATA FIRSTCALL /.TRUE./
      SAVE FIRSTCALL



*     LOAD PARAMETERS AND COVARIANCE MATRIX FROM FILE

      IF (FIRSTCALL) THEN

        OPEN (UNIT=91,FILE='F2ALLM07.inp1')
        OPEN (UNIT=92,FILE='F2ALLM07.inp2')

        do i=1,23
          READ(91,*) PAR(I),EPAR(I)
          DO J=1,23
            READ(92,*) COV(I,J)
          END DO
csk          write(6,'(2f10.3)') par(i),epar(i)
        END DO
        CLOSE (91)
        CLOSE (92)
csk        WRITE(6,*) 'FIT INFORMATION LOADED FROM FILE.'
      ENDIF

*     KINEMATIC VARIABLES
      NU   = (WSQR-MP*MP+QSQR)/(2*MP)
      IF (FLAG.EQ.2) THEN 
        XBJ  = QSQR/(2*MP*NU)
      ELSE
        WSQR=(1-XBJ)/XBJ*QSQR+MP**2
      ENDIF

*     CALCULATION OF FUNCTIONS NEEDED

      T  = LOG(LOG( (QSQR+PAR(4))/PAR(5) )/LOG( PAR(4)/PAR(5) ) )

      XP = (1. + (WSQR - MP*MP)/(QSQR + PAR(2)) )	 
      XR = (1. + (WSQR - MP*MP)/(QSQR + PAR(3)) )	 

      AP = PAR(6) + ( PAR(6) - PAR(7) )*( 1./(1+T**(PAR(8)))  - 1. ) 
      BP = PAR(9) + PAR(10)*T**(PAR(11))
      CP = PAR(12)+ ( PAR(12)- PAR(13))*( 1./(1+T**(PAR(14))) - 1. )

      AR = PAR(15) + PAR(16)*T**(PAR(17))
      BR = PAR(18) + PAR(19)*T**(PAR(20))
      CR = PAR(21) + PAR(22)*T**(PAR(23))

      FP = CP * XP**(-AP)*(1-XBJ)**BP
      FR = CR * XR**(-AR)*(1-XBJ)**BR

      IF (QSQR.NE.0.) THEN

*     ERROR PROPAGATION FOR DIS REGION

          A=1./(QSQR+PAR(1))*(FP+FR)

          DADP(1) =-1./(QSQR+PAR(1))*A
          DADP(2) =1./(QSQR+PAR(1))*CP*(1-XBJ)**BP*AP*XP**(-AP-1)*
     +             (WSQR-MP*MP)/(QSQR+PAR(2))**2
          DADP(3) =1./(QSQR+PAR(1))*CR*(1-XBJ)**BR*AR*XR**(-AR-1)*
     +             (WSQR-MP*MP)/(QSQR+PAR(3))**2
*
          DADP(4) =-XP**(-AP)*(1-XBJ)**BP*
     +        (PAR(12)-PAR(13))*PAR(14)*T**(PAR(14)-1)/(1+T**PAR(14))**2
          DADP(4) =DADP(4)-CP*XP**(-AP)*(1-XBJ)**BP*LOG(1./XP)
     +          *(PAR(6)-PAR(7))*PAR(8)*T**(PAR(8)-1)/(1+T**PAR(8))**2
          DADP(4) =DADP(4)+CP*XP**(-AP)*(1-XBJ)**BP*LOG(1-XBJ)*
     +             PAR(10)*PAR(11)*T**(PAR(11)-1)
          DADP(4) =DADP(4) +XR**(-AR)*(1-XBJ)**BR*PAR(22)*PAR(23)*
     +             T**(PAR(23)-1)
          DADP(4) =DADP(4) +CR*XR**(-AR)*(1-XBJ)**BR*LOG(1./XR)*
     +             PAR(16)*PAR(17)*T**(PAR(17)-1)
          DADP(4) =DADP(4)+CR*XR**(-AR)*(1-XBJ)**BR*LOG(1-XBJ)*
     +             PAR(19)*PAR(20)*T**(PAR(20)-1)

          DADP(4) =DADP(4)*1./(QSQR+PAR(1))*
     +             (1./(QSQR+PAR(4))*1./LOG((QSQR+PAR(4))/PAR(5))
     +             -1./PAR(4)/LOG(PAR(4)/PAR(5)))
*
          DADP(5) =-XP**(-AP)*(1-XBJ)**BP*
     +        (PAR(12)-PAR(13))*PAR(14)*T**(PAR(14)-1)/(1+T**PAR(14))**2
          DADP(5) =DADP(5)-CP*XP**(-AP)*(1-XBJ)**BP*LOG(1./XP)
     +            *(PAR(6)-PAR(7))*PAR(8)*T**(PAR(8)-1)/(1+T**PAR(8))**2
          DADP(5) =DADP(5)+CP*XP**(-AP)*(1-XBJ)**BP*LOG(1-XBJ)*
     +             PAR(10)*PAR(11)*T**(PAR(11)-1)
          DADP(5) =DADP(5) +XR**(-AR)*(1-XBJ)**BR*PAR(22)*PAR(23)*
     +             T**(PAR(23)-1)
          DADP(5) =DADP(5) +CR*XR**(-AR)*(1-XBJ)**BR*LOG(1./XR)*
     +             PAR(16)*PAR(17)*T**(PAR(17)-1)
          DADP(5) =DADP(5)+CR*XR**(-AR)*(1-XBJ)**BR*LOG(1-XBJ)*
     +             PAR(19)*PAR(20)*T**(PAR(20)-1)
          DADP(5) =DADP(5)*1./(QSQR+PAR(1))*1./PAR(5)*
     +             (1./LOG(PAR(4)/PAR(5))
     +             -1./LOG((QSQR+PAR(4))/PAR(5)))

          DADP(6) =1./(QSQR+PAR(1))*CP*(1-XBJ)**BP*XP**(-AP)*
     +             1./(1+T**PAR(8))*LOG(XP)
          DADP(7) =1./(QSQR+PAR(1))*CP*(1-XBJ)**BP*XP**(-AP)*
     +             (1-1./(1+T**PAR(8)))*LOG(1./XP)
          DADP(8) =-1./(QSQR+PAR(1))*CP*(1-XBJ)**BP*XP**(-AP)*
     +             ((PAR(6)-PAR(7))/(1+T**PAR(8))**2)*
     +             T**PAR(8)*LOG(T)*LOG(1./XP)
          DADP(9) =1./(QSQR+PAR(1))*CP*XP**(-AP)*(1-XBJ)**BP*LOG(1-XBJ)
          DADP(10)=1./(QSQR+PAR(1))*CP*XP**(-AP)*(1-XBJ)**BP*LOG(1-XBJ)*
     +             T**PAR(11)
          DADP(11)=1./(QSQR+PAR(1))*CP*XP**(-AP)*(1-XBJ)**BP*LOG(1-XBJ)*
     +             PAR(10)*T**PAR(11)*LOG(T)
          DADP(12)=1./(QSQR+PAR(1))*XP**(-AP)*(1-XBJ)**BP/(1+T**PAR(14))
          DADP(13)=1./(QSQR+PAR(1))*XP**(-AP)*(1-XBJ)**BP*
     +             (1-1./(1+T**PAR(14)))
          DADP(14)=-1./(QSQR+PAR(1))*XP**(-AP)*(1-XBJ)**BP*
     +             (PAR(12)-PAR(13))/(1+T**PAR(14))**2*T**PAR(14)*LOG(T)
          DADP(15)=1./(QSQR+PAR(1))*CR*(1-XBJ)**BR*XR**(-AR)*LOG(1./XR)
          DADP(16)=1./(QSQR+PAR(1))*CR*(1-XBJ)**BR*
     +              XR**(-AR)*LOG(1./XR)*T**PAR(17)
          DADP(17)=1./(QSQR+PAR(1))*CR*(1-XBJ)**BR*
     +              XR**(-AR)*LOG(1./XR)*PAR(16)*T**PAR(17)*LOG(T)
          DADP(18)=1./(QSQR+PAR(1))*CR*XR**(-AR)*(1-XBJ)**BR*
     +              LOG(1-XBJ)
          DADP(19)=1./(QSQR+PAR(1))*CR*XR**(-AR)*(1-XBJ)**BR*
     +              LOG(1-XBJ)*T**PAR(20)
          DADP(20)=1./(QSQR+PAR(1))*CR*XR**(-AR)*(1-XBJ)**BR*
     +              LOG(1-XBJ)*PAR(19)*T**PAR(20)*LOG(T)
          DADP(21)=1./(QSQR+PAR(1))*XR**(-AR)*(1-XBJ)**BR
          DADP(22)=1./(QSQR+PAR(1))*XR**(-AR)*(1-XBJ)**BR*T**PAR(23)
          DADP(23)=1./(QSQR+PAR(1))*XR**(-AR)*(1-XBJ)**BR*PAR(22)*
     +             T**PAR(23)*LOG(T)
*
          DA=0.
          DO I=1,23
            DO J=1,23
              DA=DA+COV(I,J)*DADP(I)*DADP(J)
            ENDDO
          ENDDO
          DA=SQRT(DA)

          SIGTOT = A* 
     +             4*(PI**2)*ALPHA/QSQR/(1-XBJ)*(QSQR+4*MP*MP*XBJ*XBJ)
     +             /QSQR*QSQR *HBARC*HBARC/MPBARN
          ESIGTOT= DA* 
     +             4*(PI**2)*ALPHA/QSQR/(1-XBJ)*(QSQR+4*MP*MP*XBJ*XBJ)
     +             /QSQR*QSQR *HBARC*HBARC/MPBARN

          F2     = A *QSQR
          EF2    = DA*QSQR

      ELSE

*     ERROR PROPAGATION FOR PHOTOPRODUCTION REGION (Q^2=0)
        
          A = 1./PAR(1)*
     +     (PAR(12)*(1+(WSQR-MP**2)/PAR(2))**(-PAR(6)) +
     +     PAR(21)*(1+(WSQR-MP**2)/PAR(3))**(-PAR(15)))

          DADP(1) =-1./PAR(1)*A
          DADP(2) = 1./PAR(1)*PAR(12)*PAR(6)*
     +             (1+(WSQR-MP**2)/PAR(2))**(-PAR(6)-1)*
     +             (WSQR-MP**2)/PAR(2)**2.
          DADP(3)   = 1./PAR(1)*PAR(21)*PAR(15)*
     +    (1+(WSQR-MP**2)/PAR(3))**(-PAR(15)-1)*(WSQR-MP**2)/PAR(3)**2
          DADP(6)   =-1./PAR(1)*PAR(12)*
     +    (1+(WSQR-MP**2)/PAR(2))**(-PAR(6))*LOG(1+(WSQR-MP**2)/PAR(2))
          DADP(12)  = 1./PAR(1)*
     +    (1+(WSQR-MP**2)/PAR(2))**(-PAR(6))
          DADP(15)  =-1./PAR(1)*PAR(21)*
     +    (1+(WSQR-MP**2)/PAR(3))**(-PAR(15))*LOG(1+(WSQR-MP**2)/PAR(3))
          DADP(21)  = 1./PAR(1)*
     +    (1+(WSQR-MP**2)/PAR(3))**(-PAR(15))
          DADP(4)=0.
          DADP(5)=0.
          DADP(7)=0.
          DADP(8)=0.
          DADP(9)=0.
          DADP(10)=0.
          DADP(11)=0.
          DADP(13)=0.
          DADP(14)=0.
          DADP(16)=0.
          DADP(17)=0.
          DADP(18)=0.
          DADP(19)=0.
          DADP(20)=0.
          DADP(22)=0.
          DADP(23)=0.

          DA=0.
          DO I=1,23
            DO J=1,23
               DA=DA+COV(I,J)*DADP(I)*DADP(J)
            ENDDO
          ENDDO
          DA=SQRT(DA)

          SIGTOT = A*(4*PI**2*ALPHA)*(HBARC**2/MPBARN)
          ESIGTOT=DA*(4*PI**2*ALPHA)*(HBARC**2/MPBARN)
          F2     = -999.
          EF2    = -999.

      ENDIF

      IF (FIRSTCALL) FIRSTCALL=.FALSE.

      END

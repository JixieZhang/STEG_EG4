*******************************************************************************
* RADVAR.INC : Include file of useful variables for radiative corrections 
* Which vary due to changes in integration variables.
*
* 2/93, LMS.
******************************************************************************
       REAL*8  T,                         ! Integration variable.
     >         MX2,                       ! Integration variable.
     >         TT,                        ! TT = T + MX2 - M2.
     >         TD,                        ! T - Y.
     >         U,                         ! SX - TT
     >         BZ, CZ, BZHAT, CZHAT,
     >         I(-1:2), IHAT(-1:2),
     >         A(2,6), AHAT(2,6),
     >         AA(0:2,3), AAHAT(0:2,3),
     >         TPOS, TNEG, ID, AD(6), AAD(6),
     >         FIR, TM, YD, A0(2), A0HAT(2), A0D

       COMMON /RADVAR/T, MX2, TT, TD, U, BZ, CZ, BZHAT, CZHAT,
     >                I,IHAT, A, AHAT, AA, AAHAT, TPOS, TNEG, ID,
     >                AD, AAD, FIR, TM, YD, A0, A0HAT, A0D


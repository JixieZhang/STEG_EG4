*******************************************************************************

        SUBROUTINE DSFs(X, Q2, YoniIndex, F1, F2, R, G1, G2, A1, A2, suppress)
*******************************************************************************
* This subroutine is used for APPLYING A1n to deuterium, using WM/YK folding
* Calls special version of Yoni.f
* 
******************************************************************************
       IMPLICIT NONE

	include 'instruct.inc'

 	REAL*8 a1calc(3840)
        REAL*8 X, Q2, W, G1n, G1p, G2n, G2p, F1n, F1p, F2n, F2p
	real*8  F1, F2, R, G1, G2, A1, A2, Nu, QoverNu
        CHARACTER*1 TARG, TARGP /'P'/, TARGN /'N'/
	integer YoniIndex
	integer pSFChoice, nSFChoice, pAsymChoice, nAsymChoice
	integer pIPOL, nIPOL, pIA1, nIA1
	logical suppress
		
csk	common /YoniDeut/ pSFChoice, nSFChoice, pAsymChoice, nAsymChoice, 
csk     *		pIPOL, nIPOL, pIA1, nIA1
	
	TARG = TARGP
csk	SFChoice = pSFChoice
csk	AsymChoice = pAsymChoice
csk	IPOL = pIPOL
csk	IA1 = pIA1
	call Yoni(X, Q2, TARG, YoniIndex, F1p, F2p, R, G1p, G2p, A1, A2, suppress)
	TARG = TARGN
csk	SFChoice = nSFChoice
csk	AsymChoice = nAsymChoice
csk	IPOL = nIPOL
csk	IA1 = nIA1
	call Yoni(X, Q2, TARG, YoniIndex, F1n, F2n, R, G1n, G2n, A1, A2, suppress)
	F1 = F1p + F1n
	F2 = F2p + F2n
	G1 = G1p + G1n
	G2 = G2p + G2n
	if (F1 .gt. 0. .and. Q2 .gt. 0.) then
	  QoverNu = 4.0D0*0.939D0*0.939D0*X*X/Q2
	  R = F2/F1/2.0/x*(1+QoverNu)-1.0
	  A1 = (G1 -QoverNu*G2)/F1
	  QoverNu = dsqrt(QoverNu)
	  A2 = QoverNu*(G1 + G2)/F1
	else
	  A1 = 0.0
	  A2 = 0.0
	  R = 0.0
	endif
	return
	end
	

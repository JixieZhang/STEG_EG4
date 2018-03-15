	subroutine Yoni(x, Q2, T, YoniIndex, F1, F2, R, g1, g2, A1, A2, suppress)
c
c	Implementation of Wally Melnitchouk/Yoni Kahn smearing to calculate F2d and g1d 
c	Created: S.E. Kuhn 30-Aug-2009
c
	
	implicit none
	integer I, J, K, L, YoniIndex, IG, IX, IXel, IYel, IQ2
	logical FirstCall /.TRUE./
	
	real*8 Pi /3.1415926536D0/
	real*8 Q2,Nu,x,R,F2,F1,FL,qsq,W2,DR,xstarel,yel,mstarsq,tauel,yeleq
	real*8 Rloop, F2loop, F1loop, G1loop, G2loop, A1loop, A2loop, FLloop
	real*8 TH, g1, g2, A1, A2, GEP,GEN,GMP,GMN
	real*8 xOVy,dF00sum,dF90sum,dF99sum,dF11sum,dF12sum,dF21sum,dF22sum,xthr,ythr
	real*8 gamma,y, rupper, rlower, tau, MassN /0.939D0/
	real*8 centroid(21,19), upedge(21,19), sum(19), censum(19), pol90sum(21,19), pol99sum(21,19),
     *	pol11sum(21,19), pol12sum(21,19), pol21sum(21,19), pol22sum(21,19), 
     *	f00in(21,1000), f11in(21,1000),  f12in(21,1000), f21in(21,1000), 
     *  f22in(21,1000),
     *	f90in(21,1000), f99in(21,1000),
     *		f00sum(21,1000), f11sum(21,1000), f12sum(21,1000), f21sum(21,1000), f22sum(21,1000), 
     *		f90sum(21,1000), f99sum(21,1000),
     *		tot(21), totpol11(21), yelvsx(40,31), xforyel(31), xtest
     
	data sum /0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.27, 0.35, 0.45,
     *		0.55,0.65, 0.73, 0.8, 0.85, 0.9, 0.93, 0.96, 0.98, 1.0/

	character*1 T	!Target Type
	logical dummy, suppress, higamma	!, newhallc
	save FirstCall, centroid, upedge, sum, censum,
     *		pol11sum, pol12sum, pol21sum, pol22sum, pol90sum, pol99sum,
     *		f00in, f11in,  f12in, f21in, f22in, f90in, f99in,
     *		f00sum, f11sum,  f12sum, f21sum, f22sum, f90sum, f99sum,
     *		tot, totpol11, yelvsx, xforyel

C******************* ! 	INITIALIZATION ! *****************************
	IF(FirstCall)THEN
	write(66,*) ' Yoni first called!'
	FirstCall = .FALSE.
	censum(1) = sum(1)/2.0D0
	do I = 2,19
	  censum(I) = (sum(I)+sum(I-1))/2.0D0
	enddo
	
	open(unit=62, file='f99.ful', status='OLD')
	open(unit=63, file='f90.ful', status='OLD')
	open(unit=64, file='f00.ful', status='OLD')
	open(unit=65, file='f11.ful', status='OLD')
	open(unit=67, file='f12.ful', status='OLD')
	open(unit=68, file='f21.ful', status='OLD')
	open(unit=69, file='f22.ful', status='OLD')
	open(unit=70, file='yelas.dat', status='OLD')
	open(unit=66, file='Yoni.out', status='NEW')
	
	do I = 1,21
	  tot(I) = 0.0
	  totpol11(I) = 0.0
	  do J = 1, 1000
	    read(62,*,end=999) gamma, y, f99in(I,J)
	      if((abs(gamma*10.0D0 - 9.0D0 - I).gt.0.01D0).or.(abs(y*5.0D2-J).gt.0.01))then
	        write(66,*) ' Wrong 99 input: ',I,gamma,J,y
	      endif
	      if (y .lt. 0.3 .or. (y/gamma) .ge. 1.5) f99in(I,J) = 0.0D0
	      
	    read(63,*,end=999) gamma, y, f90in(I,J)
	      if((abs(gamma*10.0D0 - 9.0D0 - I).gt.0.01D0).or.(abs(y*5.0D2-J).gt.0.01))then
	        write(66,*) ' Wrong 90 input: ',I,gamma,J,y
	      endif
	      if (y .lt. 0.3 .or. (y/gamma) .ge. 1.5) f90in(I,J) = 0.0D0
	  
	    read(64,*,end=999) gamma, y, f00in(I,J)
c	        write(6,*) '  00 input: ',I,gamma,J,y, f00in(I,J)
	      if((abs(gamma*10.0D0 - 9.0D0 - I).gt.0.01D0).or.(abs(y*5.0D2-J).gt.0.01))then
	        write(66,*) ' Wrong 00 input: ',I,gamma,J,y
	      endif
	      if (y .lt. 0.3 .or. (y/gamma) .ge. 1.5) f00in(I,J) = 0.0D0
	      tot(I) = tot(I) + f00in(I,J)*0.002
	    
	    read(65,*,end=999) gamma, y, f11in(I,J)
	      if((abs(gamma*10.0D0 - 9.0D0 - I).gt.0.01D0).or.(abs(y*5.0D2-J).gt.0.01))then
	        write(66,*) ' Wrong 11 input: ',I,gamma,J,y
	      endif
	      if (y .lt. 0.3 .or. (y/gamma) .ge. 1.5) f11in(I,J) = 0.0D0
	      totpol11(I) = totpol11(I) + f11in(I,J)*0.002
	    
	    read(67,*,end=999) gamma, y, f12in(I,J)
	      if((abs(gamma*10.0D0 - 9.0D0 - I).gt.0.01D0).or.(abs(y*5.0D2-J).gt.0.01))then
	        write(66,*) ' Wrong 12 input: ',I,gamma,J,y
	      endif
	      if (y .lt. 0.3 .or. (y/gamma) .ge. 1.5) f12in(I,J) = 0.0D0

	    read(68,*,end=999) gamma, y, f21in(I,J)
	      if((abs(gamma*10.0D0 - 9.0D0 - I).gt.0.01D0).or.(abs(y*5.0D2-J).gt.0.01))then
	        write(66,*) ' Wrong 21 input: ',I,gamma,J,y
	      endif
	      if (y .lt. 0.3 .or. (y/gamma) .ge. 1.5) f21in(I,J) = 0.0D0

	    read(69,*,end=999) gamma, y, f22in(I,J)
	      if((abs(gamma*10.0D0 - 9.0D0 - I).gt.0.01D0).or.(abs(y*5.0D2-J).gt.0.01))then
	        write(66,*) ' Wrong 22 input: ',I,gamma,J,y
	      endif
	      if (y .lt. 0.3 .or. (y/gamma) .ge. 1.5) f22in(I,J) = 0.0D0

	  enddo
	  write(66,*) ' Pol. Integral for I=',I,' is ',totpol11(I)
	  
	  f00sum(I,1)=f00in(I,1)*0.001/tot(I)
	  f90sum(I,1)=f90in(I,1)*0.001
	  f99sum(I,1)=f99in(I,1)*0.001
	  f11sum(I,1)=f11in(I,1)*0.001
	  f12sum(I,1)=f12in(I,1)*0.001
	  f21sum(I,1)=f21in(I,1)*0.001
	  f22sum(I,1)=f22in(I,1)*0.001

	  do J=2,1000
	    f00sum(I,J) = f00sum(I,J-1)+(f00in(I,J-1)+f00in(I,J))*0.001/tot(I)
	    f99sum(I,J) = f99sum(I,J-1)+(f99in(I,J-1)+f99in(I,J))*0.001
	    f12sum(I,J) = f90sum(I,J-1)+(f90in(I,J-1)+f90in(I,J))*0.001
	    f11sum(I,J) = f11sum(I,J-1)+(f11in(I,J-1)+f11in(I,J))*0.001
	    f12sum(I,J) = f12sum(I,J-1)+(f12in(I,J-1)+f12in(I,J))*0.001
	    f21sum(I,J) = f21sum(I,J-1)+(f21in(I,J-1)+f21in(I,J))*0.001
	    f22sum(I,J) = f22sum(I,J-1)+(f22in(I,J-1)+f22in(I,J))*0.001
	  enddo
	  ! Running total
	  
c	  write(6,*) 'Integral for gamma=',((I+9.0)/10.0),': f00 ',f00sum(I,1000),
c     *		' f11 ',f11sum(I,1000)
      	enddo
	
	close(2)
	close(3)	
	close(4)
	close(5)
	close(7)
	close(8)
	close(9)
	
	do I = 1,21
	  K = 2
	  do J = 1,19
11	    continue
	    if(censum(J).gt.f00sum(I,K))then
	      K = K+1
	      if(K.gt.999)then
	      	write(66,*) 'Cant find centroid ',I,J,K,censum(J),f00sum(I,K)
		return
	      endif
	      goto 11
	    else
	      rupper = f00sum(I,K)-censum(J)
	      rlower = censum(J) - f00sum(I,K-1)
	      centroid(I,J) = ((K-1)*rlower + (K-2)*rupper)*0.002/(rupper+rlower)
22	      continue
	      if(J.lt.19)then
	        if(sum(J).gt.f00sum(I,K))then
	          K = K+1
	          if(K.gt.999)then
	      	    write(66,*) ' Cant find up end! ',I,J,K,sum(J),f00sum(I,K)
		    return
		  endif
		  goto 22
	        else
	      	  rupper = f00sum(I,K)-sum(J)
	      	  rlower = sum(J) - f00sum(I,K-1)
		  upedge(I,J) = ((K-1)*rlower + (K-2)*rupper)*0.002/(rupper+rlower)
		  pol90sum(I,J) = (f90sum(I,K)*rlower + f90sum(I,K-1)*rupper)/(rupper+rlower)
		  pol99sum(I,J) = (f99sum(I,K)*rlower + f99sum(I,K-1)*rupper)/(rupper+rlower)
		  pol11sum(I,J) = (f11sum(I,K)*rlower + f11sum(I,K-1)*rupper)/(rupper+rlower)
		  pol12sum(I,J) = (f12sum(I,K)*rlower + f12sum(I,K-1)*rupper)/(rupper+rlower)
		  pol21sum(I,J) = (f21sum(I,K)*rlower + f21sum(I,K-1)*rupper)/(rupper+rlower)
		  pol22sum(I,J) = (f22sum(I,K)*rlower + f22sum(I,K-1)*rupper)/(rupper+rlower)
	        endif
	      else
	        upedge(I,J) = 2.0D0
	        pol90sum(I,J) = f90sum(I,1000)
	        pol99sum(I,J) = f99sum(I,1000)
	        pol11sum(I,J) = f11sum(I,1000)
	        pol12sum(I,J) = f12sum(I,1000)
	        pol21sum(I,J) = f21sum(I,1000)
	        pol22sum(I,J) = f22sum(I,1000)
	      endif
	      ! Total weight integrated up to upper end of y-interval 1...19
	    endif
	  enddo
	enddo
	
	read(70,*,end=999) (xforyel(I), I=1,31)
	do J = 1, 40
	  read(70,*,end=999) (yelvsx(J,I), I=1,31)
	enddo
	close(10)
	
	ENDIF
	
C******************* ! 	 END INITIALIZATION ! *****************************

	F1 = 0.0D0
	F2 = 0.0D0
	FL = 0.0D0
	g1 = 0.0D0
	g2 = 0.0D0
	A1 = 1.0D0
	A2 = 0.0D0
	R = 0.0

	IX = x/0.002D0
	if(x .le. 0.0D0 .or. IX .gt. 999) return
	tau = Q2/4.0D0/MassN/MassN/x/x !
	gamma = dsqrt(1.0D0 + 1/tau) !
	J = (gamma - 0.85D0)*10.0D0
	if(J .lt. 1) return
	higamma = (J .gt. 21)
	if(higamma) J = 21
	
	IQ2 = 27+13.0*(dlog(Q2/0.91883882)/dlog(10.0D0))
	if (IQ2 .lt. 1) IQ2 = 1
	if (IQ2 .gt. 40) IQ2 = 40
	xtest = (x - 1.0/30.0)*15.0 + 2.0
	IXel = xtest
	if(IXel .lt. 2) then
	  IXel = 1
	  xtest = x*30.0 + 1.0
	else if(IXel .gt. 30) then
	  IXel = 30
	  xtest = IXel + 0.5
	endif
	rupper = xtest - IXel
	rlower = 1.0 - rupper
	yel = yelvsx(IQ2,IXel)*rlower + yelvsx(IQ2,(IXel+1))*rupper

	if(YoniIndex .gt. 2) goto 1789
	if(yel .lt. 0.0001)goto 1789
	xstarel = x/yel
	if(xstarel .lt. 0.01)goto 1789
	mstarsq = 0.939*0.939 + (1.0 - 1.0/xstarel)*Q2
	if(mstarsq .lt. 0.1D0) goto 1789
	tauel = Q2/4.0D0/mstarsq/xstarel/xstarel
cccc	tauel = tau
	if(higamma) then
	  yeleq = (yel - 1.0D0)*3.0D0/gamma + 1.0D0
	  IYel = yeleq/0.002D0
	  rupper = (IYel+1) - yeleq/0.002D0
	else
	  IYel = yel/0.002D0
	  rupper = (IYel+1) - yel/0.002D0
	endif
	if(IYel .lt. 1 .or. IYel .gt. 999) goto 1789
	rlower = 1.0D0-rupper
	
	IG = 22
	call NewFORM(IG,Q2,GEP,GEN,GMP,GMN)
	F2 = (rupper*f00in(J,IYel)+rlower*f00in(J,IYel+1))
	FLloop = (rupper*f99in(J,IYel)+rlower*f99in(J,IYel+1))
	F2loop = (rupper*f90in(J,IYel)+rlower*f90in(J,IYel+1))
	if(higamma) then
	  F2 = F2*3.0D0/gamma
	  FLloop = FLloop*3.0D0/gamma*yel*yel/yeleq/yeleq
	  F2loop = F2loop*gamma/3.0D0
	endif
	if(T .eq. 'P')then
	  FL = yel*(F2loop*(GEP*GEP + tauel*GMP*GMP)/(1.0D0+tauel) +
     *		FLloop*GEP*GEP/tauel)
	  F2 = yel*F2*(GEP*GEP + tauel*GMP*GMP)/(1.0D0+tauel) !
	else
	  F2 = yel*F2*(GEN*GEN + tauel*GMN*GMN)/(1.0D0+tauel) !
	  FL = yel*(F2loop*(GEN*GEN + tauel*GMN*GMN)/(1.0D0+tauel) +
     *		FLloop*GEN*GEN/tauel)
	endif
	IG = 21
	call NewFORM(IG,Q2,GEP,GEN,GMP,GMN)
	G1loop = (rupper*f11in(J,IYel)+rlower*f11in(J,IYel+1))
	G2loop = (rupper*f12in(J,IYel)+rlower*f12in(J,IYel+1))
	if(higamma) then
	  G1loop = G1loop*3.0D0/gamma
c	  G2loop = G2loop
	endif
	if(T .eq. 'P')then
	  g1 = 0.5D0/xstarel/(1.0D0+tauel)*(G1loop*(GEP*GMP + tauel*GMP*GMP) 
     *		+ G2loop*tauel*(GEP*GMP - GMP*GMP))
	else
	  g1 = 0.5D0/xstarel/(1.0D0+tauel)*(G1loop*(GEN*GMN + tauel*GMN*GMN) 
     *		+ G2loop*tauel*(GEN*GMN - GMN*GMN))	  
	endif
	G1loop = (rupper*f21in(J,IYel)+rlower*f21in(J,IYel+1))
	G2loop = (rupper*f22in(J,IYel)+rlower*f22in(J,IYel+1))
	if(higamma) then
	  G1loop = G1loop*9.0D0/gamma/gamma
	  G2loop = G2loop*3.0D0/gamma
	endif
	if(T .eq. 'P')then
	  g2 = 0.5D0/xstarel/(1.0D0+tauel)*(G1loop*(GEP*GMP + tauel*GMP*GMP)
     *		 + G2loop*tauel*(GEP*GMP - GMP*GMP))
	else
	  g2 = 0.5D0/xstarel/(1.0D0+tauel)*(G1loop*(GEN*GMN + tauel*GMN*GMN) 
     *		+ G2loop*tauel*(GEN*GMN - GMN*GMN))	  
	endif

1789	continue	
	if(YoniIndex .lt. 2 .or. higamma) goto 451 ! Only quasi-elastic contribution
	
	xthr = Q2/(1.08D0*1.08D0 - MassN*MassN + Q2)
	ythr = x/xthr
	
	do I = 1 , 19
	  y = centroid(J,I)
	  if(y.le.ythr) goto 119
	  xOVy = x/y
	  if ((xOVy .lt. 0.001).or.(xOVy .gt. 0.98)) goto 119
	  if(I.eq.1)then
	    dF00sum = sum(I)
	    dF90sum = pol90sum(J,I)
	    dF99sum = pol99sum(J,I)
	    dF11sum = pol11sum(J,I)
	    dF12sum = pol12sum(J,I)
	    dF21sum = pol21sum(J,I)
	    dF22sum = pol22sum(J,I)
	  else
	    dF00sum = sum(I) - sum(I-1)
	    dF90sum = pol90sum(J,I) - pol90sum(J,I-1)
	    dF99sum = pol99sum(J,I) - pol99sum(J,I-1)
	    dF11sum = pol11sum(J,I) - pol11sum(J,I-1)
	    dF12sum = pol12sum(J,I) - pol12sum(J,I-1)
	    dF21sum = pol21sum(J,I) - pol21sum(J,I-1)
	    dF22sum = pol22sum(J,I) - pol22sum(J,I-1)
	  endif
	  call F1F2new (xOVy, Q2, T, F1loop, F2loop, Rloop, suppress)
C	R is here only used as placeholder for sign/sigp
	  TH = 25.0 D00 ! Arbitrary, not really needed
	  call G1G2new(xOVy, Q2, TH, T, G1loop, G2loop, A1loop, A2loop)
	  FLloop = gamma*gamma*F2loop - 2.0D0*xOVy*F1loop
	  if (FLloop .le. 0.0) FLloop = 0.0
	  F2 = F2 + dF00sum*F2loop
	  FL = FL + dF90sum*F2loop + dF99sum*FLloop
	  g1 = g1 + dF11sum*G1loop/y + dF12sum*G2loop/y
	  g2 = g2 + dF21sum*G1loop/y + dF22sum*G2loop/y
119	  continue
	enddo
	
451	continue
	F1 = (gamma*gamma*F2 - FL)/2.0D0/x
	if (F1 .gt. 0.D0 .and. x .gt. 0.D0 .and. Q2 .gt. 0.D0) then
	  R = FL/F1/2.0D0/x
	  A1 = (g1 - (gamma*gamma - 1.0D0)*g2)/F1
	  A2 = dsqrt(gamma*gamma - 1.0D0)*(g1 + g2)/F1
	else
	  F1 = 0.0D0
	  R = 0.0D0
	  A1 = 0.0D0
	  A2 = 0.0D0
	endif

	return
999	write(66,*) ' Unexpected initialization error'
	return
	end

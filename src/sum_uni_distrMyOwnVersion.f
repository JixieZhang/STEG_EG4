      subroutine sum_uni_distr(xmin,xmax,dxx,nbins,fmax,nrm,xx,ireg)
      implicit none
      include 'binning.inc'
      include 'regen.inc'
      include 'csmap.inc'       !kp: 4/13/2012
      include 'nThnPnSum.inc'   !kp: 4/13/2012
      integer iBin,fBin,Diff,mBin,dbgr !kp: 4/13/3012: initial, final, difference and middle of the bins.

      real xmin(2),xmax(2),dxx(2),nrm,xx(2)
      real fmax(n_th_el,n_p_el)
      real norm,rran,weight,fmax_hypercube
      real norm_ix(n_th_el)
      real rannumkp
      integer nbins(2),ireg(2),ix,iy,iv,iw,sBin
      common /vcache/norm_ix

 
      nrm=fprobsumm(num1Dbins)  !kp: sum of the cs of all the bins in the map, so the last in the fprobsumm array.



      rannumkp=rran()
      weight=rannumkp*nrm
c      print*,'rannum,weight,nrm',rannumkp,weight,nrm
c      weight=rran()*nrm
cc      weight=ranlux(rannumkp,1)*nrm !kp:
      iBin=1
      fBin=num1Dbins
c      Diff=fBin-iBin+1            !kp: fBin-iBin+1
      Diff=fBin-iBin            !kp: fBin-iBin+1
      mBin=iBin+int(1.0*Diff/2)
c      do while(Diff>2)          !kp: Event Sampling 2 (came back again 2/13/13)
      do while(Diff>1)  !till 8/8/12 !kp: Normal sampling
         if(weight.gt.fprobsumm(mBin)) then
            iBin=mBin
         else
            fBin=mBin
         endif
c         Diff=fBin-iBin+1 
         Diff=fBin-iBin 
         mBin=iBin+int(1.0*Diff/2)
c         if(dbgr.lt.100) print*,iBin,fBin,Diff
         dbgr=thBins(fBin)
      enddo
      sBin=fBin                 !kp: Selected or Sampled bin
c      sBin=iBin                !kp: This one gave a shifted delta peak and also gave a tiny peak below 0.9 in the W spectrum.
      ireg(1)=thBins(sBin)
      ireg(2)=pBins(sBin)
c      print*,'ireg(1),ireg(1)',ireg(1),ireg(2)



 1111 continue
      xx(1)=xmin(1)+(float(ireg(1)-1)+rran())*dxx(1)
      xx(2)=xmin(2)+(float(ireg(2)-1)+rran())*dxx(2)
      return
      end

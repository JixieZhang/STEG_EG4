      subroutine read_map(norm,acc_tot,flag_acc)
      implicit none
      include 'binning.inc'
      include 'csmap.inc'
      include 'nThnPnSum.inc'   !kp: 4/13/2012
      real norm,acc_tot,fprobSumTmp
      integer i_th_el,i_p_el,i,jj
      logical flag_acc

c     kp: fprobSum(i_th_el,i_p_el) is the sum of fprob for all the bins counting from (0,0)th to the current 
c     kp: bin (i.e., (i_th_el,i_p_el)th bin) bin. (This is crucial for my own version of STEG event sampling)
      fprobSumTmp=0.e0
      jj=0                      !kp: 4/13/2012

      print*,'num1Dbins: ',num1Dbins
c Cross section maxima
c      open(32,file='cs_max.dat',status='old')
      open(32,file='cs_map.dat',status='old')     !kp: 4/2/2012
      do i_th_el=1,n_th_el
         do i_p_el=1,n_p_el
            read(32,36) cs_max(i_th_el,i_p_el)
c            print*,i_th_el,i_p_el
         enddo
      enddo
      read(32,36) norm
      close(32)
      
      print*,'readmap.f: dbg 9/18/12'


      if(flag_acc) then
         print*,'debug1'
c Acceptance of Fiducial Cuts
         open(33,file='acc_fc.dat',status='old')
         do i_th_el=1,n_th_el
            do i_p_el=1,n_p_el
               read(33,36) acc_map(i_th_el,i_p_el)
c Product of two distributions
               fprob(i_th_el,i_p_el)=cs_max(i_th_el,i_p_el)
     <              *acc_map(i_th_el,i_p_el)
               
               fprobSumTmp=fprobSumTmp+fprob(i_th_el,i_p_el) !kp: 4/13/12
c               fprobSum(i_th_el,i_p_el)=fprobSumTmp !kp: 4/13/12
               jj=jj+1
               thBins(jj)=i_th_el
               pBins(jj)=i_p_el
               fprobsumm(jj)=fprobSumTmp
            enddo
         enddo
         read(33,36) acc_tot
         close(33)
      else
         print*,'debug2'
         do i_th_el=1,n_th_el
            do i_p_el=1,n_p_el
               fprob(i_th_el,i_p_el)=cs_max(i_th_el,i_p_el)
c               print*,i_th_el,i_p_el,fprob(i_th_el,i_p_el)
               fprobSumTmp=fprobSumTmp+fprob(i_th_el,i_p_el) !kp: 4/13/12
c               fprobSum(i_th_el,i_p_el)=fprobSumTmp !kp: 4/13/12
               jj=jj+1
               thBins(jj)=i_th_el
               pBins(jj)=i_p_el
               fprobsumm(jj)=fprobSumTmp
            enddo
         enddo
      endif

      print*,'debug3'

36    format(1pe11.4)
      return
      end

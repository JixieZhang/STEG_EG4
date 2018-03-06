      subroutine write_map(norm)
      implicit none
      include 'binning.inc'
      include 'csmap.inc'
      real norm
      integer i_th_el,i_p_el



      open(67,file='cs_max_remapped.dat',status='unknown')
      do i_p_el=1,n_p_el
         do i_th_el=1,n_th_el
            if(i_th_el.eq.n_th_el) then
               write(67,36) cs_max(i_th_el,i_p_el)
            else
               write(67,35) cs_max(i_th_el,i_p_el)
            endif
         enddo
      enddo
      write(67,36) norm
      close(67)


 35   format(1pe11.4$)
 36   format(1pe11.4)
      return
      end

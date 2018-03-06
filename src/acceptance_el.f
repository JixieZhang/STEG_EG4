      subroutine acceptance_el(Ebeam,p,theta,t_current,delta_phi)
      implicit none
      real Ebeam,p,theta,t_current,accvb1,pm,thm,
     &delta_phi
      integer part_type
      data part_type/0/
      pm=p
      thm=theta
      delta_phi=accvb1(part_type,pm,thm,t_current)
      if(delta_phi.gt.30.e0) delta_phi=30.e0
      if(delta_phi.lt.0.e0) delta_phi=0.e0
      return
      end

      real function accvb1(part_type,p,theta,t_current)
      IMPLICIT NONE
      REAL p,theta,phi,t_current,acc
      REAL t_max
      REAL phi0_el,phi0_nh,phi0_ph
      REAL theta0_el,theta0_nh,theta0_ph
      REAL thetas_el,thetas_nh,thetas_ph
      REAL p_shift,cel_ex,pel_ex
      REAL ch_ex,theta_cut
      REAL theta_min,theta_max,delta_phi,zexp
      REAL d2r,pnorm,fcut,th_min,th_max
      INTEGER part_type,electron,pos_hadron,neg_hadron
c Enlarged cuts
      data phi0_el/30./
      data theta0_el/10.6/
      data thetas_el/15./
      data theta_max/52./
      data p_shift/0.25/
      data pel_ex/0.333/
      data cel_ex/0.25/
c Standard cuts
c      data phi0_el/30./
c      data theta0_el/12.2/
c      data thetas_el/21.5/
c      data theta_max/50./
c      data p_shift/0.15/
c      data pel_ex/0.416667/
c      data cel_ex/0.25/
c Init
      data t_max/3375./
      data theta_cut/75./
      data d2r/0.0174532925/
      
      Acc=0.0e+0
      
        theta_min = theta0_el+thetas_el/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.theta_max)then
          zexp = cel_ex*(p*t_max/t_current)**pel_ex
          delta_phi = phi0_el*sin((theta-theta_min)*d2r)**zexp
          Acc=delta_phi
        endif
c	print*,theta_min,delta_phi
      
      accvb1 = acc
      return
      end

      subroutine newphi(phi,phinew,sec)
      IMPLICIT NONE
      real phinew, phi
      integer sec
        phinew = phi
        sec=1
      if (phi.gt.330.e0) then
        phinew = phi-360.e0
        sec=1
      elseif (phi.gt.30.e0.and.phi.le.90.e0) then
        phinew = phi-60.e0
        sec=2
      elseif (phi.gt.90.e0.and.phi.le.150.e0) then
        phinew = phi-120.e0
        sec=3
      elseif (phi.gt.150.e0.and.phi.le.210.e0) then
        phinew = phi-180.e0
        sec=4
      elseif (phi.gt.210.e0.and.phi.le.270.e0) then
        phinew = phi-240.e0
        sec=5
      elseif (phi.gt.270.e0.and.phi.le.330.e0) then
        phinew = phi-300.e0
        sec=6
      endif
      return
      end

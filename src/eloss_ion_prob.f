*     Calculate most probable ionization energy loss by outgoing electron in target materials      
*     input:
*       p_el in GeV, theta_el in rad
*       istick=1|2  (1 for 0.5cm, 2 for1 cm target), NH3_Pf 
*     output: 
*       p_el after energy loss
      subroutine eloss_ion_prob(p_el,theta_el,istick,NH3_Pf)

      Implicit none
      
      integer istick
      real*8 p_el,theta_el,NH3_Pf

      real*8 ttarg(2),tiw_0,tiw_pf(2),tfw_0,tfw_pf(2)
      real*8 Ztarg,Ziw,Zfw
      real*8 t_al_kapton,t_lhe_0,t_lhe_pf(2)
      real*8 tfw_50al_blw,tfw_102al_25ka,tfw_si_mylar_cerex
      
      real*8 xc,R1,dR
      data xc/108.83d0/,R1/86.54d0/,dR/0.028d0/

c     these are temp variables, I do not like Implicit, so I put them here      
      real*8 tft_nh3,tfw_lhe,axi,tgt2,dx1,dx2,dt,tfw_al,xi(3),DeltaEb(3)
      
c--------------------------------------------------------------------------------
c     EG4 target: NH3 0.5/1 cm *Pf, LHe4 2.19-0.5/1.0+0.5/1.0*(1-Pf) cm, Al 71+50 um, Kapton 25 um
c     inlet window:
c     Al 121e-4 cm/X0=8.9 cm = 0.00136
c     Kapton 25e-4 cm/X0=28.6 cm = 0.00009
c     LHe4 (2.19-0.5*Pf)/2;(2.19-1.0*Pf)/2 cm/X0=650.5 cm = 0.00168-0.000385*Pf;0.00168-0.00077*Pf r.l.
      data tiw_0/0.00313d0/
      data tiw_pf/0.000385d0,0.00077d0/
      data Ziw/2.d0/
c     NH3 0.5;1*Pf cm/X0=47.17 cm = 0.01060*Pf;0.02120*Pf
      data ttarg/0.01060d0,0.02120d0/
      data Ztarg/4.d0/
c     exit window at 0 deg.
c     LHe4 (2.19-0.5*Pf)/2;(2.19-1.0*Pf)/2 cm/X0=650.5 cm = 0.00168-0.000385*Pf;0.00168-0.00077*Pf r.l.
c     Kapton 25e-4 cm/X0=28.6 cm = 0.00009
c     Al 71e-4+25e-4/X0=8.9 cm+6e-4(MYLAR) cm/X0=28.7 cm = 0.00107865+0.0000209=0.0011
c     20+10 layers of MYLAR+CEREX 10*(0.88e-3/39.95+2*1.e-3/45.4)=0.00022+0.00044=0.00066 at theta>0.166676 rad +2x at theta>0.1192777
      data tfw_0/0.00287d0/
      data tfw_pf/0.000385d0,0.00077d0/
      data Zfw/2.d0/
      data tfw_si_mylar_cerex/0.00066d0/
      data tfw_50al_blw/0.00056d0/ ! Al 50 um beamline exit window, theta <5 deg.
      data tfw_102al_25ka/0.001236d0/ ! Al 71+6+25=102 um and Kapton 25 um
c     t(NH3)=t(N)=0.7 cm/(X0=37.99g/cm^2/rho=15*nA+61.28g/cm^2/rho=3*1*nA)=0.0162 ==> tN/tH=7.5 ==> t(H)=0.0162/(1+7.5)=0.0019
c     Z(NH3)=(7*1.3646038+1*1.35707144)/(1.3646038+1.35707144)=4.008
c     t(He4)>>t(Al+Kapton)
      data t_al_kapton/0.00145d0/,t_lhe_0/0.00336d0/,t_lhe_pf/0.00077d0,0.00154d0/


      real*8 mp,me,mpi,alpha,pi,mn,igev2mub,ialpha,cl,rel,Na
      real*8 ZAl,AAl,X0Al,rhoAl,ZHe,AHe,X0LHe,rhoLHe,ZNH3,ANH3,X0NH3,rhoNH3
      DATA mp/0.938272029e0/,me/0.510998918e-03/,mpi/0.13957018e0/,
     &     alpha/7.2970e-03/,pi/3.14159265359e0/,mn/0.93956536e0/,
     &     igev2mub/389.379e0/,ialpha/137.035e0/,cl/29.9792458e0/,
     &     rel/2.8179403267e-13/,Na/6.022140857e23/,
     &     ZAl/13.e0/,          ! number of electrons in molecule A (used to get electron density=Z*Na/A*rho)
     &     AAl/26.98e0/,        ! molar mass in [g/mol]
     &     X0Al/24.011e0/,      ! X0 in [g/cm^2]
     &     rhoAl/2.7e0/,        ! density in [g/cm^3]
     &     ZHe/2.e0/,           ! number of electrons in molecule A (used to get electron density=Z*Na/A*rho)
     &     AHe/4.003e0/,        ! molar mass in [g/mol]
     &     X0LHe/94.322e0/,     ! X0 in [g/cm^2]
     &     rhoLHe/0.145e0/,     ! density in [g/cm^3]
     &     ZNH3/10.e0/,         ! number of electrons in molecule A (used to get electron density=Z*Na/A*rho)
     &     ANH3/18.034e0/,      ! molar mass in [g/mol]
     &     X0NH3/43.255e0/,     ! X0 in [g/cm^2]
     &     rhoNH3/0.917e0/      ! density in [g/cm^3]

c      print*,'**eloss_ion_prob(p_el,theta_el,istick,NH3_Pf) input: ',p_el,theta_el,istick,NH3_Pf
      
c     Most probable ionization energy loss by outgoing electron in target materials
      tft_nh3=real(ttarg(istick)*NH3_Pf/2.d0)/cos(theta_el)
      tfw_lhe=real((t_lhe_0-t_lhe_pf(istick)*NH3_Pf)/2.d0)/cos(theta_el)
      if(abs(theta_el).lt.0.0872664626) then ! <5 deg.
         tfw_al=real(tfw_102al_25ka+tfw_50al_blw)/cos(theta_el)
      else                      ! >5 deg.
         tgt2=1.d0+(dtan(dble(theta_el)))**2
         dx1=dsqrt(xc**2-tgt2*(xc**2-R1**2))
         dx2=dsqrt(xc**2-tgt2*(xc**2-(R1+dR)**2))
         dt=dabs(dx2-dx1)/dsqrt(tgt2)
         tfw_al=real(tfw_102al_25ka)/cos(theta_el)+real(dt/8.9d0)
         if(abs(theta_el).gt.0.166676) tfw_al=tfw_al+real(tfw_si_mylar_cerex)/cos(theta_el) ! SI Mylar/Cerex x10
         if(abs(theta_el).gt.0.1192777) tfw_al=tfw_al+real(2.d0*tfw_si_mylar_cerex)/cos(theta_el) ! SI Mylar/Cerex x20
      endif
      
      axi=2.e0*pi*Na*me*rel**2  ! in [GeV*cm^2/mol]

      xi(1)=axi*(ZNH3/ANH3)*X0NH3*tft_nh3
      DeltaEb(1)=xi(1)*(log(alpha**2*X0NH3*tft_nh3/(rel*rhoNH3))+0.198e0)
      xi(2)=axi*(ZHe/AHe)*X0LHe*tfw_lhe
      DeltaEb(2)=xi(2)*(log(alpha**2*X0LHe*tfw_lhe/(rel*rhoLHe))+0.198e0)
      xi(3)=axi*(ZAl/AAl)*X0Al*tfw_al
      DeltaEb(3)=xi(3)*(log(alpha**2*X0Al*tfw_al/(rel*rhoAl))+0.198e0)
      
c      print*,'Most probable ionization energy loss (NH3,LHe,Al): ',DeltaEb(1),DeltaEb(2),DeltaEb(3),
c     &', thickness in X0 (NH3,LHe,Al): ',tft_nh3,tfw_lhe,tfw_al
     
      p_el=p_el-DeltaEb(1)-DeltaEb(2)-DeltaEb(3) ! correct outgoing electron energy, including ionization loss

      return
      end

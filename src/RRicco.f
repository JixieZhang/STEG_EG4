csk   This subroutine can be called from "H2MODEL" INSTEAD
csk   of "R1998" to use G. Ricco et al.'s function "ratio" instead of
csk   L. Whitlow's and S. Rock's parametrization for R
csk   Use at your own risk!!! Sebastian Kuhn Jun-13-2001
      
      subroutine RRicco(X,Q2,R)
csk
      REAL*8 X, Q2, R, ratio
      integer bb
      
      bb=1
      R = ratio(X,Q2,bb)
      return
      end

      function ratio(xx,qq,bb)
csk   function supplied by Mikhail Osipenko 11-Jun-2001
csk   Based on G. Ricco et al, Nucl. Phys. B555, 306(1999)
csk   Modified to prevent nonsensical behavior below Q2=0.25
csk   and to adapt to models.f conventions

      implicit none
      real*8 ratio, xx, qq
      real*8 x,qsq,w,w3,wm,x3,tet,tet3,lgg,R3,br,Rdis,dRd,R,dR,M,M2
      integer bb
      M=0.939                   ! csk
      M2=M**2
      wm=2.5
      w3=wm
      w=sqrt(M2+qq*(1./xx-1.))
csk
      if (qq .ge. 0.25D0) then
         x = xx
         qsq = qq
      else                      ! keep R the same below Q2 = 0.25
         qsq = 0.25
         x=qsq/(qsq - M2 + w*w)
      endif
csk
      x3=qsq/(qsq - M2 + w3*w3)
      tet=1.+12.*qsq/(qsq+1.)*.015625/(.015625+x**2)
      tet3=1.+12.*qsq/(qsq+1.)*.015625/(.015625+x3**2)
      lgg=log(qsq/.04)
      R3=.041*tet3/lgg+.592/qsq-.331/(.09+qsq**2)
      br=qsq*R3/(1.-x3)**3
      Rdis=.041*tet/lgg+.592/qsq-.331/(.09+qsq**2)
      dRd=.006*tet/lgg+.01/qsq+.01/(.09+qsq**2)
      R=br*(1.-x)**3/qsq
      dR=.08
      if(w.gt.w3) then
         R=Rdis
         dR=dRd
      endif
      if(R.lt.0.) R=0.0
c-----SLAC value
c     ratio=0.18
      if(bb.eq.1) then
         ratio=R
      else
         ratio=dR
      endif
      return
      end



      subroutine pot_spec(v1, v2, xmu12, npos, delr, rdeb)
!
!     *******************************************************************  
!
!	Calcul des potentiels et moments dipolaires
!	doit avoir les expressions analytiques 
!	(si numériques : doit avoir un 'read')
!
!     *******************************************************************
!
!
      integer npos
      double precision v1(npos),v2(npos),xmu12(npos)
      double precision delr, rdeb
!
      double precision z0, s, req, x1, x2, y, r
      integer l
!
      z0=.1026277d0
      s=.72d0
      req=2.d0
      x1=1.d0
      x2=-1.11d0
      y=-.055d0
!
      do l = 1, npos
         r = rdeb+(l-1)*delr
         v1(l) = z0*(dexp(-2.d0*s*(r-req)) - 2.d0*x1*dexp(-s*(r-req)))
         v2(l) = z0*(dexp(-2.d0*s*(r-req)) - 2.d0*x2*dexp(-s*(r-req)))
         xmu12(l) = 1.07d0 + (.396d0/(s*y))*(1.d0-dexp(-s*y*(r-req)))
         if((r.gt.12.d0).and.(xmu12(l).gt.0.5d0*r))then
	    xmu12(l)=0.5d0*r
         endif
      end do
!
      return
      end

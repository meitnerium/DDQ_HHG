module spo  


contains  
	  subroutine splitop (cw1, cw2, zetdt, v1, v2, xmu12, npos,champ, delr, xmu, delt)

!     *******************************************************************  
!
!	Calcul du propagateur d'évolution temporelle
!
!     *******************************************************************


      integer npos

      double precision delr,delk,kmax,kmin, champ, xmu, delt,nposreel,x(npos)
      double complex cw1(npos),cw2(npos), zetdt(npos)
      double precision v1(npos),v2(npos),xmu12(npos)

      double complex cun, cim, cnul, cwtemp, cphase1, cphase2
      double complex ctemp1, ctemp2
      double precision xmue, vp1, vp2, delta
      double precision thet,norme,normedeb,k(npos),pi
      integer ll,i


	pi=3.141593654d0
 	delk=2d0*pi/(npos*delr)
 do ll = 1, npos
         if(ll.le.npos/2)then
            k(ll) = (ll-1) * delk
         else
            k(ll) = -(npos-ll+1) * delk
         endif
 enddo
	nposreel=npos*1d0
      cun=dcmplx(1.0d0,0.d0)
      cim=dcmplx(0.d0,1.0d0)
      cnul=dcmplx(0.d0,0.d0)
      do ll = 1, npos
         xmue = xmu12(ll) * champ
	 delta = (v2(ll) - v1(ll))**2 + (2.d0*xmue)**2
	 delta = dsqrt(delta)
            thet = 0.5d0*datan((2.d0*xmue)/(v2(ll)-v1(ll)))

	 vp1 = (v2(ll) + v1(ll) - delta)*0.5d0
	 vp2 = (v2(ll) + v1(ll) + delta)*0.5d0
         cwtemp  = dcos(thet)*cw1(ll)-dsin(thet)*cw2(ll)
         cw2(ll) = dsin(thet)*cw1(ll)+dcos(thet)*cw2(ll)
         cw1(ll) = cwtemp
         cphase1 = cdexp(-cim*vp1*delt/2.d0)
         cphase2 = cdexp(-cim*vp2*delt/2.d0)
         cw1(ll) = cw1(ll)*cphase1
         cw2(ll) = cw2(ll)*cphase2
         cwtemp  = dcos(thet)*cw1(ll)+dsin(thet)*cw2(ll)
         cw2(ll) =-dsin(thet)*cw1(ll)+dcos(thet)*cw2(ll)
         cw1(ll) = cwtemp
      end do

	call zexptdt(zetdt, npos, delk, xmu, delt)

	call fourier(cw1,npos,-1,delr,delk)
	call fourier(cw2,npos,-1,delr,delk)
	


!     *******************************************************************  
!	multiplication vecteur par matrice (z pour complexe)
!     *******************************************************************

      call ZVEM(npos, zetdt(1), 1, cw1(1), 1, cw1(1), 1)
      call ZVEM(npos, zetdt(1), 1, cw2(1), 1, cw2(1), 1)

	call fourier(cw1,npos,1,delr,delk)
	call fourier(cw2,npos,1,delr,delk)
	cw1=cw1/npos
	cw2=cw2/npos

      do ll = 1, npos
         xmue = xmu12(ll) * champ
	 delta = (v2(ll) - v1(ll))**2 + (2.d0*xmue)**2
	 delta = dsqrt(delta)
	 if(dabs(xmue).gt.1.d-16)then
	    thet = 0.5d0*datan((2.d0*xmue)/(v2(ll)-v1(ll)))
         else
            thet = 0.d0
         end if
	 vp1 = (v2(ll) + v1(ll) - delta)*0.5d0
	 vp2 = (v2(ll) + v1(ll) + delta)*0.5d0
         cwtemp  = dcos(thet)*cw1(ll)-dsin(thet)*cw2(ll)
         cw2(ll) = dsin(thet)*cw1(ll)+dcos(thet)*cw2(ll)
         cw1(ll) = cwtemp
         cphase1 = cdexp(-cim*vp1*delt/2.d0)
         cphase2 = cdexp(-cim*vp2*delt/2.d0)
         cw1(ll) = cw1(ll)*cphase1
         cw2(ll) = cw2(ll)*cphase2
         cwtemp  = dcos(thet)*cw1(ll)+dsin(thet)*cw2(ll)
         cw2(ll) =-dsin(thet)*cw1(ll)+dcos(thet)*cw2(ll)
         cw1(ll) = cwtemp
      end do
      return
      end subroutine





!     *******************************************************************  
!
!	Calcul du vecteur d'energie cinetique en representation impulsion
!
!     *******************************************************************

      subroutine zexptdt(etdt, npos, xk1, xmu, delt)
      integer npos
      double complex etdt(npos)
      double precision xmu, delt
      double complex  cim
      double precision pi, xk1, xk(npos), arg
      integer ll

      pi = 3.141592654d0
      cim = dcmplx(0.d0, 1.d0)
      do ll = 1, npos
         if(ll.le.npos/2)then
            xk(ll) = (ll-1) * xk1
         else
            xk(ll) = -(npos-ll+1) * xk1
         endif
	 arg = ((xk(ll)*xk(ll))/(2.d0*xmu)) * delt
	 etdt(ll) = cdexp(-cim*arg)
      end do



      return
      end subroutine

!******************************************************
!	Transformee de Fourier
!******************************************************
	subroutine fourier(F,N,isig,dx,dk)
implicit none
integer :: N, i, j,isig,m,mmax,istep
double precision :: wr,wi,wpr,wpi,wtemp,theta,tempr,tempi
double precision :: G(2*N),dx,dk,normei,normef
double complex :: F(N)
do i=1,N
	G(2*i-1)=dreal(F(i))
	G(2*i)=dimag(F(i))
enddo
!select case (isig)
!	case (1)
!		 call simpson(N,dx,cdabs(F)**2, normei)
!	case (-1)
!		 call simpson(N,dk,cdabs(F)**2, normei)
!end select
j=1
do i=1,2*N,2
	if (j.gt.i) then
		tempr=(G(j))
		tempi=(G(j+1))
		G(j)=G(i)
		G(j+1)=G(i+1)
		G(i)=tempr
		G(i+1)=tempi
	endif
	m=N
        do while ((M.ge.2).and.(j.gt.m)) 
		j=j-m
		m=m/2
	enddo
	j=j+m
enddo
mmax=2
        do while (2*N.gt.mmax)
		istep=2*mmax
		theta=6.28318530717959d0/(isig*mmax)
		wpr=-2d0*dsin(0.5d0*theta)**2
		wpi=dsin(theta)
		wr=1d0
		wi=0d0
		do m=1,mmax,2
			do i=m,2*N,istep
				j=i+mmax
				tempr=sngl(wr)*G(j)-sngl(wi)*G(j+1)
				tempi=sngl(wr)*G(j+1)+sngl(wi)*G(j)
				G(j)=G(i)-tempr
				G(j+1)=G(i+1)-tempi
				G(i)=G(i)+tempr
				G(i+1)=G(i+1)+tempi
			enddo
			wtemp=wr
			wr=wr*wpr-wi*wpi+wr
			wi=wi*wpr+wtemp*wpi+wi
		enddo
		mmax=istep
	enddo
do i=1,N
	tempr=G(2*i-1)
	tempi=G(2*i)
	F(i)=dcmplx(tempr,tempi)
	
enddo

return
end subroutine




!***************************************************************
! 		Integration Simpson
!***************************************************************

SUBROUTINE simpson (N,H,FI,S)
!
! Subroutine for integration over f(x) with the Simpson rule.  FI:
! integrand f(x); H: interval; S: integral.  Copyright (c) Tao Pang 1997.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I
  double precision, INTENT (IN) :: H
  double precision :: S0,S1,S2
  double precision, INTENT (OUT) :: S
  double precision, INTENT (IN), DIMENSION (N) :: FI
!
  S  = 0.0
  S0 = 0.0
  S1 = 0.0
  S2 = 0.0
  DO I = 2, N-1, 2
    S1 = S1+FI(I-1)
    S0 = S0+FI(I)
    S2 = S2+FI(I+1)
  END DO
  S = H*(S1+4.0*S0+S2)/3.0
!
! If N is even, add the last slice separately
!
  IF (MOD(N,2).EQ.0) S = S &
     +H*(5.0*FI(N)+8.0*FI(N-1)-FI(N-2))/12.0
END SUBROUTINE simpson


!***************************************************************
! 		Calcul norme d'une fonction complexe
!***************************************************************
	subroutine simps(func, vint, delti, npl)

      integer j, npl
      double complex func(npl)
      double precision  vint, delti
      
      vint=0
      do j = 1, npl
         vint=vint+delti*(cdabs(func(j))**2) 
      end do
	vint=sqrt(vint)
      return
      end subroutine

	end module spo



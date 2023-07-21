	module spo     
	  subroutine splitop (cw1, cw2, zetdt, v1, v2, xmu12, npos, champ, delr, xmu, delt)
	implicit none
!     *******************************************************************  
!
!	Calcul du propagateur devolution temporelle
!
!     *******************************************************************


      integer npos

      double precision delr, champ, xmu, delt,nposreel
      double complex cw1(npos),cw2(npos), zetdt(npos)
      double precision v1(npos),v2(npos),xmu12(npos)

      double complex cun, cim, cnul, cwtemp, cphase1, cphase2
      double complex ctemp1, ctemp2
      double precision xmue, vp1, vp2, delta
      double precision xk1, xkn, arg, thet
      integer ll
      write(*,*)'qgqr'
	nposreel=npos*1d0
      cun=dcmplx(1.0d0,0.d0)
      cim=dcmplx(0.d0,1.0d0)
      cnul=dcmplx(0.d0,0.d0)
      do ll = 1, npos
         xmue = xmu12(ll) * champ
	 delta = (v2(ll) - v1(ll))**2 + (2.d0*xmue)**2
	 delta = dsqrt(delta)
            thet = .5d0*datan((2.d0*xmue)/(v2(ll)-v1(ll)))

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

	call zexptdt(zetdt, npos, delr, xmu, delt)

		Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 1, npos )	
		Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_INPLACE)
		Status = DftiSetValue( My_Desc1_Handle, DFTI_FORWARD_SCALE, 1d0/sqrt(nposreel))
		Status = DftiCommitDescriptor( My_Desc1_Handle )
		Status = DftiComputeForward( My_Desc1_Handle, cw1)
		Status = DftiComputeForward( My_Desc1_Handle, cw2)
		Status = DftiFreeDescriptor(My_Desc1_Handle)

!     *******************************************************************  
!	multiplication vecteur par matrice (z pour complexe)
!     *******************************************************************

      call ZVEM(npos, zetdt(1), 1, cw1(1), 1, cw1(1), 1)
      call ZVEM(npos, zetdt(1), 1, cw2(1), 1, cw2(1), 1)

		Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 1, npos )	
		Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_INPLACE)
		Status = DftiSetValue( My_Desc1_Handle, DFTI_FORWARD_SCALE, 1d0/sqrt(nposreel))
		Status = DftiCommitDescriptor( My_Desc1_Handle )
		Status = DftiComputeBackward( My_Desc1_Handle, cw1)
		Status = DftiComputeBackward( My_Desc1_Handle, cw2)
		Status = DftiFreeDescriptor(My_Desc1_Handle)

      do ll = 1, npos
         xmue = xmu12(ll) * champ
	 delta = (v2(ll) - v1(ll))**2 + (2.d0*xmue)**2
	 delta = dsqrt(delta)
	 if(dabs(xmue).gt.1.d-16)then
	    thet = .5d0*datan((2.d0*xmue)/(v2(ll)-v1(ll)))
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
      end





!     *******************************************************************  
!
!	Calcul du vecteur d'energie cinetique en representation impulsion
!
!     *******************************************************************

      subroutine zexptdt(etdt, npos, delr, xmu, delt)
      integer npos
      double complex etdt(npos)
      double precision xmu, delr, delt
      double complex carg, cim
      double precision pi, xk1, xkn, arg
      integer ll

      pi = 3.141592654d0
      cim = dcmplx(0.d0, 1.d0)
      xk1 = 2.d0*pi/(delr*npos)
      do ll = 1, npos
         if(ll.le.npos/2)then
            xkn = (ll-1) * xk1
         else
            xkn = -(npos-ll+1) * xk1
         endif
	 arg = ((xkn*xkn)/(2.d0*xmu)) * delt
	 etdt(ll) = cdexp(-cim*arg)
      end do
      return
      end
	end module spo








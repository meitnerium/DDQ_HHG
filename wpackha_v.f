
	      program wpack
c     Program developped by H.Abou-Rachid and Catherine Lefebvre June 2003
c     propagation de paquets d'ondes par split-operator
c     **************************************************************************
c
c	Declarations des vecteurs, variables et parametres
c
c     **************************************************************************
c
      character*24 fichier1, fichier2, fichier3, fichier4, fichier5
      logical writepot /.false./
      logical Etatpro
      integer idim, ndim
      parameter (idim=4096,ndim=4*idim)
      double complex chi1(idim), chi2(idim),psik1(ndim),psik2(ndim)
      double complex chi11(idim),chi22(idim)
      double complex zetdt(idim)
      double precision v1(idim),v2(idim),xmu12(idim)
      double precision morse, r0cut, scut, vibfunc(0:18,idim)
      double precision delr, rdeb, rmax, delt, tf, seuil,tw,theta
      double precision t0, omega, intens, freq,gama,phase,deltuv
      integer npos, ntemps, npbig
c
      double complex zwork1(ndim), zwork2(ndim), 
     &zwork3(ndim), zwork4(ndim)
      double complex cun, cim, cnul, zcutA(ndim), zcutI(ndim)
      double precision r, xmu, pi, rc0, p0, alpha, t
      double precision int1tf, int2tf, int3tf
      double precision int1t0, int2t0, int3t0
      double precision projreal, projimag, dissprob, champ,env
	double precision dissprobv
      double precision w1(idim), w2(idim), proj(idim)
      double precision xnorm1, xnorm2, work1(ndim), work2(ndim)
      double precision worka(2*idim), workb(2*ndim)
      double precision tablea(2*idim+30),tableb(2*ndim+30)
      double precision xnormk1, xnormk2, work3(ndim), work4(ndim)
      double precision dk, cte, xk, spk
      double precision dtper, timeper, evbyau, normedeb, period
      double precision delta, xmue, thet, vp1, vp2
      double complex cwtemp1, cwtemp2
      integer k, m
      integer i, j, l, nu, nudep
      external time
      integer time, cput1, cput2



      
c
c 	*********************************************************************************
c
c	Variables susceptibles d'etre changees
c
c	**********************************************************************************

      
      namelist /iofile/
     & fichier1, fichier2, fichier3, fichier4, fichier5,
     & delt, tf, t0, omega, intens, freq, phase,
     & npos, npbig, rdeb, rmax, seuil, nudep,gama,tw,theta,
     & r0cut, scut, xmu, Etatpro,deltuv
c
      cput1 = time()
c
c     *****************************************************************************************
c       
c         Donnees par defaut (se référer au input puisque la commande 'read(5,iofile)' 
c	  se trouve apres cette liste
c
c     *****************************************************************************************
      fichier1='fch1.dat'
      fichier2='fch2.dat'
      fichier3='fch3.dat'
      fichier4='fch4.dat'
      fichier5='fch5.dat'
 
c      write(*,*) 'start'
      t0 = 0.d0
      omega = 1.d0
      intens = 5.33d-2
      freq = .1381962d0
      phase=1.d0
      gama=0.5d0
      npos = 1024
      npbig = 4096
      rdeb = 2.5d0
      rmax = 60.d0
      delt = 1.778528d-01
      tf = 400.d0
      tw=500.d0
      theta=90000000000.0d0
      deltuv=90000000000.d0
c     ******************************************************************
c	Changer ces valeurs pour passer du systeme H2+ à D2+
c     ******************************************************************
      nu = 0
      xmu = 918.074d0
c     ******************************************************************
      read(5,iofile)
c
      cun = dcmplx(1.0d0,0.d0)
      cim = dcmplx(0.d0,1.0d0)
      cnul = dcmplx(0.d0,0.d0)
      pi = 3.141592654d0
      evbyau = 27.212d0
      rc0 = 1.3989d0
      p0 = 0.0d0
      alpha = 9.2085226d0
      delr = (rmax-rdeb)/(npos-1)
      dtper = 2.d0*pi/omega/1024.d0
      period = 2.d0*pi/freq
c
      open(99,file='champ.dat',status='unknown')
c
      ntemps = idnint(tf/delt+0.5d0)
c
c     *******************************************************************  
c
c	Initialisation des vecteurs FFT
c
c     *******************************************************************
c


      call zzfft(0,npos,1.0d0,chi11(1),chi11(1),tablea(1),worka(1),0)
      call zzfft(0,npbig, 1.0d0,chi22(1),chi22(1),tableb(1),workb(1),0)


c
c
c     *******************************************************************  
c
c	Calcul de la fonction d'onde initiale
c
c     *******************************************************************
c


      call eval(chi1(1), chi2(1), psik1(1), psik2(1), delr, rdeb,
     &          p0, rc0, alpha, npos, npbig)



c
c     *******************************************************************  
c
c	Calcul des potentiels et moments dipolaires
c	(doit avoir les expressions analytiques des potentiels)
c
c     *******************************************************************
c



      call pot_spec(v1(1), v2(1), xmu12(1), npos, delr, rdeb)



c
c     *******************************************************************  
c
c	Calcul de la fonction de coupure - decoupage de la fct d'onde en 2 parties
c	(dans le fichier subroutine 'airesint')
c
c     *******************************************************************
c



      call calczcut(zcutA(1), zcutI(1), rdeb, delr, npbig,
     &              r0cut, scut)



c
c     *******************************************************************  
c
c	Calcul du vecteur d'energie cinetique radial
c	(dans le fichier subroutine 'splitop')
c
c     *******************************************************************
c


      call zexptdt(zetdt(1), npos, delr, xmu, delt)


c
c     *******************************************************************  
c
c	Initialisation - Calcul des fonctions vibrationnelles du potentiel de Morse
c	(function 'morse' en relation avec le polynome de Laguerre)
c
c     *******************************************************************
c

      do nu = 0, 18
         r = rdeb - delr
         do l = 1, npos
            r = r + delr
            vibfunc(nu,l) = morse(1.026202d-1,0.72d0,xmu,2.d0,r,nu)
         end do
      end do
c
      if (Etatpro) then
         do l = 1, npos
            chi1(l) = vibfunc(nudep,l)
         end do
      end if
c

c
c     *******************************************************************  
c
c	Normalisation 
c
c     *******************************************************************
c


 
      do l = 1, npos
         work1(l) = cdabs(chi1(l))**2
      end do

c
c     *******************************************************************  
c
c	Methode Trapezoide (voir fin du programme)
c
c     *******************************************************************
c
 
      call simps(work1(1), normedeb, delr, npos)
c
c     *******************************************************************  
c
c	Renormalisation de la fonction d'onde initiale
c
c     *******************************************************************
c
      do l = 1, npos
         chi1(l) = chi1(l)/dsqrt(normedeb)
         work1(l) = cdabs(chi1(l))**2
      end do
      call simps(work1(1), normedeb, delr, npos)
c

c
c     *******************************************************************  
c
c	Initialisation des integrales pour l'analyse asymptotique
c
c     *******************************************************************
c
 
      t = 0.d0
      int1tf = 0.d0
      int2tf = 0.d0
      int3tf = 0.d0
c
      timeper = 0.d0
      do while(timeper.lt.tf)
         timeper = timeper + dtper
         
         call airesint (int1tf, int2tf, int3tf, t0, omega,
     &  phase, intens, freq,gama, timeper, dtper,tw,theta,deltuv)
      end do
       
      call airesint(int1tf, int2tf, int3tf, t0, omega,
     &  phase, intens, freq,gama, tf, tf-timeper,tw,theta,deltuv)
c
      int1t0 = 0.d0
      int2t0 = 0.d0
      int3t0 = 0.d0
      t = 0.d0
      timeper = 0.d0
      m = 0
      k = 0


      
c
c     *******************************************************************  
c	!!!!!!!!!!
c	Ouverture de la boucle sur le temps - Debut de la propagation
c	!!!!!!!!!!
c     *******************************************************************
c


      do i = 1, ntemps
         t = t + delt
c
         do while(timeper.lt.t)
            timeper = timeper + dtper

c
c     *******************************************************************  
c
c	Calcul des aires temporelles du champ
c
c     *******************************************************************
c
       
          
         call airesint(int1t0, int2t0, int3t0, t0, omega,
     &     phase, intens, freq,gama, timeper, dtper,tw,theta,deltuv)
         end do
         
         call airesint(int1t0, int2t0, int3t0, t0, omega,
     &     phase, intens, freq,gama, t, t-timeper,tw,theta,deltuv)
         timeper = t
c
         
         champ = env (t0, omega,
     &    intens, freq,gama, phase,t - 0.5d0*delt,tw,theta,deltuv)
c
         write(99,*)t,champ
c
c     *******************************************************************  
c
c	Calcul du propagateur d'evolution temporelle
c
c     *******************************************************************
c
c
         call splitop(chi1(1), chi2(1), zetdt(1), v1(1), v2(1),
     &                xmu12(1), npos, champ, delr, xmu, delt,
     &                 tablea(1), worka(1))
c

c
c     *******************************************************************  
c
c	Analyse asymptotique en cours de propagation
c
c     *******************************************************************
c


         if( (cdabs(chi1(npos)).gt.seuil).or.
     &                    (cdabs(chi2(npos)).gt.seuil) )then
            call asympt(t, tf, psik1(1), psik2(1), chi1(1), chi2(1),
     &                  zcutA(1), npbig, npos, xmu,
     &                  int1tf, int2tf, int3tf,
     &                  int1t0, int2t0, int3t0, delr,
     &                  zwork1(1), zwork2(1), zwork3(1),
     &                  worka(1), workb(1), tablea(1), tableb(1))


c
c     *******************************************************************  
c
c	Recuperation de ce qui reste du paquet d'onde
c
c     *******************************************************************
c

            call ZVEM(npos,chi1(1),1,zcutI(1),1,chi1(1),1)
            call ZVEM(npos,chi2(1),1,zcutI(1),1,chi2(1),1)
         endif

c
c     *******************************************************************  
c
c	Analyse des resultats a chaque tranche de temps
c
c     *******************************************************************
c
c
         if(mod(i,10).eq.0)then
            dissprob = 0.d0
            do nu = 0, 18
               do j = 1, npos
                  proj(j) = vibfunc(nu,j)*dreal(chi1(j))
               end do
               call simps(proj(1),projreal,delr,npos)
               do j = 1, npos
                  proj(j) = vibfunc(nu,j)*dimag(chi1(j))
               end do
               call simps(proj(1),projimag,delr,npos)

c
c     *******************************************************************  
c
c	Calcul de la probabilite (P_lie)
c	Somme sur tous les etats vibrationnels pour tous les r
c
c     *******************************************************************
c
		dissprobv=projreal**2 + projimag**2
		write(30+nu,1003) t, dissprobv, 1-dissprobv
               dissprob = dissprob + projreal**2 + projimag**2
            end do
            write(29,1002) t, dissprob
         endif

c
c     *******************************************************************  
c
c	Cette partie est en commentaire
c
c     *******************************************************************
c
c
c         if (((t.gt.0.d0).and.(t.le.0.d0+delt)).or.
c    &((t.gt.1.d0*period/4.d0).and.(t.le.1.d0*period/4.d0+delt)).or.
c     &((t.gt.2.d0*period/4.d0).and.(t.le.2.d0*period/4.d0+delt)).or.
c     &((t.gt.3.d0*period/4.d0).and.(t.le.3.d0*period/4.d0+delt)).or.
c     &((t.gt.4.d0*period/4.d0).and.(t.le.4.d0*period/4.d0+delt)).or.
c
c     &((t.gt.76.d0*period/4.d0).and.(t.le.76.d0*period/4.d0+delt)).or.
c     &((t.gt.77.d0*period/4.d0).and.(t.le.77.d0*period/4.d0+delt)).or.
c     &((t.gt.78.d0*period/4.d0).and.(t.le.78.d0*period/4.d0+delt)).or.
c     &((t.gt.79.d0*period/4.d0).and.(t.le.79.d0*period/4.d0+delt)).or.
c     &((t.gt.80.d0*period/4.d0).and.(t.le.80.d0*period/4.d0+delt)))then
c            m = m + 1
c            write(30+m,1002) t, t/period
c            r = rdeb - delr
c
c     *******************************************************************  
c
c	Evolution du p.o. a differente periode du cycle optique
c
c     *******************************************************************
c
c            do j = 1, npos
c               r = r + delr
c               xmue = xmu12(j) * champ
c	       delta = (v2(j) - v1(j))**2 + (2.d0*xmue)**2
c	       delta = dsqrt(delta)
c	       vp1 = (v2(j) + v1(j) - delta)*0.5d0
c	       vp2 = (v2(j) + v1(j) + delta)*0.5d0
c               write(30+m,1007) r, dsqrt(dreal(chi1(j))**2+
c     &  dimag(chi1(j))**2),
c     & dsqrt(dreal(chi2(j))**2+ dimag(chi2(j))**2), vp1, vp2
c            end do
c         end if
c
c
c     *******************************************************************  
c
c	Energie cinetique a differente periode du sycle optique
c
c     *******************************************************************
c
c         if (((t.gt.period).and.(t.le.period+delt)).or.
c     &((t.gt.2.d0*period).and.(t.le.2.d0*period+delt)).or.
c     &((t.gt.3.d0*period).and.(t.le.3.d0*period+delt)).or.
c     &((t.gt.4.d0*period).and.(t.le.4.d0*period+delt)).or.
c     &((t.gt.5.d0*period).and.(t.le.5.d0*period+delt)).or.
c     &((t.gt.6.d0*period).and.(t.le.6.d0*period+delt)).or.
c     &((t.gt.7.d0*period).and.(t.le.7.d0*period+delt)).or.
c     &((t.gt.8.d0*period).and.(t.le.8.d0*period+delt)).or.
c     &((t.gt.9.d0*period).and.(t.le.9.d0*period+delt)).or.
c     &((t.gt.10.d0*period).and.(t.le.10.d0*period+delt)).or.
c     &((t.gt.11.d0*period).and.(t.le.11.d0*period+delt)).or.
c     &((t.gt.12.d0*period).and.(t.le.12.d0*period+delt)).or.
c     &((t.gt.13.d0*period).and.(t.le.13.d0*period+delt)).or.
c    &((t.gt.14.d0*period).and.(t.le.14.d0*period+delt)).or.
c    &((t.gt.15.d0*period).and.(t.le.15.d0*period+delt)).or.
c     &((t.gt.16.d0*period).and.(t.le.16.d0*period+delt)).or.
c     &((t.gt.17.d0*period).and.(t.le.17.d0*period+delt)).or.
c     &((t.gt.18.d0*period).and.(t.le.18.d0*period+delt)).or.
c     &((t.gt.19.d0*period).and.(t.le.19.d0*period+delt)).or.
c     &((t.gt.20.d0*period).and.(t.le.20.d0*period+delt)))then
c            k = k + 1
c            write(50+k,1002) t, t/period
c            dk = 2.d0*pi/(npbig*delr)
c            cte = 0.5d0/xmu
c
c     *******************************************************************  
c
c	Energie cinetique a temps courant
c
c     *******************************************************************
c
c            do j = npbig/2 +2, npbig
c               xk = dk*(j-1 - npbig)
c               work1(j) = cdabs(psik1(j))**2
c               work2(j) = cdabs(psik2(j))**2
c               work1(j) = work1(j)*(-xmu/(xk*evbyau))
c               work2(j) = work2(j)*(-xmu/(xk*evbyau))
c               spk = work1(j) + work2(j)
c               xk = -xk*xk*0.5/xmu
c               write(50+k,1004) xk*evbyau, work1(j), work2(j), spk
c            end do
c
c            do j = 1, npbig/2 +1
c               xk = dk*(j-1)
c               work1(j) = cdabs(psik1(j))**2
c               work2(j) = cdabs(psik2(j))**2
c               if (xk.eq.0.d0) then
c                  work1(j) = 0.d0
c                  work2(j) = 0.d0
c               else
c                  work1(j) = work1(j)*(xmu/(xk*evbyau))
c                  work2(j) = work2(j)*(xmu/(xk*evbyau))
c               end if
c               spk = work1(j) + work2(j)
c               xk = xk*xk*0.5/xmu
c               write(50+k,1004) xk*evbyau, work1(j), work2(j), spk
c            end do
c         end if
c
c         write(9,*) 100.d0*float(i)/ntemps,'% done'
c         backspace(9)
      end do
c
c     *******************************************************************  
c	!!!!!!!!
c	Fin de la boucle sur le temps - Fin de la propagation !!!!
c	!!!!!!!!
c     *******************************************************************
c
c
c     *******************************************************************  
c
c	Analyse asymptotique a la fin de la propagation
c
c     *******************************************************************
c
      call asympt(tf, tf, psik1(1), psik2(1), chi1(1), chi2(1),
     &            zcutA(1), npbig, npos, xmu,
     &            int1tf, int2tf, int3tf,
     &            int1t0, int2t0, int3t0, delr,
     &            zwork1(1), zwork2(1), zwork3(1),
     &            worka(1), workb(1), tablea(1), tableb(1))
      call ZVEM(npos,chi1(1),1,zcutI(1),1,chi1(1),1)
      call ZVEM(npos,chi2(1),1,zcutI(1),1,chi2(1),1)
c
c      open(1,file=fichier1,status='unknown')
      open(2,file=fichier2,status='unknown')
      r = rdeb - delr
c
c     *******************************************************************  
c
c	P.o. au temps final
c
c     *******************************************************************
c
      do l = 1, npos
	 r = r + delr
	 w1(l) = cdabs(chi1(l))**2
	 w2(l) = cdabs(chi2(l))**2
c	 write(1,1004) dreal(chi1(l)),dimag(chi1(l)),
c     $              dreal(chi2(l)),dimag(chi2(l))
	 write(2,1003) r, w1(l), w2(l)
      end do
c      close(1,status='keep')
      close(2,status='keep')
c
      open(1,file=fichier3,status='unknown')
      open(2,file=fichier4,status='unknown')
      open(3,file=fichier5,status='unknown')
      dk = 2.d0*pi/(npbig*delr)
      cte = 0.5d0/xmu
c
c
c     *******************************************************************  
c
c	Energie cinetique au temps final (partie negative)
c
c     *******************************************************************
c
      do i = npbig/2 +2, npbig
         xk = dk*(i-1 - npbig)
         work1(i) = cdabs(psik1(i))**2
         work2(i) = cdabs(psik2(i))**2
         work3(i) = work1(i)
         work4(i) = work2(i)
         work1(i) = work1(i)*(-xmu/(xk*evbyau))
         work2(i) = work2(i)*(-xmu/(xk*evbyau))
         spk = work1(i) + work2(i)
         xk = -xk*xk*0.5/xmu
         write(1,1002) xk*evbyau, work1(i)
         write(2,1002) xk*evbyau, work2(i)
         write(3,1002) xk*evbyau, spk
      end do
c
c
c     *******************************************************************  
c
c	Energie cinetique au temps final (partie positive)
c
c     *******************************************************************
c
      do i = 1,npbig/2 +1
         xk = dk*(i-1)
         work1(i) = cdabs(psik1(i))**2
         work2(i) = cdabs(psik2(i))**2
         work3(i) = work1(i)
         work4(i) = work2(i)
         if (xk.eq.0.d0) then
            work1(i) = 0.d0
            work2(i) = 0.d0
         else
            work1(i) = work1(i)*(xmu/(xk*evbyau))
            work2(i) = work2(i)*(xmu/(xk*evbyau))
         end if
         spk = work1(i) + work2(i)
         xk = xk*xk*0.5/xmu
         write(1,1002) xk*evbyau, work1(i)
         write(2,1002) xk*evbyau, work2(i)
         write(3,1002) xk*evbyau, spk
      end do
c
      close(1)
      close(2)
      close(3)
      close(99)
c
      call simps(w1(1), xnorm1, delr, npos)
      call simps(w2(1), xnorm2, delr, npos)
      write(*,*)
      write(*,*) 'Niveau vibrationnel de depart'
      write(*,*) nu
      write(*,*)
      write(*,*) 'NORME de depart'
      write(*,*) normedeb
      write(*,*)
      write(*,*) 'NORMES population restante'
      write(*,*) '    SIGMA g           SIGMA u           TOTAL'
      write(*,1003) xnorm1, xnorm2, xnorm1 + xnorm2
      call simps(work3(1), xnormk1, dk, npbig)
      call simps(work4(1), xnormk2, dk, npbig)
      write(*,*)
      write(*,*) 'NORMES population dissociee'
      write(*,*) '    SIGMA g           SIGMA u           TOTAL'
      write(*,1003) xnormk1, xnormk2, xnormk1 + xnormk2
      write(*,*)
      write(*,*)
      write(*,*) '1 - NORMES population restante:', 1 - xnorm1 - xnorm2
      write(*,*)
      write(*,*) 'Erreur sur normes:',
     &           1 - xnormk1 - xnormk2 - xnorm1 - xnorm2
c

1002  format(2(e16.8e3,2x))
1003  format(3(e16.8e3,2x))
1004  format(4(e16.8e3,2x))
1005  format(5(e16.8e3,2x))
1006  format(6(e16.8e3,2x))
1007  format(7(e16.8e3,2x))
      write(*,*) 'end'
c
c
c     *******************************************************************  
c
c	Écriture du temps de calcul
c
c     *******************************************************************
c

      cput2 = time() - cput1
c
      write(*,*) idint(cput2/3600.d0),' hr',
     &           idint((cput2 - mod(cput2,60) -
     &                  3600*idint(cput2/3600.d0))/60.d0),' min',
     &           mod(cput2,60),' sec'
c
      stop
      end
c
c    ******************************************************************************************
c    ******************************************************************************************
      subroutine simps(func, vint, delti, npl)
c
      integer j, npl
      double precision func(npl), vint, delti
      double precision s
c
      s = -.5d0 * (func(1) + func(npl))
      do j = 1, npl
         s = s + func(j)
      end do
      vint = delti * s
      return
      end

      function dargz(c)
      implicit double complex (a-c)
      implicit double precision (d-h,o-z)
      dargz=datan2(dimag(c),dreal(c))
      return
      end
      

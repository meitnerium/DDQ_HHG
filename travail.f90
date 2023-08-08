program travail

use spo

implicit none

integer :: npos,ntps, n,v,i,j,l,ph,m, maxpot1,minpot2,ideb,ii
!npos pour la grille de position et d'impulsion, ntps pour la grille de temps

parameter (npos=512, ntps=24000,ideb=85)
double precision, allocatable :: pot(:,:), x(:), ep(:,:)
double precision, allocatable :: xmu12(:),HHG(:),HHG_simple_v(:)
double precision :: xmin, xmax, requ, diss,massreduite, morse,lieprobv, temp_omega, domega
double precision :: work1(npos),table1(npos),normedeb, champ(ntps),e0,rc0
double precision :: tablea(npos),worka(npos),projreal, projimag, projreal_HHG, projimag_HHG,lieprob 
double complex :: chi1(npos), chi2(npos),cun,cim,cnul,zetdt(npos),ctemp(npos),chilie(npos),chi1init(npos)
double precision :: alpha,delr,p0,rdeb,proj(npos) ,proj_R(npos), proj_I(npos), xmue,auto_correl(ntps), tmpreal
double precision :: t(ntps),delt,t0,tf,pi,omega,freq,phase, dtper, norme,periode,delta,vp1(npos),vp2(npos),sigma,tmax,rmoyen(ntps),rmoyenlie(ntps),rclapet1(ntps),rclapet2(ntps)
real,dimension(npos-ideb) :: vp1reel,vp2reel
! 
!
double complex, allocatable ::xmu_chi1(:,:),xmu_chi2(:,:),proj_xmu_chi(:)
double precision, allocatable :: proj_HHG(:,:)    ! pour calcul du spectre HHG
!
!***********************************************************************
!         Valeurs des paramètres, allocation des variables
!***********************************************************************

 cun = dcmplx(1.0d0,0.d0)
 cim = dcmplx(0.d0,1.0d0)
 cnul = dcmplx(0.d0,0.d0)
 pi = 3.141592654d0


requ=2.d0 !valeur de r à l'equilibre = 2*a0 (a0=1 en u.a) 
diss=2.7925d0/27.2d0 !potentiel de dissociation de H2+
massreduite=918.0762887608628d0 !=masse proton/2 en u.a
v=18 !donne 19 niveaux vibrationnels en partant de v=0
xmin=2.d-3 
xmax=30.d0
rdeb = xmin !paramètre pour la construction de chi1 et chi2
rc0 = 1.3989d0 !position du paquet d'onde à t0
p0 = 0.0d0 !impulsion du paquet d'onde à t0
alpha = 13.019d0 !paramètre pour chi1 et chi2
t0=0d0
e0 = 3.77d-2!racine carrée de l'intensité du champ

freq = 4.3d-3 !fréquence du champ
!freq=freq/1.25d0

periode=2*pi/freq
!tf=13*periode
delt=periode/(1023) ! grille de 1024 temps utilisee pour calcul du spectre HHG. Ceci definit ici delt
tf=t0+delt*ntps
!delt=(tf-t0)/(ntps) 
write(*,*) 'tf= ',tf,'delt =',delt

allocate (x(npos))
allocate (xmu12(npos)) !moment dipolaire de H2+
allocate (pot(2,npos)) ! 2 potentiels
allocate (ep(0:v,npos)) !les etats propres vont de 0 à v
!
allocate (proj_xmu_chi(1024))
allocate (proj_HHG(1024,npos))
allocate (xmu_chi1(1024,npos))! xmu_chi1(t_n,R_K) sur  grille de 1024 temps utilisee pour calcul du spectre HHG
allocate (xmu_chi2(1024,npos))! xmu_chi1(t_n,R_K) sur  grille de 1024 temps utilisee pour calcul du spectre HHG
allocate (HHG(1024))
allocate (HHG_simple_v(1024))
!

!***********************************************************************
! Mise en place de la grille des positions, potentiels de Morse et
! etats propres associés
!***********************************************************************
delr=(xmax-xmin)/(npos-1)
	do i=1,npos
	 x(i)=xmin+(i-1)*delr !création de l'axe des positions
	enddo
 
 call pot_spec(pot(1,:),pot(2,:),xmu12, npos,delr,xmin) !construction des 2 potentiels de H2+ et du moment dipolaire
 
 open(unit=1,name="potentiel.dat",status='replace') 
	do i=1,npos
	  write(1,*)x(i),pot(1,i),pot(2,i) 
	enddo
 close(1)

 do n=0,v

	do i=1,npos
 	  ep(n,i)=morse(diss,0.72d0,massreduite,requ,x(i),n) !construction des n etats vibrationnels sur la grille
	work1(i)=(dabs(ep(n,i)))**2
        enddo

	  call simpson(npos,delr,work1,norme)
	  ep(n,:)=ep(n,:)/sqrt(norme) !normalisation des etats vibrationnels

	open(unit=(10+n),status='replace')
	do i=1,npos
		write(10+n,*)x(i),ep(n,i)!ecriture des etats propres d'indice n dans fort.10+n
	enddo	
	close(10+n)
 enddo



! *******************************************************************  	
!    		Construction de la grille temporelle	
! *******************************************************************
     
      do i = 1, ntps
         t(i) = t0 + (i-1)*delt
      enddo

!***********************************************************************
!            Calcul de la fonction d'onde initiale
!***********************************************************************

!do ph=0,1
!ph=0
	m=0 
	phase=ph*pi/2d0 !calcul pour 2 phases, 0 et pi/2
 call eval(chi1, chi2, delr, rdeb, p0, rc0, alpha, npos)
 !do l=1,npos
!	chi1(l)=dcmplx(ep(0,l),0d0)
!	chi2(l)=dcmplx(0d0,0d0)
! enddo
 do j=1,npos
	chi1init(j)=chi1(j)
 enddo


!***********************************************************************
! 	      	Normalisation de la fonction d'onde initiale
!***********************************************************************
	do l=1,npos
	worka(l)=(cdabs(chi1(l)))**2
	enddo

   call simpson( npos,delr,worka, normedeb)
 open(unit=4,name='chi1init.dat',status='replace')
 open(unit=5,name='chi2init.dat',status='replace')
   do l = 1, npos
         chi1(l) = chi1(l)/dsqrt(normedeb)
	write(4,*)x(l),dreal(chi1(l)),dimag(chi1(l))
	write(5,*)x(l),dreal(chi2(l)),dimag(chi2(l))
   enddo
 close(4)
 close(5)


!********************************************************************
! 	   		Construction du champ
!********************************************************************
 sigma=periode/(2d0*2.3548d0)
 tmax=2d0*periode
  open(unit=2+ph,status='replace')	
 	do j=1,ntps	
		!champ(j)=e0*dexp(-(tmax-t(j))**2/(2*sigma**2))*dcos(freq*t(j)+phase)
		champ(j)=e0*dcos(freq*t(j)+phase)
	  write(2+ph,*)t(j), champ(j)
	enddo
   close (2+ph)

!********************************************************************
!	Application du split operator-Ouverture de la boucle temps 
!********************************************************************

 do i=1,ntps 
 write(*,*)'t = ',t(i),' u.a',' phase = ',ph
	do j=1,npos
		worka(j)=x(j)*(cdabs(chi1(j)))**2 !position moyenne du paquet d'onde de l'état fondamental
		tablea(j)=(cdabs(chi1(j)))**2 !densité de probabilité du paquet d'onde de l'état fondamental
	enddo
	call simpson(npos,delr,worka,normedeb)
	call simpson(npos,delr,tablea,norme)
	rmoyen(i)=normedeb/norme
	write(7+ph,*)t(i),rmoyen(i)!r moyen du paquet d'onde total

	do j=1,npos
		chilie(j)=dcmplx(0d0,0d0)
		do n=0,v
			ctemp(j)=ep(n,j)*chi1(j)
			call simpson(npos,delr,dreal(ctemp),normedeb)
			call simpson(npos,delr,dimag(ctemp),norme)
			chilie(j)=chilie(j)+ep(n,j)*(normedeb+cim*norme)
		enddo
		worka(j)=x(j)*(cdabs(chilie(j)))**2
		tablea(j)=(cdabs(chilie(j)))**2
	enddo
	
	call simpson(npos,delr,worka,normedeb)
	call simpson(npos,delr,tablea,norme)
	rmoyenlie(i)=normedeb/norme
	write(30+ph,*)t(i),rmoyenlie(i)!rmoyen de la partie liee du paquet d'onde
		

	         !if (((t(i).gt.(0.d0)).and.(t(i).le.(0.d0+delt))).or.((t(i).gt.(1.d0*periode/4.d0)).and.(t(i).le.(1.d0*periode/4.d0+delt))).or.((t(i).gt.(2.d0*periode/4.d0)).and.(t(i).le.(2.d0*periode/4.d0+delt))).or.((t(i).gt.(3.d0*periode/4.d0)).and.(t(i).le.(3.d0*periode/4.d0+delt))).or.((t(i).gt.(4.d0*periode/4.d0)).and.(t(i).le.(4.d0*periode/4.d0+delt))).or.((t(i).gt.(12.d0*periode/4.d0)).and.(t(i).le.(12.d0*periode/4.d0+delt))).or.((t(i).gt.(13.d0*periode/4.d0)).and.(t(i).le.(13.d0*periode/4.d0+delt))).or.((t(i).gt.(14.d0*periode/4.d0)).and.(t(i).le.(14.d0*periode/4.d0+delt))).or.((t(i).gt.(15.d0*periode/4.d0)).and.(t(i).le.(15.d0*periode/4.d0+delt))).or.((t(i).gt.(16.d0*periode/4.d0)).and.(t(i).le.(16.d0*periode/4.d0+delt)))) then
	
   !write(*,*) "fort.",2000+i+(ntps+2000)*ph," : wavepacket at time ",t(i)," with phase ",ph
!	open(unit=2000+i+(ntps+2000)*ph)
 !  write(*,*) "fort.",2*ntps+5000+i+(ntps+2000)*ph," : V(t) at time ",t(i)," with phase ",ph
!	open(unit=2*ntps+5000+i+(ntps+2000)*ph)
	do j=1,npos
		xmue=xmu12(j)*champ(i)
		delta=(pot(2,j)-pot(1,j))**2+(2d0*xmue)**2
		delta=dsqrt(delta)	
		vp1(j)=(pot(2,j)+pot(1,j)-delta)*0.5d0
		vp2(j)=(pot(2,j)+pot(1,j)+delta)*0.5d0
		if (j.ge.ideb) then
			vp1reel(j-ideb)=vp1(j)
			vp2reel(j-ideb)=vp2(j)
		endif
		!write(2000+i+(ntps+2000)*ph,*) x(j), dsqrt((dreal(chi1(j))**2+dimag(chi1(j))**2)),dsqrt((dreal(chi2(j))**2+ dimag(chi2(j))**2))
!construction des potentiels E- et E+
		!write(2*ntps+5000+i+(ntps+2000)*ph,*)x(j), vp1(j), vp2(j)
	enddo
	!close(2000+i+(ntps+2000)*ph)
!	close(2*ntps+5000+i+(ntps+2000)*ph)
			!endif

 if ((t(i).gt.(11071.75d0)).and.(t(i).le.(11072.5d0)))then
	do j=1,npos
		write(51,*)x(j),dreal(chi1(j)),dimag(chi1(j))
	enddo

 endif
 if ((t(i).gt.(22143.75d0)).and.(t(i).le.(22144.5d0))) then
	do j=1,npos
		write(50,*)x(j),dreal(chi1(j)),dimag(chi1(j))
	enddo
 endif

 do j=1,npos
	ctemp=chi1init(j)*chi1(j)
	worka(j)=(cdabs(ctemp(j)))**2 !a revoir (calcul d'autocorrelation)
 enddo
 call simpson (npos,delr,worka,auto_correl(i))
 write(52,*)t(i),auto_correl(i)

 maxpot1=maxloc(vp1reel,1)
 minpot2=minloc(vp2reel,1)
 rclapet1(i)=x(ideb+maxpot1)
 rclapet2(i)=x(ideb+minpot2)
 write(200000+ph,*)t(i),rclapet1(i),rclapet2(i)

	call splitop(chi1, chi2, zetdt,pot(1,:),pot(2,:),xmu12, npos, champ(i), delr, massreduite, delt)

!!!!!!!!!!!!!!
       if(i.le.1024) then
           do j=1,npos
            xmu_chi1(i,j)=xmu12(j)*chi2(j) !initialiser xmu_chi1(i,:)
            xmu_chi2(i,j)=xmu12(j)*chi1(j)  !initialiser xmu_chi2(i,:)
            ctemp(j)= dconjg(chi1(j))*xmu_chi1(i,j)+dconjg(chi2(j))*xmu_chi2(i,j) 
!            proj_R(j)=dreal(ctemp(j) )
           enddo
               call simpson(npos,delr,dreal(ctemp),projreal_HHG)
                call simpson(npos,delr,dimag(ctemp),projimag_HHG)
               proj_xmu_chi(i)=dcmplx(projreal_HHG,projimag_HHG)
!!!!!!!!!!!!!!! BOUCLE de t(i) à t(1)=t0 
         do ii=i,1,-1
	   call splitop(xmu_chi1(i,:), xmu_chi1(i,:), zetdt,pot(1,:),pot(2,:),xmu12, npos, champ(ii), delr, massreduite, -delt)! (back-) propagation de xmu_chi1 et  xmu_chi1
         enddo
       endif
!!!!!!!!!!!!!!

!********************************************************************
!	    Calcul de probabilité de dissociation
!********************************************************************
	 lieprob = 0.d0

            do n = 0, v
	       do j = 1, npos
                  proj(j) = ep(n,j)*dreal(chi1(j))
	       end do
               call simpson(npos,delr,proj,projreal)
               do j = 1, npos
                  proj(j) = ep(n,j)*dimag(chi1(j))
	       end do
	   call simpson(npos,delr,proj,projimag)
	   lieprobv=(projreal**2 + projimag**2)
		write(100+n+20*ph,*) t(i), lieprobv, 1-lieprobv
	        lieprob = lieprob + lieprobv
		if (i.eq.1) then		
			write(49,*)n,lieprobv
		endif
            end do
            write(1000+ph,*) t(i), lieprob

 enddo


!********************************************************************
!	    Calcul du spectre `HHG`: A TRAVAILLER 
!********************************************************************
domega= 2.d0*pi/(1024*delt)
open(unit=1010,name="HHG.dat",status='replace') 
call fourier(proj_xmu_chi,1024,1,delt,domega)
!! ATTENTION: Il faut replacer les elements du vecteur-resultat de cette FT
HHG_simple_v=cdabs(proj_xmu_chi)**2
do j = 1, npos	
     call fourier(xmu_chi1(:,j),1024,1,delt,domega)	
     call fourier(xmu_chi2(:,j),1024,1,delt,domega)
	proj_HHG(:,j)=cdabs(xmu_chi1(:,j))**2+cdabs(xmu_chi2(:,j))**2
enddo
do i=1,1024
 !
         if(i.le.512)then
            temp_omega  = (i-1) * domega! ( (i-1)+514 )* domega
         else
             temp_omega  = (-(512-i+1)+1024) * domega !( -(512-i)+1 )* domega
         endif
	   call simpson(npos,delr,proj_HHG(i,:),HHG(i))  !!! ATTENTION: Il faut replacer les elements du vecteur-resultat des FT s
 !	    	Écriture  du spectre `HHG` (version longue et version simplifiée) dans un fichier.
	  write(1010,*)temp_omega ,HHG(i),HHG_simple_v(i)
 enddo
  close(1010)

end program travail


!*******************************************************************  
!
!	Calcul du paquet d'onde initial
!
!*******************************************************************


     subroutine eval(cw1, cw2, delr, rdeb,p0, rc0, alpha, npos)

      double complex cw1(npos), cw2(npos)
      double precision delr, rdeb, p0, rc0, alpha
      integer npos
      double precision pi, r
      double complex cnul, cim, cpoi, cval, arg
      integer l

      cnul = dcmplx(0.d0,0.d0)
      cim = dcmplx(0.d0,1.d0)
      pi = 3.141592654d0
      cpoi = cdsqrt(cdsqrt(dcmplx(2.d0*alpha/pi,0.d0)))
      r = rdeb-delr
      do l = 1, npos
         r = r + delr
         arg = dcmplx(-alpha*(r-rc0)**2, p0*(r-rc0))
         cval = cpoi*cdexp(arg)
         cw1(l) = cval
         cw2(l) = cnul

      end do
 
      return
      end

!***************************************************************
! 		Calcul norme d'une fonction complexe
!***************************************************************
	subroutine simps(func, vint, delti, npl)

      integer j, npl
      double complex func(npl)
      double precision  vint, delti
      
      vint=0d0
      do j = 1, npl-1
         vint=vint+delti*sqrt(cdabs(func(j))**2) 
      end do
      return
      end

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


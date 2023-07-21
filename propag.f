!****************************************************************************
!Application de la fonction de propagation sur une fonction psi par 
!multiplication du produit d'exponentielles issu du split-operator.
!****************************************************************************	

 subroutine propag(psi, elaser,xmu,deltat,r,p,t)
 integer :: i,l,c	
 double precision :: deltat, npts,h,massreduite
 double precision :: psi, r,p, t,xmu, elaser,opcin
 double precision, dimension(1:2,1:2) :: ident,oppot

 h=1
 massreduite=918.0762887608628d0
 oppot(1,1)=
 oppot(2,2)=
 oppot(1,2)=xmu*elaser
 oppot(2,1)=xmu*elaser
 opcin=-h**2/(2*mass) 
 do l=1,2
	do c=1,2
	 if l=c then ident(l,c)=1
	 else ident(l,c)=0
	enddo
 enddo


 return	
 end

 


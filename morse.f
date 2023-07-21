      double precision function morse(diss,smalla,xmu,requ,r,nu)
!
      double precision diss,smalla,xmu
      double precision biga,bigc,enu,alpha
      integer nu,m
      double precision x,arg,r,requ
      double precision laguerrel,factrl
      double precision norm
	
         if(r.lt.64.d0)then
         biga = (dsqrt(2.d0*xmu))/smalla 
	 bigc = biga*dsqrt(diss) 
         enu = -((bigc-nu-.5d0)**2)/biga**2 
         alpha = bigc-nu-0.5d0   
	 arg=dexp(-smalla*(r-requ)) 
         x=2.d0*bigc*arg 
         m = 2*(idint(bigc)-nu)-1 
	 morse = laguerrel(m,nu,x)
	 morse = morse * (x**(idint(bigc)))
         morse = morse / (x**nu)
         morse = morse / dsqrt(x)
         morse = morse * dexp(-x/2.d0)
	
         norm = smalla * factrl(nu) * (2.d0*bigc - 2.d0*nu - 1.d0)
         norm = norm / factrl(m+nu)
         norm = dsqrt(norm)
         morse = morse*norm
      else
         morse = 0.d0
      endif
      return
      end

      double precision function laguerrel(a,n,x)
      integer a,n,j
      double precision x,lagtmp
      double precision factrl
	
      lagtmp=0.d0 
      do 22 j=0,n
      lagtmp=lagtmp+(1.d0/factrl(j))*((-x)**j)*factrl(n+a)/
     $(factrl(n-j)*factrl(j+a))
   22 continue
      laguerrel=lagtmp
      return
      end

      double precision function factrl(n)
      integer n,j
      double precision a(100)
	
      data ntop,a(1)/0,1./
      if(n.lt.0) then
         pause 'negative factorial'
      else if(n.le.ntop) then
         factrl=a(n+1)
      else if(n.le.100) then
         do 11 j=ntop+1,n
            a(j+1)=j*a(j)
   11    continue
         ntop=n
         factrl=a(n+1)
      else
         factrl = dsqrt(6.283185307d0*n)*(n*dexp(-1.d0))**n
      endif
      return
      end




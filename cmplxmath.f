      subroutine  ZSCALE(n,za,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex za,zx(*)
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = za * zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = za * zx(i)
   30 continue
      return
      end

      subroutine ZAXPY(n,za,zx,incx,zy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),za
      integer i,incx,incy,ix,iy,n
      if(n.le.0)return
      if (cdabs(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end

      subroutine  ZCOPY(n,zx,incx,zy,incy)
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 4/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zx(i)
   30 continue
      return
      end

      subroutine ZVEA(n,zx,incx,zy,incy,zz,incz)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),zz(*)
      integer i,incx,incy,incz,ix,iy,iz,n
      if(n.le.0)return
      if (incx.eq.1.and.incy.eq.1.and.incz.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      iz = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      if(incz.lt.0)iz = (-n+1)*incz + 1
      do 10 i = 1,n
        zz(iz) = zx(ix) + zy(iy)
        ix = ix + incx
        iy = iy + incy
        iz = iz + incz
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zz(i) = zx(i) + zy(i)
   30 continue
      return
      end

      subroutine ZVES(n,zx,incx,zy,incy,zz,incz)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),zz(*)
      integer i,incx,incy,incz,ix,iy,iz,n
      if(n.le.0)return
      if (incx.eq.1.and.incy.eq.1.and.incz.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      iz = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      if(incz.lt.0)iz = (-n+1)*incz + 1
      do 10 i = 1,n
        zz(iz) = zx(ix) - zy(iy)
        ix = ix + incx
        iy = iy + incy
        iz = iz + incz
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zz(i) = zx(i) - zy(i)
   30 continue
      return
      end

      subroutine ZVEM(n,zx,incx,zy,incy,zz,incz)
c
c     constant times a vector plus a vector.
c     jack dongarra, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*),zy(*),zz(*)
      integer i,incx,incy,incz,ix,iy,iz,n
      if(n.le.0)return
      if (incx.eq.1.and.incy.eq.1.and.incz.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      iz = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      if(incz.lt.0)iz = (-n+1)*incz + 1
      do 10 i = 1,n
        zz(iz) = zx(ix) * zy(iy)
        ix = ix + incx
        iy = iy + incy
        iz = iz + incz
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zz(i) = zx(i) * zy(i)
   30 continue
      return
      end

program test
implicit none

integer :: i,n
parameter(n=8000)
double precision :: g(n)
do i=1,n
	g(i)=dexp(-((i*6000d0/8000d0-4d0*3.141592654d0/2.15d-2)**2)/(2d0*(2d0*3.141592654d0/2.15d-2/2.3548d0)**2))
write(56,*)i*6000d0/8000d0,g(i)
enddo





end program test



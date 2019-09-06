program splitfft
implicit none

complex, dimension(256):: psi
real, dimension(256):: v,x
integer i,j,h
real a,x0,k0,m,xmin,dx,dt
complex:: iota =cmplx(0.,1.)
character*5 p
real:: pi= acos(-1.0)

a=20.0
x0=-0.5
k0=50.0
m=14500.0
xmin=-2.0
dx=0.02
dt=0.1

! defining x grid ---
x(1)=xmin
do i=2,256
	x(i) = x(i-1) + dx
enddo

! defining potential--
do i=1,256
	if(x(i) .ge. 0.0 .and. x(i) .le. 0.5) then 
                v(i) = 0.1
        else
                v(i) = 0.0
        endif
enddo

! initial psi values---
open(1,file='40001')
do h=1,256
	psi(h) = sqrt(sqrt(2.*a/pi)) * exp(-a*(x(h)-x0)**2) * exp(iota*(k0*(x(h)-x0)))
	write(1,*) x(h), (abs(psi(h))**2)
enddo

do j=2,5000
	do h=1,256
		psi(h) = exp(-iota*v(h)*dt/2.) * psi(h)
	enddo

	call fft(psi,256,+1)

	do h=1,256
		psi(h) = exp(-iota*dt*(k(h)**2.)/(2.*m)) * psi(h)
	enddo

	call fft(psi,256,-1)
	psi = psi/256.

	do h=1,256
		psi(h) = exp(-iota*v(h)*dt/2.) * psi(h)
	enddo

	if(mod(j,100).eq.1) then
		write(p,'(I5)')40000+j
		open(unit=j, file=p )
		do h=1,256
			write(j,*) x(h), (abs(psi(h))**2)
		enddo
	endif
enddo

contains

real function k(h1)
integer h1
if (h1 .le. 128)then  !N/2 - 1 as N = 256 points of x axis
k = 2*pi*(h1-1)/(256.*dx)
else
k = 2*pi*(h1-1-256)/(256.*dx)
endif
return
end function k

end
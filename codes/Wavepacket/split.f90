program tdseft
implicit none

complex(kind = 8), dimension(256,5000):: psi
real(kind = 8), dimension(256):: v,x
complex(kind = 8), dimension(256):: p1,p2
integer i,j,h
real(kind = 8) a,x0,k0,m,xmin,dx,dt
complex(kind = 8):: iota= cmplx(0.,1.)
character*5 p
real:: pi= acos(-1.0)



open(10, file='psi')
open(20, file='psi2')

a=20.0
x0=-0.5
k0=50.0
m=14500.0
xmin=-2.0
dx=0.02
dt=0.1

! defining x and v grid ---
x(1)=xmin
do i=2,256
	x(i) = x(i-1) + dx
enddo

do i=1,256
	if(x(i) .ge. 0.0 .and. x(i) .le. 0.5) then 
                v(i) = 0.1
        else
                v(i) = 0.0
        endif
enddo

! initial value of psi ---
do i=1,256
	psi(i,1) = sqrt(sqrt(2*a/pi))*exp(iota*k0*(x(i)-x0))*exp(-a*((x(i)-x0)**2))
enddo

! psi value for rest time ---
do j=1,4999
	! fourier transformation ---
	do h=1,256
		p1(h) = cmplx(0.,0.)
		do i=1,256
			p1(h) = p1(h) + exp(-iota*v(i)*dt/2.)*psi(i,j)*exp(-iota*k(h)*x(i))/sqrt(256.)
		enddo
	enddo

	! inverse fourier transformation ---
	do i=1,256
		p2(i) = cmplx(0.,0.)
		do h=1,256
			p2(i) = p2(i) + exp(-iota*(k(h)**2)*dt/(2.*m))*p1(h)*exp(iota*k(h)*x(i))/sqrt(256.)
		enddo
	enddo

	do i=1,256
		psi(i,j+1) = exp(-iota*v(i)*dt/2.) * p2(i)
	enddo
	
	write(*,*) j
enddo

! writing values --
do j=1,5000,100
	write(p,'(I5)')30000+j
	open(unit=j, file=p ) 
	do i=1,256
		write(j,*) x(i),(abs(psi(i,j))**2)
	enddo
	close(j)
enddo

contains


real(kind=8) function k(h1)
integer h1
if (h1 .le. 128)then  !N/2 - 1 as N = 256 points of x axis
k = 2*pi*(h1-1)/(256.*dx)
else
k = 2*pi*(h1-1-256)/(256.*dx)
endif
return
end function k

end
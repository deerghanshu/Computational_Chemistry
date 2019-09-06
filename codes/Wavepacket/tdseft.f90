program tdseft
implicit none

complex, dimension(256,5000):: psi,dp
complex, dimension(256):: p1,p2
real, dimension(256):: v,x,z
integer i,j,h
real a,x0,k0,m,xmin,dx,dt
character*5 p
real:: pi = acos(-1.0)
complex:: iota=cmplx(0.0,1.0)

a=20.0
x0=-0.5
k0=50.0
m=14500.0
xmin=-2.0
dx=0.02
dt=0.1

x(1)=xmin
do i=2,256
	x(i) = x(i-1) + dx
enddo

do i=1,256
	if(x(i).ge.0) then 
                v(i) = 0.1
        else
                v(i) = 0.
        endif
enddo


! initial value of psi ---
do i=1,256
	psi(i,1) = sqrt(sqrt(2*a/pi))*exp(iota*k0*(x(i)-x0))*exp(-a*((x(i)-x0)**2))
enddo

! psi value for rest time ---
do j=1,5000
	! fourier transformation ---
	do h=1,256
		p1(h) = cmplx(0.,0.)
		do i=1,256
			p1(h) = p1(h) + psi(i,j)*exp(-iota*k(h)*x(i))/sqrt(256.)
		enddo
	enddo

	! inverse fourier transformation ---
	do i=1,256
		p2(i) = cmplx(0.,0.)
		do h=1,256
			p2(i) = p2(i) + (-k(h)**2)*p1(h)*exp(iota*k(h)*x(i))/sqrt(256.)
		enddo
	enddo

	do i=1,256
		dp(i,j) = p2(i)
	enddo
	
	! psi value for time t=2 ----
	if(j.eq.1) then
		do i=1,256
			psi(i,2) = psi(i,1) + dt*iota*( dp(i,1)/(2.*m) - v(i)*psi(i,1) )
		enddo
	endif

	! psi value for rest time values--- 
	if(j.ge.2) then
		do i=1,256
			psi(i,j+1) = psi(i,j-1) - (2.*dt)*iota*( -dp(i,j)/(2.*m) + v(i)*psi(i,j) )
		enddo
	endif
	
	write(*,*) j
enddo




do j=1,5000,100
	write(p,'(I5)')20000+j
	open(unit=j, file=p ) 
	do i=1,256
		write(j,*) x(i),(real(psi(i,j))**2+aimag(psi(i,j))**2)
	enddo
	close(j)
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

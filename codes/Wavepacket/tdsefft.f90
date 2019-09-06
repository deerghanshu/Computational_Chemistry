program tdseft
implicit none

complex, dimension(256,5000):: psi
complex, dimension(256):: ps
real, dimension(256):: v,x
integer i,j,h
real a,x0,k0,m,xmin,dx,dt
real:: pi = acos(-1.0)
complex:: iota=cmplx(0.0,1.0)
character*5 p

a=20.0
x0=-0.5
k0=50.0
m=14500.0
xmin=-2.0
dx=0.02
dt=0.1

! defining x grid
x(1)=xmin
do i=2,256
	x(i) = x(i-1) + dx
enddo

! defining potential
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
do j=1,5000
	do i = 1,256
		ps(i) = psi(i,j)
	enddo
	! fourier transformation ---
	call fft(ps,256,+1)
	do h=1,256
		ps(h) = -ps(h)*(k(h)**2)
	enddo
	! inverse fourier transformation ---
	call fft(ps,256,-1)
	ps = ps/256.

	! psi value for time t=2 ----
	if(j.eq.1) then
	do i=1,256
		psi(i,2) = psi(i,1) + dt*iota*( ps(i)/(2.*m) - v(i)*psi(i,1) )
	enddo
	endif

	! psi value for rest time values--- 
	if(j.ge.2) then
	do i=1,256
		psi(i,j+1) = psi(i,j-1) - (2.*dt)*iota*( -ps(i)/(2.*m) + v(i)*psi(i,j) )
	enddo
	endif
		
	if(mod(j,100).eq.1) then
		write(p,'(I5)')50000+j
		open(unit=j, file=p )
		do h=1,256
			write(j,*) x(h), (abs(psi(h,j))**2)
		enddo
	endif

	write(*,*) j
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

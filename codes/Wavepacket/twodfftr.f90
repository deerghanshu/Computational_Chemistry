program twodfft
implicit none

complex, dimension(128,128):: psi
complex, dimension(128):: psix,psiy,psi1,psi2
real, dimension(128,128)::v
real, dimension(128):: x,y,vx,vy
integer i,j,h,t
real a,x0,k0,m,xmin,dx,dt,dy,y0,ymin
complex:: iota =cmplx(0.,1.)
character*5 p
real:: pi= acos(-1.0)

a=20.0
x0=-0.5
y0=-0.5
k0=50.0
m=14500.0
xmin=-2.0
ymin=-2.0
dx=0.04
dy=0.04
dt=0.1

open (10, file = 'a')
! defining x and y grid ---
x(1)=xmin
do i=2,128
	x(i) = x(i-1) + dx
enddo
y(1)=xmin
do j=2,128
	y(j) = y(j-1) + dy
enddo

! defining potential--
do i=1,128
	if(x(i) .ge. 0.0 .and. x(i) .le. 0.5) then 
                vx(i) = 0.8
        else
                vx(i) = 0.0
        endif
enddo
do j=1,128
	if(y(j) .ge. 0.0 .and. y(j) .le. 0.5) then 
                vy(j) = 0.8
        else
                vy(j) = 0.0
        endif
enddo
do i=1,128
	do j=1,128
		v(i,j) = vx(i)*vy(j)
		write(10,*) x(i),y(j), v(i,j) 
	enddo
enddo

! initial psi values---
open(0,file='60000')
do i=1,128
	psix(i) = sqrt(sqrt(2.*a/pi)) * exp(-a*(x(i)-x0)**2) * exp(iota*(k0*(x(i)-x0)))
enddo
do j=1,128
	psiy(j) = sqrt(sqrt(2.*a/pi)) * exp(-a*(y(j)-y0)**2) * exp(iota*(k0*(y(j)-x0)))	
enddo
do i=1,128
	do j=1,128
		psi(i,j) = psix(i)*psiy(j)
		write(0,*) x(i),y(j), (abs(psi(i,j))**2)
	enddo
enddo



do t=2,5000
	! psi definition---
	do j=1,128
		do i=1,128
			psi1(i) = exp(-iota*v(i,j)*dt/2.) * psi(i,j)
		enddo	

		call fft(psi1,128,+1)

		do i=1,128
			psi1(i) = exp(-iota*dt*(k(i)**2.)/(2.*m)) * psi1(i)
		enddo

		call fft(psi1,128,-1)
		psi1 = psi1/128.

		do i=1,128
			psi(i,j) = exp(-iota*v(i,j)*dt/2.) * psi1(i)
		enddo
	enddo

	do i=1,128	
		do j=1,128
			psi2(j) = exp(-iota*v(i,j)*dt/2.) * psi(i,j)
		enddo

		call fft(psi2,128,+1)

		do j=1,128
			psi2(j) = exp(-iota*dt*(k(j)**2.)/(2.*m)) * psi2(j)
		enddo
	
		call fft(psi2,128,-1)
		psi2 = psi2/128.
	
		do j=1,128
			psi(i,j) = exp(-iota*v(i,j)*dt/2.) * psi2(j)
		enddo
	enddo

	if(mod(t,100).eq.0) then
		write(p,'(I5)')60000+t
		open(unit=t, file=p )
		do i=1,128
			do j=1,128
				write(t,*) x(i), y(j), (abs(psi(i,j))**2)
			enddo	
				write(t,*)
		enddo
	endif
enddo

contains

real function k(h1)
integer h1
if (h1 .le. 64)then  
k = 2*pi*(h1-1)/(128.*dx)
else
k = 2*pi*(h1-1-128)/(128.*dx)
endif
return
end function k

end

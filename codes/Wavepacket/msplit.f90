program tdse
implicit none
complex,dimension(256,5000) :: psi
complex,dimension(256) :: tmp
real,dimension(256) :: x,V
real,dimension(5000) :: t
complex :: z
complex,dimension(256) :: f_k,if_k
integer :: i,j,k
real :: dx,dt,hc,alpha,x0,p0,m,k0

!---constans---------
dx = 0.02
dt = 0.1
m = 14500.0
alpha = 20.
x0 = -0.50
p0 = 50.0
hc = 1.
k0 = p0/hc

z = cmplx(0.0,1.0) !-- this is iota------

!-----x axis-----------
x(1) = -2.0
do i=1,255
x(i+1) = x(i) + dx
enddo

!---- t axis -------------------
t(1) = 0.
do i=1,4999
t(i+1) = t(i) + dt
enddo

!---- potential defi --- V(i)---------
do i=1,256
if(x(i) .lt. 0.0)then
V(i) = 0.0
else if(x(i) .gt. 0.0 .and. x(i) .lt. 0.5)then
V(i) = 0.1
else
V(i) = 0.
endif
enddo


!-----psi matrix first column i.e., at t=0 psi(x)------
open(unit=10,file="fo")
do i=1,256
psi(i,1) = psi_(x(i),x0,k0,alpha)
enddo

do j=1,4999

do i=1,256
psi(i,j)=psi(i,j)*exp(-z*v(i)*dt/2.)
tmp(i)=psi(i,j)
enddo

call fft(tmp,256,1)

do i=1,256
tmp(i)=tmp(i)*exp( -z*(km(i)**2)*dt/(2.*m))
enddo

call fft(tmp,256,-1)
tmp=tmp/256.

do i=1,256
psi(i,j+1)=tmp(i)*exp(-z*v(i)*dt/2.)
enddo
     



        
write(*,*) j
enddo


do i=1,256
write(10,*) abs(psi(i,1))**2,(abs(psi(i,j*100))**2,j=1,50),x(i),v(i)
enddo


contains
complex(kind=8) function psi_(x,x0,k0,alpha)
implicit none
real(kind=8) :: x,k0,x0,alpha,pi
complex(kind=8) :: z !--iota----
z = cmplx(0.,1.)
pi = atan(1.)*4.0
psi_ = sqrt(sqrt(2*alpha/pi))*exp(z*k0*(x-x0))*exp(-alpha*((x-x0)**2))
return
end function psi_


real(kind=8) function km(m)
implicit none
integer :: m,k
real(kind=8) :: pi,L
pi = atan(1.)*4.0
L = 256*0.02 ! N*dx
k = m-1
if (k .le. 127)then  !N/2 - 1 as N = 256 points of x axis
km = 2*pi*k/L
else
km = 2*pi*(k-256)/L
endif
return
end function km




end program tdse

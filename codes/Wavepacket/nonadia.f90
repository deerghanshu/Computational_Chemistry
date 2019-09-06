program nonadia
implicit none

complex*16, dimension(2048,2):: psi
complex*16, dimension(2048):: p11,p12
complex*16, dimension(2,2):: va,vb
complex*16, dimension(2,2,2048):: v
complex*16, dimension(2):: fv
real*8, dimension(2048):: x,v1ad,v2ad,v11,v12,v22,f,d
real*8:: m, beta, sigma, x0, xmin, xmax, n, v1, v2, v3, pi, v4, b1, b2, vasym, v5, b3, vlower, x1, x2 ,x3 ,x4 ,b4, dt, &
		 dx, dxmask, xmask,p0
integer i,j
complex*16:: iota = cmplx(0.0,1.0)
character*5 p

open(10, file = 'pot')
open(20, file = 'damp')
open(30, file = 'v')
open(40, file = 'f')

m = 3474.057
pi = acos(-1.0)
v1 = 4.0167971782296E-2
v2 = 4.79833373E-3
v3 = 9.8998917754E-1
v4 = 1.122019E-2
v5 = 7.9781762366E-1
vlower = 0.0
vasym = 3.61196179E-1
sigma = 0.3
beta = 1.0/(4.0*(sigma**2))
b1 = 5.5
b2 = 4.9818195151
b3 = 2.3471780470
b4 = 1.0487590725
x0 = 9.0
x1 = -4.364721325998E-2
x2 = 5.0012635420E-2
x3 = -7.6042693477E-1
x4 = 8.1790045179E-1
xmin = -45.0
xmax = 45.0
dx = (xmax - xmin)/2048.0
dt = 8.0

do i = 1,2048
	x(i) = xmin + (i-1)*dx
enddo

do i = 1,2048
	f(i) = (1.0 - tanh(b4*(x(i)-x4)))/2.0
	write(40,*) x(i),f(i)
enddo

do i = 1,2048
	v1ad(i) = v1*exp(b1*(x(i)-x1))/(1.0+exp(b1*(x(i)-x1)))**2.0 + v2*exp(b1*(x(i)-x1))/(1.0+exp(b1*(x(i)-x1)))
	v2ad(i) = vasym - v3*exp(b2*(x(i)-x2))/(1.0+exp(b2*(x(i)-x2)))**2.0 - v4*exp(b2*(x(i)-x2))/(1.0+exp(b2*(x(i)-x2))) - v5* &
	exp(b3*(x(i)-x3))/(1.0+exp(b3*(x(i)-x3)))**2.0 - vlower
enddo

do i = 1,2048
	v11(i) = (1.0-f(i))*v1ad(i) + f(i)*v2ad(i)
	v22(i) = (1.0-f(i))*v2ad(i) + f(i)*v1ad(i)
	v12(i) = - sqrt(f(i)*(1.0-f(i)))*(v2ad(i) - v1ad(i))
	write(10,*) x(i),v11(i),v12(i),v22(i)
enddo

do i=1,2048
	if(x(i) .ge. 0 ) then
		xmask = 10
	else
		xmask = -30
	endif
	if(abs(x(i)).ge. abs(xmask)) then 
		dxmask = xmax - abs(xmask)
		d(i) = sin(pi*(abs(xmask)+ dxmask - abs(x(i)))/(2*dxmask))
	else
		d(i) = 1.0
	endif
	write(20,*) x(i), d(i)
enddo


p0 = -sqrt(2*m*abs(0.029-v11(1230)))
open(0, file= '70000')
do i = 1,2048
	psi(i,1) = sqrt(sqrt(1./(2.*pi*(sigma**2)))) * exp(-beta*((x(i)-x0)**2.)) * exp(iota*p0*(x(i)-x0))
	psi(i,2) = cmplx(0.0,0.0)
	write(0,*) x(i), abs(psi(i,1))**2, abs(psi(i,2))**2
enddo

do j=1,999
if (mod(j,100) .eq. 0) then
write(p,'(I5)')70000+j*8
open(j, file=p )
endif
enddo

do i=1,2048
	va(1,1) = exp(-iota*v11(i)*dt/4.)
	va(2,2) = exp(-iota*v22(i)*dt/4.)
	va(1,2) = cmplx(0.0,0.0)
	va(2,1) = cmplx(0.0,0.0)
	
	vb(1,1) = cmplx(cos(v12(i)*dt/2.),0.0)
	vb(1,2) = cmplx(0.0,-sin(v12(i)*dt/2.))
	vb(2,1) = cmplx(0.0,-sin(v12(i)*dt/2.))
	vb(2,2) = cmplx(cos(v12(i)*dt/2.),0.0)
	
	v(:,:,i) = matmul(vb,va)
	v(:,:,i) = matmul(va,v(:,:,i))
	write(30,*) x(i),abs(v(:,:,i))**2
enddo

do j=1,999
	write(*,*) j
	do i=1,2048
		fv = psi(i,:)
		fv = matmul(v(:,:,i),fv)
		p11(i) = fv(1)
		p12(i) = fv(2)
		!print*, abs(p11(i))**2, abs(p12(i))**2	
	enddo
	!print*, abs(p11)**2
	call fft(p11,2048,+1)
	call fft(p12,2048,+1)
	!print*, abs(p11)**2, abs(p12)**2
	do i=1,2048
		p11(i) = exp(-iota*dt*(k(i)**2.)/(2.*m)) * p11(i)
		p12(i) = exp(-iota*dt*(k(i)**2.)/(2.*m)) * p12(i)
	enddo
	!print*, abs(p11)**2, abs(p12)**2
	call fft(p11,2048,-1)
	call fft(p12,2048,-1)
	!print*, abs(p11)**2, abs(p12)**2
	p11 = p11/2048.
	p12 = p12/2048.
	
	do i=1,2048
		fv(1) = p11(i)
		fv(2) = p12(i)
		fv = matmul(v(:,:,i),fv)
		psi(i,:) = fv*d(i)	
		if(mod(j,100) .eq. 0) then
			write(j,*) x(i), abs(psi(i,1))**2, abs(psi(i,2))**2
		endif
	enddo
enddo


contains

real*8 function k(h1)
implicit none
integer h1,m
real*8 ::  dx
m = h1-1
dx = 90.0/2048.0
if (m .le. 1023)then
k = 2*pi*(m)/(2048.*dx)
else
k = 2*pi*(m-2048)/(2048.*dx)
endif
return
end function k

end


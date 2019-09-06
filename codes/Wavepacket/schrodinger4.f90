program schrodinger
implicit none
real a,x0,p0,m,xmin,dx,dt,t,k0,n
complex iota
character(len=90) :: filename
integer j,q,i
real, dimension(256) :: x,v,k,k2
complex, dimension(256) :: psi0,psi1,psi2,psi
real , parameter :: pi = acos(-1.0)
iota = (0.0,1.0)
t = 5000
dt = 0.1
dx = 0.02
xmin = -2.0
m = 14500.0
p0 = 50.0
x0 = -0.5
k0 = p0
a = 20.0
n= 256.0
open(unit=10,file="psi-0001")
do j=1,256
	x(j) = xmin + ((j-1)*dx)	
	if(x(j).LT.0) then
		v(j) = 0
	else
		v(j) = 1
	end if
end do

do i =1,256
	if(i.LE.128) then
		k(i) = 2*pi*(i-1)/(n*dx)
		k2(i) = k(i)**2
	else
		k(i) = 2*pi*(i-1-n)/(n*dx)
		k2(i) = k(i)**2
	end if	
end do

do j =1,256
	psi0(j) = ((2*a/pi)**0.25) * exp(iota*p0*(x(j)-x0)) * exp((-a)*((x(j)-x0)**2))
	psi(j) = psi0(j)
end do

call dft(psi,x,k,k2)
psi = psi/(2.0*m)

do j=1,256
	psi1(j) = psi0(j) + iota*dt*(psi(j) - v(j)*psi0(j))
	write(10,*)x(j),v(j),realpart(psi0(j))**2 + imagpart(psi0(j))**2,realpart(psi1(j))**2 + imagpart(psi1(j))**2
end do

do j =1,256
	psi(j) = psi1(j)
end do

do q = 3,5000
	write(filename,30)q
	30 format('psi-',i4.4)
	open(unit=20,file=filename)
	call dft(psi,x,k,k2)	
	psi = psi/(2.0*m)
	do j=1,256
		psi2(j) = psi0(j) + (2*iota*dt*(psi(j) - v(j)*psi1(j)))
		write(20,*)x(j),v(j),realpart(psi2(j))**2 + imagpart(psi2(j))**2
	end do
	do j = 1,256
		psi(j) = psi2(j)
	end do
	do j = 1,256
		psi0(j) = psi1(j)
	end do
	do j = 1,256
		psi1(j) = psi2(j)
	end do	
end do
end

subroutine dft(psi,x,k,k2)
implicit none
integer i,j
real,dimension(256) :: x,k,k2
real pi
complex eye
complex,dimension(256) :: phi,psi
eye = (0.0,1.0)
pi = acos(-1.0)
do i = 1,256
	phi(i) = (0.0,0.0)
	do j = 1,256
		phi(i) = phi(i) + psi(j)*exp(-eye*k(i)*x(j))
	end do
end do
phi = phi/sqrt(256.0)

do i =1,256
phi(i) = -k2(i)*phi(i)
end do

do j=1,256
	psi(j) = (0.0,0.0)
	do i = 1,256
		psi(j) = psi(j) + phi(i)*exp(eye*k(i)*x(j))
	end do
end do
psi = psi/sqrt(256.0)
return
end

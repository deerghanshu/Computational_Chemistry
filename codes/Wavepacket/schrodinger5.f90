program schrodinger
implicit none
real a,x0,p0,m,xmin,dx,dt,t,k0,n
complex iota
character(len=90) :: filename
integer j,m1,q
real, dimension(256) :: x,v,k,k2
complex, dimension(256) :: psi0,psi1,psi11,psisquare0,psisquare1,psi2,test
real , parameter :: pi = acos(-1.0)

iota = (0.0,1.0)
t = 5000
dt = 0.1
dx = 0.02
xmin = -2.0
m = 14500.0
p0 = 20.0
x0 = -0.5
k0 = p0
a = 20.0
n= 256.0
open(unit=10,file="output3")
do j=1,256
	x(j) = xmin + ((j-1)*dx)
	if(j.LE.128) then
		k(j) = 2*pi*j/(n*dx)
		k2(j) = k(j)**2
	else
		k(j) = 2*pi*(j-n)/(n*dx)
		k2(j) = k(j)**2
	end if	
	if(x(j).LT.0) then
		v(j) = 0
	else
		v(j) = 1
	end if
	psi0(j) = initial(x(j))
	psisquare0(j) = real(psi0(j))**2 + aimag(psi0(j))**2
	test(j) = psi0(j)
end do
call dft(psi0,x,k,k2)
do j=1,256
	psi11(j) = psi0(j) + (iota*dt*((psi1(j)/(2*m)) - (v(j)*psi0(j))))
	psisquare1(j) = real(psi11(j))**2 + aimag(psi11(j))**2
	write(10,*)x(j),v(j),real(psisquare0(j)),real(psisquare1(j))
end do

do q=3,70
	call dft(psi11,x,k,k2)
	write(filename,30)q
	30 format('psi-',i2)
	open(unit=20,file=filename)	
	do j=1,256
		psi2(j) = psi0(j) + (2*iota*dt*((psi1(j)/(2.0*m)) - (v(j)*psi11(j))))
		write(20,*)x(j),v(j),real(real(psi2(j))**2 + aimag(psi2(j))**2)
	end do
	psi0 = psi11
	psi11 = psi2
end do

contains
complex function initial(x)
implicit none
real x
initial = ((2*a/pi)**0.25) * exp(iota*k0*(x-x0)) * exp(-a*((x-x0)**2))
return
end function
end

subroutine dft(psi0,x,k,k2)
implicit none
integer i,j
real x,k,k2,pi
complex ftmp,psi2,phi,psi0,psi1
complex eye
dimension phi(256),k(256),k2(256),x(256),ftmp(256),psi0(256),psi1(256)
eye = (0.0,1.0)
pi = acos(-1.0)
do i = 1,256
	phi(i) = 0.0
	do j = 1,256
		phi(i) = phi(i) + psi0(j)*exp(-eye*k(i)*x(j))
	end do
end do
phi = phi/sqrt(256.0)

do i =1,256
	phi(i) = -k2(i)*phi(i)
end do

do j=1,256
	psi1(j) = 0.0
	do i = 1,256
		psi1(j) = psi1(j) + phi(i)*exp(eye*k(i)*x(j))
	end do
end do
psi1 = psi1/sqrt(256.0)
return
end

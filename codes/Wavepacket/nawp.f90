program nawp
implicit none
real*8 a,x0,m,xmin,dx,dt,t,k0,xmax,sigma,beta,x1,x2,x3,x4,b1,b2,vasym,b3,v1,v2,v3,v4,v5,b4,vlower,dxmask,p0,v1adp0,v2adp0,v11p0
complex iota
character(len=90) :: filename
integer j,q,i,n
real*8, dimension(2048) :: x,v,k,k2,d,v11,v12,v22,v1ad,v2ad
complex*32, dimension(2048) :: timestamp,timestamp2
complex*32, dimension(2,2048) :: psi0,psi1,psi2,psit,psi3,psi4,psi5,psi6
real, parameter :: pi = acos(-1.0)
iota = (0.0,1.0)
t = 8000
dt = 0.8
xmin = -45.0
xmax = 45.0
n = 2048
dx = (xmax-xmin)/(n-1)
m = 3474.057
sigma = 0.3
beta = 1/(4.0*(sigma**2))
x0 = 9.0
v1 = 0.040167971782296
b1 = 5.5
x1 = -0.04364721325998 
v2 = 0.00479833373
vasym = 0.361196179
v3 = 0.98998917754
b2 = 4.9818195151
x2 = 0.050012635420
v4 = 0.01122019
v5 = 0.79781762366
b3 = 2.3471780470
x3 = -0.76042693477	
vlower = 0.0

!x-grid calculation
open(unit=10,file="psi-0001")
do j=1,2048
	x(j) = xmin + ((j-1)*dx)
end do

!k-grid calculation
do j =1,2048
	if(j.LE.1024) then
		k(j) = 2*pi*(j-1)/(n*dx)
		k2(j) = k(j)**2
	else
		k(j) = 2*pi*(j-1-n)/(n*dx)
		k2(j) = k(j)**2
	end if	
end do

do j=1,2048
	v1ad(j) = (v1*exp(b1*(x(j)-x1)))/(1+exp(b1*(x(j)-x1)))**2 + (v2*exp(b1*(x(j)-x1)))/(1+exp(b1*(x(j)-x1)))
	v2ad(j) = vasym - ((v3*exp(b2*(x(j)-x2)))/((1+exp(b2*(x(j)-x2)))**2)) - ((v4*exp(b2*(x(j)-x2)))/(1+exp(b2*(x(j)-x2))))
	v2ad(j) = v2ad(j) - ((v5*exp(b3*(x(j)-x3)))/((1+exp(b3*(x(j)-x3)))**2)) - vlower
end do

!Potentials and coupling
do j =1,2048
	v11(j) = ((1-f(x(j)))*v1ad(j)) + (f(x(j))*v2ad(j))
	v22(j) = (f(x(j))*v1ad(j)) + ((1-f(x(j)))*v2ad(j))
	v12(j) = -sqrt(f(x(j))*(1-f(x(j)))) * (v2ad(j) - v1ad(j))	
end do

!Damping function
do j =1,2048
	if(x(j).LE.-30) then
		dxmask = -45.0 + 30.0
		d(j) = sin(pi*(-30 + dxmask - x(j))/(2.0*dxmask))
	else if(x(j).GE.10) then
		dxmask = 45 - 10
		d(j) = sin(pi*(10 + dxmask - x(j))/(2.0*dxmask))
	else
		d(j) = 1.0
	end if
end do

!Calculating p0
v1adp0 = ((v1*exp(b1*(x0-x1)))/((1+exp(b1*(x0-x1)))**2)) + ((v2*exp(b1*(x0-x1)))/(1+exp(b1*(x0-x1))))
v2adp0 = vasym - ((v3*exp(b2*(x0-x2)))/((1+exp(b2*(x0-x2)))**2)) - ((v4*exp(b2*(x0-x2)))/(1+exp(b2*(x0-x2))))
v2adp0 = v2adp0 - ((v5*exp(b3*(x0-x3)))/((1+exp(b3*(x0-x3)))**2))
v11p0 = ((1-f(x0))*v1adp0) + (f(x0)*v2adp0)
p0 = -sqrt(2.0*m*(0.029-v11p0))


!Initial wavepacket
do j =1,2048
	psi0(1,j) = ((1/(2*pi*(sigma**2)))**0.25) * exp(-beta*((x(j)-x0)**2)) * exp(iota*p0*(x(j)-x0))
	psi0(2,j) = (0.0,0.0)
	write(10,*)x(j),abs(psi0(1,j))**2,abs(psi0(2,j))**2
end do

!Wavepacket Propagation
do q = 2,8000
	write(filename,30)q
	30 format('psi-',i4.4)
	open(unit=20,file=filename)
	do j =1,2048
		psi1(1,j) = exp(-iota*v11(j)*dt/4.0) * psi0(1,j)
		psi1(2,j) = exp(-iota*v22(j)*dt/4.0) * psi0(2,j)
	end do
	do j = 1,2048
		psi2(1,j) = (cos(v12(j)*dt/2.0)*psi1(1,j)) + (-iota*sin(v12(j)*dt/2.0)*psi1(2,j))
		psi2(2,j) = (-iota*sin(v12(j)*dt/2.0)*psi1(1,j)) + (cos(v12(j)*dt/2.0)*psi1(2,j))
	end do
	do j =1,2048
		psi3(1,j) = exp(-iota*v11(j)*dt/4.0) * psi2(1,j)
		psi3(2,j) = exp(-iota*v22(j)*dt/4.0) * psi2(2,j)
	end do

	do j =1,2048
		timestamp(j) = psi3(1,j)
	end do
	call fft(timestamp,n,1)
	do j = 1,2048
		timestamp(j) = timestamp(j) * exp(-iota*k2(j)*dt/(2.0*m))
	end do
	call fft(timestamp,n,-1)
	do j = 1,2048
		psi6(1,j) = timestamp(j)/n
	end do

	do j =1,2048
		timestamp(j) = psi3(2,j)
	end do
	call fft(timestamp,n,1)
	do j = 1,2048
		timestamp(j) = timestamp(j) * exp(-iota*k2(j)*dt/(2.0*m))
	end do
	call fft(timestamp,n,-1)
	do j = 1,2048
		psi6(2,j) = timestamp(j)/n
	end do

	do j =1,2048
		psi4(1,j) = exp(-iota*v11(j)*dt/4.0) * psi6(1,j)
		psi4(2,j) = exp(-iota*v22(j)*dt/4.0) * psi6(2,j)
	end do
	do j = 1,2048
		psi2(1,j) = (cos(v12(j)*dt/2.0)*psi4(1,j)) + (-iota*sin(v12(j)*dt/2.0)*psi4(2,j))
		psi2(2,j) = (cos(v12(j)*dt/2.0)*psi4(2,j)) + (-iota*sin(v12(j)*dt/2.0)*psi4(1,j))
	end do
	do j =1,2048
		psi5(1,j) = exp(-iota*v11(j)*dt/4.0) * psi2(1,j) 
		psi5(2,j) = exp(-iota*v22(j)*dt/4.0) * psi2(2,j) 
	end do

	do j = 1,2048
		psit(1,j) = psi5(1,j) *d(j)
		psit(2,j) = psi5(2,j) *d(j)
		write(20,*)x(j),abs(psit(1,j))**2,abs(psit(2,j))**2
	end do

	do j=1,2048		
		psi0(1,j) = psi5(1,j)
		psi0(2,j) = psi5(2,j)
	end do
end do

contains
real*8 function f(x)
implicit none
real*8 x,b4,x4
b4 = 1.0487590725
x4 = 0.81790045179
f = 0.5 * (1 - tanh(b4*(x-x4)))
return
end function
end


!FFT subroutine
subroutine fft(psi,n,isign)
        implicit real*8(a-h,o-z)
        parameter(npts=10,nptx=2*2**npts)
        complex*32 s,v,w,psi(n),cstore
        dimension cstore(nptx)
        complex conjg
        data ntbl/0/
        if(n.gt.ntbl)then
        ntbl=n
        pi=4.d0*atan(1.0d0)
        j=1
        icnt=0
10      s=pi*(0.d0,1.d0)/float(j)
        do 20 k=0,j-1
        icnt=icnt+1
20      cstore(icnt)=exp(s*float(k))
        j=j+j
        if(j.lt.n)go to 10
        end if
        j=1
        do 30 i=1,n
        if(i.le.j)then
        v=psi(j)
        psi(j)=psi(i)
        psi(i)=v
        end if
        m=n/2
25      continue
        if(j.gt.m)then
        j=j-m
        m=m/2
        if(m.ge.1)go to 25
        else
        j=j+m
        end if
30      continue
        j=1
        icnt=0
40       jj=j+j
         do 50 k=1,j
        icnt=icnt+1
        w=cstore(icnt)
        if(isign.lt.0)w=conjg(w)
        do 50 i=k,n,jj
        v=w*psi(i+j)
        psi(i+j)=psi(i)-v
50      psi(i)=psi(i)+v
        j=jj
        if(j.lt.n) go to 40

        return
end

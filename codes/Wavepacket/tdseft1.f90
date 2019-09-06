program tdseft
implicit none

complex(kind = 8), dimension(256):: psi
real(kind = 8), dimension(256):: v,x,z
real(kind = 8), dimension(5000):: t
integer i,j,h
real(kind = 8) a,x0,k0,m,xmin,dx,dt
complex:: iota=cmplx(0.,1.)

open(10, file='psi')

a=20.0
x0=-0.5
k0=20.0
m=14500.0
xmin=-2.0
dx=0.02
dt=0.1

! defining v and x grid ---
x(1)=xmin
V(1)=0
do i=2,256
	x(i) = x(i-1) + dx
	if(x(i).ge.0) then 
                v(i) = 0.1
        else
                v(i) = 0.0
        endif
enddo

psi0(1) = sqrt(sqrt(2.*a/3.14)) * exp(-a*(x(1)-x0)**2) * exp(iota*(k0*(x(1)-x0))) 

! defining  k grid ---
do h =1,256
        if(h .le. 128) then
		k(h) = (2*3.1415)*(h-1)/(256.*dx)
	else
		k(h) = (2*3.1415)*(h-1-256)/(256.*dx)
	endif
	write(10,*) x(h),k(h)
enddo

do h=1,256
	phi(i) = cmplx(0,0)
	do i=1,20
		phi(h) = phi(h) + (1./16)*psi(i,j)*exp(-iota*k(h)*x(i)) 
	enddo
enddo

psi1(1) = psi0(1) + dt*iota*( dp(1,1)/(2*m) - v(i)*psi0(1) )
do j=3,5000
	psi(1,j) = psi(1,j-2) - (2*dt)*iota * ( -dp(1,j-1)/(2*m) + v(1)*psi(1,j-1) )
enddo

do i=1,256
	psi0(i) = sqrt(sqrt(2.*a/3.14)) * exp(-a*(x(i)-x0)**2) * exp(iota*(k0*(x(i)-x0)))
	psi1(i) = psi0(i) + dt*iota*( (1/(2*m))*dp(i,1) - v(i)*psi0(i) )
enddo


do j=3,5000
	do i=2,20
		psi(i,j) = psi(i,j-2) - (2*dt)*iota * ( -1/(2*m)*dp(i,j-1) + v(i)*psi(i,j-1) )
	enddo
enddo

do j=1,5000,100
	write(p,'(I5)')20000+j
	open(unit=j, file=p ) 
	do i=1,20
		write(j,*) x(i),(real(psi(i,j))**2+aimag(psi(i,j))**2)
	enddo
	close(j)
enddo

contains

complex(kind = 16) function dp(i1,j1)
integer i1,j1,m1,h
complex(kind = 16) dd
complex, dimension(20) :: phi,k
do m1=1,20
	if(m1.le.9) then
		k(m1) = (6.28/5.12)*m1
	else
		k(m1) =  (6.28/5.12)*(m1-50)
	endif
enddo


dd = cmplx(0,0)
do m1=1,20
	dd = dd - (k(m1)**2)*(1./16)*phi(m1)*exp(iota*k(m1)*x(i1))
enddo
dp=dd

end function dp

end

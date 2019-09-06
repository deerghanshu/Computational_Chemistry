program tdseft
implicit none
real a,x0,m,xmin,dx,dt,k0
character*5 :: p
integer i,j,m1
real, dimension(256) :: x,v,k
complex, dimension(256:5000) :: psi0,psi1,psi2,psi
dt = 0.1
dx = 0.02
xmin = -2.0
m = 14500.0
k0 = 50.0
x0 = -0.5
a = 20.0


do i=1,256
        x(i) = xmin + ((i-1)*dx)
        if(x(i).lt.0) then
                v(i) = 0
        else
                v(i) = 1
        end if
enddo

do m1 =1,256
        if(m1.le.127) then
                k(m1) = 6.28*(m1-1)/5.12
        else
                k(m1) = 6.28*(m1-1-n)/5.12
        end if
enddo

do i =1,256
        psi0(i) = 1.89 * exp(iota*k0*(x(i)-x0)) *exp((-a)*((x(i)-x0)**2))
	psi(i)= psi0(i)
end do

call dft(psi,x,k)

do i=1,256
        psi1(i) = psi0(i) + cmplx(0.,1.)*dt*(psi(j) - v(j)*psi0(j))
        write(10,*)x(j),v(j),realpart(psi0(j))**2 +imagpart(psi0(j))**2,realpart(psi1(j))**2 + imagpart(psi1(j))**2
	psi(j) = psi1(j)
end do

do q = 3,5000
        write(filename,30)q
        30 format('psi-',i4.4)
        open(unit=20,file=filename)
        call dft(psi,x,k)
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

subroutine dft(psi,x,k)
implicit none
integer i,j
real,dimension(256) :: x,k
complex,dimension(256) :: phi,psi
do m1 = 1,256
        phi(m1) = cmplx(0.0,0.0)
        do i = 1,256
                phi(m1) = phi(m1) + (1./16)*psi(i)*exp(-cmplx(0.,1.)*k(m1)*x(i))
        end do
end do

do i=1,256
        psi(i) = cmplx(0.0,0.0)
        do m1 = 1,256
                psi(i) = psi(i) + (1./16)*(-k(m1)**2)*phi(m1)*exp(cmplx(0.,1.)*k(m1)*x(i))
        enddo
enddo
return
end

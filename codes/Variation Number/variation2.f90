program code
implicit none
real*8,allocatable::s(:,:),h(:,:),b(:,:),a(:,:),c(:,:)
complex*16,allocatable::psi0(:,:),psi(:,:)
real*8,allocatable::x(:),v(:),sigma(:),work(:)
integer::m,n,i,j,ifail,lwork
real*8,allocatable::k(:)
real::dx,PI
complex*16 :: iota
m=5
n=5
lwork=64*n
allocate(k(m),sigma(m),work(lwork),psi(m,121),x(121),v(121),psi0(m,121),a(m,m),s(m,m),h(m,m),b(m,m),c(m,m))
PI = 4.0*datan(1.d0)
!potential
!if(abs(x(i)).gt.1) then
!v(i)=infinity
!else 
!v(i)=0
!end if

!x assignment
dx=2.0/120.0
x(1)=-1.0
do i=1,120
x(i+1)=x(i)+dx
end do
print*,"x=",x(121)


do i = 1,n
  k(i) = (i-1)*PI/2.
enddo

!psi0 definition
iota = cmplx(0.0,1.0)

do j=1,n !state
do i=1,121
psi0(j,i) = 0.5*exp(iota*k(j)*x(i))
enddo
enddo
print*,"psi0=",psi0(1,80)

!matrix s calculation

do i = 1,n
  do j = 1, n
     if (i .eq. j) then
      s(i,j) = 1.0
      h(i,j) =(k(j)**2.) - 0.5
     else
      s(i,j) = 0.0
      h(i,j) = (-sin(k(i) - k(j))/(2.*(k(i) - k(j))))
     endif
  enddo
enddo

call  dsyev('v','u',n,s,n,sigma,work,lwork,ifail)

do i=1,m
	do j=1,n
	a(i,j)=s(i,j)/sqrt(sigma(j))
	end do
end do
b=matmul(transpose(a),matmul(h,a))

call  dsyev('v','u',n,b,n,sigma,work,lwork,ifail)
 c=matmul(a,b)
 c=transpose(c)

do i=1,121
psi(:,i)=matmul(c,psi0(:,i))
write(10,*)x(i),abs(psi(1,i))
write(11,*)x(i),abs(psi(2,i))
write(12,*)x(i),abs(psi(3,i))
write(13,*)x(i),abs(psi(4,i))

end do

write(*,*)"sigma=",sigma

end program
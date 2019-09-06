program deepwell
implicit none

integer n,i,j,iplusj,nrot,ifail,lwork

real*8 smatrix, hmatrix, diag, x, res, pi, eigvec
allocatable::smatrix(:,:),hmatrix(:,:),diag(:), eigvec(:,:)

real*8 k,amat,a1mat, hmat,work
allocatable:: k(:),amat(:,:),a1mat(:,:),hmat(:,:),work(:)
real*8, dimension(0:120) :: xaxis
complex*16, dimension(:),allocatable :: psi
complex*16,dimension(:,:),allocatable :: basis
complex*16 :: iota
pi = 4.0*datan(1.d0)
write (*,*) 'Give nr. of basis states'
read (*, *) n

allocate(smatrix(0:n-1,0:n-1),hmatrix(0:n-1,0:n-1),diag(0:n-1),eigvec(0:n-1,0:n-1))
allocate(amat(0:n-1,0:n-1),a1mat(0:n-1,0:n-1),hmat(0:n-1,0:n-1),work(64*n),k(0:n-1))
allocate(psi(0:n-1))
allocate(basis(0:120,0:n-1))

do i = 0,n-1
  k(i) = (i)*pi/2.
enddo

do i = 0,n-1
  do j = 0, n-1

     if (i .eq. j) then
      smatrix(i,j) = 1.0
      hmatrix(i,j) =(k(j)**2.) - 0.5
     else
      smatrix(i,j) = 0.0
      hmatrix(i,j) = (-sin(k(i) - k(j))/(2.*(k(i) - k(j))))
     endif
  enddo
enddo

!diagonalize S-mat first
lwork=64*n
call dsyev('v','u',n,smatrix,n,diag,work,lwork,ifail)

do i = 0, n-1
 do j = 0, n-1
  amat(i,j) = smatrix(i,j)/sqrt(diag(j))
 enddo
enddo

hmat=matmul(transpose(amat),matmul(hmatrix,amat))

!diagonalize h-mat
call dsyev('v','u',n,hmat,n,diag,work,lwork,ifail)

!carry out Av=c
a1mat=matmul(amat,hmat)

write (*,*) 'Variational     Exact'
do i=0, n-1
  write (*,'(5F12.4)') diag(i), (k(i)*k(i) - 0.5)
enddo


!------------calculating psi-------------------
xaxis(0) = -1.0
do i=0,119
xaxis(i+1) = xaxis(i) + (2.0/121.)
enddo

iota = cmplx(0.0,1.0)
do i=0,120
do j=0,n-1 !state
basis(i,j) = 0.5*exp(iota*k(j)*xaxis(i))
enddo
enddo
open(unit=1230,file="psi")
do i=0,120
psi = matmul(transpose(a1mat),basis(i,:))
write(1230,*) xaxis(i), (abs(psi(j))**2,j=0,n-1)
enddo


end
program deepwell
implicit none

integer n,i,j,iplusj,nrot,ifail,lwork

real*8 smatrix, hmatrix, diag, x, res, pi, eigvec,L,k
allocatable::smatrix(:,:),hmatrix(:,:),diag(:), eigvec(:,:),k(:)

real*8 amat,a1mat, hmat,work
allocatable:: amat(:,:),a1mat(:,:),hmat(:,:),work(:)


write (*,*) 'Give nr. of basis states'
read (*, *) n

allocate(smatrix(0:n-1,0:n-1),hmatrix(0:n-1,0:n-1),diag(0:n-1),eigvec(0:n-1,0:n-1),k(0:n-1))
allocate(amat(0:n-1,0:n-1),a1mat(0:n-1,0:n-1),hmat(0:n-1,0:n-1),work(64*n))

pi = 4.0*datan(1.d0)
l=2.d0

do i=0,n-1
	k(I) = (i)*pi/L
enddo

do i=0,n-1
  do j=0,n-1
    IF (i .eq. j) then
      smatrix(I,J) = 1.d0
      hmatrix(I,J) = k(i)**2.d0 -1.d0/L
          ELSE
      smatrix(i,j) = 0.d0
      hmatrix(i,j) =  - sin(k(i) - k(j))/(l*(k(i)-k(j)))
    ENDIF
  ENDDO
ENDDO

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

! Output the variational eigenvalues to the screen  together with the exact ones
! exact = n^2*pi^2/length^2

write (6,*) 'Variational     Exact'
do i=0, n-1
  write (6,'(5F12.4)') diag(i), k(i)**2 - 0.5
enddo

end

program deepwell
implicit none

integer n,i,j,k,iplusj,nrot,ifail,lwork

real*8 smatrix, hmatrix, diag, x, res, pi, eigvec
allocatable::smatrix(:,:),hmatrix(:,:),diag(:), eigvec(:,:)

real*8 amat,a1mat, hmat,work
allocatable:: amat(:,:),a1mat(:,:),hmat(:,:),work(:)

write (*,*) 'Give nr. of basis states'
read (*, *) n

allocate(smatrix(0:n-1,0:n-1),hmatrix(0:n-1,0:n-1),diag(0:n-1),eigvec(0:n-1,0:n-1))
allocate(amat(0:n-1,0:n-1),a1mat(0:n-1,0:n-1),hmat(0:n-1,0:n-1),work(64*n))

do i=0,n-1
  do j=0,n-1
    iplusj = i+j
    IF (mod(iplusj,2) .eq. 0) then
      smatrix(I,J) = 2.d0/(iplusj+5.d0) - 4.d0/(iplusj+3.d0) + 2.d0/(iplusj+1.d0)
      hmatrix(I,J) = 8*(I+J+2.d0*I*J-1)/(iplusj+3.d0)/(iplusj+1.d0)/(iplusj-1.d0)
          ELSE
      smatrix(i,j) = 0.d0
      hmatrix(i,j) = 0.d0
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
pi = 4.0*datan(1.d0)
write (6,*) 'Variational     Exact'
do i=0, n-1
  write (6,'(5F12.4)') diag(i), (i+1)*(i+1)*pi*pi/4.d0
enddo

end


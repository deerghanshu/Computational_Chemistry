program chevy
implicit none
real,parameter::pi = 3.14159265
integer::i,j,k,n
real::fun,pn,inp,term
real,dimension(:),allocatable::c,x

print*, "Enter the value of n :"
read*,n

n = n+1

allocate(c(n),x(n))

!print*,cos(pi),cos(pi/8),cos(pi/2),cos(pi/3)

j = 0
do i=1,n
x(i) = cos(((2*j + 1)*pi)/(2*n))  !2*n because i did n = n+1 before
j = j+1
enddo

print*,"The roots are :"
write(*,*)(x(i),i=1,n)
write(*,*)

do j=1,n
	c(j) = 0
	do k = 1,n
		if (j.eq.1)then
			c(j) = c(j) + (exp(x(k)))/(n)
		else
			c(j) = c(j) + ((2*exp(x(k)))*cos(n*acos(x(k))))/n
		endif
	enddo
enddo

print*,"The coeffecients are : "
write(*,*)(c(i),i = 1,n)
end program chevy

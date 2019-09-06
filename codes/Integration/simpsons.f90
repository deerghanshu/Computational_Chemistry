PROGRAM simpsons
real, dimension(:),allocatable :: x,y
real a,b,h,s
integer n,i
write(*,*)'To find the value of integral of a function using Simpsons rule. '
write(*,*)
write(*,*)'Enter the limits of integration:'
read(*,*)a,b
write(*,*)'Enter the number of subintervals(an even no.) :'
read(*,*)n
allocate(x(n),y(n))
h=(b-a)/n
x(1)=a
y(1)=f(x(1))
DO i=2,n+1
	x(i)=x(i-1)+h
	y(i)=f(x(i))
end do
s=y(1)+y(n+1)
DO i=2,n,2
	s=s+4*y(i)
end do
do i=3,n-1,2
	s=s+2*y(i)
end do
s=h*s/3
print*, 'The value of the integral is: ',s
STOP
END

real function f(x)
f=exp(x)
return
END

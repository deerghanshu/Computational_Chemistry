program composite_trapezoidal
implicit none
real a,b,sum,integral
integer i,n
a=0
b=1
print*,'Enter n: '
read*,n
sum = 0
do i=1,n-1
	sum = sum + func(a+(i*((b-a)/n)))
end do
integral = ((b-a)/n)*((func(a)+func(b))/2 + sum)
print*,'The integral is ',integral
contains
real function func(x)
implicit none
real x
func = exp(x)
return 
end function
end

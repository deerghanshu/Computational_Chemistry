program euler
implicit none
real x0,y0,h
integer i,n
real, dimension(:),allocatable :: x,y
print*,'Enter value of h'
read*,h
n=(1/h)+1
allocate(x(n),y(n))
x(0) = 0
y(0) = 1
open(unit=10,file="output")
do i=0,n-1
	x(i+1) = x(i) + h
	y(i+1) = y(i) + (h*func(x(i),y(i)))
	write(10,*) x(i),y(i)
end do
contains 
real function func(x,y)
real x,y
func = y
return 
end function
end

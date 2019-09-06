program runge_kutta
implicit none
real h,s1,s2,s3,s4,n
integer i
real,dimension(0:99) :: x,y
print*, 'Enter value of h'
read*,h
print*, 'Enter value of x'
read*,n
i=0
y(0) = 0
x(0) = 0
open(unit=10,file='output3')
do while(x(i).LE.n)
	x(i+1)=x(i) +  h
	s1 = h*func(x(i),y(i))
	s2 = h*func(x(i) + h/2, y(i) + s1/2)
	s3 = h*func (x(i) + h/2, y(i) + s2/2)
	s4 = h*func (x(i) + h, y(i) + s3)
	y(i+1) = y(i) + (s1+(2*s2)+ (4*s3)+ s4)/6
	write(10,*)x(i),y(i)
	i=i+1
end do

contains
real function func(x,y)
implicit none
real x,y
func = x**2 + y**2
end function
end


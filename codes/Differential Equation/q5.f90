program euler
implicit none
real,dimension(5):: x,y
real,parameter:: h=0.25
real sl,sr
real s1,s2,s3,s4
integer i,j


x(1)=0.0
y(1)=0.5
do i=1,4
x(i+1)=x(i)+h
enddo

do i=1,4
y(i+1)=y(i)+h*f(x(i),y(i))
enddo
print*,'For Euler',y(5)

do i=1,4
sl = f(x(i),y(i))
sr = f(x(i+1),y(i)+h*sl)
y(i+1)=y(i)+h*(sl+sr)/2
enddo
print*,'For RK2  ',y(5)

do i=1,4
s1=h*f(x(i),y(i))
s2=h*f(x(i)+h/2,y(i)+s1/2)
s3=h*f(x(i)+h/2,y(i)+s2/2)
s4=h*f(x(i)+h,y(i)+s3)
y(i+1)=y(i)+(s1+2*s2+2*s3+s4)/6
enddo
print*,'For RK4  ',y(5)


contains

real function f(a,b)
real a,b
f=a-b
return
end function f

end
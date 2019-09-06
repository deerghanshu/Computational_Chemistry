program sho_runge
implicit none
real m,k,w,period,h,t,x1,x2,x3,x4,p1,p2,p3,p4,pos0,vel0,x,p,k1pos,m1vel
integer i
real, parameter :: pi = 3.1416
m = 1.0
k = 1.0
w = 1.0
period = 2.*pi
h = 0.02*period
x=1.0
p=0.0
pos0 = x
vel0 = p
open(unit=10,file='rk4')
do i =1,200
	t = h * float(i)
	x1 = p
	p1 = -x
	x2 = p + (p1*h/2)
	p2 = -(x+(x1*h/2))
	x3 = p + (p2*h/2)
	p3 = -(x+(x2*h/2))
	x4 = p + (h*p3)		
	p4 = -(x+(h*x3))
	k1pos = pos0 + (h*((x1+(2*x2)+(2*x3)+x4)/6))
	m1vel = vel0 + (h*((p1+(2*p2)+(2*p3)+p4)/6))
	x = k1pos
	p = m1vel
	write(10,*)x,p,t/period,pos0**2 + vel0**2
	pos0 = k1pos
	vel0 = m1vel	
enddo
end

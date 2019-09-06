program harmonic_oscillator
implicit none
real m,k,w,period,h,t,x1,x2,p1,p2,pos0,vel0,x,p,k1pos,m1vel
integer i
real, parameter :: pi = acos(-1.0)
m = 1.0
k = 1.0
w = 1.0
period = 2.0*pi
h = 0.02*period
open(unit=10,file='rk2')
x=1.0
p=0.0
pos0 = x
vel0 = p
do i =1,200
	t = h * float(i)
	x1 = p
	p1 = -x
	x2 = p + (h*p1)
	p2 = -(x+(h*x1))
	k1pos = pos0 + ((h/2)*(x1+x2))
	m1vel = vel0 + ((h/2)*(p1+p2))
	x = k1pos
	p = m1vel
	write(10,*)x,p,t/period,pos0**2 + vel0**2
	pos0 = k1pos
	vel0 = m1vel
	
enddo
end

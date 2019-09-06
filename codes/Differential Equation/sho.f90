program sho
implicit none
real dt,tp
real,dimension(201):: t,p1,x1,E1,p2,x2,E2,p3,x3,E3
real s1,s2,s3,s4,sl,sr
integer i

open(unit=10, file='p.dat')
open(unit=20, file='E.dat')
open(unit=30, file='x.dat')

tp=2*3.14
t(1) = 0.0
p1(1) = 0.0
x1(1) = 1.0
E1(1) = 0.5
p2(1) = 0.0
x2(1) = 1.0
E2(1) = 0.5
p3(1) = 0.0
x3(1) = 1.0
E3(1) = 0.5

dt = 0.02*tp

do i=1,200
	t(i+1)=t(i)+dt
enddo

do i=1,200
	x1(i+1)=x1(i)+dt*fx(p1(i))
	p1(i+1)=p1(i)+dt*fp(x1(i))
	E1(i+1)=(x1(i+1)**2+p1(i+1)**2)*0.5
enddo

do i=1,200
	sl=fx(p2(i))
	sr=fx(p2(i)+dt*sl)
	x2(i+1)=x2(i)+dt*0.5*(sl+sr)

	sl=fp(x2(i))
	sr=fp(x2(i)+dt*sl)
	p2(i+1)=p2(i)+dt*0.5*(sl+sr)

	E2(i+1)=(x2(i+1)**2+p2(i+1)**2)*0.5
enddo

do i=1,200
	s1=dt*fx(p3(i))
	s2=dt*fx(p3(i)+s1/2)
	s3=dt*fx(p3(i)+s2/2)
	s4=dt*fx(p3(i)+s3)
	x3(i+1)=x3(i)+(s1+2*s2+2*s3+s4)/6
	
	s1=dt*fp(x3(i))
	s2=dt*fp(x3(i)+s1/2)
	s3=dt*fp(x3(i)+s2/2)
	s4=dt*fp(x3(i)+s3)
	p3(i+1)=p3(i)+(s1+2*s2+2*s3+s4)/6

	E3(i+1)=(x3(i+1)**2+p3(i+1)**2)*0.5
	
enddo
	
do i=1,201
	write(10,*) x1(i),p1(i),x2(i),p2(i),x3(i),p3(i)
	write(20,*) (t(i)/tp),2*E1(i),2*E2(i),2*E3(i)
	write(30,*) (t(i)/tp),x1(i),x2(i),x3(i)
enddo

contains

real function fx(p)
real p
fx=p
return
end function fx

real function fp(x)
real x
fp=-x
return
end function fp

end
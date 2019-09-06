program kinetics
implicit none
real kab,kba,kbc,kcb,kcd
real dt
real,dimension(600001):: t,a,b,c,d
real s1,s2,s3,s4
integer i
kab=1
kba=3
kbc=4.2
kcb=7.3
kcd=0.4
dt=0.0001

open(10, file='a.dat')
open(20, file='b.dat')
open(30, file='c.dat')
open(40, file='d.dat')

t(1)=0.0
a(1)=5.0
b(1)=0.0
c(1)=0.0
d(1)=0.0

do i=1,600000
t(i+1)=t(i)+dt
enddo

do i=1,600000
s1=dt*fa(a(i),b(i))
s2=dt*fa(a(i)+s1/2,b(i)+s1/2)
s3=dt*fa(a(i)+s2/2,b(i)+s2/2)
s4=dt*fa(a(i)+s3,b(i)+s3)
a(i+1)=a(i)+(s1+2*s2+2*s3+s4)/6

s1=dt*fb(a(i),b(i),c(i))
s2=dt*fb(a(i)+s1/2,b(i)+s1/2,c(i)+s1/2)
s3=dt*fb(a(i)+s2/2,b(i)+s2/2,c(i)+s2/2)
s4=dt*fb(a(i)+s3,b(i)+s3,c(i)+s3)
b(i+1)=b(i)+(s1+2*s2+2*s3+s4)/6

s1=dt*fc(b(i),c(i))
s2=dt*fc(b(i)+s1/2,c(i)+s1/2)
s3=dt*fc(b(i)+s2/2,c(i)+s2/2)
s4=dt*fc(b(i)+s3,c(i)+s3)
c(i+1)=c(i)+(s1+2*s2+2*s3+s4)/6

s1=dt*fd(c(i))
s2=dt*fd(c(i)+s1/2)
s3=dt*fd(c(i)+s2/2)
s4=dt*fd(c(i)+s3)
d(i+1)=d(i)+(s1+2*s2+2*s3+s4)/6
enddo

do i=1,600001
write(10,*) t(i),a(i)
write(20,*) t(i),b(i)
write(30,*) t(i),c(i)
write(40,*) t(i),d(i)
enddo

contains

real function fa(a1,b1)
real a1,b1
fa=-kab*a1+kba*b1
return
end function fa

real function fb(a1,b1,c1)
real a1,b1,c1
fb=-kbc*b1+kcb*(c1**2)+kab*a1-kba*b1
return
end function fb

real function fc(b1,c1)
real c1,b1
fc=2*kbc*b1-2*kcb*(c1**2)-kcd*c1
return
end function fc

real function fd(c1)
real c1
fd=kcd*c1
return
end function fd

end

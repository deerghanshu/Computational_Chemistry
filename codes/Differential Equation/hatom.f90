program hatom
implicit none
real dr
real,dimension(10000):: r,p3,x3
real,dimension(21):: E
real s1,s2,s3,s4
integer i,j
 character*12 file1,file2

do j=1,21
E(j)=-0.6+(j-1)*0.01
enddo



do j=1,21

	write(file1,'("rad.",I3)')100+j
	write(file2,'("rad.",I3)')200+j
	open(unit=j*2-1, file=file1)
	open(unit=j*2, file=file2)

	r(1)=0.0005
	x3(1)=0.000001
	p3(1)=-1000.0
	dr=0.0005

do i=1,9999
	r(i+1)=r(i)+dr
enddo

do i=1,9999
	
	s1=dr*fp(j,r(i),p3(i),x3(i))
	s2=dr*fp(j,r(i)+dr/2,p3(i)+s1/2,x3(i)+s1/2)
	s3=dr*fp(j,r(i)+dr/2,p3(i)+s2/2,x3(i)+s2/2)
	s4=dr*fp(j,r(i)+dr,p3(i)+s3,x3(i)+s3)
	p3(i+1)=p3(i)+(s1+2*s2+2*s3+s4)/6

	s1=dr*fx(p3(i))
	s2=dr*fx(p3(i)+s1/2)
	s3=dr*fx(p3(i)+s2/2)
	s4=dr*fx(p3(i)+s3)
	x3(i+1)=x3(i)+(s1+2*s2+2*s3+s4)/6
	
	
enddo
	
	do i=1,10000
	write(j*2-1,*) r(i),x3(i)
	write(j*2,*) r(i),(abs(r(i)*x3(i)))**2
	enddo

enddo

contains

real function fx(p)
real p
fx=p
return
end function fx

real function fp(k,r,p,x)
real r,p,x
integer k
fp=-2*p/r-2*x*(E(k)+1/r)
return
end function fp

end
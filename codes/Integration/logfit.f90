program fitting
implicit none
integer i
real a,b,sumx,sumy,sumxy,sumx2
real,dimension(5) :: y,x,x2,y2,y3
sumx= 0
sumy=0
sumxy =0
sumx2 = 0
do i =1,5
	x(i) = i
end do
y(1) = 0.5
y(2) = 2
y(3) = 4.5
y(4) = 8 
y(5) = 12.5
open(unit=20,file='output2')
do i=1,5
	x2(i) = log(x(i))
	y2(i) = log(y(i))
	sumx= sumx + x2(i)
	sumy = sumy + y2(i)
	sumxy = sumxy + (x2(i)*y2(i))
	sumx2 = sumx2 + (x2(i)**2)
end do
b = ((5*sumxy) - (sumx*sumy))/((5*sumx2)-(sumx**2))
a = exp((sumy - (b*sumx))/5)
do i=1,5
	y3(i) = a*(x(i)**b)
	write(20,*)x(i),y(i),y3(i)	
end do
end

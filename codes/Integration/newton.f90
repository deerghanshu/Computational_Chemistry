program newton
implicit none
real num,s,p
integer i,j,n
real, dimension(:),allocatable :: x
real, dimension(:,:),allocatable :: d
print*,'enter n'
read*,n
allocate(x(n),d(n,n))
print*,'enter num'
read*,num
do i=1,n
	read*,x(i),d(i,1)
end do
do j=2,n
	do i = 1,n-j+1
		d(i,j) = (d(i+1,j-1) - d(i,j-1))/(x(i+j-1) - x(i))
	end do
end do
s = d(1,1)
do i=2,n
	p=1
	do j = 1,i-1
		p = p*(num-x(j))
	end do
	s = s+p*d(1,i)
end do
print*,'root of ',num,' is ',	s
end

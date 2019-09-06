program newton_ralphson
implicit none
real a,b
integer i
i=1
print*, 'Enter the starting no.:' 
read*, a

do
if(f(a).eq.0) then
	print*, 'Root is',a,'.'
	exit
else
	b = a - f(a)/fd(a)
	print*, i,b,f(b)
	if(abs(b/a-1).lt.0.000001) then
		print*, 'Root is',b,'.'
		exit
	endif
	i=i+1
	a=b
endif
	
enddo
contains

real function f(x)
	real x
	f =x**2 - 3*x + 2
	return
end function f

real function fd(x)
	real x
	fd =2*x - 3
	return
end function fd

end
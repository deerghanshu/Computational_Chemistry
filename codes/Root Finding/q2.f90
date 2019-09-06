program false_position
implicit none
real a,b,c
integer i
i=1
print*, 'Enter two no. between which root is to be found :' 
read*, a,b

do
if(f(a).eq.0) then
	print*, 'Root is',a,'.'
	exit

elseif(f(b).eq.0) then
	print*, 'Root is',b,'.'
	exit

elseif(f(a)*f(b).lt.0) then
	c=a-((b-a)*f(a))/(f(b)-f(a))
	if(f(a)*f(c).lt.0) then
		b=c
	else 
		a=c
	endif
	print*, i,c,f(c)
	i=i+1
else 
	print*, 'No or two roots between',a,'and',b,'.'
	exit

endif
	
enddo
contains

real function f(x)
	real x
	f =x**2 + x - 2
	return
end function f

end
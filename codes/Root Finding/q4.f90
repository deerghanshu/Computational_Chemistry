program secant
implicit none
real a,b,c
integer i
i=1
print*, 'Enter two starting nos.:' 
read*, a,b

do
if(f(a).eq.0) then
	print*, 'Root is',a,'.'
	exit

elseif(f(b).eq.0) then
	print*, 'Root is',b,'.'
	exit

else
	c = b - f(b)*((b-a))/(f(b)-f(a))
	print*, i,c,f(c)
	if(abs(c/b-1).lt.0.000001) then
		print*, 'Root is',c,'.'
		exit
	endif
	i=i+1
	a=b	
	b=c
endif

enddo
contains

real function f(x)
	real x
	f =x**2 - 4*x - 10
	return
end function f

end
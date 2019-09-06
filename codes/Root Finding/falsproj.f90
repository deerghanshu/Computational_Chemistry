program falseproj
implicit none
real x1,x2,x3
print*,'Enter values of starting points'
read*,x1,x2
if(func(x1)*func(x2).GT.0) then
	print*,'No roots exist'
	Stop
end if
20 x3 = x1-(((x2-x1)*func(x1))/(func(x2)-func(x1)))
if (func(x1)*func(x3).LT.0) then
	x2 = x3
else
	x1 = x3
end if
if(abs(x2-x1).LT.(0.000001)) then
	print*,'The root is ',x3
else
	goto 20
end if
contains
real function func(x)
implicit none
real x
func = (x**2)+ x - 2
return
end function
end


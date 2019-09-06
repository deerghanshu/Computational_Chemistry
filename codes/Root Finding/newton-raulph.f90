program newton
implicit none
real x0,x1,f0,fd0
x0 = 0
3 x1 = x0 - (f(x0)/((2*(x0))- 3))
if(abs((x1-x0)/x1).LT.(0.000001)) then
        print*,'The root is', x0
        stop
else
        x0 = x1
        goto 3
end if

contains
real function f(x)
real x
f = (x*x)-(3*x) +2
return
end function f
end

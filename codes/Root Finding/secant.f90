program secant
implicit none
real x1,x2,f1,f2,x3,e
x1=4
x2=2
f1 = f(x1)
f2 = f(x2)
10 x3 = ((f2*x1)-(f1*x2))/(f2-f1)
if (abs((x3-x2)/x3).LT.(0.000001)) then
        print*,'The root is' ,x3
        stop
else
        x1=x2
        f1=f2
        x2=x3
        f2=f(x3)
        goto 10
end if
contains
real function f(x)
real x
f = (x**2) - (4*x) - 10
return
end function
end

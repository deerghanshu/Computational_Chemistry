PROGRAM bisection
implicit none
real a,b,x1,x2,f1,x0,f0,f2,c,root,s,count
print*,'Enter numbers between which the root is to be found:'
read*,a,b
x1 = a
x2 = b
f1 = func(x1)
f2 = func(x2)
if(f1*f2.gt.0) then
        s=0
        print*,'No roots are there'
	stop
endif
count =1
do
        x0=(x1+x2)/2.0
        f0=func(x0)
        if(f0.eq.0) then
                s=1
                root=x0
                print*,'The root is',root
                return
        endif

        if (f1*f0.lt.0) then
                x2=x0
        else
                x1=x0
        endif
        print*,count,x0

        if (abs((x2-x1)/x2).lt.(0.000001)) then
                s=2
                root = (x1+x2)/2.0
                print*,'The root is', root
                return
        else
                count = count + 1
        endif
end do
contains
real function func(x)
implicit none
real x
func = (x**2)+x-2
return
end function func
end

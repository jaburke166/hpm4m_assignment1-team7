program integral2d
implicit none
integer comm,ierr,rank,nproc
real(kind=8) a,b,real_int,dx,dy,integral,error
real(kind=8) xl,xu,yl,yu,start,finish
integer i,j,n

a=20.d0
b=20.d0
n=100000        !subdivision per side
call cpu_time(start)

real_int=(a+b-a*cos(b**2)-b*cos(a**2))/2
dx=a/n
dy=b/n
integral=0.d0

do i=1,n
  xl=dx*(i-1)
  xu=dx*i
  do j=1,n
    yl=dy*(j-1)
    yu=dy*j
    integral=integral+fun(xl,yl)+fun(xl,yu)+fun(xu,yl)+fun(xu,yu)
  enddo
enddo
integral=integral*dx*dy/4
error=abs(integral-real_int)
call cpu_time(finish)
print '(a,f10.6)','Real value:    ',real_int
print '(a,f10.6)','Approximation: ',integral
print '(a,e10.2)','Error:         ',error
print *,'Time:      ',finish-start
contains
  real function fun(x,y)
  real(kind=8), intent(in):: x,y
  fun=x*sin(x**2)+y*sin(y**2)
  end function fun
end program integral2d

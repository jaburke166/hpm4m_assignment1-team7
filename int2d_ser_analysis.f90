program integral2d
implicit none
real(kind=8):: a,b,real_int,dx,dy,integral,error
real(kind=8):: xl,xu,yl,yu
integer(kind=8):: n,i,j,k,kmax


a=20.d0
b=20.d0
kmax=17

real_int=(a+b-a*cos(b**2)-b*cos(a**2))/2

do k=1,kmax
n=2**k
integral=int(a,b,n)
error=abs(integral-real_int)
print '(i6,a,e10.3)',n,' subdivisions, error:',error
enddo




contains
  real function int(a,b,n)
  implicit none      
  real(kind=8),intent(in):: a,b
  real(kind=8):: dx,dy,xu,xl,yu,yl,integral
  integer(kind=8),intent(in):: n
  integer(kind=8):: i,j
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
  int=integral*dx*dy/4
  end function int





  real function fun(x,y)
  implicit none
  real(kind=8), intent(in):: x,y
  fun=x*sin(x**2)+y*sin(y**2)
  end function fun
end program integral2d

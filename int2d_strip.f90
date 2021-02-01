program integral2d_parallel_strips
implicit none
include 'mpif.h'
integer comm,ierr,rank,nproc
real(kind=8) a,b,real_int,dx,dy,intproc,integral,error
real(kind=8) xl,xu,yl,yu,start,finish
integer i,j
integer(kind=8)::n,rowl,rowu,col4proc

comm = MPI_COMM_WORLD
call MPI_INIT(ierr)
call MPI_COMM_RANK(comm, rank, ierr)
call MPI_COMM_SIZE(comm, nproc, ierr)

a=20.d0
b=20.d0
n=100000               !subdivision per side
col4proc=n/nproc      !columns for each processor
n=col4proc*nproc      !make n a multiple of nproc

if (rank==0) then
  real_int=(a+b-a*cos(b**2)-b*cos(a**2))/2
  call cpu_time(start)
endif

dx=a/n
dy=b/n
intproc=0.d0

rowl=rank*col4proc+1
rowu=(rank+1)*col4proc

do i=rowl,rowu
  xl=dx*(i-1)
  xu=dx*i
  do j=1,n
    yl=dy*(j-1)
    yu=dy*j
    intproc=intproc+fun(xl,yl)+fun(xl,yu)+fun(xu,yl)+fun(xu,yu)
  enddo
enddo

call mpi_reduce(intproc,integral,1,mpi_real8,mpi_sum,0,comm,ierr)

if (rank==0) then
  integral=integral*dx*dy/4
  error=abs(integral-real_int)
  print *,'Processors: ',nproc
  print '(a,i10)','Subdivisions:  ',n
  print '(a,f10.6)','Real value:    ',real_int
  print '(a,f10.6)','Approximation: ',integral
  print '(a,e10.2)','Error:         ',error
  call cpu_time(finish)
  print '(a,f10.3)','Time:          ',finish-start
endif
call mpi_finalize(ierr)

contains
  real function fun(x,y)
  real(kind=8), intent(in):: x,y
  fun=x*sin(x**2)+y*sin(y**2)
  end function fun
end program integral2d_parallel_strips

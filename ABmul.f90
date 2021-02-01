program ABmul
implicit none
include "mpif.h"
integer :: comm, rank, nproc, ierr
integer,parameter :: N=5
integer:: B(1:N,1:N),i,j,Arow(1:N),c(1:N),D(1:N,1:N)
integer,allocatable:: cc(:,:),Dt(:,:)
integer:: row4proc
real*8:: start,finish

comm = MPI_COMM_WORLD
call MPI_INIT(ierr)
call MPI_COMM_RANK(comm, rank, ierr)
call MPI_COMM_SIZE(comm, nproc, ierr)



if (nproc>=N) then



!build B and broadcast
if (rank==0) then
  call cpu_time(start)
  do i=1,N
      B(i,:)=(/((j+i)*(N-j+1),j=1,N)/)
  enddo
endif
call MPI_bcast(B,N**2,MPI_int,0,comm,ierr)

!build row of A
if (rank<N) then
do j=1,N
  Arow(j)=(N-j+rank+2)*(rank+1)
enddo

!multiply row by B
c=matmul(Arow,B)
endif

!gather and transpose D=A*B
call mpi_gather(c,N,mpi_int,D,N,mpi_int,0,comm,ierr)
D=transpose(D)

!display D
if (rank==0) then
  do i=1,N
    print *, (D(i,j),j=1,N)
  enddo
endif



else 



!how many rows of A per processor (the last one will have junk rows)        
row4proc=N/nproc+1
allocate(cc(1:row4proc,1:N))
allocate(Dt(1:N,1:row4proc*nproc))

!build B and broadcast
if (rank==0) then
  do i=1,N
      B(i,:)=(/((j+i)*(N-j+1),j=1,N)/)
  enddo
endif
call MPI_bcast(B,N**2,MPI_int,0,comm,ierr)

!build row of A
do i=1,row4proc
  do j=1,N
    Arow(j)=(N-j+rank*row4proc+i+1)*(rank*row4proc+i)
  enddo
  cc(i,:)=matmul(Arow,B)
enddo
print *,'Proc ',rank,'. Rows: '
do i=1,row4proc
print *,cc(i,:)
enddo
call mpi_gather(transpose(cc),row4proc*N,mpi_int,Dt,row4proc*N,mpi_int,0,comm,ierr)

!display D
if (rank==0) then 
  D=Dt(1:N,1:N)
  D=transpose(D)
  call cpu_time(finish)
  print *,'Time = ',finish-start
  print *,'A * B = '
  do i=1,N
    print *, (D(i,j),j=1,N)
  enddo
endif

deallocate(cc,Dt)

endif

call MPI_FINALIZE(ierr)
end program ABmul

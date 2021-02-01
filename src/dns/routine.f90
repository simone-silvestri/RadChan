
subroutine globalArr(input,output,nx,ny,nz,R_coord,C_coord,myid,scaling)

  implicit none
  include "mpif.h" 
  include "param.txt" 
  integer ii
  integer realsize  

  integer :: nx, ny, nz 

  real(4), dimension(nx,ny,nz)                          :: output
  real(8), dimension(0:nx-1,0:ny/p_row+1,0:nz/p_col+1)  :: input
  Real(4), dimension(:,:,:), allocatable                :: Temp0
  Real(8)                                               :: scaling,stime2
  Integer               :: ierr,  num_procs, myid, newtype, resizedsd, resizedrv

  Integer, Dimension(2) :: R_coord, C_coord
  Integer, Dimension(3) :: sizes, subsizes, starts     
  integer,dimension(:),allocatable :: displacement, recvcnt
  integer,dimension(0:p_row*p_col-1) :: sending,receivers
  integer(kind=mpi_address_kind) :: lb, extent 

  num_procs = p_row*p_col

!--------------=======--------------
!----------Obtain displacement-----
!--------------=======--------------  

  stime2 = mpi_wtime()
  allocate(Temp0(1:nx,R_coord(1):R_coord(2),C_coord(1):C_coord(2)))!
  allocate(displacement(num_procs),recvcnt(num_procs)) 
 
  Do k=C_coord(1),C_coord(2)
  Do j=R_coord(1),R_coord(2)
  Do i=1,nx
       Temp0(i,j,k)=real(input(i-1,j-R_coord(1)+1,k-C_coord(1)+1)*scaling,4) 
  End Do; End Do; End Do

  Do i=1,num_procs
     ii = i-1
     displacement(i)= (ii/p_col)*nx*(ny/p_row) + mod(ii,p_col)*nx*ny*(nz/p_col)
  Enddo  

!--------------=======--------------
! ---  Create the same block type ---
!--------------=======--------------  

  sizes(1) = nx
  sizes(2) = ny
  sizes(3) = nz

  subsizes(1) = nx
  subsizes(2) = R_coord(2)-R_coord(1)+1
  subsizes(3) = C_coord(2)-C_coord(1)+1

  starts(1) = 0  ! 0-based index
  starts(2) = 0
  starts(3) = 0
  recvcnt(:)= 1
  call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,         &
        MPI_ORDER_FORTRAN, MPI_REAL4, newtype, ierr)

  call MPI_Type_size(MPI_REAL4, realsize, ierr)
  extent = 1*realsize
  
  lb = 0
  call MPI_Type_create_resized(newtype, lb, extent, resizedrv, ierr)
  call MPI_Type_commit(resizedrv, ierr)

  Call MPI_AllGatherV(Temp0(1,R_coord(1),C_coord(1)), subsizes(1)*subsizes(2)*subsizes(3), MPI_REAL4,  &
                 output, recvcnt,displacement, resizedrv, MPI_COMM_WORLD, ierr)
  call MPI_TYPE_FREE(resizedrv,ierr)

  deallocate(Temp0)
  
  Call MPI_Barrier(MPI_COMM_WORLD, ierr) 

end subroutine globalArr



subroutine gatherArr(input,output,nx,ny,nz,R_coord,C_coord,myid,gpu_flag,scaling)

  implicit none
  include "mpif.h" 
  include "param.txt" 
  integer ii
  integer realsize  

  integer :: nx, ny, nz 

  real(4), dimension(nx,ny,nz)                          :: output
  real(8), dimension(0:nx-1,0:ny/p_row+1,0:nz/p_col+1)  :: input
  Real(4), dimension(:,:,:), allocatable                :: Temp0
  Real(8)                                               :: scaling,stime
  Integer               :: ierr,  num_procs, myid, newtype, resizedsd, resizedrv, gpu_flag, rcv

  Integer, Dimension(2) :: R_coord, C_coord
  Integer, Dimension(3) :: sizes, subsizes, starts     
  integer,dimension(:),allocatable :: displacement, recvcnt
  integer,dimension(0:p_row*p_col-1) :: sending,receivers
  integer(kind=mpi_address_kind) :: lb, extent 

  num_procs = p_row*p_col

!--------------=======--------------
!----------Obtain displacement-----
!--------------=======--------------  

  stime = mpi_wtime()
  allocate(Temp0(1:nx,R_coord(1):R_coord(2),C_coord(1):C_coord(2)))!
  allocate(displacement(num_procs),recvcnt(num_procs)) 
 
  Do k=C_coord(1),C_coord(2)
  Do j=R_coord(1),R_coord(2)
  Do i=1,nx
       Temp0(i,j,k)=real(input(i-1,j-R_coord(1)+1,k-C_coord(1)+1)*scaling,4) 
  End Do; End Do; End Do

  Do i=1,num_procs
     ii = i-1
     displacement(i)= (ii/p_col)*nx*(ny/p_row) + mod(ii,p_col)*nx*ny*(nz/p_col)
  Enddo  

!--------------=======--------------
! ---  Create the same block type ---
!--------------=======--------------  

  sizes(1) = nx
  sizes(2) = ny
  sizes(3) = nz

  subsizes(1) = nx
  subsizes(2) = R_coord(2)-R_coord(1)+1
  subsizes(3) = C_coord(2)-C_coord(1)+1

  starts(1) = 0  ! 0-based index
  starts(2) = 0
  starts(3) = 0
  recvcnt(:)= 1
  call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,         &
        MPI_ORDER_FORTRAN, MPI_REAL4, newtype, ierr)

  call MPI_Type_size(MPI_REAL4, realsize, ierr)
  extent = 1*realsize
  
  lb = 0
  call MPI_Type_create_resized(newtype, lb, extent, resizedrv, ierr)
  call MPI_Type_commit(resizedrv, ierr)

  sending = 0
  if(gpu_flag==1) sending(myid) = 1
  call mpi_allreduce(sending,receivers,num_procs,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
 
  do rcv = 0,num_procs-1
   if(receivers(rcv)==1) then
    Call MPI_GatherV(Temp0(1,R_coord(1),C_coord(1)), subsizes(1)*subsizes(2)*subsizes(3), MPI_REAL4,  &
                 output, recvcnt,displacement, resizedrv, rcv, MPI_COMM_WORLD, ierr)
   endif
  enddo
  call MPI_TYPE_FREE(resizedrv,ierr)

  deallocate(Temp0)
  
  Call MPI_Barrier(MPI_COMM_WORLD, ierr) 

  if(myid==0) write(*,*) 'Communication time = ',mpi_wtime()-stime

end subroutine gatherArr


subroutine bcastArr(input,output,nx,ny,nz,R_coord,C_coord,myid,communicator,gpu_flag,scaling)

  implicit none
  include "mpif.h" 
  include "param.txt" 
  integer ii
  integer realsize  

  integer :: nx, ny, nz 

  real(4), dimension(0:nx-1,0:ny+1,0:nz+1)              :: output
  real(8), dimension(0:nx-1,0:ny/p_row+1,0:nz/p_col+1)  :: input
  Real(4), dimension(:,:,:), allocatable                :: Temp0,buffer
  Real(8)                                               :: scaling,stime
  Integer               :: ierr,  num_procs, myid, newtype, resizedsd, resizedrv, communicator, gpu_flag, bdim 

  Integer, Dimension(2) :: R_coord, C_coord
  Integer, Dimension(3) :: sizes, subsizes, starts     
  integer,dimension(:),allocatable :: displacement, recvcnt
  integer(8) :: totsize
  integer(kind=mpi_address_kind) :: lb, extent 

  num_procs = p_row*p_col

!--------------=======--------------
!----------Obtain displacement-----
!--------------=======--------------  

  stime = mpi_wtime()
  allocate(Temp0(1:nx,R_coord(1):R_coord(2),C_coord(1):C_coord(2)))!
  allocate(displacement(num_procs),recvcnt(num_procs)) 
 
  Do k=C_coord(1),C_coord(2)
  Do j=R_coord(1),R_coord(2)
  Do i=1,nx
       Temp0(i,j,k)=real(input(i-1,j-R_coord(1)+1,k-C_coord(1)+1)*scaling,4) 
  End Do; End Do; End Do

  Do i=1,num_procs
     ii = i-1
     displacement(i)= (ii/p_col)*nx*(ny/p_row) + mod(ii,p_col)*nx*ny*(nz/p_col)
  Enddo  

!--------------=======--------------
! ---  Create the same block type ---
!--------------=======--------------  


  sizes(1) = nx
  sizes(2) = ny
  sizes(3) = nz

  subsizes(1) = nx
  subsizes(2) = R_coord(2)-R_coord(1)+1
  subsizes(3) = C_coord(2)-C_coord(1)+1

  starts(1) = 0  ! 0-based index
  starts(2) = 0
  starts(3) = 0
  recvcnt(:)= 1
  call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,         &
        MPI_ORDER_FORTRAN, MPI_REAL4, newtype, ierr)

  call MPI_Type_size(MPI_REAL4, realsize, ierr)
  extent = 1*realsize
  
  lb = 0

  call MPI_Type_create_resized(newtype, lb, extent, resizedrv, ierr)
  call MPI_Type_commit(resizedrv, ierr)


  Call MPI_GatherV(Temp0(1,R_coord(1),C_coord(1)), subsizes(1)*subsizes(2)*subsizes(3), MPI_REAL4,  &
                 output(0:nx-1,1:ny,1:nz), recvcnt,displacement, resizedrv, 0, MPI_COMM_WORLD, ierr)

  call mpi_barrier(communicator,ierr)

  do i=1,p_row
   if(mod(nx,i)==0) then
     bdim = nx/i
     if(bdim*ny*nz<=2e8) exit 
   endif
  enddo
  if(myid==0) write(*,*) 'dimension of buff ',bdim

  if(gpu_flag==1) then
   allocate(buffer(bdim,ny,nz))
   do i=1,nx-bdim+1,bdim
     buffer = output(i-1:i+bdim-2,1:ny,1:nz)
     call MPI_BCAST(buffer,bdim*ny*nz,MPI_REAL4,0,communicator,ierr)

     output(i-1:i+bdim-2,1:ny,1:nz) = buffer
     call mpi_barrier(communicator,ierr)
   enddo
   deallocate(buffer)
  endif

  call MPI_TYPE_FREE(resizedrv,ierr)
  
  deallocate(Temp0)
  
  Call MPI_Barrier(MPI_COMM_WORLD, ierr) 
  if(myid==0) write(*,*) 'Communication time = ',mpi_wtime()-stime
end subroutine bcastArr

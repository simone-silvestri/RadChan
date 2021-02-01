
      subroutine read_kP(kP,Height,Tnb,nrank)
      implicit none
      include "param.txt"
      include "common.txt"
      real kP(nTemp),Tnb(nTemp),dummy(nTemp),dummy2(nTemp),Height
      integer nrank
      open(unit=1,file="tables/planck-mean.txt")
      do i=1,nTemp       
        read(1,*) Tnb(i),kP(i),dummy(i),dummy2(i)
      enddo
      close(1)
      kP = kP*Height
      end

      subroutine fix_kappa(kP,Tnb,T1,Tin,Knew)
      implicit none
      include "param.txt"
      include "common.txt"
      real kP(nTemp),Tnb(nTemp)
      real T1(0:i1,0:j1,0:k1)
      real ttemp(0:i1,0:j1,0:k1)
      real Knew (0:i1,jmax,kmax)
      real Tin
      ttemp = T1 * Tin
      do i=0,i1
       do j=1,jmax
        do k=1,kmax
         call linear_int(kP,Tnb,ttemp(i,j,k),Knew(i,j,k),nTemp)
        enddo
       enddo
      enddo
      end

      subroutine linear_int(y,x,x0,y0,n)
      implicit none
      include "param.txt"
      real y(n),x(n),x0,y0
      integer t,n
      t = int(x0 - x(1)) / int(x(2) - x(1)) + 1
      y0 = (y(t+1) - y(t)) / (x(t+1) - x(t)) * (x0 - x(t)) + y(t)
      end

      subroutine fix_temp_old(Th, TMC, Tin, startc, endc, ranks, comms, gpu_flag)
      implicit none
      include "mpif.h"
      include "param.txt"
      real   ::  Th(0:i1,0:j1,0:k1),stime
      real*4, dimension(0:i1,0:jmax_tot+1,0:kmax_tot+1) :: TMC
      real*4, dimension(0:i1,0:j1,0:k1) :: Tbuf
      real Tin
      integer istat(mpi_status_size), startc(3), endc(3),ierr, ranks(5), comms(5),gpu_flag,rk
      integer bffsnd(4),bffrcv(4)

      stime = mpi_wtime()
      Tbuf = real(Th,4)*Tin
      if(ranks(1).eq.0) then
       do i=0,i1
        do j=startc(2),endc(2)
         do k=startc(3),endc(3)
          TMC(i,j,k) = Tbuf(i,j-startc(2)+1,k-startc(3)+1)
         enddo
        enddo
       enddo
      endif

      bffrcv(1) = startc(2)
      bffrcv(2) = endc(2)
      bffrcv(3) = startc(3)
      bffrcv(4) = endc(3)


      do rk=1,p_row*p_col-1      
         if(ranks(1).eq.rk) then
           call mpi_send( Tbuf, (imax+2)*(jmax+2)*(kmax+2),MPI_REAL4, 0, 
     +                   0, comms(1), ierr ) 
           call mpi_send( bffrcv, 4 , MPI_INTEGER, 0, 
     +                     1, comms(1), ierr )
         endif
         if(ranks(1).eq.0) then 
           call mpi_recv( Tbuf, (imax+2)*(jmax+2)*(kmax+2),MPI_REAL4, rk, 
     +                     0, comms(1), istat, ierr ) 
           call mpi_recv( bffrcv, 4 , MPI_INTEGER, rk, 
     +                     1, comms(1), istat,  ierr )
           do i=0,i1
            do j=bffrcv(1),bffrcv(2)
             do k=bffrcv(3),bffrcv(4)
              TMC(i,j,k) = Tbuf(i,j-bffrcv(1)+1,k-bffrcv(3)+1)
             enddo
            enddo
           enddo
         endif
      enddo

      call mpi_bcast(TMC,(imax+2)*(jmax_tot+2)*(kmax_tot+2), MPI_REAL4, 0, comms(5), ierr)

      end



      subroutine fix_temp(Th, TMC, Tin, startc, endc, ranks, comms, gpu_flag)
      implicit none
      include "mpif.h"
      include "param.txt"
      real   ::  Th(0:i1,0:j1,0:k1),stime
      real*4, dimension(0:i1,0:jmax_tot+1,0:kmax_tot+1) :: TMC
      real*4, dimension(0:i1,0:j1,0:k1) :: Tbuf
      real Tin
      integer istat(mpi_status_size), startc(3),endc(3), ierr, ranks(5), comms(5),gpu_flag,rk
      integer bffrcv(4),rcv,t1,t2
      character(len=50) :: filename

      stime = mpi_wtime()
      Tbuf = real(Th*Tin,4) 
      if(gpu_flag.eq.1) then
       do i=0,i1
        do j=startc(2),endc(2)
         do k=startc(3),endc(3)
          TMC(i,j,k) = Tbuf(i,j-startc(2)+1,k-startc(3)+1)
         enddo
        enddo
       enddo
      endif

      bffrcv(1) = startc(2)
      bffrcv(2) = endc(2)
      bffrcv(3) = startc(3)
      bffrcv(4) = endc(3)

      do t1=0,p_col/p_row-1
       do t2=0,p_row-1
        rcv=t1*p_row + t2*(p_col+1)
        do rk=0,p_row*p_col-1      
         if(ranks(1).eq.rk.and.ranks(1).ne.rcv) then
           call mpi_send( Tbuf, (imax+2)*(jmax+2)*(kmax+2),MPI_REAL4, rcv, 
     +                   0, comms(1), ierr ) 
           call mpi_send( bffrcv, 4 , MPI_INTEGER, rcv, 
     +                     1, comms(1), ierr )
         endif
         if(ranks(1).eq.rcv.and.ranks(1).ne.rk) then 
           call mpi_recv( Tbuf, (imax+2)*(jmax+2)*(kmax+2),MPI_REAL4, rk, 
     +                     0, comms(1), istat, ierr ) 
           call mpi_recv( bffrcv, 4 , MPI_INTEGER, rk, 
     +                     1, comms(1), istat,  ierr )
           do i=0,i1
            do j=bffrcv(1),bffrcv(2)
             do k=bffrcv(3),bffrcv(4)
              TMC(i,j,k) = Tbuf(i,j-bffrcv(1)+1,k-bffrcv(3)+1)
             enddo
            enddo
           enddo
 	   Tbuf = real(Th*Tin,4)
           bffrcv(1) = startc(2)
           bffrcv(2) = endc(2)
           bffrcv(3) = startc(3)
           bffrcv(4) = endc(3)
         endif
        enddo
       enddo
      enddo

!      if(gpu_flag.eq.1) then
!       write(filename,'(A4,I1.1,A4)') 'temp',ranks(1),'.cpu'
!       open(1,file=filename)
!       do i=0,i1
!        do j=1,jmax_tot
!         do k=1,kmax
!          write(1,*) i,j,k,TMC(i,j,k),startc(3)
!         enddo
!        enddo
!       enddo
!       close(1)
!      endif

      if(ranks(1).eq.0) write(*,*) 'received temperature to 0 in time: ',mpi_wtime()-stime
      end


      subroutine fix_rad(Qnew, Varnew, qr, vr, Tin, Height, MPI_COMM_K, startc, endc, gpu_flag)
      implicit none
      include "mpif.h"
      include "param.txt"
      real,   dimension(0:i1,jmax, kmax)        :: Qnew, Varnew
      real*4, dimension(0:i1,0:jmax_tot+1,0:k1) :: qr,vr,Qbuf
      real Tin,tavg(0:i1),Height
      integer istat(mpi_status_size), startc(3), endc(3),ierr
      integer MPI_COMM_K, gpu_flag 
      character(len=50) :: filename,filename2
  
!      if(gpu_flag.eq.1) then
!       write(filename,'(A3,I5.5)') 'rad',startc(3)
!       open(1,file=filename)
!       write(filename2,'(A8,I5.5)') 'comprad',startc(3)
!       open(2,file=filename2)
!       tavg = 0
!       do i=0,i1
!        do j=1,jmax_tot
!         do k=1,kmax
!          tavg(i)=tavg(i)+real(qr(i,j,k)/jmax_tot/kmax,8)
!          write(2,*) i,j,k,qr(i,j,k)
!         enddo
!        enddo
!        write(1,*) i, tavg(i)
!       enddo
!       close(1)
!       close(2)
!      endif
         
 
      Qbuf = 0
      Qnew = 0
      call mpi_allreduce(qr, Qbuf, (imax+2)*(jmax_tot+2)*(kmax+2), MPI_REAL4,
     +                    MPI_SUM, MPI_COMM_K, ierr)
      do i=0,i1
       do j=startc(2),endc(2)
        do k=1,kmax
         Qnew(i,j-startc(2)+1,k) = real(Qbuf(i,j,k),8)/stefan/(Tin**4.)
        enddo
       enddo
      enddo
      do j=startc(2),endc(2)
       do k=1,kmax
        Qnew(0 ,j-startc(2)+1,k) = real(Qbuf(0 ,j,k),8)/stefan/(Tin**4.)/Height
        Qnew(i1,j-startc(2)+1,k) = real(Qbuf(i1,j,k),8)/stefan/(Tin**4.)/Height
       enddo
      enddo
      call mpi_barrier(MPI_COMM_K,ierr)
      Qbuf   = 0
      Varnew = 0
      call mpi_allreduce(vr, Qbuf, (imax+2)*(jmax_tot+2)*(kmax+2),MPI_REAL4,
     +                    MPI_SUM, MPI_COMM_K, ierr)

      do i=0,i1
       do j=startc(2),endc(2)
        do k=1,kmax
         Varnew(i,j-startc(2)+1,k) = real(Qbuf(i,j,k),8)/stefan/(Tin**4.)
        enddo
       enddo
      enddo
      call mpi_barrier(MPI_COMM_K,ierr)

      end

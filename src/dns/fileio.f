

c********************************************************************
c     read table from file
c********************************************************************
      subroutine readTable(rank)
      implicit none
      include 'param.txt'
      include 'mpif.h'
      if (EOSmode.eq.0) call readTableIG(rank)
      if (EOSmode.eq.1) call readTableRG(rank)
      end

c********************************************************************
c     read table from file
c********************************************************************
      subroutine readTableIG(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr

      if (rank.eq.0) then
        do i=1,nTab
          tempTab(i)   = ((1.0*i-1)/(nTab-1.0)-0.1)*3.0 + 1.0
          rhoTab(i)    = 1.0/tempTab(i)
          muTab(i)     = 1.0
          lamTab(i)    = 1.0
          cpTab(i)     = 1.0
          enthTab(i)   = tempTab(i) - 1.0
          lamocpTab(i) = lamTab(i)/cpTab(i)
        enddo
      endif

      call MPI_BCAST(tempTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rhoTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(muTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cpTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(enthTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamocpTab, nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)

      end

c********************************************************************
c     read table from file
c********************************************************************
      subroutine readTableRG(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr

      if (rank.eq.0) then
	 open(27,file='properties/h2o_500_1800K.dat')
         do i=1,nTab
            read (27,*) tempTab(i),rhoTab(i),muTab(i),lamTab(i),cpTab(i),enthTab(i)
          lamocpTab(i) = lamTab(i)/cpTab(i)
        enddo
         close(27)
      endif



      call MPI_BCAST(tempTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rhoTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(muTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cpTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(enthTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamocpTab, nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)


      call spline(enthTab, rhoTab,    nTab, rho2Tab)
      call spline(enthTab, muTab,     nTab, mu2Tab)
      call spline(enthTab, lamTab,    nTab, lam2Tab)
      call spline(enthTab, cpTab,     nTab, cp2Tab)
      call spline(enthTab, lamocpTab, nTab, lamocp2Tab)
      call spline(enthTab, tempTab,   nTab, temp2Tab)

      end

      subroutine outputProf2(U1,V1,W1,nrank,istap)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real*8 U1(0:i1,jmax,kmax), V1(0:i1,jmax,kmax),W1(0:i1,jmax,kmax) 
      integer istap,ierr,nrank
      real*8 um(0:imax+1)
      real*8 vm(0:imax+1)
      real*8 wm(0:imax+1)
      real*8 umm(0:imax+1)
      real*8 vmm(0:imax+1)
      real*8 uwm(0:imax+1)
      real*8 wmm(0:imax+1)
      real*8 ur(0:imax+1)
      real*8 vr(0:imax+1)
      real*8 uw(0:imax+1)
      real*8 wr(0:imax+1)
      real*8 uv(0:imax+1)
      umm = 0
      vmm = 0
      wmm = 0
      do k=1,kmax
         do j=1,jmax
            do i=0,imax+1
               umm(i)=umm(i)+u1(i,j,k)
               vmm(i)=vmm(i)+v1(i,j,k)
               wmm(i)=wmm(i)+w1(i,j,k)
            enddo
         enddo
      enddo
      call mpi_allreduce(umm,um,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(vmm,vm,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(wmm,wm,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      um = um /(jmax*kmax*px)
      vm = vm /(jmax*kmax*px)
      wm = wm /(jmax*kmax*px)
      umm = 0
      vmm = 0
      wmm = 0
      uwm = 0
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               umm(i)=umm(i)+(u1(i,j,k)-um(i))**2
               vmm(i)=vmm(i)+(v1(i,j,k)-vm(i))**2
               wmm(i)=wmm(i)+(w1(i,j,k)-wm(i))**2
               uwm(i)=uwm(i)+(w1(i,j,k)+w1(i,j,k-1)-2*wm(i))*
     &              (u1(i,j,k)+u1(i-1,j,k)-2*um(i))/4
            enddo
         enddo
      enddo
      call mpi_allreduce(umm,ur,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(vmm,vr,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(wmm,wr,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(uwm,uw,imax+2,mpi_real8,
     &     mpi_sum,mpi_comm_world,ierr)
      ur = ur /(jmax*kmax*px)
      vr = vr /(jmax*kmax*px)
      wr = wr /(jmax*kmax*px)
      uw = uw /(jmax*kmax*px)

      if(nrank.eq.0)then
      open(7,file = 'prof_ins2')
      do i=1,imax
         write(7,'(20E18.6)') Rp(i),um(i),vm(i),wm(i),sqrt(ur(i))
     &        ,sqrt(vr(i)),sqrt(wr(i)),uw(i),
     &        -(wm(i+1)-wm(i-1))/(ru(i)-ru(i-1))/(2*Re)+uw(i)
      enddo
      close(7)
      endif
      end

      subroutine outputProf(U1,V1,W1,nrank,istap)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real*8 U1(0:i1,0:j1,0:k1), V1(0:i1,0:j1,0:k1),W1(0:i1,0:j1,0:k1) 
      integer istap,ierr,nrank
      real*8 um(0:imax+1)
      real*8 vm(0:imax+1)
      real*8 wm(0:imax+1)
      real*8 umm(0:imax+1)
      real*8 vmm(0:imax+1)
      real*8 uwm(0:imax+1)
      real*8 wmm(0:imax+1)
      real*8 ur(0:imax+1)
      real*8 vr(0:imax+1)
      real*8 uw(0:imax+1)
      real*8 wr(0:imax+1)
      real*8 uv(0:imax+1)
      umm = 0
      vmm = 0
      wmm = 0
      do k=1,kmax
         do j=1,jmax
            do i=0,imax+1
               umm(i)=umm(i)+u1(i,j,k)
               vmm(i)=vmm(i)+v1(i,j,k)
               wmm(i)=wmm(i)+w1(i,j,k)
            enddo
         enddo
      enddo
      call mpi_allreduce(umm,um,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(vmm,vm,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(wmm,wm,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      um = um /(jmax*kmax*px)
      vm = vm /(jmax*kmax*px)
      wm = wm /(jmax*kmax*px)
      umm = 0
      vmm = 0
      wmm = 0
      uwm = 0
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               umm(i)=umm(i)+(u1(i,j,k)-um(i))**2
               vmm(i)=vmm(i)+(v1(i,j,k)-vm(i))**2
               wmm(i)=wmm(i)+(w1(i,j,k)-wm(i))**2
               uwm(i)=uwm(i)+(w1(i,j,k)+w1(i,j,k-1)-2*wm(i))*
     &              (u1(i,j,k)+u1(i-1,j,k)-2*um(i))/4
            enddo
         enddo
      enddo
      call mpi_allreduce(umm,ur,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(vmm,vr,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(wmm,wr,imax+2,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(uwm,uw,imax+2,mpi_real8,
     &     mpi_sum,mpi_comm_world,ierr)
      ur = ur /(jmax*kmax*px)
      vr = vr /(jmax*kmax*px)
      wr = wr /(jmax*kmax*px)
      uw = uw /(jmax*kmax*px)

      open(7,file = 'prof_ins')
      if(nrank.eq.0)then
      do i=1,imax
         write(7,'(20E18.6)') Rp(i),um(i),vm(i),wm(i),sqrt(ur(i))
     &        ,sqrt(vr(i)),sqrt(wr(i)),uw(i),
     &        -(wm(i+1)-wm(i-1))/(ru(i)-ru(i-1))/(2*Re)+uw(i)
      enddo
      close(7)
      endif
      end


c

      subroutine output2d(U1,V1,W1,C1,R1,P1,Q2,M1,L1,nrank,istap)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*5 cha
      integer istap,jstart,nrank
      real*8 ,dimension(0:i1,0:j1,0:k1)::U1,V1,W1,P1,C1,R1
      real*8 ,dimension(0:i1,jmax,kmax)::Q2,M1,L1
      if (nrank.lt.p_col)then
      write(cha,'(I5.5)')nrank
      open(15,file='plane/tecp.'//cha)
      if (nrank.eq.0) then
        write(15,*) 'VARIABLES ="X", "Y", "U", "V", "W", "C", "P", "R", "MU", "Lcp"'
        write(15,*) 'ZONE I=  ', imax,' J=  ',kmax_tot,' F=POINT '
      endif
      do k=1,kmax
      do i=0,i1
          j =jmax/2 
              write(15,'(11E20.10)')(k+nrank*(kmax))*dz,rp(i),(u1(i,j,k)+u1(i-1,j,k))/2.,
     1          (v1(i,j,k)+v1(i,j-1,k))/2.,(w1(i,j,k)+w1(i,j,k-1))/2.,C1(i,j,k),R1(i,j,k),P1(i,j,k),
     1           Q2(i,j,k),M1(i,j,k),L1(i,j,k)
         enddo
      enddo
      close(15)
      endif
      end
      subroutine pack_cbuffer_name(chars,buffer,buffer_count,reset)
      implicit none
      integer, dimension(*), intent(inout) :: buffer
      character(len=*), intent(in) :: chars
      integer, intent(inout) :: buffer_count
      integer :: i, n, reset
        if (reset.eq.0) buffer_count = 0
        n = len_trim(chars)
        do i = 1,n
          buffer(buffer_count+i) = ICHAR(chars(i:i))
        end do
        buffer(buffer_count+n+1) = 0
        buffer_count = buffer_count+n+1
      end
     

       subroutine output3d(U1,V1,W1,P1,nrank,istap)
       implicit none

       include 'param.txt'
       include 'common.txt'
       include 'mpif.h'

      integer*4 buffer(1000), tec_version, tec_byte, tec_prec,tec_nvalues_3D,tec_nvalues_2D, tec_format, tec_color, tec_repeat
      real*4 tec_marker1, tec_marker2
      integer bufferCount,nrank
      real*8 U1(0:i1,0:j1,0:k1),V1(0:i1,0:j1,0:k1),W1(0:i1,0:j1,0:k1),P1(0:i1,0:j1,0:k1)
      integer istap, rank_j, rank_k
      real*8 theta
      real*8 Unode,Vnode,Wnode,Pnode,Ux,Uy
      character*5 cha


      rank_j=int(nrank/p_col)
      rank_k=mod(nrank,p_col)

      tec_byte = 1             ! byte order integer value
      tec_marker1 = 299.0
      tec_marker2 = 357.0
      tec_format = 1           ! 1..structured data
      tec_color = -1
      tec_repeat = 0
      tec_prec = 2             ! 1..single, 2..double precission

      tec_nvalues_3D = 7          ! number of variables to be written



!      if (nrank.lt.p_col)then
      write(cha,'(I5.5)')nrank
      open(45,file='3D/3D_plot.'//cha,form='unformatted',access='stream')

!      if (nrank.eq.0) then
        write(45) '#!TDV75 '
        write(45) tec_byte
        call pack_cbuffer_name('DNS pipe',buffer,bufferCount, 0)
        write(45) buffer(1:bufferCount)
        write(45) tec_nvalues_3D
          call pack_cbuffer_name('x',buffer,bufferCount, 0)
          call pack_cbuffer_name('y',buffer,bufferCount, 1)
          call pack_cbuffer_name('z',buffer,bufferCount, 1)
          call pack_cbuffer_name('Vel x',buffer,bufferCount, 1)
          call pack_cbuffer_name('Vel y',buffer,bufferCount, 1)
          call pack_cbuffer_name('Vel z',buffer,bufferCount, 1)
          call pack_cbuffer_name('P',buffer,bufferCount, 1)
          write(45) buffer(1:bufferCount)

        ! begin tecplot binary header section
        write(45) tec_marker1
          call pack_cbuffer_name('flowfield',buffer,bufferCount, 0)
          write(45) buffer(1:bufferCount)
          write(45) tec_format
        write(45) tec_color
        write(45) imax+2, jmax+2, kmax+2
        write(45) tec_marker2
!       end tecplot binary header section
!       begin tecplot data
        write(45) tec_marker1
          write(45) tec_repeat
          do i=1,tec_nvalues_3D
            write(45) tec_prec   ! precision of the data
          enddo
 !     endif

        do k=0,k1
         do j=0,j1
         theta = (j*0.7+jmax*rank_j)*dphi
            do i=0,i1

               Unode =0.0
               Vnode =0.
               Wnode =W1(i,j,k)

               Ux = Unode
               Uy = Unode
               Pnode = 0

!               Pnode = 1./8.*(P1(i  ,j  ,k)+P1(i  ,j  ,k+1)
!     &                       +P1(i  ,j+1,k)+P1(i  ,j+1,k+1)
!     &                       +P1(i+1,j  ,k)+P1(i+1,j  ,k+1)
!     &                       +P1(i+1,j+1,k)+P1(i+1,j+1,k+1))

               write(45)(ru(i)+0.2)*cos(theta),(ru(i)+0.2)*sin(theta),(k+1.3*rank_k*kmax)*dz,Ux,Uy,Wnode,Pnode
            enddo
         enddo
      enddo
      close(45)

!      endif
        end

      subroutine load_save2(U1,V1,W1,C1)!,R2,Wx,Wy,Wz,dold,rank,istart,counter_post)
      use decomp_2d
      use decomp_2d_io
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*4 cha
      character*5 cha2
      integer ierr,istart,ini,counter_post
      real*8,  dimension(0:i1,jmax,kmax) :: U1,V1,W1,C1
      real*8 delta,Yplus
      call decomp_2d_read_one(1,U1  ,'Restart/U.dns')
      call decomp_2d_read_one(1,V1  ,'Restart/V.dns')
      call decomp_2d_read_one(1,W1  ,'Restart/W.dns')
      call decomp_2d_read_one(1,C1  ,'Restart/C.dns')
      end

      subroutine load_save1(ini,U1,V1,W1,C1,R2,Wx,Wy,Wz,dold,rank,istart,counter_post)
      use decomp_2d
      use decomp_2d_io
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*4 cha
      character*5 cha2
      integer ierr,istart,ini,counter_post
      real*8,  dimension(0:i1,jmax,kmax) :: Wx,Wy,Wz,dold
      real*8,  dimension(0:i1,0:j1,0:k1) :: U1,V1,W1,C1,R2
      real*8 delta,Yplus
!      write(cha,'(I4.4)')counter_post
      if (ini.eq.0) then
!      if (rank.eq.0)then
      open(11,file='Restart/Init_stop')
      read(11,*)istart,counter_post
      close(11)
!      endif
!      call MPI_BCAST(istart,1,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call decomp_2d_read_one(1,wx  ,'Restart/wx.dns')
      call decomp_2d_read_one(1,wy  ,'Restart/wy.dns')
      call decomp_2d_read_one(1,wz  ,'Restart/wz.dns')
      call decomp_2d_read_one(1,dold,'Restart/dold.dns')
       write(cha2,'(I5.5)')rank
         open(13,file='Restart/OLD_var.'//cha2,form='unformatted')
           read(13) U1,V1,W1,C1,R2
         close(13)
      endif
   
      if (ini.eq.1) then
      if (rank.eq.0)then
      open(11,file='Restart/Init_stop')
      write(11,*)istart,counter_post
      close(11)
      endif
      call decomp_2d_write_one(1,wx  ,'Restart/wx.dns')
      call decomp_2d_write_one(1,wy  ,'Restart/wy.dns')
      call decomp_2d_write_one(1,wz  ,'Restart/wz.dns')
      call decomp_2d_write_one(1,dold,'Restart/dold.dns')
      write(cha2,'(I5.5)')rank
        open(13,file='Restart/OLD_var.'//cha2,form='unformatted')
          Write(13) U1,V1,W1,C1,R2
        close(13)
      endif

      if (ini.eq.-1) then
      delta=1.
      U1 =0.
      V1 =0.
      C1 =0.
      Wx=0.
      Wy=0.
      Wz=0.
      dold=0.
            
!       do k=0,k1
!          do j=0,j1
!             do i=1,imax
!       Yplus = (0.5-Rp(i))*Re
!       if (Yplus .gt. 11)then
!          W1(i,j,k)=(1./0.41)*log(Yplus)+5.3
!       else
!          W1(i,j,k)=Yplus
!       endif
!               u1(i,j,k)=u1(i,j,k)+delta*(2.5*log(Re/2.)+5.)*cos(float(i*j*k+rank*(k+2)))
!!        vnew(i,j,k)=vnew(i,j,k)+0.25*(2.5*log(Re_i/2.)+5.)*sin(float(i*j*k))
!               w1(i,j,k)=W1(i,j,k)+delta*(2.5*log(Re/2.)+5.)*sin(float(i*j*k+rank*(k+2)))
!             enddo
!          enddo
!       enddo
       U1=0. 
       V1=0.
       do k=0,k1
         do j=0,j1
           do i=1,imax
             W1(i,j,k)=(Re)*rp(i)*(1.-0.5*rp(i))
           enddo
         enddo
       enddo
      endif
      end
      subroutine load_save(ini,Unew,Vnew,Wnew,Cnew,R2,Wx,Wy,Wz,dold,LAST,FIRST,rank,istart,counter_post)
      use decomp_2d
      use decomp_2d_io
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*4 cha
      character*5 cha2
      integer ierr,istart,ini,counter_post
      real*8,  dimension(0:i1,jmax,kmax) :: Unew,Vnew,Wnew,Cnew,R2,Wx,Wy,Wz,dold
      real*8,  dimension(0:i1,jmax_tot,6) :: LAST,FIRST
      real*8 delta,Yplus
      if (ini.eq.0) then
       open(21,file='Restart/FIRST_LAST',form='unformatted')
       read(21) FIRST,LAST
       close(21)
      open (11,file='Restart/Init_stop')
      read(11,*)istart,counter_post
      close(11)
      call decomp_2d_read_one(1,Unew  ,'Restart/Unew.dns')
      call decomp_2d_read_one(1,Vnew  ,'Restart/Vnew.dns')
      call decomp_2d_read_one(1,Wnew  ,'Restart/Wnew.dns')
      call decomp_2d_read_one(1,Cnew  ,'Restart/Cnew.dns')
      call decomp_2d_read_one(1,dold  ,'Restart/Dold.dns')
      call decomp_2d_read_one(1,wx    ,'Restart/Wx.dns')
      call decomp_2d_read_one(1,wy    ,'Restart/Wy.dns')
      call decomp_2d_read_one(1,wz    ,'Restart/Wz.dns')
      call decomp_2d_read_one(1,R2    ,'Restart/R2.dns')
      endif
      if (ini.eq.1) then
      if (rank.eq.0)then
        open(21,file='Restart/FIRST_LAST',form='unformatted')
        WRITE(21) FIRST,LAST
        close(21)
      open (11,file='Restart/Init_stop')
      write(11,*)istart,counter_post
      close(11)
      endif
      call decomp_2d_write_one(1,Unew  ,'Restart/Unew.dns')
      call decomp_2d_write_one(1,Vnew  ,'Restart/Vnew.dns')
      call decomp_2d_write_one(1,Wnew  ,'Restart/Wnew.dns')
      call decomp_2d_write_one(1,Cnew  ,'Restart/Cnew.dns')
      call decomp_2d_write_one(1,dold  ,'Restart/Dold.dns')
      call decomp_2d_write_one(1,wx    ,'Restart/Wx.dns')
      call decomp_2d_write_one(1,wy    ,'Restart/Wy.dns')
      call decomp_2d_write_one(1,wz    ,'Restart/Wz.dns')
      call decomp_2d_write_one(1,R2    ,'Restart/R2.dns')
      endif
      end

      subroutine POST_Video(U1,V1,W1,C1,Q1,counter_video,rank)
      use decomp_2d
      use decomp_2d_io
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*5 cha
      integer ierr,counter_video
      real*8,  dimension(0:i1,jmax,kmax) :: U1,V1,W1,C1,Q1
      write(cha,'(I5.5)')counter_video
      call decomp_2d_write_one(1,U1  ,'video/U.'//cha)
      call decomp_2d_write_one(1,V1  ,'video/V.'//cha)
      call decomp_2d_write_one(1,W1  ,'video/W.'//cha)
      call decomp_2d_write_one(1,C1  ,'video/C.'//cha)
      call decomp_2d_write_one(1,Q1  ,'video/Q.'//cha)
      end



      subroutine POST_Data(U1,V1,W1,C1,P1,Q1,counter_post,rank)
      use decomp_2d
      use decomp_2d_io
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*5 cha
      integer ierr,counter_post
      real*8,  dimension(0:i1,jmax,kmax) :: U1,V1,W1,C1,P1,Q1
      write(cha,'(I5.5)')counter_post
      call decomp_2d_write_one(1,U1  ,'DATA/U.'//cha)
      call decomp_2d_write_one(1,V1  ,'DATA/V.'//cha)
      call decomp_2d_write_one(1,W1  ,'DATA/W.'//cha)
      call decomp_2d_write_one(1,C1  ,'DATA/C.'//cha)
      call decomp_2d_write_one(1,P1  ,'DATA/P.'//cha)
      call decomp_2d_write_one(1,Q1  ,'DATA/Q.'//cha)
      end


      subroutine cmpinf(Wnew,T1,Bulk,Cbulk,Stress,istap,rank)
      implicit none
c
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr,istap
      real*8 waver(imax),waver2(imax),taver(imax)
      real*8, dimension(0:i1,jmax,kmax) :: Wnew
      real*8, dimension(0:i1,0:j1,0:k1) :: T1 
      real*8  Cbulk,Bulk,Stress,temp
      Waver = 0.0
      taver = 0.0
      do k=1,kmax
         do j=1,jmax
            do  i=1,imax
               Waver(i) = Waver(i) + Wnew(i,j,k)
               taver(i) = taver(i) + T1(i,j,k)
            enddo
         enddo
      enddo
      call mpi_allreduce(waver,waver2,imax, mpi_real8,mpi_sum,mpi_comm_world,ierr)
      waver = waver2/(jmax*kmax*px)
      call mpi_allreduce(taver,waver2,imax, mpi_real8,mpi_sum,mpi_comm_world,ierr)
      taver = waver2/(jmax*kmax*px)
      Bulk = 0.0
      Cbulk = 0.0
      do i=1,imax
         Bulk  = Bulk  + Waver(i) * dr(i) * 0.5
         Cbulk = Cbulk + taver(i) * dr(i) * 0.5
      enddo
      Stress =  Waver(imax) /(0.5-rp(imax))/Re
      return
      end

      subroutine outputX_h(W1,C1,R1,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*2 cha
      integer ierr,istap,ktabhi,ktablo
      real*8,  dimension(0:i1,0:j1,0:k1) :: W1,C1,R1

      real*8  massflow(kmax),enthflow(kmax),enth_b(kmax),Twall(kmax),T_bulk(kmax),Qp(kmax),
     &      tmp(kmax),massfl,enth,laminter,tempinter,cpinter,Gb(kmax),rhob(kmax),
     &     addedHeat,addedHeatTot,subheat,subheatTot,w_T

      twall    = 0.0
      Qp       = 0.0
      massflow = 0.0
      enthflow = 0.0
      W_T=0.


      do k=1,kmax
         do j=1,jmax
            ktabhi = 0
            ktablo = 0
            call splint(enthTab,tempTab, temp2Tab,nTab,0.5*(c1(i1,j,k)+c1(imax,j,k)),tempinter,ktabhi,ktablo)
            twall(k) = twall(k) + tempinter
            do i=1,imax
               massfl = 0.5*w1(i,j,k)*(r1(i,j,k)+r1(i,j,k+1))*Rp(i)*dphi*dr(i)
               massflow(k) = massflow(k) + massfl

               enthflow(k) = enthflow(k) +massfl*(C1(i,j,k)+C1(i,j,k+1))/2.

            enddo
         enddo
      enddo
      twall=twall/real(jmax)
      write(cha,'(I2.2)')rank
      open(9,file = 'H_b/profX.'//cha)
      do k=1,kmax
      write(9,*) massflow(k),enthflow(k),Twall(k)
      enddo
      close(9)
      return
      end


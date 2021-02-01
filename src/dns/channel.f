      use decomp_2d
      implicit none

      include      'param.txt'
      include      'common.txt'
      include      'mpif.h'
      character*3  pha
      character*5  inflow
      character*(MPI_MAX_PROCESSOR_NAME) node_name
      logical    periodic_bc(3)
      integer    ploc,ierr,istart,iavg,select_init,counter_post,acc,istep_fic
      integer    result_proc,color,node_size,gpu_flag,counter_video
      integer    rank_k    ,MPI_COMM_K                        
      integer    rank_node ,MPI_COMM_NODE
      integer    rank_intra,MPI_COMM_INTRA                        
      integer    rank_gpu  ,MPI_COMM_GPU                        
      integer    ranks(5)  ,comms(5), coordj(2), size_core
      integer*8  hash

      real*8, allocatable, dimension(:,:,:) :: Uh,Vh,Wh,Ph,Ch,Cuh,Cvh,Cwh,dummy2
      real*8, dimension(0:i1,0:j1,0:k1)     :: ekm,ekh,cpp,con,R1,R2,R3,RS,Rh,RSS,dummy
      real*8, dimension(0:i1,0:j1,0:k1)     :: U1,V1,W1,C1,T1
      real*8, dimension(0:i1,jmax,kmax)     :: Unew,Vnew,Wnew,Cnew,dudt,dvdt,dwdt,dcdt,p,Wx,Wy,Wz,dold,US,VS,WS,CS,R22
      real*8, dimension(0:i1,jmax,kmax)     :: Qnew,Gnew,Knew,Varnew
      real*8, dimension(0:i1,jmax,kmax)     :: Cu,Cv,Cw,dudt_t,dvdt_t,dwdt_t
      real*8, dimension(0:i1,jmax_tot,6)    :: LAST,LAST_tot,FIRST,FIRST_tot
      real*4, allocatable, dimension(:,:,:) :: TMC
      real*4, allocatable, dimension(:,:,:) :: qr,vr 
      real*8, dimension(nTemp)          :: kPlanck, Tnb
      real*8, dimension(0:i1,jmax_tot,nInflow_sub*ratioI,3) :: U_bound_tot
      real*8  bulk,cbulk,stress,stime,stime2,stime3,time1,time2,timer,time3,A1,B1,Tin,Height,underrelax

      call mpi_init(ierr)
      call cpu_time(time1)
      periodic_bc(1) = .false.
      periodic_bc(2) = .true.
      periodic_bc(3) = .false.
      if (periodic.eq.1) periodic_bc(3) = .true.
      call MPI_COMM_SIZE(MPI_COMM_WORLD, size_core, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierr)
      if (nrank.eq.0) write(*,*) 'n_proc = ', size_core

      call decomp_2d_init(imax_tot+2,jmax_tot,kmax_tot,p_row,p_col,periodic_bc)
      if (nrank.eq.0) then
      write(6,*) xsize(1),xsize(2),xsize(3),xsize(1)*xsize(2)*xsize(3)
      write(6,*) ysize(1),ysize(2),ysize(3),ysize(1)*ysize(2)*ysize(3)
      write(6,*) zsize(1),zsize(2),zsize(3),zsize(1)*zsize(2)*zsize(3)
      endif

      allocate(TMC(0:i1,0:jmax_tot+1,0:kmax_tot+1))
      allocate(qr(0:i1,0:jmax_tot+1,0:k1)) 
      allocate(vr(0:i1,0:jmax_tot+1,0:k1)) 
 

!!!!!!!!!!!!!!!!!!!IDENTIFYING PROCESSING NODES, DEFINING NEW !COMMUNICATOR!!!!!!!!!!!!!!!!!!!!!!!!!!!

      color = xstart(3)
      call mpi_comm_split(MPI_COMM_WORLD, color, nrank, MPI_COMM_K,ierr)
      call mpi_comm_rank(MPI_COMM_K , rank_k, ierr)

      call mpi_get_processor_name(node_name, result_proc, ierr)
      call name_to_hash(hash, node_name, MPI_MAX_PROCESSOR_NAME)
      call mpi_comm_split(MPI_COMM_WORLD, hash, nrank, MPI_COMM_NODE, ierr)
      call mpi_comm_rank(MPI_COMM_NODE , rank_node, ierr)
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_barrier(MPI_COMM_NODE, ierr)

      color = rank_node
      call mpi_comm_split(MPI_COMM_WORLD, color, color, MPI_COMM_INTRA,ierr)
      call mpi_comm_rank(MPI_COMM_INTRA , rank_intra, ierr)

      gpu_flag = 0
      color = rank_intra/(p_col/p_row)

      if(rank_node.eq.color) gpu_flag = 1

      !! communicate who communicates with the GPU
      call mpi_comm_split(MPI_COMM_WORLD, gpu_flag, nrank, MPI_COMM_GPU,ierr)
      call mpi_comm_rank(MPI_COMM_GPU , rank_intra, ierr)


      ranks(1) = nrank
      ranks(2) = rank_k
      ranks(3) = rank_node
      ranks(4) = rank_intra
      ranks(5) = rank_gpu

      comms(1) = MPI_COMM_WORLD
      comms(2) = MPI_COMM_K
      comms(3) = MPI_COMM_NODE
      comms(4) = MPI_COMM_INTRA
      comms(5) = MPI_COMM_GPU

      coordj(1) = xstart(2)
      coordj(2) = xend(2)
      coordk(1) = xstart(3)
      coordk(2) = xend(3)

      call mpi_barrier(MPI_COMM_WORLD, ierr)


!!!!!!!!!!!!!!!!!!!STARTING CODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      call readTable(nrank)
      Height= 0.5
      call read_kP(kPlanck,Height,Tnb,nrank)
      FIRST=0.
      LAST=0.
      dtmax = 2.e-5*ratioI
      Tin   = T_0 
      call mkgrid(nrank)
      dt = dtmax
      underrelax = 0.0      

!!!!!!!!!!!!!!!!!!!!START O
      call load_save(0,Unew,Vnew,Wnew,Cnew,R22,Wx,Wy,Wz,dold,LAST,FIRST,nrank,istart,counter_post)
      counter_video=100
      call update_halo(Unew,Uh,1)
      call update_halo(Vnew,Vh,1)
      call update_halo(Wnew,Wh,1)
      call update_halo(Cnew,Ch,1)
      call update_halo(R22,dummy2,1)
      call insert(U1,V1,W1,C1,R2,Uh,Vh,Wh,Ch,dummy2,nrank)
      if (xstart(3).eq.(p_col-1)*kmax+1) then
         do j=1,jmax
          do  i=0,i1
          U1(i,j,k1) = LAST(i,j+jmax*(nrank/p_col),1)! U1(i,j,kmax)!
          V1(i,j,k1) = LAST(i,j+jmax*(nrank/p_col),2)! V1(i,j,kmax)!
          W1(i,j,k1) = LAST(i,j+jmax*(nrank/p_col),3)! W1(i,j,kmax)!
          C1(i,j,k1) = LAST(i,j+jmax*(nrank/p_col),4)! C1(i,j,kmax)!
          R2(i,j,k1) = LAST(i,j+jmax*(nrank/p_col),5)! R2(i,j,kmax)!
         enddo
       enddo
       endif
       if (xstart(3).eq.1) then
        do j=1,jmax
          do  i=0,i1
          U1(i,j,0) = FIRST(i,j+jmax*(nrank/p_col),1)! U1(i,j,1) !
          V1(i,j,0) = FIRST(i,j+jmax*(nrank/p_col),2)! V1(i,j,1) !
          W1(i,j,0) = FIRST(i,j+jmax*(nrank/p_col),3)! W1(i,j,1) !
          C1(i,j,0) = FIRST(i,j+jmax*(nrank/p_col),4)! C1(i,j,1) !
          R2(i,j,0) = FIRST(i,j+jmax*(nrank/p_col),5)! R2(i,j,1) !
         enddo
       enddo
       endif

       call state(C1,T1,R1,ekm,ekh)
       istart =istart+1
       R3=R2
       R2=R1
       if(nrank.eq.0) write(*,*) 'Calling the output'
       call output2d(U1,V1,W1,C1,R1,U1,Qnew,Gnew,Knew,nrank,istep)

!!!!!!!!!!!!!!!!!!!!START OF MAIN LOOP!!!!!!!!!!!!!!!!!!!!
!       istart=1
       do istep=istart,nstep
        if(nrank.eq.0) write(*,*) 'Starting the time step'
        underrelax = 1.0 !min(istep/5000.0,1.0)
        time = time + dt
        if(mod(istep,10).eq.0) underrelax = min(1.0,underrelax+0.001)
        isave=mod(istep-1,nInflow_sub)+1
        if (isave.eq.1.and.xstart(3).eq.1) then
          istep_fic = istep*ratioI
          acc=int((istep_fic-1.-int(istep_fic/(nInflowProf*(nInflow_sub*ratioI)))
     &            *nInflowProf*(nInflow_sub*ratioI))/(nInflow_sub*ratioI))+1
          write(inflow,'(I5.5)') acc
          open(72,file='../Inflow/Inflow.'//inflow,
     &                                  form='unformatted')
            read(72) U_bound_tot
          close(72)
          do j=1,jmax
            do i=0,i1
              do k=1,nInflow_sub
                u_bound(i,j,k)=U_bound_tot(i,j+jmax*(nrank/p_col),k*ratioI,1)
                v_bound(i,j,k)=U_bound_tot(i,j+jmax*(nrank/p_col),k*ratioI,2)
                w_bound(i,j,k)=U_bound_tot(i,j+jmax*(nrank/p_col),k*ratioI,3)
              enddo
            enddo
          enddo
        endif

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REMEMBER TO CHANGE
!        if(istep.eq.istart) then
!          call MPI_BCAST(u_bound, (imax+2)*jmax*nInflow_Sub ,mpi_real8, 0 ,MPI_COMM_WORLD,ierr)
!          call MPI_BCAST(v_bound, (imax+2)*jmax*nInflow_Sub ,mpi_real8, 0 ,MPI_COMM_WORLD,ierr)
!          call MPI_BCAST(w_bound, (imax+2)*jmax*nInflow_Sub ,mpi_real8, 0 ,MPI_COMM_WORLD,ierr)
!          do k=1,kmax
!            do j=1,jmax
!              do i=0,i1
!                Unew(i,j,k) = u_bound(i,j,1)
!                Vnew(i,j,k) = v_bound(i,j,1)
!                Wnew(i,j,k) = w_bound(i,j,1)
!                Cnew(i,j,k) = 0.0
!                R22 (i,j,k) = 1.0
!              enddo
!            enddo
!          enddo
!          call update_halo(Unew,Uh,1)
!          call update_halo(Vnew,Vh,1)
!          call update_halo(Wnew,Wh,1)
!          call update_halo(Cnew,Ch,1)
!          call update_halo(R22,dummy2,1)
!          call insert(U1,V1,W1,C1,R2,Uh,Vh,Wh,Ch,dummy2,nrank)
!          do j=1,jmax
!            do  i=0,i1
!              LAST(i,j+jmax*(nrank/p_col),1)=U1(i,j,k1)
!              LAST(i,j+jmax*(nrank/p_col),2)=V1(i,j,k1)
!              LAST(i,j+jmax*(nrank/p_col),3)=W1(i,j,k1)
!              LAST(i,j+jmax*(nrank/p_col),4)=C1(i,j,k1)
!              LAST(i,j+jmax*(nrank/p_col),5)=R2(i,j,k1)
!            enddo
!          enddo
!          do j=1,jmax
!            do  i=0,i1
!              FIRST(i,j+jmax*(nrank/p_col),1)=U1(i,j,0)
!              FIRST(i,j+jmax*(nrank/p_col),2)=V1(i,j,0)
!              FIRST(i,j+jmax*(nrank/p_col),3)=W1(i,j,0)
!              FIRST(i,j+jmax*(nrank/p_col),4)=C1(i,j,0)
!              FIRST(i,j+jmax*(nrank/p_col),5)=R2(i,j,0)
!            enddo
!          enddo
!          call state(C1,T1,R1,ekm,ekh)
!          R2=R1
!          R3=R2
!        endif      
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! REMBER TO CHANGE
 
        A1=1.5
        B1=-0.5
        if(istep.eq.1)then
          A1=1.
          B1=0.
        endif
        RS=2.5*R1-2.*R2+0.5*R3

!!!!!!!!!!!!!!!!!!!KICK-OFF RADIATIVE CALCULATION!!!!!!!!!!!!!!!!!!!!

      if(istep.eq.istart) then
        call bcastArr(T1,TMC,imax+2,jmax_tot,kmax_tot,coordj,coordk,nrank,comms(5),gpu_flag,Tin)
        if(gpu_flag.eq.1) then
          call MC_GPU(TMC, xstart(3), nrank)
        endif
      endif

!!!!!!!!!!!!!!!!!RADIATIVE LOOP!!!!!!!!!!!!!!!!!!!!!!!! I HAVE TO USE C1 !!!!!!
      stime2 = MPI_WTIME()
      if(mod(istep,radStep)==0.or.istep==istart) then
        call bcastArr(T1,TMC,imax+2,jmax_tot,kmax_tot,coordj,coordk,nrank,comms(5),gpu_flag,Tin)
        if(gpu_flag.eq.1) then
          call GET_RESULTS(qr, vr, nrank)
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(gpu_flag.eq.1) then
          call MC_GPU(TMC, xstart(3), nrank)                     !!!!! ALSO HAVE TO CALCULATE RADIATIVE FLUX TO THE BOUNDARY !!!!!
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call fix_rad(Qnew,Varnew,qr,vr,Tin,Height,MPI_COMM_K,xstart,xend,gpu_flag)
        call fix_kappa(kPlanck,Tnb,T1,Tin,Knew)
        do i=0,i1
          do j=1,jmax
            do k=1,kmax
              Qnew(i,j,k) = Qnew(i,j,k)*underrelax
              Gnew(i,j,k) = 4*T1(i,j,k)**4. - Qnew(i,j,k)/Knew(i,j,k)
            enddo
          enddo
        enddo
        call mpi_barrier(MPI_COMM_WORLD,ierr)
      endif
      stime3 = MPI_WTIME()

      call limiter(Cu,Cv,Cw,C1,U1,V1,W1)
      call update_halo(CU,CUh,1)
      call update_halo(CV,CVh,1)
      call update_halo(CW,CWh,1)
      if(xstart(3).eq.1)then
        CUh(:,:,0)=0.
        CVh(:,:,0)=0.
        CWh(:,:,0)=0.
      endif
      call advance_h(dcdt,C1,R1,C1,R1,U1,V1,W1,Cuh,Cvh,Cwh,dold,RS,Qnew,ekh,A1,B1,1,nrank)
      call bound_i_c(dcdt,Qnew,ekh,nrank,istart,xstart(3))
      call update_halo(dcdt,Ch,1)
      call bound_c(Ch,W1,C1,nrank)
      call advance_m(dudt,dvdt,dwdt,U1,V1,W1,R1,U1,V1,W1,R1,Wx,Wy,Wz,ekm,A1,B1,1,nrank)
      call bound_i_u(dudt,dvdt,dwdt,nrank)
      dudt_t=dudt
      dvdt_t=dvdt
      dwdt_t=dwdt
      call update_halo(dudt,Uh,1)
      call update_halo(dvdt,Vh,1)
      call update_halo(dwdt,Wh,1)

      if (periodic.ne.1) call bound_m(Uh,Vh,Wh,W1,RS,R1,R2,nrank)
      call fillps(p,Uh,Vh,Wh,RS,R1,R2,nrank)
      call SOLVEpois(p,Ru,Rp,DPHI,dz,nrank)
      call update_halo(p,Ph,1)
      call correc(Unew,Vnew,Wnew,Uh,Vh,Wh,Ph,RS,nrank)
      call bound_i_u(Unew,Vnew,Wnew,nrank)

      call update_halo(Unew,Uh,1)
      call update_halo(Vnew,Vh,1)
      call update_halo(Wnew,Wh,1)


      if (periodic.ne.1) call bound_v(Uh,Vh,Wh,U1,V1,W1,nrank)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call state(Ch,T1,Rh,dummy,dummy)
      call limiter(Cu,Cv,Cw,Ch,Uh,Vh,Wh)
      call update_halo(CU,CUh,1)
      call update_halo(CV,CVh,1)
      call update_halo(CW,CWh,1)
      if(xstart(3).eq.1)then
        CUh(:,:,0)=0.
        CVh(:,:,0)=0.
        CWh(:,:,0)=0.
      endif
      RSS=RS
      RS=R1+0.5*Rh-0.5*R2
      call advance_h(dcdt,C1,R1,Ch,RSS,Uh,Vh,Wh,Cuh,Cvh,Cwh,dold,RS,Qnew,ekh,0.5,0.5,0,nrank)
      call bound_i_c(dcdt,Qnew,ekh,nrank,istart,xstart(3))
      Cnew=dcdt
      call update_halo(dcdt,Ch,1)
      call bound_c(Ch,W1,C1,nrank)
      dudt=dudt_t
      dvdt=dvdt_t
      dwdt=dwdt_t
      call update_halo(dudt,Uh,1)
      call update_halo(dvdt,Vh,1)
      call update_halo(dwdt,Wh,1)
      if (periodic.ne.1) call bound_m(Uh,Vh,Wh,W1,RS,R1,R2,nrank)
      call fillps(p,Uh,Vh,Wh,RS,R1,R2,nrank)
      call SOLVEpois(p,Ru,Rp,DPHI,dz,nrank)
      call update_halo(p,Ph,1)
      call correc(Unew,Vnew,Wnew,Uh,Vh,Wh,Ph,RS,nrank)
      call bound_i_u(Unew,Vnew,Wnew,nrank)

      call update_halo(Unew,Uh,1)
      call update_halo(Vnew,Vh,1)
      call update_halo(Wnew,Wh,1)
      if (periodic.ne.1) call bound_v(Uh,Vh,Wh,U1,V1,W1,nrank)
      if (mod(istep,divchk) .eq.0)then
        call chkdiv(Uh,Vh,Wh,RS,R1,R2,nrank)
      endif
      R3=R2
      R2=R1
      call state(Ch,T1,R1,ekm,ekh)
      U1=Uh
      V1=Vh
      W1=Wh
      C1=Ch
      if (mod(istep,500).eq.0)then
        counter_post=counter_post+1
        call POST_Data(Unew,Vnew,Wnew,Cnew,p,Qnew,counter_post,nrank)
      endif
!      if (mod(istep,50).eq.0.and.(istep-istart).lt.25000) then
!        counter_video=counter_video+1
!        call POST_Video(Unew,Vnew,Wnew,Cnew,Qnew,counter_video,nrank)
!      endif


      if (mod(istep,1000).eq.0)then
        FIRST=0.
        LAST=0.
        LAST_tot=0.
        FIRST_tot=0.
        if (coordk(2).eq.kmax_tot) then  
          do j=1,jmax
            do  i=0,i1
              LAST(i,j+jmax*(nrank/p_col),1)=U1(i,j,k1)
              LAST(i,j+jmax*(nrank/p_col),2)=V1(i,j,k1)
              LAST(i,j+jmax*(nrank/p_col),3)=W1(i,j,k1)
              LAST(i,j+jmax*(nrank/p_col),4)=C1(i,j,k1)
              LAST(i,j+jmax*(nrank/p_col),5)=R2(i,j,k1)
            enddo
          enddo
        endif
        if (coordk(1).eq.1) then  
          do j=1,jmax
            do  i=0,i1
              FIRST(i,j+jmax*(nrank/p_col),1)=U1(i,j,0)
              FIRST(i,j+jmax*(nrank/p_col),2)=V1(i,j,0)
              FIRST(i,j+jmax*(nrank/p_col),3)=W1(i,j,0)
              FIRST(i,j+jmax*(nrank/p_col),4)=0.
              FIRST(i,j+jmax*(nrank/p_col),5)=1.
            enddo
          enddo
        endif
        do k=1,kmax
          do j=1,jmax
            do i=0,i1
              R22(i,j,k)=R2(i,j,k)
            enddo
          enddo
        enddo
        call  mpi_allreduce(LAST,LAST_tot,(i1+1)*jmax_tot*6,mpi_real8,mpi_sum,mpi_comm_world,ierr)
        call  mpi_allreduce(FIRST,FIRST_tot,(i1+1)*jmax_tot*6,mpi_real8,mpi_sum,mpi_comm_world,ierr)
        call load_save(1,Unew,Vnew,Wnew,Cnew,R22,Wx,Wy,Wz,dold,LAST_tot,FIRST_tot,nrank,istep,counter_post)
      endif
      if (mod(istep,out2D).eq.0) call output2d   (Uh,Vh,Wh,C1,R1,Ph,Qnew,Gnew,Knew,nrank,istep)
      dt=dtmax

      if (mod(istep,visbulk).eq. 0)call cmpinf(Wnew,T1,bulk,cbulk,stress,istep,nrank)
      if (nrank.eq.0 .and. mod(istep,visbulk).eq. 0)then
        call cpu_time(time2)
        timer=time3+time2-time1
        write(6,444) istep,time,dt,bulk,cbulk,stress,timer,stime3-stime2
444     format('step',i8,'t=',F10.4,'   dt=',F12.7,
     &           '  bulk-stress-timer=',4F20.7,' radtime=',F14.7)
      endif
      enddo  

!!!!!!!!!!!!!!!!!!!!!END OF MAIL LOOO!!!!!!!!!!!!!!!!!!!!!!
      call decomp_2d_finalize
      call mpi_finalize(ierr)
      stop
      end


      subroutine advance_h(dcdt,C11,R11,C1,R1,U1,V1,W1,Cu1,Cv1,Cw1,dold,RS,Qh,ekh,ABC1,ABC2,setold,rank)
      use decomp_2d
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer  setold,jm,jp
      real*8 ,  dimension(0:i1,0:j1,0:k1)::U1,V1,W1,C1,Cu1,Cv1,Cw1,R1,ekh,RS
      real*8 ,  dimension(0:i1,0:j1,0:k1)::U11,V11,W11,R11,C11
      real*8 ,  dimension(0:i1,jmax,kmax)::dcdt,dold,Qh
      real*8 ,  dimension(0:i1,jmax,kmax)::dnew,Crank,rhs,a,b,c
      real*8 ,  dimension(0:imx,jmax_tot,kmax)::rhst,at,bt,ct
      real*8 t1,t2,t3,t4,ABC_c,dif,adv,ABc1,ABc2
      real*8     eps,r2,phi2
      eps=1.0e-16
      ABC_c=0.5
      dif=1.0
      adv=1.0
      dnew=0.

      do k=0,kmax
        do j=0,jmax
          do i=0,imax
             U11(i,j,k)=U1(i,j,k)*0.5*(R1(i,j,k)+R1(i+1,j,k))
             V11(i,j,k)=V1(i,j,k)*0.5*(R1(i,j,k)+R1(i,j+1,k))
             W11(i,j,k)=W1(i,j,k)*0.5*(R1(i,j,k)+R1(i,j,k+1))
          enddo
        enddo
      enddo
      if(coordk(1).eq.1.and.k.eq.0) W11(:,:,0)=W1(:,:,0)
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            dnew(i,j,k) = - (
     1              (U11(i,j,k)*cu1(i,j,k)- U11(i-1,j,k)*cu1(i-1,j,k))/dr(i)
     &              +
     2              (V11(i,j,k)*cv1(i,j,k)- V11(i,j-1,k)*cv1(i,j-1,k))/dphi
     &              +
     3              (W11(i,j,k)*cw1(i,j,k)- W11(i,j,k-1)*cw1(i,j,k-1))/dz
     &              ) 
              enddo
           enddo
      enddo

      call diffcE(dnew,C1,Qh,ekh,rank,dif)

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
              rhs(i,j,k)= C11(i,j,k)*R11(i,j,k)+dt*(ABc1*dnew(i,j,k)+ABc2*dold(i,j,k)-Qh(i,j,k)/(Re*Prq*Pl)) 
            enddo
         enddo
      enddo
      if (setold.eq.1) dold = dnew

       do k=1,kmax
         do j=1,jmax
            do i=1,imax
               dcdt(i,j,k)=rhs(i,j,k)/RS(i,j,k)
           enddo
         enddo
      enddo
      return
      end

      subroutine advance_m(dudt,dvdt,dwdt,U11,V11,W11,R11,U1,V1,W1,R1,Wx,Wy,Wz,ekm,ABC1,ABC2,setold,rank)
      use decomp_2d
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer  setold
      real*8 ,  dimension(0:i1,0:j1,0:k1)::U11,V11,W11,R11,U1,V1,W1,R1,ekm
      real*8 ,  dimension(0:i1,jmax,kmax)::dudt,dvdt,dwdt,Wx,Wy,Wz
      real*8 ,  dimension(0:i1,jmax,kmax)::dnew,Crank,rhs,a,b,c
      real*8 ,  dimension(0:imx,jmax_tot,kmax)::rhst,at,bt,ct
      real*8 t1,t2,t3,t4,dpdz,ABC_c,dif,adv,ABc1,ABc2,ekma,ekmc,rhoa,rhob,rhoc
      
      ABC_c=0.5
      dif=1.0
      adv=1.0
      dnew=0.
      call advecu(dnew,U1,V1,W1,R1,adv,rank)
      call diffuE(dnew,U1,V1,W1,ekm,dif,rank)
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
              rhs(i,j,k)= U11(i,j,k)*(R11(i,j,k)+R11(i+1,j,k))*0.5+dt*(ABc1*dnew(i,j,k)+ABc2*wx(i,j,k))
            enddo
         enddo
      enddo

      if (setold.eq.1) wx = dnew

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            dUdt(i,j,k)=rhs(i,j,k)
          enddo
        enddo
      enddo
c********************************************************************
c     CALCULATE advection, diffusion and Force in t-direction
c     at the old(n-1) and new(n) timelevels
c********************************************************************
      dnew=0.
      call advecv(dnew,U1,V1,W1,R1,adv,rank)
      call diffvE(dnew,U1,V1,W1,ekm,dif,rank)
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
           rhs(i,j,k)= V11(i,j,k)*(R11(i,j,k)+R11(i,j+1,k))*0.5+dt*(ABc1*dnew(i,j,k)+ABc2*wy(i,j,k)) !+ABC_c*crank(i,j,k))
          enddo
        enddo
      enddo

      if (setold.eq.1) wy = dnew

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            dVdt(i,j,k)=rhs(i,j,k)
          enddo
        enddo
      enddo

c********************************************************************
c     CALCULATE advection, diffusion and Force in z-direction
c     at the old(n-1) and new(n) timelevels
c********************************************************************
      dnew=0.
      call advecw(dnew,U1,V1,W1,R1,adv,rank)
      call diffwE(dnew,U1,V1,W1,ekm,dif,rank)
      dpdz=0.
      if (periodic.eq.1)dpdz=1.
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            rhs(i,j,k)= W11(i,j,k)*(R11(i,j,k)+R11(i,j,k+1))*0.5+dt*(ABc1*dnew(i,j,k)+ABc2*wz(i,j,k)+dpdz)
          enddo
        enddo
      enddo

      if (setold.eq.1) wz = dnew

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            dWdt(i,j,k)=rhs(i,j,k)
          enddo
        enddo
      enddo
      return
      end

      subroutine bound_i_c(C1,Qh,ekh,rank,istart,cstart)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real*8 ,  dimension(0:i1,jmax,kmax)::C1,Qh
      real*8 ,  dimension(0:i1,0:j1,0:k1)::ekh
      integer ierr,istart,cstart
      real*8 Q,cb1,cb2,facttot,factor

      factor = istep*1.0/10000
      facttot = min(factor,1.0)
      do k=1,kmax
        do j=1,jmax
          if (cstart+k-1.lt.K_start_heat)then
            cb1=C1(1   ,j,k)
            cb2=C1(imax,j,k)
          else
            cb1=-0.15 !*(1-facttot)-0.2*facttot 
            cb2=-0.15 !*(1-facttot)-0.2*facttot
          endif
            C1(0,j,k)  =(C1(1,j,k)   -cb1)/(rp(1)   -ru(0))   *(rp(0)-ru(0))    +cb1
            C1(i1,j,k) =(C1(imax,j,k)-cb2)/(rp(imax)-ru(imax))*(rp(i1)-ru(imax))+cb2 
        enddo
      enddo
      return
      end

      subroutine bound_i_u(U1,V1,W1,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real*8 , intent(inout), dimension(0:i1,jmax,kmax)::U1,V1,W1
      integer ierr

      do k=1,kmax
        do j=1,jmax
           U1(0,j,k)    =   0.0
           U1(imax,j,k) =   0.0
           U1(i1,j,k)   = - U1(imax-1,j,k)
           V1(0,j,k)    = - V1(1,j,k)
           V1(i1,j,k)   = - V1(imax,j,k)
           W1(0,j,k)    = - W1(1,j,k)
           W1(i1,j,k)   = - W1(imax,j,k)
        enddo
      enddo
      return
      end

      subroutine insert(U1,V1,W1,C1,R2,Uh,Vh,Wh,Ch,dummy2,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real*8 , dimension(0:i1,0:j1,0:k1)::U1,V1,W1,C1,R2,Uh,Vh,Wh,Ch,dummy2
      U1=Uh
      V1=Vh
      W1=Wh
      C1=Ch
      R2=dummy2

      if (coordk(1).eq.1)then
        do j=0,j1
          do i=0,i1
            U1(i,j,0)=U1(i,j,1)
            V1(i,j,0)=V1(i,j,1)
            W1(i,j,0)=W1(i,j,1)
            C1(i,j,0)=0.
          enddo
        enddo
      endif
      if (coordk(2).eq.0)then
        do j=0,j1
          do i=0,i1
            U1(i,j,k1)=U1(i,j,kmax)
            V1(i,j,k1)=V1(i,j,kmax)
            W1(i,j,k1)=W1(i,j,kmax)
            C1(i,j,k1)=C1(i,j,kmax)
          enddo
        enddo
      endif

      return
      end

      subroutine bound_v(U1,V1,W1,U2,V2,W2,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real*8 , intent(inout), dimension(0:i1,0:j1,0:k1)::U1,V1,W1,U2,V2,W2

      real*8 Ub,Ub_tot
      integer ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!Inlet!!!!!!!!!!!!!!!!!!!
      if (coordk(1).eq.1) then
        do j=1,jmax
          do i=0,i1
            U1(i,j,0)=U_bound(i,j,isave)
            V1(i,j,0)=V_bound(i,j,isave)
            W1(i,j,0)=W_bound(i,j,isave)
          enddo
        enddo
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      Ub=0.
      if (coordk(2).eq.kmax_tot)then
        do j=1,jmax
          do i=1,imax
            Ub = Ub + W2(i,j,kmax)*(ru(i)-ru(i-1))*dphi
          enddo
        enddo
      endif
      call mpi_allreduce(Ub,Ub_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr) 
      Ub_tot = Ub_tot/(2.0*2.0*LoY)
      
      if (coordk(2).eq.kmax_tot)then
        do i=1,imax
          do j=1,jmax
            U1(i,j,k1) = U2(i,j,k1) - dt*Ub_tot*(U2(i,j,k1)-U2(i,j,kmax))/dz
            V1(i,j,k1) = V2(i,j,k1) - dt*Ub_tot*(V2(i,j,k1)-V2(i,j,kmax))/dz
            W1(i,j,k1) = W2(i,j,k1) - dt*Ub_tot*(W2(i,j,k1)-W2(i,j,kmax))/dz
          enddo
        enddo
      endif

      return
      end

      subroutine bound_c(C1,W2,C2,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real*8 ,  dimension(0:i1,0:j1,0:k1)::C1,C2,W2
      real*8 Ub,Ub_tot
      integer ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!Inlet!!!!!!!!!!!!!!!!!!!
      if (coordk(1).eq.1)then
        C1(:,:,0)=0.
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Ub=0.
      if (coordk(2).eq.kmax_tot)then
        do j=1,jmax
          do i=1,imax
            Ub = Ub + W2(i,j,kmax)*(ru(i)-ru(i-1))*dphi
          enddo
        enddo
      endif
      call mpi_allreduce(Ub,Ub_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      Ub_tot = Ub_tot/(2*2.0*LoY)

      if (coordk(2).eq.kmax_tot)then
        do i=0,i1
          do j=0,j1
            C1(i,j,k1) = C2(i,j,k1) - dt*Ub_tot*(C2(i,j,k1)-C2(i,j,kmax))/dz
          enddo
        enddo
      endif

      return
      end

      subroutine bound_m(U1,V1,W1,W2,RS,R1,R2,rank)
      implicit none
      include      'param.txt'
      include      'common.txt'
      include      'mpif.h'
      real*8 , intent(inout), dimension(0:i1,0:j1,0:k1)::U1,V1,W1,W2,RS,R1,R2
      real*8 flux,mass_balance,flux_tot,wfunc,wfunc_tot,Ub
      real*8 Wr(imax),Wr_tot(imax)
      integer ierr
      flux=0.
      wfunc=0.
      Wr=0.
      if (coordk(1).eq.1)then
        do j=1,jmax
          do i=0,i1
            U1(i,j,0)=U_bound(i,j,isave)
            V1(i,j,0)=V_bound(i,j,isave)
            W1(i,j,0)=W_bound(i,j,isave)
          enddo
        enddo
      endif


      if (coordk(2).eq.kmax_tot)then
        do i=1,imax
          do j=1,jmax
            Wr(i)=Wr(i)+W2(i,j,kmax)/jmax
          enddo
        enddo
      endif
      call mpi_allreduce(Wr,Wr_tot,imax,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      Wr_tot=Wr_tot/real(p_row)
      if (coordk(2).eq.kmax_tot) then
        do j=1,jmax
          do i=1,imax
            wfunc = wfunc + (ru(i)-ru(i-1))*dphi*Wr_tot(i)
          enddo
        enddo
      endif
      call mpi_allreduce(wfunc,wfunc_tot,1,mpi_real8,mpi_sum, mpi_comm_world,ierr)

      Ub=wfunc_tot/(2*2.0*LoY)
      if (coordk(2).eq.kmax_tot) then
        do j=0,j1
          do i=1,imax
            W1(i,j,kmax) = W2(i,j,kmax) - dt*Ub*(W2(i,j,kmax)-W2(i,j,kmax-1))/dz
            W1(i,j,kmax) = W1(i,j,kmax)*0.5*(RS(i,j,kmax)+RS(i,j,k1))
          enddo
        enddo
      endif
      if (coordk(2).eq.kmax_tot) then
        do j=1,jmax
          do i=1,imax
            flux=flux-W1(i,j,kmax)*(ru(i)-ru(i-1))*dphi
          enddo
        enddo
      endif
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            if (fDM.eq.1)then
              flux = flux - (3.*RS(i,j,k)-4.*R1(i,j,k)+R2(i,j,k))/(2.*dt)*dr(i)*dphi*dz
            else
              flux = flux - (RS(i,j,k)-R1(i,j,k))/(dt)*dr(i)*dphi*dz
            endif
          enddo
        enddo
      enddo
      if (coordk(1).eq.1) then
        do j=1,jmax
          do i=1,imax
            flux=flux+W1(i,j,0)*(ru(i)-ru(i-1))*dphi
          enddo
        enddo
      endif
      call mpi_allreduce(flux,flux_tot,1,mpi_real8,mpi_sum, mpi_comm_world,ierr)

      if (coordk(2).eq.kmax_tot)then
        mass_balance=flux_tot/wfunc_tot
        do j=1,jmax
          do i=1,imax
            W1(i,j,kmax) = W1(i,j,kmax)+mass_balance*Wr_tot(i)
          enddo
        enddo
      endif
      return
      end


      subroutine diffcE(putout,putin,Qh,ekh,rank,dif)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 ,  dimension(0:i1,0:j1,0:k1)::ekh,putin
      real*8 ,  dimension(0:i1,jmax,kmax)::putout,Qh
      integer   im,ip,jm,jp,km,kp
      real*8     dif,ekh_b,Q
c
      do k=1,kmax
        kp=k+1
        km=k-1
        do j=1,jmax
          jp=j+1
          jm=j-1
          do i=1,imax
            ip=i+1
            im=i-1
            ekh_b=ekh(i,j,k)+ekh(i,j,km)
            if(coordk(1).eq.1.and.k.eq.1) ekh_b=2./(Re*Pr)
            putout(i,j,k) = putout(i,j,k) +  0.5*(
     1             ( (ekh(i,j,k)+ekh(ip,j,k))*(putin(ip,j,k)-putin(i,j,k) )/(Rp(ip)-Rp(i))-
     1               (ekh(i,j,k)+ekh(im,j,k))*(putin(i,j,k) -putin(im,j,k))/(Rp(i)-Rp(im)) )/dr(i)
     &               +
     2         dif*( (ekh(i,j,k)+ekh(i,jp,k))*(putin(i,jp,k)-putin(i,j,k) ) -
     2               (ekh(i,j,k)+ekh(i,jm,k))*(putin(i,j,k) -putin(i,jm,k)))
     2             / (dphi*dphi)
     &              +
     3             ( (ekh(i,j,k)+ekh(i,j,kp))*(putin(i,j,kp)-putin(i,j,k) )-
     3               (ekh_b)                 *(putin(i,j,k) -putin(i,j,km)) )/(dz*dz)
     &              )

           enddo
        enddo
      enddo
!      do k=1,kmax
!       if (mod(rank,p_col).eq.0.and.k.lt.K_start_heat) then
!        Q=0.
!       else
!        Q=Qwall
!       endif
!       kp=k+1
!       km=k-1
!       do j=1,jmax
!         jp=j+1
!         jm=j-1
!         i=1
!         ip=i+1
!         im=i-1
!         ekh_b=ekh(i,j,k)+ekh(i,j,km)
!         if(mod(rank,p_col).eq.0.and.k.eq.1) ekh_b=2./(Re*Pr)
!         putout(i,j,k) = putout(i,j,k) +  0.5* (
!     1           (+Qh(0,j,k)/(Re*Prq*Pl)+
!     1           (ekh(ip,j,k)+ekh(i,j,k))*(putin(ip,j,k) -putin(i,j,k))/(Rp(ip)-Rp(i)))/(dr(i))
!     &               +
!     2      dif*((ekh(i,j,k)+ekh(i,jp,k))*(putin(i,jp,k)-putin(i,j,k)) -
!     2           (ekh(i,j,k)+ekh(i,jm,k))*(putin(i,j,k) -putin(i,jm,k)))
!     2         / (dphi*dphi)
!     &           +
!     3           ((ekh(i,j,k)+ekh(i,j,kp))*(putin(i,j,kp)-putin(i,j,k))-
!     3           (ekh_b)*(putin(i,j,k) -putin(i,j,km)))/(dz*dz)
!     &           )
!
!         enddo
!      enddo
!      do k=1,kmax
!       if (mod(rank,p_col).eq.0.and.k.lt.K_start_heat) then
!        Q=0.
!       else
!        Q=Qwall
!       endif
!       kp=k+1
!       km=k-1
!       do j=1,jmax
!         jp=j+1
!         jm=j-1
!         i=imax
!         ip=i+1
!         im=i-1
!         ekh_b=ekh(i,j,k)+ekh(i,j,km)
!         if(mod(rank,p_col).eq.0.and.k.eq.1) ekh_b=2./(Re*Pr)
!         putout(i,j,k) = putout(i,j,k) +  0.5* (
!     1           (Q/(Re*Pr)-Qh(i1,j,k)/(Re*Prq*Pl)-
!     1           (ekh(i,j,k)+ekh(im,j,k))*(putin(i,j,k) -putin(im,j,k))/(Rp(i)-Rp(im)))/(dr(i))
!     &               +
!     2      dif*((ekh(i,j,k)+ekh(i,jp,k))*(putin(i,jp,k)-putin(i,j,k)) -
!     2           (ekh(i,j,k)+ekh(i,jm,k))*(putin(i,j,k) -putin(i,jm,k)))
!     2         / (dphi*dphi)
!     &           +
!     3           ((ekh(i,j,k)+ekh(i,j,kp))*(putin(i,j,kp)-putin(i,j,k))-
!     3           (ekh_b)*(putin(i,j,k) -putin(i,j,km)))/(dz*dz)
!     &           )
!
!         enddo
!      enddo
      return
      end


      subroutine crank_c(putout,putin,ekh)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer   im,ip,jm,jp,km,kp

      real*8 ,  dimension(0:i1,0:j1,0:k1)::ekh,putin
      real*8 ,  dimension(0:i1,jmax,kmax)::putout
      do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do i=1,imax
               ip=i+1
               im=i-1

               putout(i,j,k) = 0.5*(
     1           ( (ekh(i,j,k)+ekh(i,jp,k)) * (putin(i,jp,k)-putin(i,j,k) ) -
     2        (ekh(i,j,k)+ekh(i,jm,k)) * (putin(i,j,k) -putin(i,jm,k))  )
     2         / ( Rp(i)*dphi*Rp(i)*dphi)
     &             )

            enddo
         enddo
      enddo
      return
      end

      subroutine diffuE(putout,Uvel,Vvel,Wvel,ekm,dif,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
c
c*****************************************************************
c
c     diffu calculates the diffusion of u-velocity, which is
c     the velocity in the radial direction.
c
c
c     In formula:  (4 terms)
c
c     1  d                  1  d                     d
c     - -- (r Sigma(r r)) + - ---- (Sigma(phi r)) + -- (Sigma(z r)) -
c     r dr                  r dphi                  dz
c
c
c     1
c     - Sigma(phi phi)
c     r
c
c     r   : direction  ---> explicit (subroutine diffu)
c     phi : direction  ---> implicit (subroutine predic)
c     z   : direction  ---> explicit (subroutine diffu)
c
c     on input :
c
c     putout            : advection part
c     Uvel,Vvel,Wvel    : contain velocities at n-1
c     ekm               : diffusion coefficients (for velocity) in
c     center points of the grid cells
c     dr,dphi,dz        : grid spacing in r, phi and z-direction
c     i1,j1,k1          : parameters for array-dimensions
c     ib,imax,jb,jmax,kb,kmax : range of gridpoints for which the
c     diffusion part has to be calculated
c     Ru,Rp             : radial positions of the U-velocity
c     component and the pressure location
c     respectively
c
c     on output :
c
c     putout            : advection and diffusion part
c     other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp
      real*8   eppo,epmo,epop,epom,drp,dzi,divUim,divUip,divUi,dif
      real*8 , dimension(0:i1,jmax,kmax)::putout
      real*8 , intent(in), dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,ekm

c
      dzi =1./dz
      do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do i=1,imax
               ip=i+1
               im=i-1
               drp = Rp(ip)-Rp(i)

               eppo = 0.25*(ekm(i,j,k)+ekm(ip,j,k) + ekm(ip,jp,k) + ekm(i,jp,k))
               epmo = 0.25*(ekm(i,j,k)+ekm(ip,j,k) + ekm(ip,jm,k) + ekm(i,jm,k))
               epop = 0.25*(ekm(i,j,k)+ekm(ip,j,k) + ekm(ip,j,kp) + ekm(i,j,kp))
               epom = 0.25*(ekm(i,j,k)+ekm(ip,j,k) + ekm(ip,j,km) + ekm(i,j,km))
               if(coordk(1).eq.1.and.k.eq.1)epom=1./Re
               divUim = (   Uvel(i,j,k) -        Uvel(im,j,k))/dr(i)
     &              + (     Vvel(i,j,k) -        Vvel(i,jm,k))/dphi 
     &              + (     Wvel(i,j,k) -        Wvel(i,j,km))/dz

               divUip = (   Uvel(ip,j,k) -  Uvel(i,j,k))/dr(ip)
     &              + (     Vvel(ip,j,k) -  Vvel(ip,jm,k))/dphi
     &              + (     Wvel(ip,j,k) -  Wvel(ip,j, km))/dz

               divUi = ((Uvel(ip,j,k)+Uvel(i,j,k)) - (Uvel(i,j,k)  +Uvel(im,j,k)))
     &                 /(2.*(Rp(ip)-Rp(i)))
     &              +  ((Vvel(ip,j,k)+Vvel(i,j,k)) - (Vvel(ip,jm,k)+Vvel(i,jm,k)))
     &                 /(2.*dphi)
     &              +  ((Wvel(ip,j,k)+Wvel(i,j,k)) - (Wvel(ip,j,km)+Wvel(i,j,km)))
     &                 /(2.*dz)


               putout(i,j,k) = putout(i,j,k) +
     1              (ekm(ip,j,k)*((Uvel(ip,j,k)-Uvel(i,j,k) )/
     1              dr(ip)-1./3.*divUip)-ekm(i,j,k) *((Uvel(i,j,k)
     1              -Uvel(im,j,k))/dr(i) -1./3.*divUim))/(0.5*(drp))
     &              +
     2              (eppo * ((Vvel(ip,j,k) - Vvel(i,j,k))/
     2              ( Rp(ip) - Rp(i) )
     2              +dif*(Uvel(i,jp,k) - Uvel(i,j,k))  / (dphi)
     2              )             -
     2               epmo * ((Vvel(ip,jm,k) - Vvel(i,jm,k))/
     2              (drp )
     2              +dif*(Uvel(i,j,k)  - Uvel(i,jm,k)) / (dphi)
     2              ) ) / (dphi)
     &              +
     3              (epop*((Uvel(i,j,kp)  - Uvel(i,j,k) )*dzi
     3              + (Wvel(ip,j,k)  - Wvel(i,j,k) ) / (Rp(ip) - Rp(i))
     3              )             -
     3               epom*((Uvel(i,j,k)   - Uvel(i,j,km))*dzi
     3              + (Wvel(ip,j,km) - Wvel(i,j,km)) / (Rp(ip) - Rp(i))
     3              ) ) * dzi
            enddo
         enddo
      enddo
      return
      end

      subroutine crank_u(putout,Uvel,Vvel,Wvel,ekm,Rho,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
c
c*****************************************************************
c
c     diffu calculates the diffusion of u-velocity, which is
c     the velocity in the radial direction.
c
c
c     In formula:  (4 terms)
c
c     1  d                  1  d                     d
c     - -- (r Sigma(r r)) + - ---- (Sigma(phi r)) + -- (Sigma(z r)) -
c     r dr                  r dphi                  dz
c
c
c     1
c     - Sigma(phi phi)
c     r
c
c     r   : direction  ---> explicit (subroutine diffu)
c     phi : direction  ---> implicit (subroutine predic)
c     z   : direction  ---> explicit (subroutine diffu)
c
c     on input :
c
c     putout            : advection part
c     Uvel,Vvel,Wvel    : contain velocities at n-1
c     ekm               : diffusion coefficients (for velocity) in
c     center points of the grid cells
c     dr,dphi,dz        : grid spacing in r, phi and z-direction
c     i1,j1,k1          : parameters for array-dimensions
c     ib,imax,jb,jmax,kb,kmax : range of gridpoints for which the
c     diffusion part has to be calculated
c     Ru,Rp             : radial positions of the U-velocity
c     component and the pressure location
c     respectively
c
c     on output :
c
c     putout            : advection and diffusion part
c     other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp
      real*8     putout(0:i1,j1-1,k1-1),eppo,epmo,rhojp,rhojm
       real*8 , intent(in), dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,ekm,Rho

      do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do i=1,imax
               ip=i+1
               im=i-1
               eppo = 0.25*(ekm(i,j,k)+ekm(ip,j,k) + ekm(ip,jp,k) + ekm(i,jp,k))
               epmo = 0.25*(ekm(i,j,k)+ekm(ip,j,k) + ekm(ip,jm,k) + ekm(i,jm,k))

               putout(i,j,k) = 
     2              ( eppo * ((Uvel(i,jp,k)  - Uvel(i,j,k) ) / ( Ru(i) * dphi )
     2              )             -
     2              epmo * ((Uvel(i,j,k)   - Uvel(i,jm,k)) / ( Ru(i) * dphi )
     2              ) ) / ( Ru(i) * dphi )
           enddo
         enddo
      enddo
!       deallocate(Uvel,Vvel,Wvel,ekm)
      return
      end

      subroutine diffvE(putout,Uvel,Vvel,Wvel,ekm,dif,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
c
c*****************************************************************
c
c     diffv calculates the diffusion of v-velocity, which is
c     the velocity in the tangential direction.
c
c
c     In formula:  (4 terms)
c
c     1  d                    1  d                       d
c     - -- (r Sigma(r phi)) + - ---- (Sigma(phi phi)) + -- (Sigma(z phi)) +
c     r dr                    r dphi                    dz
c
c
c     1
c     - Sigma(r phi)
c     r
c
c     r   : direction  ---> explicit (subroutine diffv)
c     phi : direction  ---> implicit (subroutine predic)
c     z   : direction  ---> explicit (subroutine diffv)
c
c     on input :
c
c     putout            : advection part
c     Uvel,Vvel,Wvel    : contain velocities at n-1
c     ekm               : diffusion coefficients (for velocity) in
c     center points of the grid cells
c     dr,dphi,dz        : grid spacing in r, phi and z-direction
c     i1,j1,k1          : parameters for array-dimensions
c     ib,imax,jb,jmax,kb,kmax : range of gridpoints for which the
c     diffusion part has to be calculated
c     Ru,Rp             : radial positions of the U-velocity
c     component and the pressure location
c     respectively
c
c     on output :
c
c     putout            : advection and diffusion part
c     other parameters  : all unchanged
c
c     Note :
c
c     The 'source' term [ Sigma (r phi) ] / r has been
c     incorporated into the radial derivative to avoid
c     interpolation problems at the centerline.
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp
      real*8   eppo,empo,eopp,eopm,divUjm,divUjp,dif
      real*8 , dimension(0:i1,jmax,kmax)::putout
      real*8 , intent(in), dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,ekm
c     -------------------------------------------start k-loop
      do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do i=1,imax
               ip=i+1
               im=i-1

               eppo = 0.25*(ekm(i,j,k)+ekm(ip,j,k) + ekm(ip,jp,k) + ekm(i,jp,k) )
               empo = 0.25*(ekm(i,j,k)+ekm(im,j,k) + ekm(i,jp,k)  + ekm(im,jp,k))
               eopp = 0.25*(ekm(i,j,k)+ekm(i,j,kp) + ekm(i,jp,k)  + ekm(i,jp,kp))
               eopm = 0.25*(ekm(i,j,k)+ekm(i,j,km) + ekm(i,jp,k)  + ekm(i,jp,km))
               if(coordk(1).eq.1.and.k.eq.1)eopm=1./Re
               divUjm = (  Uvel(i,j,k) -  Uvel(im,j,k))/dr(i)
     &              + (    Vvel(i,j,k) -  Vvel(i,jm,k))/dphi 
     &              + (    Wvel(i,j,k) -  Wvel(i,j,km))/dz

               divUjp = (  Uvel(i,jp,k)-  Uvel(im,jp,k ))/dr(i)
     &              + (    Vvel(i,jp,k)-  Vvel(i, j, k ))/dphi 
     &              + (    Wvel(i,jp,k)-  Wvel(i, jp,km))/dz


               putout(i,j,k) = putout(i,j,k) +
     1              (eppo*((Vvel(ip,j,k) - Vvel(i,j,k))/
     1              (Rp(ip) - Rp(i))
     1              + (Uvel(i,jp,k)  - Uvel(i,j,k) ) / (dphi)
     1              ) -
     1               empo*((Vvel(i,j,k)  - Vvel(im,j,k))/
     1              (Rp(i) - Rp(im))
     1              + (Uvel(im,jp,k) - Uvel(im,j,k)) / (dphi)
     1              ) ) / (dr(i))
     &              +
     2              (2*ekm(i,jp,k)*(
     2              dif* (Vvel(i,jp,k) - Vvel(i,j,k)  ) / dphi
     2              - 1./3.*divUjp)    -
     2               2*ekm(i,j,k) *(
     2              dif* (Vvel(i,j,k)  - Vvel(i,jm,k) ) / dphi
     2              - 1./3.*divUjm))/(dphi)
     &              +
     3              (eopp*((Vvel(i,j,kp)  - Vvel(i,j,k) ) / dz
     3              +      (Wvel(i,jp,k)  - Wvel(i,j,k) ) / dphi
     3              ) -
     3               eopm*((Vvel(i,j,k)   - Vvel(i,j,km)) / dz
     3              +      (Wvel(i,jp,km) - Wvel(i,j,km)) / dphi
     3              )  ) / dz
            enddo
         enddo
      enddo
!      deallocate(Uvel,Vvel,Wvel,ekm)
      return
      end

      subroutine crank_v(putout,Uvel,Vvel,Wvel,ekm,rho,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
c
c*****************************************************************
c
c     diffv calculates the diffusion of v-velocity, which is
c     the velocity in the tangential direction.
c
c
c     In formula:  (4 terms)
c
c     1  d                    1  d                       d
c     - -- (r Sigma(r phi)) + - ---- (Sigma(phi phi)) + -- (Sigma(z phi)) +
c     r dr                    r dphi                    dz
c
c
c     1
c     - Sigma(r phi)
c     r
c
c     r   : direction  ---> explicit (subroutine diffv)
c     phi : direction  ---> implicit (subroutine predic)
c     z   : direction  ---> explicit (subroutine diffv)
c
c     on input :
c
c     putout            : advection part
c     Uvel,Vvel,Wvel    : contain velocities at n-1
c     ekm               : diffusion coefficients (for velocity) in
c     center points of the grid cells
c     dr,dphi,dz        : grid spacing in r, phi and z-direction
c     i1,j1,k1          : parameters for array-dimensions
c     ib,imax,jb,jmax,kb,kmax : range of gridpoints for which the
c     diffusion part has to be calculated
c     Ru,Rp             : radial positions of the U-velocity
c     component and the pressure location
c     respectively
c
c     on output :
c
c     putout            : advection and diffusion part
c     other parameters  : all unchanged
c
c     Note :
c
c     The 'source' term [ Sigma (r phi) ] / r has been
c     incorporated into the radial derivative to avoid
c     interpolation problems at the centerline.
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp
      real*8 , dimension(0:i1,jmax,kmax)::putout
      real*8 , intent(in), dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,ekm,rho
       do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do i=1,imax
               ip=i+1
               im=i-1


               putout(i,j,k) = 
     2              ( ekm(i,jp,k) * ((Vvel(i,jp,k) - Vvel(i,j,k)  ) / dphi)    -
     2              ekm(i,j,k)  * ((Vvel(i,j,k)  - Vvel(i,jm,k) ) / dphi))/(0.5*Rp(i)*Rp(i)*dphi)
           enddo
         enddo
      enddo
!      deallocate(Uvel,Vvel,Wvel,ekm)
      return
      end

      subroutine diffwE(putout,Uvel,Vvel,Wvel,ekm,dif,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer  im,ip,jm,jp,km,kp
      real*8   epop,emop,eopp,eomp,divUkm,divUkp,dif
      real*8 , dimension(0:i1,jmax,kmax)::putout
      real*8 , intent(in), dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,ekm
c
      do k=1,kmax
        kp=k+1
        km=k-1
        do j=1,jmax
          jp=j+1
          jm=j-1
          do i=1,imax
            ip=i+1
            im=i-1
            epop = 0.25*(ekm(i,j,k)+ekm(i,j,kp) + ekm(ip,j,k) + ekm(ip,j,kp) )
            emop = 0.25*(ekm(i,j,k)+ekm(i,j,kp) + ekm(im,j,k) + ekm(im,j,kp) )
            eopp = 0.25*(ekm(i,j,k)+ekm(i,j,kp) + ekm(i,jp,k) + ekm(i,jp,kp) )
            eomp = 0.25*(ekm(i,j,k)+ekm(i,j,kp) + ekm(i,jm,k) + ekm(i,jm,kp) )

            divUkm = (      Uvel(i,j,k) -   Uvel(im,j,k))/dr(i)
     &             + (      Vvel(i,j,k) -   Vvel(i,jm,k))/dphi 
     &             + (      Wvel(i,j,k) -   Wvel(i,j,km))/dz

            divUkp = (      Uvel(i,j,kp)-   Uvel(im,j, kp))/dr(i)
     &             + (      Vvel(i,j,kp)-   Vvel(i, jm,kp))/dphi 
     &             + (      Wvel(i,j,kp)-   Wvel(i, j, k ))/dz


            putout(i,j,k) =  putout(i,j,k)+
     1             (epop*((Uvel(i,j,kp)  - Uvel(i,j,k) ) / dz
     1              +     (Wvel(ip,j,k)  - Wvel(i,j,k))  / (Rp(ip)-Rp(i))
     1                           ) -
     1              emop*((Uvel(im,j,kp) - Uvel(im,j,k)) / dz
     1              +     (Wvel(i,j,k)   - Wvel(im,j,k)) / (Rp(i)-Rp(im))
     1                          ) ) / (dr(i))
     &              +
     2             (eopp*((Vvel(i,j,kp)  - Vvel(i,j,k) ) / dz
     2              +dif* (Wvel(i,jp,k)  - Wvel(i,j,k) ) / dphi 
     2              ) -
     2              eomp*((Vvel(i,jm,kp) - Vvel(i,jm,k)) / dz
     2              +dif* (Wvel(i,j,k)   - Wvel(i,jm,k)) / dphi 
     2                          ) ) / (dphi)
     &              +
     3              (2.*ekm(i,j,kp)*((Wvel(i,j,kp)-Wvel(i,j,k))/dz - 1./3.*divUkp)-
     3               2.*ekm(i,j,k) *((Wvel(i,j,k)-Wvel(i,j,km))/dz - 1./3.*divUkm))/dz
            enddo
         enddo
      enddo
      return
      end

      subroutine crank_w(putout,Uvel,Vvel,Wvel,ekm,Rho,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
c
c*****************************************************************
c
c     diffw calculates the diffusion of w-velocity, which is
c     the velocity in the axial direction.
c
c
c     In formula:  (3 terms)
c
c     1  d                  1  d                     d
c     - -- (r Sigma(r z)) + - ---- (Sigma(phi z)) + -- (Sigma(z z))
c     r dr                  r dphi                  dz
c
c     r   : direction  ---> explicit (subroutine diffw)
c     phi : direction  ---> implicit (subroutine predic)
c     z   : direction  ---> explicit (subroutine diffw)
c
c     on input :
c
c     putout            : advection part
c     Uvel,Vvel,Wvel    : contain velocities at n-1
c     ekm               : diffusion coefficients (for velocity) in
c     center points of the grid cells
c     dr,dphi,dz        : grid spacing in r, phi and z-direction
c     i1,j1,k1          : parameters for array-dimensions
c     ib,imax,jb,jmax,kb,kmax : range of gridpoints for which the
c     diffusion part has to be calculated
c     Ru,Rp             : radial positions of the U-velocity
c     component and the pressure location
c     respectively
c
c     on output :
c
c     putout            : advection and diffusion part
c     other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp
      real*8   eopp,eomp,rhojp,rhojm
      real*8 , dimension(0:i1,jmax,kmax)::putout
      real*8 , intent(in), dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,ekm,rho
c
      do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do i=1,imax
               ip=i+1
               im=i-1
               eopp = 0.25*(ekm(i,j,k)+ekm(i,j,kp) + ekm(i,jp,k) + ekm(i,jp,kp) )
               eomp = 0.25*(ekm(i,j,k)+ekm(i,j,kp) + ekm(i,jm,k) + ekm(i,jm,kp) )



               putout(i,j,k) =  
     2              (  eopp * ((Wvel(i,jp,k)- Wvel(i,j,k)  ) /( Rp(i) * dphi )
     2              ) -
     2              eomp * ((Wvel(i,j,k) - Wvel(i,jm,k) ) /( Rp(i) * dphi )
     2              ) ) / ( Rp(i) * dphi )
           enddo
         enddo
      enddo
!      deallocate(Uvel,Vvel,Wvel,ekm)
      return
      end

      subroutine limiter(Cu,Cv,Cw,C1,U1,V1,W1)
      implicit none
      include 'param.txt'
      integer  im,ip,jm,jp,km,kp
      real*8,   dimension(0:i1,0:j1,0:k1) :: U1,V1,W1,C1,dcu,dcv,dcw
      real*8,  dimension(0:i1,jmax,kmax) :: Cu,Cv,Cw
      real*8     eps,r1,phi1,r2,phi2,r3,phi3,fak


      do k=0,kmax
         do j=0,jmax
            do i=0,imax
               dcu(i,j,k) = C1(i+1,j,k)-C1(i,j,k)
               dcv(i,j,k) = C1(i,j+1,k)-C1(i,j,k)
               dcw(i,j,k) = C1(i,j,k+1)-C1(i,j,k)
            enddo
         enddo
      enddo
      eps=1.0e-16
      fak=1./3.
      do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do i=1,imax
               ip=i+1
               im=i-1
              if(U1(i,j,k).ge.(0.0))then
                 r1=(dcu(i,j,k)+eps)/(dcu(im,j,k)+eps)
                 phi1=max(0.,min(2.*r1,min(fak*(1.+2.*r1),2.)))
                 Cu(i,j,k)=C1(i,j,k)+0.5*phi1*(dcu(im,j,k))
              else
                 r1=(dcu(i,j,k)+eps)/(dcu(ip,j,k)+eps)
                 phi1=max(0.,min(2.*r1,min(fak*(1.+2.*r1),2.)))
                 Cu(i,j,k)=C1(ip,j,k)+0.5*phi1*(-dcu(ip,j,k))
              endif
 
              if(V1(i,j,k).ge.(0.0))then
                 r2=(dcv(i,j,k)+eps)/(dcv(i,jm,k)+eps)
                 phi2=max(0.,min(2.*r2,min(fak*(1.+2.*r2),2.)))
                 Cv(i,j,k)=C1(i,j,k)+0.5*phi2*(dcv(i,jm,k))
              else
                 r2=(dcv(i,j,k)+eps)/(dcv(i,jp,k)+eps)
                 phi2=max(0.,min(2.*r2,min(fak*(1.+2.*r2),2.)))
                 Cv(i,j,k)=C1(i,jp,k)+0.5*phi2*(-dcv(i,jp,k))
              endif
 
               if(W1(i,j,k).ge.(0.0))then
                  r3=(dcw(i,j,k)+eps)/(dcw(i,j,km)+eps)
                  phi3=max(0.,min(2.*r3,min(fak*(1.+2.*r3),2.)))
                  Cw(i,j,k)=C1(i,j,k)+0.5*phi3*(dcw(i,j,km))
               else
                  r3=(dcw(i,j,k)+eps)/(dcw(i,j,kp)+eps)
                  phi3=max(0.,min(2.*r3,min(fak*(1.+2.*r3),2.)))
                  Cw(i,j,k)=C1(i,j,kp)+0.5*phi3*(-dcw(i,j,kp))
               endif
            enddo
         enddo
      enddo
      return
      end

      subroutine advecu(putout,Uvel,Vvel,Wvel,rho,adv,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer  im,ip,jm,jp,km,kp
      real*8 rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,adv
      real*8 , dimension(0:i1,jmax,kmax)::putout
      real*8 , dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,rho

      do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do  i=1,imax
               ip=i+1
               im=i-1
               rhojp=0.25*(rho(i,j,k)+rho(i,jp,k)+rho(ip,j,k)+rho(ip,jp,k))
               rhojm=0.25*(rho(i,j,k)+rho(i,jm,k)+rho(ip,j,k)+rho(ip,jm,k))
               rhokp=0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
               rhokm=0.25*(rho(i,j,k)+rho(i,j,km)+rho(ip,j,k)+rho(ip,j,km))
               if(coordk(1).eq.1.and.k.eq.1)rhokm=1.
               putout(i,j,k) = 0.0
               putout(i,j,k) = - 0.25 * (
     1              ((Uvel(i,j,k)+Uvel(ip,j,k))*(Uvel(i,j,k)+Uvel(ip,j,k))
     1              *rho(ip,j,k) -
     1               (Uvel(im,j,k)+Uvel(i,j,k))*(Uvel(i,j,k)+Uvel(im,j,k))
     1              *rho(i ,j,k)  )
     1              / (Rp(ip)-Rp(i))
     &              +
     2              adv*((Vvel(i,j,k) +Vvel(ip,j,k))*(Uvel(i,j,k)+Uvel(i,jp,k))*rhojp -
     2              (Vvel(i,jm,k)+Vvel(ip,jm,k))*(Uvel(i,j,k)+Uvel(i,jm,k))*rhojm  )
     2              / dphi 
     &              +
     3              ((Wvel(i,j,k)+Wvel(ip,j,k))*(Uvel(i,j,k)+Uvel(i,j,kp))*rhokp -
     3              (Wvel(i,j,km)+Wvel(ip,j,km))*(Uvel(i,j,k)+Uvel(i,j,km))*rhokm  )
     3              / dz
     &              )
            enddo
         enddo
      enddo

      return
      end

      subroutine advecv(putout,Uvel,Vvel,Wvel,rho,adv,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
c
c********************************************************************
c
c     advecv calculates the advection of the v-velocity, which is
c     the velocity in the tangential direction.
c
c     In formula:
c
c     1 d(ruv)     1 d(vv)     d(vw)     uv
c     - (  - ------  +  - -----  +  -----  +  --  )
c     r   dr       r  dphi      dz        r
c
c     on input :
c
c     putout            : "empty" (initialised to zero)
c     Uvel,Vvel,Wvel    : contain velocities at former timestep
c     Vtmp              : contains velocity at oldest timestep
c     dr,dphi,dz        : grid spacing in r, phi and z-direction
c     i1,j1,k1          : parameters for array-dimensions
c     ib,imax,jb,jmax,kb,kmax : range of gridpoints for which the
c     advection has to be calculated
c     Ru,Rp             : radial positions of the U-velocity
c     component and the pressure location
c     respectively
c
c     on output :
c
c     putout            : advection part
c     other parameters  : all unchanged
c
c********************************************************************
      integer  im,ip,jm,jp,km,kp
      real*8 , dimension(0:i1,jmax,kmax)::putout
      real*8 rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,adv
      real*8 , intent(in), dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,rho

c
      do k=1,kmax
        kp=k+1
        km=k-1
        do j=1,jmax
          jp=j+1
          jm=j-1
          do  i=1,imax
            ip=i+1
            im=i-1

            rhoip =0.25*(rho(i,j,k)+rho(ip,j,k)+rho(i,jp,k)+rho(ip,jp,k))
            rhoim =0.25*(rho(i,j,k)+rho(im,j,k)+rho(i,jp,k)+rho(im,jp,k))
            rhojp =rho(i,jp,k)
            rhojm =rho(i,j,k )

            rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
            rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(i,jp,k)+rho(i,jp,km))
            if(coordk(1).eq.1.and.k.eq.1)rhokm=1.
            putout(i,j,k) = 0.0
            putout(i,j,k) = - 0.25 * (
     1              (rhoip*
     1              (Uvel(i,j,k) +Uvel(i,jp,k) )*(Vvel(i,j,k)+Vvel(ip,j,k)) -
     1               rhoim*
     1              (Uvel(im,j,k)+Uvel(im,jp,k))*(Vvel(i,j,k)+Vvel(im,j,k))  )
     1              / (dr(i))
     &              +
     2              adv*((Vvel(i,j,k) +Vvel(i,jp,k))*(Vvel(i,j,k)+Vvel(i,jp,k))*rhojp -
     2                   (Vvel(i,jm,k)+Vvel(i,j,k) )*(Vvel(i,j,k)+Vvel(i,jm,k))*rhojm )
     2              / (dphi)
     &              +
     3              ((Wvel(i,j,k) +Wvel(i,jp,k) )*(Vvel(i,j,k)+Vvel(i,j,kp))*rhokp-
     3               (Wvel(i,j,km)+Wvel(i,jp,km))*(Vvel(i,j,k)+Vvel(i,j,km))*rhokm)
     3              / (dz)
     &              )
            enddo
         enddo
      enddo

      return
      end

      subroutine advecw(putout,Uvel,Vvel,Wvel,rho,adv,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
c********************************************************************
c
c     advecw calculates the advection of the w-velocity, which is
c     the velocity in the axial direction.
c
c     In formula:
c
c     1 d(ruw)     1 d(wv)     d(ww)
c     - (  - ------  +  - -----  +  -----  )
c     r   dr       r  dphi      dz
c
c     on input :
c
c     putout            : "empty" (initialised to zero)
c     Uvel,Vvel,Wvel    : contain velocities at former timestep
c     Wtmp              : contains velocity at oldest timestep
c     dr,dphi,dz        : grid spacing in r, phi and z-direction
c     i1,j1,k1          : parameters for array-dimensions
c     ib,imax,jb,jmax,kb,kmax : range of gridpoints for which the
c     advection has to be calculated
c     Ru,Rp             : radial positions of the U-velocity
c     component and the pressure location
c     respectively
c
c     on output :
c
c     putout            : advection part
c     other parameters  : all unchanged
c
c********************************************************************
      integer   im,ip,jm,jp,km,kp
      real*8 , dimension(0:i1,jmax,kmax)::putout
      real*8 rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,advcecw_w,adv
      real*8 , intent(in), dimension(0:i1,0:j1,0:k1)::Uvel,Vvel,Wvel,rho

        do k=1,kmax
         kp=k+1
         km=k-1
         do j=1,jmax
            jp=j+1
            jm=j-1
            do i=1,imax
               ip=i+1
               im=i-1

               rhoip=0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
               rhoim=0.25*(rho(i,j,k)+rho(i,j,kp)+rho(im,j,k)+rho(im,j,kp))
               rhojp=0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
               rhojm=0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jm,k)+rho(i,jm,kp))

               putout(i,j,k) = 0.0
               putout(i,j,k) = - 0.25 * (
     1              ((Uvel(i,j,k) +Uvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(ip,j,k))
     1              *rhoip -
     1               (Uvel(im,j,k)+Uvel(im,j,kp))*(Wvel(i,j,k)+Wvel(im,j,k))
     1              *rhoim )
     1              / (dr(i))
     &              +
     2              adv*((Vvel(i,j,k) +Vvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(i,jp,k))*rhojp -
     2                   (Vvel(i,jm,k)+Vvel(i,jm,kp))*(Wvel(i,j,k)+Wvel(i,jm,k))*rhojm  )
     2              / (dphi)
     &              +
     3                ((Wvel(i,j,k) +Wvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(i,j,kp))*rho(i,j,kp) -
     3                (Wvel(i,j,km)+Wvel(i,j,k)  )*(Wvel(i,j,k)+Wvel(i,j,km))*rho(i,j,k )  ) / dz
     3                                   )
            enddo
         enddo
      enddo

      return
      end

      subroutine name_to_hash(hash, node, length)
      implicit none
      integer    length,i
      integer*8  hash
      character*(length) node
      hash = 0
      do i=1,len(trim(node))
        hash = hash + (ichar(trim(node(i:i)))+1)*10**i
      enddo
      end subroutine


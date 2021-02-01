      subroutine chkdiv(U1,V1,W1,RS,R1,R2,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr,ll
      real*8 ,  dimension(0:i1,0:j1,0:k1)::U1,V1,W1,RS,R1,R2
      real*8   div,div1,divmax,divbar,divmax_tot,divbar_tot,rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
      divbar = 0.0
      divmax = 0.0
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               rhoip = 0.5*(RS(i,j,k)+RS(i+1,j,k))
               rhoim = 0.5*(RS(i,j,k)+RS(i-1,j,k))
               rhojp = 0.5*(RS(i,j,k)+RS(i,j+1,k))
               rhojm = 0.5*(RS(i,j,k)+RS(i,j-1,k))
               rhokp = 0.5*(RS(i,j,k)+RS(i,j,k+1))
               rhokm = 0.5*(RS(i,j,k)+RS(i,j,k-1))
               if (coordk(1).eq.1.and.k.eq.1)rhokm =1.
               if (fDM.eq.1)then
                  div =
     &                 (U1(i,j,k)*rhoip-U1(i-1,j,k)*rhoim)*dphi*dz   +
     &                 (V1(i,j,k)*rhojp-V1(i,j-1,k)*rhojm)*dr(i)*dz+
     &                 (W1(i,j,k)*rhokp-W1(i,j,k-1)*rhokm)*dphi*dr(i)+
     &                 ((3*RS(i,j,k)-4*R1(i,j,k)+R2(i,j,k))/(2.*dt))
     &                 *dr(i)*dphi*dz
              else
                 div =
     &                 (U1(i,j,k)*rhoip-U1(i-1,j,k)*rhoim)*dphi*dz +
     &                 (V1(i,j,k)*rhojp-V1(i,j-1,k)*rhojm)*dr(i)*dz+
     &                 (W1(i,j,k)*rhokp-W1(i,j,k-1)*rhokm)*dphi*dr(i)+
     &                 ((RS(i,j,k)-R1(i,j,k))/(dt))
     &                 *dr(i)*dphi*dz
              endif

               divbar = divbar+div
               div1    = abs(div)
               divmax = max(divmax,div1)
            enddo
         enddo
      enddo
      call mpi_allreduce(divbar,divbar_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(divmax,divmax_tot,1,mpi_real8,mpi_max,mpi_comm_world,ierr)


      if (rank.eq.0) write(*,100) divbar_tot,divmax_tot
 100  format('Mass_loss/gain:Tot= ',e15.6,' Max= ',e15.6)
      return
      end



       subroutine chkdt(Unew,Vnew,Wnew,ekm,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr
      real*8  dr2,dz2,df,df2,kcoeff,tmp1,tmp2,tmp3,Courant,dtmp
      real*8 , dimension(0:i1,jmax,kmax)::Unew,Vnew,Wnew,ekm
      dt = dtmax
      Courant = 0.15
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
               df = Rp(i)*dphi
               df2 = df * df
               dr2 = dr(i) * dr(i)
               dz2 = dz    * dz
               kcoeff =ekm(i,j,k)
               tmp1 = ( abs(Unew(i,j,k)) / ( Rp(i+1)-Rp(i) ) ) +
     &             ( 0*abs(Vnew(i,j,k)) /         df        ) +
     &              ( abs(Wnew(i,j,k)) /         dz        )
               tmp2 = (1.0/dr2 + 1.0/dz2+.0/df2)
               tmp3 = 1.0 / ( 1.0 * tmp2 * kcoeff + tmp1 )
               tmp3 = Courant*tmp3
               dt = min( dt , tmp3 )
            enddo
         enddo
      enddo
      dtmp = dt
      call mpi_allreduce(dtmp,dt,1,mpi_real8,mpi_min,mpi_comm_world,ierr)
      return
      end

      subroutine correc(Unew,Vnew,Wnew,U1,V1,W1,P1,RS,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 , dimension(0:i1,0:j1,0:k1)::U1,V1,W1,P1,RS
      real*8 , dimension(0:i1,jmax,kmax)::Unew,Vnew,Wnew

      if(coordk(2).eq.kmax_tot.and.periodic.ne.1)then
        do i=0,i1
          do j=0,j1
            p1(i,j,k1)=p1(i,j,kmax)
          enddo
        enddo
      endif

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            U1(i,j,k)=U1(i,j,k)-dt*(p1(i+1,j,k)-p1(i,j,k))/(Rp(i+1)-Rp(i))
          enddo
        enddo
      enddo
      do k=1,kmax
        do i=1,imax
          do j=1,jmax
             V1(i,j,k)=V1(i,j,k)-dt*(p1(i,j+1,k)-p1(i,j,k))/(dphi)
          enddo
        enddo
      enddo
      do j=1,jmax
        do i=1,imax
          do k=1,kmax
             W1(i,j,k)=W1(i,j,k)-dt*(p1(i,j,k+1)-p1(i,j,k))/dz
          enddo
        enddo
      enddo

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            Unew(i,j,k)=U1(i,j,k)/(0.5*(RS(i,j,k)+RS(i+1,j,k)))
            Vnew(i,j,k)=V1(i,j,k)/(0.5*(RS(i,j,k)+RS(i,j+1,k)))
            Wnew(i,j,k)=W1(i,j,k)/(0.5*(RS(i,j,k)+RS(i,j,k+1)))
            if(coordk(1).eq.1.and.k.eq.1) Wnew(i,j,k)=W1(i,j,k)
          enddo
        enddo
      enddo
      return
      end

      subroutine fillps(p,U1,V1,W1,RS,R1,R2,rank)
      implicit none
c
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr

      real*8 sumps,sumps_tot
      real*8 , dimension(0:i1,0:j1,0:k1)::U1,V1,W1,RS,R1,R2
      real*8 , dimension(0:i1,jmax,kmax)::p
      sumps=0.
      p=0.

      if (fDM.eq.1)then
        do  k=1,kmax
          do j=1,jmax
            do i=1,imax
              p(i,j,k)  = (
     1                   ( U1(i,j,k) - U1(i-1,j,k) ) /
     1                 (dr(i))
     2                 + ( V1(i,j,k) - V1(i,j-1,k) ) /
     2                 (dphi )
     3                 + ( W1(i,j,k) - W1(i,j,k-1) ) /
     3                 ( dz         ) )/dt
     &                 + (3.*RS(i,j,k)-4.*R1(i,j,k)+R2(i,j,k))/(2.*dt*dt)

!              sumps = sumps + p(i,j,k)*Rp(i)*dr(i)*dphi*dz

            enddo
          enddo
        enddo
      else
        do  k=1,kmax
          do j=1,jmax
            do i=1,imax
              p(i,j,k)  = (
     1                   ( U1(i,j,k) - U1(i-1,j,k) ) /
     1                 (dr(i))
     2                 + ( V1(i,j,k) - V1(i,j-1,k) ) /
     2                 (dphi)
     3                 + ( W1(i,j,k) - W1(i,j,k-1) ) /
     3                 ( dz         ) )/dt
     &                 +      (RS(i,j,k)-R1(i,j,k))/(dt*dt)

!                  sumps = sumps + p(i,j,k)*Rp(i)*dr(i)*dphi*dz

              enddo
            enddo
          enddo
        endif


!      open(1,file='test/wf')
!      do j=1,jmax
!      do i=0,i1
!      write(1,'(3I5,2F40.20)')i,j,k,W1(i,j,1),W1(i,j,0)
!      enddo
!      enddo
!      close(1)
!      open(1,file='test/pf')
!      do k=1,1
!      do j=1,jmax
!      do i=0,i1
!      write(1,'(3I5,F10.5)')i,j,k,p(i,j,k)
!      enddo
!      enddo
!      enddo
!      close(1)
        
        !      call mpi_allreduce(sumps,sumps_tot,1,mpi_real8,
!     &     mpi_sum,mpi_comm_world,ierr)

!      if (rank.eq.(px/2)) write(*,*) 'sumps_tot: ',sumps_tot
!       deallocate(U1,V1,W1,R1)
      return
      end


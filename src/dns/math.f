      SUBROUTINE SOLVEpois(rhs,Ru,Rp,dphi,dz,myid)
      use decomp_2d
      implicit none
      include 'param.txt'
      include 'mpif.h'
      integer     rank_j,rank_k,ierr,myid
      real*8      Ru(0:i1),Rp(0:i1),DPHI
      real*8      dz,dzi,pi,d(imax,jmax,kmax),bbb,z,norm
      real*8      a(imax),b(imax),c(imax),dd(imax)
      real*8      zrt(kmax_tot),yrt(jmax_tot)
      real*8      fk(kmax_tot),dfk(kmax_tot),fj(jmax_tot),dfj(jmax_tot)
      real*8      wj(4*jmax_tot+15),wk(4*kmax_tot+15),bb(imax)
      real*8      pj (0:imx,jmax_tot,kmax)
      real*8      pk(0:imx,jmax_tot/p_col,kmax_tot)
      real*8      rhs(0:i1,jmax,kmax)

c     generate tridiagonal systems

      pi = 4.*atan(1.)

      do i=1,imax
         a(i)= 1.0/((Rp(I)-Rp(I-1))*(Ru(I)-Ru(I-1)))
         b(i)=-(1.0/(Rp(I+1)-Rp(I))+1.0/(Rp(I)-Rp(I-1)))/
     &        (Ru(I)-Ru(I-1))
         c(i)= 1.0/((Rp(I+1)-Rp(I))*(Ru(I)-Ru(I-1)))
      end do
      b(1)=   -(1.0/(Rp(2)-Rp(1)))/((Ru(1)-Ru(0)))
      b(imax)=b(imax)+c(imax)
      c(imax)=0.
      a(1)=0.

      do i=1,imax
         dd(i) = 1.0 
      enddo
       dzi = 1./dz


      if (periodic.eq.1) then
         zrt(1)=0.
         do k=3,kmax_tot,2
            zrt(k-1)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax_tot)))**2
            zrt(k)=zrt(k-1)
         enddo
         zrt(kmax_tot)=-4.*dzi*dzi
      else
         do k=1,kmax_tot
            zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax_tot)))**2
         enddo
      endif


c     J --> direction      (yrt)
      yrt(1)=0.
      do j=3,jmax_tot,2
         yrt(j-1)=(-4./(dphi*dphi))*(sin(float((j-1))*pi/(2.*jmax_tot)))**2
         yrt(j  )=yrt(j-1)
      enddo 
      yrt(jmax_tot)= -4./(dphi*dphi)
  



      if (periodic.eq.1) then
         call vrffti(kmax_tot,wk)
      else
         call vcosqi(kmax_tot,wk)
      endif
      call vrffti(jmax_tot,wj)

      call transpose_x_to_y(rhs,pj)

!-----------------------------------------------------
c     J --> direction
      do k=1,kmax
         do i=0,imx
            do j=1,jmax_tot
               fj(j)=pj(i,j,k)
            enddo
            call vrfftf(1,jmax_tot,fj,dfj,1,wj)
            do j=1,jmax_tot
               pj(i,j,k)=fj(j)
            enddo
         enddo
      enddo

      call transpose_y_to_z(pj,pk)
!-----------------------------------------------------
!     K --> direction
      do j=1,jmax_tot/p_col
         do i=0,imx
            do k=1,kmax_tot
               fk(k)=pk(i,j,k)
            enddo
      if (periodic.eq.1) then
       call vrfftf(1,kmax_tot,fk,dfk,1,wk)
      else
       call vcosqb(1,kmax_tot,fk,dfk,1,wk)
      endif
            do k=1,kmax_tot
               pk(i,j,k)=fk(k)
            enddo
         enddo
      enddo
      call transpose_z_to_y(pk,pj)
      call transpose_y_to_x(pj,rhs)


      rank_j=int(myid/p_col)
      rank_k=mod(myid,p_col)

       do k=1,kmax
         do j=1,jmax
            !bbb        = b(1)+yrt(j+jmax*rank_j)*dd(1)+zrt(k+kmax*rank_k)
            bbb        = b(1)+yrt(j+xstart(2)-1)*dd(1)+zrt(k+xstart(3)-1)
            z          = 1./bbb
            d(1,j,k)   = c(1)*z
            rhs(1,j,k) = rhs(1,j,k)*z
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            do i=2,imax-1
               !bb(i)      = b(i)+yrt(j+jmax*rank_j)*dd(i)+zrt(k+kmax*rank_k)
               bb(i)      = b(i)+yrt(j+xstart(2)-1)*dd(i)+zrt(k+xstart(3)-1)
               z          = 1./(bb(i)-a(i)*d(i-1,j,k))
               d(i,j,k)   = c(i)*z
               rhs(i,j,k) = (rhs(i,j,k)-a(i)*rhs(i-1,j,k))*z
            enddo
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            !z = b(imax)+yrt(j+jmax*rank_j)*dd(imax)+zrt(k+kmax*rank_k)-a(imax)*d(imax-1,j,k)
            z = b(imax)+yrt(j+xstart(2)-1)*dd(imax)+zrt(k+xstart(3)-1)-a(imax)*d(imax-1,j,k)
            if (z .ne. 0.) then
               rhs(imax,j,k) = (rhs(imax,j,k)-a(imax)*rhs(imax-1,j,k))/z
            else 
               rhs(imax,j,k) = 0.
            endif
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            do  i=imax-1,1,-1
               rhs(i,j,k) = rhs(i,j,k)-d(i,j,k)*rhs(i+1,j,k)
            enddo
         enddo
      enddo

      call transpose_x_to_y(rhs,pj)

!     BACKWARD FFT ---> J direction
      do k=1,kmax
         do i=0,imx
            do j=1,jmax_tot
               fj(j)=pj(i,j,k)
            enddo
        call vrfftb(1,jmax_tot,fj,dfj,1,wj)
            do j=1,jmax_tot
               pj(i,j,k)=fj(j)
            enddo
         enddo
      enddo

      call transpose_y_to_z(pj,pk)
!-----------------------------------------------------
!     K --> direction
      do j=1,jmax_tot/p_col
         do i=0,imx
            do k=1,kmax_tot
               fk(k)=pk(i,j,k)
            enddo
      if (periodic.eq.1) then
          call vrfftb(1,kmax_tot,fk,dfk,1,wk)
      else
          call vcosqf(1,kmax_tot,fk,dfk,1,wk)
      endif
            do k=1,kmax_tot
               pk(i,j,k)=fk(k)
            enddo
         enddo
      enddo

      call transpose_z_to_y(pk,pj)
      call transpose_y_to_x(pj,rhs)
      rhs(i1,:,:)=rhs(imax,:,:)
      rhs(0,:,:) =rhs(1,:,:)
      if (rank.eq.0)norm=rhs(1,1,1)
         call MPI_BCAST(norm,1,mpi_real8,0,MPI_COMM_WORLD,ierr)
         do  k=1,kmax
            do j=1,jmax
               do i=1,imax
                rhs(i,j,k)=rhs(i,j,k)-norm
               enddo
            enddo
         enddo


      return
      end

      subroutine matrix(jmax,A,B,C,RHS)
******************************************************
****  Bendiks Jan Boersma (boersma@ymp.sara.nl)  ****
******************************************************
      implicit none
      integer jmax,j
      integer jmaxm1,jmaxm2
      real*8    A(jmax),B(jmax),C(jmax),RHS(jmax),
     &     D(jmax),tmpvar
      jmaxm1 = jmax - 1
      jmaxm2 = jmax - 2
         tmpvar = C(jmax)
         D(1)    = A(1)
      do  j=1,jmax-3
            B(j+1)    =  B  (j+1)  - (A(j+1)/B(j))  * C  (j)
            RHS(j+1)  =  RHS(j+1)  - (A(j+1)/B(j))  * RHS(j)
            D(j+1)    =              - (A(j+1)/B(j))  * D  (j)
            RHS(jmax) =  RHS(jmax) - (tmpvar/B(j)) * RHS(j)
            B(jmax)   =  B  (jmax) - (tmpvar/B(j)) * D  (j)
            tmpvar   =              - (tmpvar/B(j)) * C  (j)
      enddo
         B(jmaxm1)   = B(jmaxm1)   - ( A(jmaxm1)/B(jmaxm2) ) *
     &        C(jmaxm2)
         C(jmaxm1)   = C(jmaxm1)   - ( A(jmaxm1)/B(jmaxm2) ) *
     &        D(jmaxm2)
         RHS(jmaxm1) = RHS(jmaxm1) - ( A(jmaxm1)/B(jmaxm2) ) *
     &        RHS(jmaxm2)
         A(jmax)   = A(jmax)   - (tmpvar/B(jmaxm2))*C  (jmaxm2)
         B(jmax)   = B(jmax)   - (tmpvar/B(jmaxm2))*D  (jmaxm2)
         RHS(jmax) = RHS(jmax) - (tmpvar/B(jmaxm2))*RHS(jmaxm2)
         B(jmax)   = B(jmax)   - (A(jmax)/B(jmaxm1))*C  (jmaxm1)
         RHS(jmax) = RHS(jmax) - (A(jmax)/B(jmaxm1))*RHS(jmaxm1)
         RHS(jmax) = RHS(jmax) / B(jmax)
         RHS(jmaxm1) = ( RHS(jmaxm1) - C(jmaxm1) * RHS(jmax) ) /
     &        B(jmaxm1)
      do j=jmax-2,1,-1
            RHS(j) = ( RHS(j) - C(j)*RHS(j+1) - D(j)*RHS(jmax) )/ B(j)
      enddo
      end

c********************************************************************
c     spline (numerical recipes)
c********************************************************************
      subroutine spline(x, y, n, y2)
      implicit none
      integer   i, k, n, nmax
      parameter  (nmax=5000)
      real*8    yp1, ypn, x(n), y(n), y2(n), p, qn, sig, un, u(nmax)

      y2(1) = 0.
      u(1)  = 0.
      do i=2, n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     &        (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo

      qn=0.
      un=0.
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

      do k=n-1, 1, -1
         y2(k)=y2(k)*y2(k+1)+u(k)
      enddo

      return
      end

c********************************************************************
c     spline (numerical recipes)
c********************************************************************
      subroutine splint(xa,ya,y2a,n,x,y,khi,klo)
      implicit none
      integer n,k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n), a,b,h

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif

      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
       stop 'bad xa input in splint'
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end

c********************************************************************
c     spline (numerical recipes)
c********************************************************************
      subroutine splint_(xa,ya,y2a,n,x,y)
      implicit none
      integer n,k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n),a,b,h

      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif


      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end
c********************************************************************
c     Newton solver for wall boundary condition
c********************************************************************
      subroutine funcNewtonSolve(enth_i1,enth_imax,ekh_imax,Q)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 enth_i1,  enth_imax,ekh_imax,Q
      if (EOSmode.eq.0) call funcNewtonSolveIG(enth_i1,enth_imax,Q)
      if (EOSmode.eq.1) call funcNewtonSolveRG(enth_i1,enth_imax,ekh_imax,Q)
      end

c********************************************************************
c     Newton solver for wall boundary condition with PIG
c********************************************************************
      subroutine funcNewtonSolveIG(enth_i1, enth_imax,Q)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 enth_i1, enth_i11, enth_imax, fxValue, fxValue1,Q
      enth_i1 = enth_imax + (Rp(i1)-Rp(imax))*Q
      end

      subroutine funcNewtonSolveRG(enth_i1, enth_imax,ekh_imax,Q)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 enth_i1, enth_i11, enth_imax,ekh_imax,Q
      enth_i1 = enth_imax + (Rp(i1)-Rp(imax))*Q/(ekh_imax*Re*pr)
      end
c********************************************************************
c     Newton solver for wall boundary condition RG
c********************************************************************
      subroutine funcNewtonSolveRG1(enth_i1, enth_imax,ekh_imax,Q)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer nIterNewton,success
      real*8 enth_i1, enth_imax, fxValue, fxValue1,ekh_imax,Q,dx_step
      dx_step=1.0e-6
      success = 1
      fxValue = 1000.0
      nIterNewton = 0

      if (enth_imax.gt.2)    enth_imax =  2.0
      if (enth_imax.lt.-0.1) enth_imax = -0.1
      enth_i1 = enth_imax + (Rp(i1)-Rp(imax))*Q/(ekh_imax*Re*pr)

      do while (abs(fxValue).gt.1.0e-7)
         call funcNewtonBC(enth_i1,        enth_imax, fxValue,Q)
         call funcNewtonBC(enth_i1+dx_step, enth_imax, fxValue1,Q)

         enth_i1 = enth_i1 - fxValue/((fxValue1-fxValue)/dx_step)

         if (nIterNewton.gt.400) then
            fxValue = 0.0
            success = 0
         endif

         nIterNewton = nIterNewton + 1
      enddo

      if (success.eq.0) then
         write (*,*) 'newton didnt converge, enthi1= ',enth_i1,', ', nIterNewton, ', ', enth_i1
         stop
      endif
      end

c********************************************************************
c     function for wall boundary condition used by Newton solver
c********************************************************************
      subroutine funcNewtonBC(enth, enthIMAX, fxValue,Q)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer tabkhi,tabklo
      real*8 enth,laminter,cpinter,enthIMAX,fxValue,Q
      tabkhi = 0.0
      tabklo = 0.0
      call splint(enthTab,lamTab,lam2Tab,nTab,0.5*(enth+enthIMAX),laminter,tabkhi,tabklo)
      call splint(enthTab,cpTab,cp2Tab,nTab,0.5*(enth+enthIMAX),cpinter,tabkhi,tabklo)
      fxValue = enth - enthIMAX - (Rp(i1)-Rp(imax))*Q/(laminter/cpinter)
      end
c********************************************************************
c     PIG equation of state
c********************************************************************
      subroutine state(enth,temp,rho,mu,lam)
      implicit none
      include 'param.txt'
      integer istap
      real*8 enth(0:i1,0:j1,0:k1)
      real*8 rho (0:i1,0:j1,0:k1)
      real*8 mu  (0:i1,0:j1,0:k1)
      real*8 lam (0:i1,0:j1,0:k1)
      real*8 temp(0:i1,0:j1,0:k1)
      if (EOSmode.eq.0) call stateIG(enth,temp,rho,mu,lam)
      if (EOSmode.eq.1) call stateRG(enth,temp,rho,mu,lam)
      end

c********************************************************************
c     PIG equation of state IG
c********************************************************************
      subroutine stateIG(enth,temp,rho,mu,lam)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer istap
      real*8 enth(0:i1,0:j1,0:k1)
      real*8 rho (0:i1,0:j1,0:k1)
      real*8 mu  (0:i1,0:j1,0:k1)
      real*8 lam (0:i1,0:j1,0:k1)
      real*8 temp(0:i1,0:j1,0:k1)
      do k=0,k1
         do j=0,j1
            do i=0,i1
               rho(i,j,k) = 1.0 !/(enth(i,j,k)+1.0)
            enddo
         enddo
      enddo
      mu  = 1./Re
      lam = 1./(Re*Pr)
      ! temp= enth + 1.0
      return
      end

c********************************************************************
c     real gas equation of state RG
c********************************************************************
      subroutine stateRG(enth,temp,rho,mu,lam)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer tabkhi,tabklo,istap

      real*8 enth(0:i1,0:j1,0:k1)
      real*8 rho (0:i1,0:j1,0:k1)
      real*8 mu  (0:i1,0:j1,0:k1)
      real*8 lam (0:i1,0:j1,0:k1)
      real*8 temp(0:i1,0:j1,0:k1)

      do k=0,k1
         do j=0,j1
            do i=0,i1
               tabkhi = 0
               tabklo = 0
               call splint(enthTab,rhoTab,   rho2Tab,   nTab,enth(i,j,k),rho(i,j,k), tabkhi,tabklo)
               call splint(enthTab,muTab,    mu2Tab,    nTab,enth(i,j,k),mu(i,j,k),  tabkhi,tabklo)
               call splint(enthTab,lamocpTab,lamocp2Tab,nTab,enth(i,j,k),lam(i,j,k), tabkhi,tabklo)
               call splint(enthTab,tempTab,  temp2Tab,  nTab,enth(i,j,k),temp(i,j,k),tabkhi,tabklo)
               mu(i,j,k)  = mu(i,j,k)/(Re) 
               lam(i,j,k) = lam(i,j,k)/(Re*Pr)
            enddo
         enddo
      enddo
      return
      end





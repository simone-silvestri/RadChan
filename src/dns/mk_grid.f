
      subroutine mkgrid(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 Lz,rp_c(0:imax),rp_s(imax),x,dx,rnorm,rmax,drtot,pi
      pi    = 4.0*atan(1.0)
      dphi  = 2.0*LoY/(jmax_tot+1)
      dz    = 2.0*LoD/(kmax_tot+1)
      drtot = 2.0/imax
      rmax  = 2.0
      ru(0)=0.
      do i=1,imax/2
          x  = 1.*i/imax
          dx = 0.5-fact_mesh*(x-0.5)**2.0
          ru(i)=ru(i-1)+dx
      enddo
      rnorm = ru(imax/2)
      do i=1,imax/2
          ru(i)=2.0/2.*ru(i)/rnorm
      enddo
      do i=imax,imax/2+1,-1
          ru(i)=2.0-ru(imax-i)
      enddo

      do i=1,imax
        rp(i)=0.5*(ru(i)+ru(i-1))
        dr(i)=ru(i)-ru(i-1)
      enddo

      rp(0) = - rp(1)
      rp(i1)=ru(imax)+(Ru(imax)-rp(imax))
      if (rank.eq.0) then
      open(11,file = 'grid.txt')
      write(11,*) Re,Ru(imax)
      do i=1,imax
         write(11,'(i5,4F12.6)') i,Ru(i),Rp(i),dr(i)
      enddo
      endif
      close(11)



      end




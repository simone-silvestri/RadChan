!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! A few utility routines to find factors of integer numbers

  subroutine findfactor(iproc, factors, nfact)
    
    implicit none

    integer, intent(IN) :: iproc
    integer, intent(OUT) :: nfact
    integer, dimension(20), intent(OUT) :: factors
    integer :: i, imax

    imax = int(sqrt(real(iproc+1)))
    nfact = 0 
    i = 1
    do while(i<=imax)
       if (mod(iproc,i)==0) then
          nfact = nfact + 1
          factors(nfact) = i
       end if
       if (nfact==20) exit
       i = i+1
    end do

    return

  end subroutine findfactor


  subroutine primefactors(num, factors, nfact)

    implicit none
  
    integer, intent(IN) :: num
    integer, intent(OUT), dimension(*)::factors
    integer, intent(INOUT) :: nfact

    integer :: i, n
    
    i = 2  
    nfact = 1
    n = num 
    do
       if (mod(n,i) == 0) then
          factors(nfact) = i
          nfact = nfact + 1
          n = n / i
       else
          i = i + 1
       end if
       if (n == 1) then
          nfact = nfact - 1
          exit
       end if
    end do
    
    return

  end subroutine primefactors
  

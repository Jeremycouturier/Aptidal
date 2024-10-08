*********************************************************
*
!  routines de sorties et integrales premieres
!  derive de integ_jacobi8 (P.Robutel)
!
! modif M. GASTINEAU : extraction de mon_cin et renommage en mon_cin2
*********************************************************
#include "realc.f"

      Subroutine mon_cin2(nplan, mpl,X,XC,C)
!------- --------------------------------------------------------*
!       moment cinetique ( canonique helio)
!----------------------------------------------------------------*
      implicit none
      integer, intent(in) :: nplan
      real(TREAL),dimension(3), intent(out) :: C
      real(TREAL),dimension(3,nplan), intent(in) :: X,XC    
      real(TREAL),dimension(0:nplan), intent(in) :: mpl    
      integer :: i
      
      C = REALC(0e0)
      do i=1,nplan
         C(1) = C(1) + (X(2,i)*XC(3,i) - X(3,i)*XC(2,i))*mpl(i)
         C(2) = C(2) + (X(3,i)*XC(1,i) - X(1,i)*XC(3,i))*mpl(i)
         C(3) = C(3) + (X(1,i)*XC(2,i) - X(2,i)*XC(1,i))*mpl(i)
      end do         
      end 
 
   
   
!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file gmsun.f 
!!  \brief  routine de calcul du gm du soleil
!!
! history : creation 16/10/2017
!***********************************************************************
#include "realc.f"
      module mod_gmsun
      
        real(TREAL),parameter :: kgauss = REALC(0.01720209895e0)  ! constante de Gauss UA,jour
        real(TREAL),parameter :: GM_Sun_IAU2015_B3 = REALC(1.3271244e20)  ! constante de GM soleil m^3/s^2
        real(TREAL),parameter :: AU_IAU2012  = REALC(1.495978707e11) ! unite astronomique UAI 2012 en metre

      contains
      
!-----------------------------------------------
!!> calcul du valeur du GM du soleil en UA/an a partir de la constante de gauss 
!----------------------------------------------------
      function calcGM_sun_gauss() result(cG)
      implicit none
      real(TREAL) :: cG !< constante du valeur du GM du soleil   en AU, an
      cG = kgauss**2*(REALC(365.25e0))**2
      write (*,*) ' cGM_Sun Gauss =', cG
      end function calcGM_sun_gauss

      
!-----------------------------------------------
!!> calcul du valeur du GM en UA/an du soleil a partir de la Table 1 de "NOMINAL VALUES FOR SELECTED SOLAR AND PLANETARY QUANTITIES: IAU 2015 RESOLUTION B3"
!----------------------------------------------------
      function calcGM_sun_iau() result(cG)
      implicit none
      real(TREAL) :: cG !< constante du valeur du GM du soleil   en AU, an
      real(TREAL) :: fact !< constante du valeur du GM du soleil   en AU, an
      fact = (REALC(86400e0)*REALC(365.25e0))**2/(AU_IAU2012**3)
      cG = GM_Sun_IAU2015_B3*fact
      write (*,*) ' cGM_Sun IAU2015 =', cG
      end function calcGM_sun_iau

!-----------------------------------------------
!!> calcul du valeur du GM du soleil en UA/an
!----------------------------------------------------
      function calcGM_sun(refvalue) result(cG)
      implicit none
      integer, intent(in) :: refvalue
      real(TREAL) :: cG !< constante de Gauss  en AU, an
        select case(refvalue)
         case (0)
            cG = calcGM_sun_iau()
         case (1)
            cG = calcGM_sun_gauss()
            
         case default
          stop 'ref_gmsun invalid'
        end select 
      end function calcGM_sun

      end module mod_gmsun

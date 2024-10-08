!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file coordpart.f 
!!  \brief  Routines de transformation des positions vitesses 
!!   entre les systemes de coordonnees heliocentrique,
!!   barycentrique et de Jacobi pour les particules
!!
! history : creation 26/08/2014
!***********************************************************************

#include "realc.f"
 
       module mod_coordpart

       contains
       
      subroutine coord_vh2vb_part(npl,mpl,Vh, npart, Vhpart, Vbpart)
!*******************************************************
!> Transformation des vitesses heliocentriques (Vhpart)
!! en vitesses barycentriques (Vbpart) pour les particules.
!! issu de coord_vh2vb.
!*******************************************************
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masse des planetes et en 0 masse de l'etoile
      real(TREAL),dimension(3,npl),intent(in):: Vh !< vitesses heliocentriques des planetes
      integer, intent(in) :: npart !< nombre de particules
      real(TREAL),dimension(3,npart),intent(in):: Vhpart !< vitesses heliocentriques des particules
      real(TREAL),dimension(3,npart),intent(out):: Vbpart !< vitesses barycentriques des particules
      real(TREAL),dimension(3) :: vb0
      real(TREAL) :: masstot
      integer i
!----------- calcul de la vitesse barycentrique du soleil
      masstot = sum(mpl)
      do i=1,3
         Vb0(i) = -dot_product(mpl(1:),Vh(i,:))/masstot
      enddo
!------ 
      do i=1,3
        Vbpart(i,:) = Vhpart(i,:)+Vb0(i)
      end do
      end subroutine coord_vh2vb_part    

      subroutine coord_vb2vh_part(npl,mpl,Vb, npart, Vbpart, Vhpart)
!*******************************************************
!> Transformation des vitesses barycentriques (Vbpart)
!! en vitesses heliocentriques (Vhpart) pour les particules.
!! issu de coord_vb2vh
!*******************************************************
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masse des planetes et en 0 masse de l'etoile
      real(TREAL),dimension(3,npl),intent(in):: Vb !< vitesses barycentriques des planetes
      integer, intent(in) :: npart !< nombre de particules
      real(TREAL),dimension(3,npart),intent(in):: Vbpart !< vitesses barycentriques des particules
      real(TREAL),dimension(3,npart),intent(out):: Vhpart !< vitesses heliocentriques des particules
      real(TREAL),dimension(3) :: vb0
      real(TREAL) :: masstot
      integer i
!----------- calcul de la vitesse barycentrique du soleil
      masstot = sum(mpl)
      do i=1,3
         Vb0(i) = -dot_product(mpl(1:),Vb(i,:))/mpl(0)
      enddo
!------ 
      do i=1,3
        Vhpart(i,:) = Vbpart(i,:)-Vb0(i)
      end do
      end subroutine coord_vb2vh_part    

      subroutine coord_b2h_part(npl,mpl,Xb, npart, Xbpart, Xhpart)
!*******************************************************
!> Transformation d'une quantite (position ou vitesse) 
!! barycentrique en heliocentrique pour les particules.
!
!*******************************************************
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masse des planetes + masse de l'etoile en 0
      real(TREAL),dimension(3,npl),intent(in):: Xb   !< quantite barycentrique des planetes
      integer, intent(in) :: npart !< nombre de particules
      real(TREAL),dimension(3,npart),intent(in):: Xbpart   !< quantite barycentrique des particules
      real(TREAL),dimension(3,npart),intent(out):: Xhpart  !< quantite heliocentrique des particules
      integer i

      do i=1,3
        Xhpart(i,:) =  dot_product(mpl(1:),Xb(i,:))/mpl(0) + Xbpart(i,:)
      end do
      end  subroutine coord_b2h_part   

      end module



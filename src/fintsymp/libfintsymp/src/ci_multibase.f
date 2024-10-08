!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ci_multibase.f 
!!  \brief representation de plusieurs types de conditions initiales
!!      e.g. : condition initiale + distance minimale
!!
!!     Cette classe t_ci_multibase est utilisee pour stockee les informations lues par t_ciread_ci_multibase.
!!
! history : creation 29/01/2018
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour stocker plusieurs types de conditions initiales
!***********************************************************************
      module mod_ci_multibase
       use mod_ci_base


!***********************************************************************
!> @class t_ci_multibase
!! classe de stockage de types de conditions initiales
!!      e.g. : condition initiale + distance minimale
!!  
!***********************************************************************
      type, extends(t_ci_base) :: t_ci_multibase

          integer :: m_ci_nb = 0 !< nombre de type de conditions initiales
          
          type(t_ci_base_ptr), dimension(:), allocatable :: m_arci !<  type de ci (1:m_ci_nb)  


       contains
                    
          final ::  ci_multibase_destructor ! destructor       

          procedure:: set_ci_nb => ci_multibase_set_ci_nb ! fixe le nombre de type de ci

          procedure:: debug => ci_multibase_debug ! debug de la structure
          
      end type t_ci_multibase  

      contains
      
!***********************************************************************
!> @brief destructor
!***********************************************************************
       subroutine ci_multibase_destructor (this)
        implicit none
        type ( t_ci_multibase ) :: this  !< dummy argument
        integer k

        if (allocated(this%m_arci)) then
              do k=1 ,this%m_ci_nb
                 if (associated(this%m_arci(k)%m_ptr)) then 
                   deallocate(this%m_arci(k)%m_ptr) 
                   nullify(this%m_arci(k)%m_ptr) 
                 endif                
              enddo
              deallocate(this%m_arci)
        endif
       end subroutine ci_multibase_destructor

!***********************************************************************
!> @brief fixe le nombre de condtions et alloue la memoire
!***********************************************************************
      subroutine ci_multibase_set_ci_nb(this, nci)
       implicit none
       class(t_ci_multibase), intent(inout):: this    !< conteneur de conditions initiales
       integer, intent(in) :: nci                 !< nombre de type de ci
       
       this%m_ci_nb = nci
       
       allocate(this%m_arci(1:nci))
       
      end  subroutine ci_multibase_set_ci_nb  
          
!***********************************************************************
!> @brief affiche le debug de la structure 
!***********************************************************************
      subroutine ci_multibase_debug(this)
       implicit none
       class(t_ci_multibase), intent(in):: this    !< conteneur de conditions initiales
       
         integer k
         write(*,*) "ci_multibase_debug - start"
         call ci_base_debug(this)
         write(*,*) "m_ci_nb", this%m_ci_nb
         write(*,*) "distance minimale :"
         do k=1, this%m_ci_nb
          write(*,*) k," : "
          call this%m_arci(k)%m_ptr%debug()
         enddo
         write(*,*) "-----------------------"
           
      end  subroutine ci_multibase_debug  


      end module mod_ci_multibase

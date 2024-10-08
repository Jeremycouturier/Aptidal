!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciwrite_multibase.f 
!!  \brief ecriture de plusieurs fichiers texte de conditions initiales quelconques basees sur un stockage par ligne
!! l'ordre dans chaque fichier doit être le même : 1ere colonne de chaque fichier identique.  
!! un systeme a des informations dans plusieurs fichiers. 
!!     
!!     Le format du fichier est decrit dans ciread_multibase.f
!!
! history : creation 29/01/2018
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture plusieurs fichiers texte de conditions initiales
!***********************************************************************
      module mod_ciwrite_multibase
       use mod_ciwrite_base
       use mod_ci_multibase
       
!***********************************************************************
!> @class t_ciwrite_multibase
!! classe d'ecriture de plusieurs fichiers texte de conditions initiales
!!  
!***********************************************************************
      type, extends(t_ciwrite_base)  :: t_ciwrite_multibase
          
          type(t_ciwrite_base_ptr), dimension(:), allocatable, private   &
     &        :: m_writers !< ecrivains

       contains
          
          
          procedure :: set_writers => ciwrite_multibase_set_writers ! fixe les ecrivains
                    
          procedure :: fopen => ciwrite_multibase_fopen ! ouverture d'un fichier de conditions initiales
          procedure :: fclose => ciwrite_multibase_fclose ! fermeture d'un fichier de conditions initiales
          procedure :: fwrite => ciwrite_multibase_fwrite ! ecriture de toutes les conditions initiales issues d'un reader

          procedure :: fwriteline => ciwrite_multibase_fwriteline ! ecriture d'une seule ligne de conditions initiales         
          
      end type t_ciwrite_multibase  

     
      contains

!***********************************************************************
!> @brief fixe les ecrivains
!***********************************************************************
       subroutine ciwrite_multibase_set_writers(this, writers)
        implicit none
        class(t_ciwrite_multibase), intent(inout):: this  !< lecteur de ci
        type(t_ciwrite_base_ptr), dimension(:),intent(in) :: writers  !< tableau des ecrivains
        allocate(this%m_writers(1:size(writers)), source=writers)
       end  subroutine ciwrite_multibase_set_writers
 

!***********************************************************************
!> @brief ouvre le fichier en ecriture
!***********************************************************************
            subroutine ciwrite_multibase_fopen(this)
             implicit none
             class(t_ciwrite_multibase), intent(inout):: this  !< ecrivain de ci
             integer k

             do k=1, size(this%m_writers)
               call this%m_writers(k)%m_ptr%fopen()
             enddo

            end  subroutine ciwrite_multibase_fopen

!***********************************************************************
!> @brief ferme le fichier
!***********************************************************************
            subroutine ciwrite_multibase_fclose(this)
             implicit none
             class(t_ciwrite_multibase), intent(inout):: this  !< ecrivain de ci
             integer k
             
             do k=1, size(this%m_writers)
               call this%m_writers(k)%m_ptr%fclose()
             enddo
             
            end  subroutine ciwrite_multibase_fclose
    
!***********************************************************************
!> @brief ouvre le fichier du lecteur, lit toutes les conditions intiiales via le lecteur puis ferme le fichier du lecteur.
!! Chaque condition initiale est ensuite ecrite par l'ecrivain en ouvrant le fichier puis le fermant a la fin.
!***********************************************************************
         subroutine ciwrite_multibase_fwrite(this, reader) 
          use mod_ciread_multibase
          implicit none
          class(t_ciwrite_multibase), intent(inout):: this  !< ecrivain de ci
          class(t_ciread_base), intent(inout):: reader  !< lecteur de ci
          integer k
          
        select type(reader)
         class is(t_ciread_multibase)
          do k=1, size(this%m_writers)
        call this%m_writers(k)%m_ptr%fwrite(reader%m_readers(k)%m_ptr)
          enddo
         
         class default
          stop 'ci in ciwrite_multibase_fwrite : bad class'
        end select 
          
         end  subroutine ciwrite_multibase_fwrite

!***********************************************************************
!> @brief ecriture d'une ligne dans chaque fichier
!! Elle suppose que le fichier ouvert est bien au debut de la ligne a ecrire.
!***********************************************************************
       subroutine ciwrite_multibase_fwriteline(this, ci) 
        implicit none
        class(t_ciwrite_multibase), intent(inout) :: this !< lecteur du fichier de ci
        class(t_ci_base), intent(in) :: ci !< condition initiale a ecrire
        
        integer k
             
        select type(ci)
         class is(t_ci_multibase)
             write(*,*) 'ciwrite_multibase_fwri', size(this%m_writers)
          do k=1, size(this%m_writers)
             write(*,*) 'ciwrite_multibase_fwriteline k', k
         call this%m_writers(k)%m_ptr%fwriteline(ci%m_arci(k)%m_ptr)
              write(*,*) 'ciwrite_multibase_fwriteline fin k', k
         enddo
         
         class default
          stop 'ci in ciwrite_multibase_fwriteline : bad class'
        end select 

             write(*,*) 'ciwrite_multibase_fwriteline fin'
       end subroutine ciwrite_multibase_fwriteline

      end module mod_ciwrite_multibase

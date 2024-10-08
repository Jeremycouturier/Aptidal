!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciread_multibase.f 
!!  \brief lecture de plusieurs fichiers texte de conditions initiales quelconques basees sur un stockage par ligne
!! l'ordre dans chaque fichier doit être le même : 1ere colonne de chaque fichier identique.  
!! un systeme a des informations dans plusieurs fichiers. 
!!
!!
! history : creation 29/01/2018
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la lecture de plusieurs fichiers texte de conditions initiales
!***********************************************************************
      module mod_ciread_multibase
       use mod_ci_multibase
       use mod_ciread_base
!***********************************************************************
!> @class t_ciread_multibase
!! classe de lecture de plusieurs fichiers texte de conditions initiales
!!  
!***********************************************************************
      type, extends(t_ciread_base) :: t_ciread_multibase
          
          type(t_ciread_base_ptr),dimension(:),allocatable :: m_readers !< lecteurs

       contains
          
          
          procedure :: set_readers => ciread_multibase_set_readers ! fixe les lecteurs
          procedure :: fopen => ciread_multibase_fopen ! ouverture d'un fichier de conditions initiales
          procedure :: fclose => ciread_multibase_fclose ! fermeture d'un fichier de conditions initiales
          procedure :: get_countline => ciread_multibase_get_countline ! retourne le nombre de conditions initiales

          procedure :: fgoto_line => ciread_multibase_fgoto_line !se positionne a la ligne  specifiee        

          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          procedure :: freadline=>ciread_multibase_freadline ! lecture d'une seule ligne de conditions initiales
          
      end type t_ciread_multibase  

    
      contains
      
!***********************************************************************
!> @brief fixe les lecteurs
!***********************************************************************
       subroutine ciread_multibase_set_readers(this, readers)
        implicit none
        class(t_ciread_multibase), intent(inout):: this  !< lecteur de ci
        type(t_ciread_base_ptr), dimension(:),intent(in) :: readers  !< tableau des lecteurs
        allocate(this%m_readers(1:size(readers)), source=readers)
       end  subroutine ciread_multibase_set_readers
 

!***********************************************************************
!> @brief ouvre le fichier en lecture
!***********************************************************************
            subroutine ciread_multibase_fopen(this)
             implicit none
             class(t_ciread_multibase), intent(inout):: this  !< lecteur de ci
             integer k
             
             do k=1, size(this%m_readers)
               call this%m_readers(k)%m_ptr%fopen()
             enddo
            end  subroutine ciread_multibase_fopen

!***********************************************************************
!> @brief ferme le fichier
!***********************************************************************
            subroutine ciread_multibase_fclose(this)
             implicit none
             class(t_ciread_multibase), intent(inout):: this  !< lecteur de ci
             integer k

             do k=1, size(this%m_readers)
               call this%m_readers(k)%m_ptr%fopen()
             enddo
            end  subroutine ciread_multibase_fclose

!***********************************************************************
!> @brief ouvre le fichier/lit toutes les lignes et ferme le fichier.
!!
!! calcule le nombre de lignes (donc le nombre de conditions initiales)
!! verifie que tous les fichiers ont le meme nombre de lignes
!!
!! @return retourne le nombre de conditions initiales (nombre de ligne)
!***********************************************************************
            function ciread_multibase_get_countline(this) result(ret)
             implicit none
             class(t_ciread_multibase), intent(inout):: this  !< lecteur de ci
             integer :: ret ! valeur de retour
             integer :: k
             
             ret = this%m_readers(1)%m_ptr%get_countline()
             do k=2, size(this%m_readers)
              if (this%m_readers(k)%m_ptr%get_countline().ne.ret) then
                write(*,*) "nombre different de lignes entre "
                write(*,*) this%m_readers(1)%m_ptr%m_name
                write(*,*) " et "
                write(*,*) this%m_readers(k)%m_ptr%m_name
                stop
              endif
             enddo

            end  function ciread_multibase_get_countline
 

!***********************************************************************
!> @brief lecture d'une ligne de chacun des fichiers
!! Elle suppose que le fichier ouvert est bien au debut de la ligne a lire.
!! Elle verifie que chaque fichier a bien la meme ci
!! @return retourne .true. en cas de succes, retourne .false. en cas d'eof
!***********************************************************************
       function ciread_multibase_freadline(this, ci) result(ret)
        implicit none
        class(t_ciread_multibase), intent(inout) :: this !< lecteur du fichier de ci
        class(t_ci_base), pointer, intent(out) :: ci !< condition initiale lue en sortie
        logical :: ret !< =true => condition lue, =false => aucune condition lue
        
        class(t_ci_multibase), pointer :: cipla
        integer :: nci
        integer k
        
        ! allocation de la condition initiale
        nci = size(this%m_readers)
        allocate(cipla)
        call cipla%set_ci_nb(nci)
        
        ! lecture
        do k=1, nci
           ret=this%m_readers(k)%m_ptr%freadline(cipla%m_arci(k)%m_ptr)
           if (ret.eqv..false.) then
              if (k.gt.1) then
                   write(*,*) 'condition initiale non concordante !!!!'
                   write(*,*) '1er fichier ligne ', k
                   write(*,*) trim(this%m_readers(1)%m_ptr%m_name)
                   write(*,*) trim(cipla%m_arci(1)%m_ptr%m_id)
                   write(*,*) 'dernier fichier'
                   write(*,*) trim(this%m_readers(k)%m_ptr%m_name)
              endif
              exit
           else
             if (k.gt.1) then
                if (cipla%m_arci(k)%m_ptr%m_id.ne.                            &
     &              cipla%m_arci(1)%m_ptr%m_id) then
                   write(*,*) 'condition initiale non concordante !!!'
                   write(*,*) '1er fichier'
                   write(*,*) trim(cipla%m_arci(1)%m_ptr%m_id)
                   write(*,*) 'dernier fichier'
                   write(*,*) trim(cipla%m_arci(k)%m_ptr%m_id)
                   ret = .false.
                   exit
                endif
             endif
           endif
           
        enddo
        
        if (ret) then
          call cipla%set_id(cipla%m_arci(1)%m_ptr%m_id)
          ci => cipla
        else
          deallocate(cipla)
          ci =>NULL()
        endif
        
       end function ciread_multibase_freadline
      
!***********************************************************************
!> @brief se positionne au debut de la ligne specifiee nbline
!!
!***********************************************************************
            subroutine ciread_multibase_fgoto_line(this, nbline)
             implicit none
             class(t_ciread_multibase), intent(inout):: this  !< lecteur de ci
             integer, intent(in):: nbline  !< nombre de lignes a eviter depuis le debut du fichier
             integer k

             do k=1, size(this%m_readers)
               call this%m_readers(k)%m_ptr%fgoto_line(nbline)
             enddo
            end  subroutine ciread_multibase_fgoto_line

      end module mod_ciread_multibase

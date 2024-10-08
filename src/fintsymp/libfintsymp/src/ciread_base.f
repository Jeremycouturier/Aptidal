!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciread_base.f 
!!  \brief lecture d'un fichier texte de conditions initiales quelconques basees sur un stockage par ligne
!!
!!
! history : creation 14/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la lecture d'un fichier texte de conditions initiales quelconques
!***********************************************************************
      module mod_ciread_base
       use mod_ci_base
!***********************************************************************
!> @class t_ciread_base
!! classe de lecture d'un fichier texte de conditions initiales quelconques
!!  
!***********************************************************************
      type, abstract :: t_ciread_base
          
          integer, private :: m_nf !< numero fortran du fichier 
          character(len=512) :: m_name !< nom du fichier

       contains
          
          
          procedure :: set_name => ciread_base_set_name ! fixe le nom du fichier et le numero de fichier
          procedure :: fopen => ciread_base_fopen ! ouverture d'un fichier de conditions initiales
          procedure :: fclose => ciread_base_fclose ! fermeture d'un fichier de conditions initiales
          procedure :: get_countline => ciread_base_get_countline ! retourne le nombre de conditions initiales

          procedure :: fgoto_line => ciread_base_fgoto_line !se positionne a la ligne  specifiee        
          procedure :: fread_strline => ciread_base_fread_strline ! lecture d'une ligne sous forme de caracteres

          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          procedure(freadline), deferred :: freadline  ! lecture d'une seule ligne de conditions initiales
          
      end type t_ciread_base  

           abstract interface

            ! interface pour la lecture d'une ligne
            function freadline(this, ci) result(ret)
             import t_ciread_base
             import t_ci_base
             class(t_ciread_base), intent(inout) :: this !< lecteur du fichier de ci
             class(t_ci_base), pointer, intent(out) :: ci !<  condition initiale lue en sortie
             logical :: ret !< =true => condition lue, =false => aucune condition lue
            end function freadline


          end interface
    
!***********************************************************************
!> @class t_ciread_base_ptr
!! type pour les tableaux de pointeurs de lecteurs
!***********************************************************************
      type :: t_ciread_base_ptr 
          class(t_ciread_base),pointer :: m_ptr => NULL() !< pointeur vers un seul lecteur
      end type t_ciread_base_ptr 


      contains
      
!***********************************************************************
!> @brief fixe le nom du fichier et le numero fortran
!***********************************************************************
            subroutine ciread_base_set_name(this, name, nf)
             implicit none
             class(t_ciread_base), intent(inout):: this  !< lecteur de ci
             character(len=*), intent(in) :: name  !< nom du fichier
             integer, intent(in) :: nf  !< numero fortran du fichier
             this%m_name = name
             this%m_nf = nf
            end  subroutine ciread_base_set_name

!***********************************************************************
!> @brief ouvre le fichier en lecture
!***********************************************************************
            subroutine ciread_base_fopen(this)
             implicit none
             class(t_ciread_base), intent(inout):: this  !< lecteur de ci
             open(unit=this%m_nf,form="formatted",file=this%m_name,      &
     &        status='old')
            end  subroutine ciread_base_fopen

!***********************************************************************
!> @brief ferme le fichier
!***********************************************************************
            subroutine ciread_base_fclose(this)
             implicit none
             class(t_ciread_base), intent(inout):: this  !< lecteur de ci
             close(this%m_nf)
            end  subroutine ciread_base_fclose

!***********************************************************************
!> @brief lecture d'une seule ligne sous forme de chaine dans le fichier
!! @return retourne .true. en cas de succes, retourne .false. en cas d'eof
!***********************************************************************
            function ciread_base_fread_strline(this, str) result(ret)
             implicit none
             class(t_ciread_base), intent(inout):: this  !< lecteur de ci
             character(len=*), intent(inout):: str  !< chaine lue
             logical :: ret ! valeur de retour
             integer status
             
             ret = .true.
             read(this%m_nf,fmt='(A)', iostat=status) str
             if (status<0) then
             ! eof
              ret = .false.
             else if (status==0) then
             ! tout est ok
              ret = .true.
             else
             ! erreur
             write(*,*) "read reports an error ", status
             stop
             endif
             
            end  function ciread_base_fread_strline
             

!***********************************************************************
!> @brief ouvre le fichier/lit toutes les lignes et ferme le fichier.
!!
!! calcule le nombre de lignes (donc le nombre de conditions initiales) 
!!
!! @return retourne le nombre de conditions initiales (nombre de ligne)
!***********************************************************************
            function ciread_base_get_countline(this) result(ret)
             implicit none
             class(t_ciread_base), intent(inout):: this  !< lecteur de ci
             integer :: ret ! valeur de retour
             character(len=4096) :: str  ! chaine lue
             
             ret = 0
             call this%fopen()
             do while(this%fread_strline(str))
              ret = ret+1
             enddo
             call this%fclose()

            end  function ciread_base_get_countline
 
!***********************************************************************
!> @brief se positionne au debut de la ligne specifiee nbline
!!
!***********************************************************************
            subroutine ciread_base_fgoto_line(this, nbline)
             implicit none
             class(t_ciread_base), intent(inout):: this  !< lecteur de ci
             integer, intent(in):: nbline  !< nombre de lignes a eviter depuis le debut du fichier
             integer :: j
             logical :: ret
             character(len=4096) :: str  ! chaine lue
             
             rewind(this%m_nf)
             do j=1, nbline-1
              ret = this%fread_strline(str)
             enddo
            end  subroutine ciread_base_fgoto_line
 
      
      end module mod_ciread_base

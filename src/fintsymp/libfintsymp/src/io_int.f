!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file io_int.f 
!!  \brief ecriture sur fichier des integrales premieres
!!
!!     cette classe t_io_int gere la sortie des integrales premieres
!!
!!
! history : creation 01/08/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture des integrales premieres dans un fichier
!***********************************************************************
      module mod_io_int
       use mod_buffer
       
!***********************************************************************
!> @class t_io_int
!! classe d'ecriture des integrales premieres dans un fichier
!!  
!***********************************************************************
      type :: t_io_int

          integer, private :: m_nf !< numero fortran du fichier 
          character(len=512), private :: m_name !< nom du fichier
          type(t_buffer_consumer),pointer, private :: m_input_buffer     &
     &        =>NULL()      !< consommateur du buffer en entree

       contains
          
          procedure :: set_name =>t_io_int_set_name ! fixe le nom du fichier et le numero de fichier
          procedure :: fcreate => t_io_int_fcreate  ! cree le fichier et l'ouvre
          procedure :: fclose  => t_io_int_fclose   ! ferme le fichier
          procedure :: getconsumer  => t_io_int_getconsumer   ! retourne le consommateur du buffer
                  
          final ::  t_io_int_destructor ! destructor       
      end type t_io_int  

     
      contains
         
!***********************************************************************
!> @brief destructeur
!***********************************************************************
            subroutine t_io_int_destructor(this)
             implicit none
             type(t_io_int), intent(inout):: this  !< dummy argument
             if (associated(this%m_input_buffer)) then
                deallocate(this%m_input_buffer)
             endif
            end  subroutine t_io_int_destructor

!***********************************************************************
!> @brief fixe le nom du fichier et le numero fortran
!***********************************************************************
            subroutine t_io_int_set_name(this, name, nf)
             implicit none
             class(t_io_int), intent(inout):: this  !< dummy argument
             character(len=*), intent(in) :: name  !< nom du fichier
             integer, intent(in) :: nf  !< numero fortran du fichier
             procedure(t_procbufferfull), pointer :: pfcnfull
             procedure(t_procgraphdot), pointer :: pfcndot
             
             this%m_name = name
             this%m_nf = nf
             allocate(this%m_input_buffer)
             pfcnfull => t_io_int_onbufferfull
             pfcndot => t_io_int_ongraphdot
             call this%m_input_buffer%set_owner(this,                    &
     &                pfcnfull, pfcnfull, pfcndot)
            end  subroutine t_io_int_set_name

!***********************************************************************
!> @brief cree le fichier et l'ouvre
!***********************************************************************
            subroutine t_io_int_fcreate(this)
             implicit none
             class(t_io_int), intent(inout):: this  !< dummy argument
             open(this%m_nf,file=this%m_name,status='unknown')
            end  subroutine t_io_int_fcreate

!***********************************************************************
!> @brief ferme le fichier
!***********************************************************************
            subroutine t_io_int_fclose(this)
             implicit none
             class(t_io_int), intent(inout):: this  !< dummy argument
             close(this%m_nf)
            end  subroutine t_io_int_fclose

!***********************************************************************
!> @brief retourne le consommateur du buffer d'entree
!***********************************************************************
            function t_io_int_getconsumer(this) result(cs)
             implicit none
             class(t_io_int), intent(inout):: this  !< dummy argument
#if __INTEL_COMPILER <= 1600
             type(t_buffer_consumer), pointer :: cs 
#else
             class(t_buffer_consumer), pointer :: cs 
#endif
             cs => this%m_input_buffer
            end  function t_io_int_getconsumer


!***********************************************************************
!> @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
            subroutine t_io_int_onbufferfull(this, userdata, poswrite)
             implicit none
             class(t_buffer), intent(inout) :: this !< buffer de sortie
             class(*), intent(inout) :: userdata    !< de type t_io_int
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: j,k
             real(TREAL), dimension(5) :: R
             
              ! write(*,*) "t_io_int_onbufferfull"
             
             select type(userdata)
              class is(t_io_int)
               do j=1, poswrite
               call this%readdata(j, R)
               write(userdata%m_nf,1001) (R(k),k=1,5) 
1001  format(5(1X,FMT_TREAL))       !(etendu)
               enddo
              class default
               stop 't_io_int_onbufferfull : bad class'
             end select 
             
             call this%empty()
            
            end subroutine t_io_int_onbufferfull

!***********************************************************************
!> @brief  fonction appellee pour generer le graph au format dot
!***********************************************************************
       subroutine t_io_int_ongraphdot(this,userdata,nodesrc,             &
     &                                nodecounter, num_file)
        implicit none
        class(t_buffer), intent(inout) :: this !< buffer ayant genere l'appel
        class(*), intent(inout) :: userdata !< donnee utilisateur
        integer, intent(in) :: nodesrc !< numero du noeud  source 
        integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
        integer, intent(in) :: num_file !<numero du fichier de sortie
        
        call buffer_dot_writelabel(num_file, nodecounter+1,"I/O int")
        call buffer_dot_writenode(num_file,nodesrc,nodecounter+1)
        nodecounter = nodecounter+1
       end subroutine t_io_int_ongraphdot
 
       end module mod_io_int

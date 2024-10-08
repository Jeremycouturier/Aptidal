!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ctrl_ener.f 
!!  \brief Controle de la variation de l'erreur relative de l'energie du systeme
!!
!!
! history : creation 09/04/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le controle de la variation de l'erreur relative de l'energie du systeme
!***********************************************************************
      module mod_ctrl_ener
       use mod_finalconsumer
       
!***********************************************************************
!> @class t_ctrl_ener
!! Controle de la variation de l'erreur relative de l'energie du systeme.\n
!! Notifie le buffer predecesseur si l'erreur relative de l'energie varie trop.\n
!! Prend comme reference la valeur de l'energie au temps initial
!! Cela provoque une erreur de code -5.\n
!! Le buffer d'entree doit contenir les integrales premieres du systeme + le temps.\n
!!  
!***********************************************************************
      type, extends(t_finalconsumer) :: t_ctrl_ener
          private
          
         real(TREAL) ::   m_relenermax !< difference maximale de l'erreur relative sur l'energie. si l'erreur relative de l'energie >m_relenermax =-1, alors on genere une erreur 
         real(TREAL) ::   m_energie0   !< valeur initiale de l'energie
         logical     ::   m_energie0_init = .false. !< = true => m_energie0 est initialise

      contains
          procedure :: set_relenermax => ctrl_ener_set_relenermax ! definit la variation maximale de l'erreur relative sur l'energie
          
      end type t_ctrl_ener  

      contains
           

!***********************************************************************
!> @brief definit la variation maximale de l'erreur relative sur l'energie
!***********************************************************************
            subroutine ctrl_ener_set_relenermax(this, relener)
             implicit none
             class(t_ctrl_ener), intent(inout):: this  !< dummy argument
             real(TREAL), intent(in) :: relener  !< variation maximale de l'erreur relative sur l'energie
             procedure(t_procbufferfull), pointer :: pfcnfull
             
             call this%set_graphdotname("CTRL energie")
             pfcnfull=>ctrl_ener_onbufferfull
             call this%setconsumer(pfcnfull, pfcnfull)
     
             this%m_relenermax = relener
             this%m_energie0_init = .false.
             
            end  subroutine ctrl_ener_set_relenermax


!***********************************************************************
!> @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
            subroutine ctrl_ener_onbufferfull(this, userdata, poswrite)
             implicit none
             class(t_buffer), intent(inout) :: this !< buffer de sortie
             class(*), intent(inout) :: userdata    !< de type t_ctrl_ener
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: j
             real(TREAL), dimension(5) :: R
             real(TREAL) relener, t, ener, e0
             type(t_arret) :: arret
                          
             select type(userdata)
              class is(t_ctrl_ener)
               e0 = userdata%m_energie0
               do j=1, poswrite
                call this%readdata(j, R)
                t = R(1)
                ener = R(2)
                !write(*,*) "t", t, ener
                if (.not.userdata%m_energie0_init) then
                 userdata%m_energie0_init = .true.
                 userdata%m_energie0 = ener
                else
                 relener = abs(ener/e0)
                 if (relener.gt.userdata%m_relenermax) then
                  ! valeur trop grande de l'energie => on notifie le predeccesseur
                  call arret%set_error(-5)
                  call arret%set_time(t)
                  call arret%set_value(relener)
                  call this%notify_successor_error(arret)
                 endif
                endif
              enddo
              class default
               stop 'ctrl_ener_onbufferfull : bad class'
             end select 
             
             call this%empty()
            
            end subroutine ctrl_ener_onbufferfull
 
      end module mod_ctrl_ener

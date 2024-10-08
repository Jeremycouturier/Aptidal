!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ctrl_distpartstar.f 
!!  \brief Controle de la distance des particules a l'etoile
!!
!!
! history : creation 17/09/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le controle de de la distance des particules a l'etoile
!***********************************************************************
      module mod_ctrl_distpartstar
       use mod_finalconsumer
       
!***********************************************************************
!> @class t_ctrl_distpartstar
!! Controle de la distance des particules a l'etoile.\n
!! Notifie le buffer predecesseur si la particule est trop proche ou trop loin.\n
!! Cela provoque une erreur de code -6 ou -7 mais cette erreur ne provoque pas l'arret des integrations.\n
!! Le buffer d'entree doit contenir les positions/vitesses heliocentriques  + le temps.\n
!!  
!***********************************************************************
      type, extends(t_finalconsumer) :: t_ctrl_distpartstar
          private
          
         real(TREAL) ::   m_distmin !< distance minimale a l'etoile. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
         real(TREAL) ::   m_distmax !< distance maximale a l'etoile. si m_distmax =-1, pas de controle, si m_distmax>=0, genere une erreur si distance>m_distmax
         integer     ::   m_plan_nb !<nombre de planetes (identique au nombre de planetes dans le buffer)
         integer     ::   m_part_nb !<nombre de particules a verifier (identique au nombre de particules dans le buffer)

      contains
          procedure :: set_minmax => ctrl_distpartstar_set_minmax ! definit la distance minimale et maximale a l'etoile
          
      end type t_ctrl_distpartstar  

      contains
           

!***********************************************************************
!> @brief definit la distance minimale ou maximale des particules a l'etoile
!***********************************************************************
            subroutine ctrl_distpartstar_set_minmax(this,nbplan,nbpart,         &
     &                      mins, maxs)
             implicit none
             class(t_ctrl_distpartstar), intent(inout):: this  !< dummy argument
             integer, intent(in) :: nbplan  !<nombre de planetes  (identique au nombre de planetes dans le buffer)
             integer, intent(in) :: nbpart  !<nombre de particules a verifier (identique au nombre de planetes dans le buffer)
             real(TREAL), intent(in) :: mins  !< distance minimale a l'etoile. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
             real(TREAL), intent(in) :: maxs  !< distance maximale a l'etoile. si m_distmax =-1, pas de controle, si m_distmax>=0, genere une erreur si distance>m_distmax
             procedure(t_procbufferfull), pointer :: pfcnfull
             
             this%m_part_nb = nbpart
             this%m_plan_nb = nbplan
             this%m_distmin = mins
             this%m_distmax = maxs
             
             call this%set_graphdotname("CTRL distpartstar min/max")
             pfcnfull => ctrl_distpartstar_onbufferfull
             call this%setconsumer(pfcnfull, pfcnfull)
             
            end  subroutine ctrl_distpartstar_set_minmax

!***********************************************************************
!> @brief verifie les distances dans le buffer
!***********************************************************************
            subroutine ctrl_distpartstar_processbuffer(ctrl, buffer,           &
     &              poswrite)
             implicit none
             class(t_buffer), intent(inout) :: buffer !< buffer de sortie
             class(t_ctrl_distpartstar), intent(inout) :: ctrl    !< de type t_ctrl_distpartstar
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: j
             integer k, ideb, ifin
             real(TREAL), dimension(1+6*ctrl%m_plan_nb+                        &
     &              6*ctrl%m_part_nb) :: R
             real(TREAL), dimension(3,ctrl%m_part_nb) :: Ph
             real(TREAL)  t, dmin2, dmax2, d
             
             integer nplan, npart
             type(t_arret) :: arret
             
               nplan = ctrl%m_plan_nb
               npart = ctrl%m_part_nb
               dmin2 = ctrl%m_distmin**2
               dmax2 = ctrl%m_distmax**2
               do j=1, poswrite
                call buffer%readdata(j, R)
                t = R(1)
                ideb = 2+6*nplan
                ifin = 1+3*npart+6*nplan
                Ph = reshape(R(ideb:ifin), shape(Ph))
                do k=1, npart
                 if (buffer%m_mce(k)%m_stop.eq.0) then
                     d=dot_product(Ph(:,k),Ph(:,k))
                     if (d.lt.dmin2) then
                      ! valeur trop petite de la distance  => on notifie le predecesseur
                      call arret%set_error(-6)
                      call arret%set_time(t)
                      call arret%set_value(sqrt(d))
                      call arret%set_body(0)
                      buffer%m_mce(k) = arret
                     else if (d.gt.dmax2) then
                      ! valeur trop grande de la distance  => on notifie le predecesseur
                      call arret%set_error(-7)
                      call arret%set_time(t)
                      call arret%set_value(sqrt(d))
                      call arret%set_body(0)
                      buffer%m_mce(k) = arret
                     endif
                 endif
                enddo
               enddo
            
            end subroutine ctrl_distpartstar_processbuffer

!***********************************************************************
!> @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
            subroutine ctrl_distpartstar_onbufferfull(this, userdata,           &
     &              poswrite)
             implicit none
             class(t_buffer), intent(inout) :: this !< buffer de sortie
             class(*), intent(inout) :: userdata    !< de type t_ctrl_distpartstar
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)

             select type(userdata)
              class is(t_ctrl_distpartstar)
               call ctrl_distpartstar_processbuffer(userdata, this,           &
     &              poswrite)

              class default
               stop 'ctrl_distpartstar_onbufferfull : bad class'
             end select 
             
             call this%empty()
            
            end subroutine ctrl_distpartstar_onbufferfull
 
      end module mod_ctrl_distpartstar

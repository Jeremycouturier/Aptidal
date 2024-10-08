!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ctrl_distpartpla.f 
!!  \brief Controle de la distance des particules aux planetes
!!
!!
! history : creation 17/09/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le controle de de la distance des particules aux planetes
!***********************************************************************
      module mod_ctrl_distpartpla
       use mod_finalconsumer
       
!***********************************************************************
!> @class t_ctrl_distpartpla
!! Controle de la distance des particules aux planetes.\n
!! Notifie le buffer predecesseur si la particule est trop proche ou trop loin.\n
!! Cela provoque une erreur de code -9 mais cette erreur ne provoque pas l'arret des integrations.\n
!! Le buffer d'entree doit contenir les positions/vitesses heliocentriques  + le temps.\n
!!  
!***********************************************************************
      type, extends(t_finalconsumer) :: t_ctrl_distpartpla
          private
          
         real(TREAL) ::   m_distmin !< distance minimale aux planetes. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
         integer     ::   m_plan_nb !<nombre de planetes (identique au nombre de planetes dans le buffer)
         integer     ::   m_part_nb !<nombre de particules a verifier (identique au nombre de particules dans le buffer)

      contains
          procedure :: set_minmax => ctrl_distpartpla_set_minmax ! definit la distance minimale et maximale a l'etoile
          
      end type t_ctrl_distpartpla  

      contains
           

!***********************************************************************
!> @brief definit la distance minimale des particules aux planetes
!***********************************************************************
            subroutine ctrl_distpartpla_set_minmax(this, nbplan, nbpart,         &
     &                      mins)
             implicit none
             class(t_ctrl_distpartpla), intent(inout):: this  !< dummy argument
             integer, intent(in) :: nbplan  !<nombre de planetes  (identique au nombre de planetes dans le buffer)
             integer, intent(in) :: nbpart  !<nombre de particules a verifier (identique au nombre de planetes dans le buffer)
             real(TREAL), intent(in) :: mins  !< distance minimale a l'etoile. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
             procedure(t_procbufferfull), pointer :: pfcnfull
             this%m_part_nb = nbpart
             this%m_plan_nb = nbplan
             this%m_distmin = mins
             
             call this%set_graphdotname("CTRL distpartpla min/max")
             pfcnfull => ctrl_distpartpla_onbufferfull
             call this%setconsumer(pfcnfull, pfcnfull)
             
            end  subroutine ctrl_distpartpla_set_minmax

!***********************************************************************
!> @brief verifie les distances dans le buffer
!***********************************************************************
            subroutine ctrl_distpartpla_processbuffer(ctrl,  buffer,           &
     &              poswrite)
             implicit none
             class(t_buffer), intent(inout) :: buffer !< buffer de sortie
             class(t_ctrl_distpartpla), intent(inout) :: ctrl    !< de type t_ctrl_distpartpla
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: j
             integer k, p, ideb, ifin
             real(TREAL), dimension(1+6*ctrl%m_plan_nb                         &
     &             +6*ctrl%m_part_nb) :: R
             real(TREAL), dimension(3,ctrl%m_plan_nb) :: Phpla
             real(TREAL), dimension(3,ctrl%m_part_nb) :: Phpart
             real(TREAL), dimension(3) :: XPP
             real(TREAL)  t, dmin2, d
             
             integer nplan, npart
             type(t_arret) :: arret
             
               nplan = ctrl%m_plan_nb
               npart = ctrl%m_part_nb
               dmin2 = ctrl%m_distmin**2
               do j=1, poswrite
                call buffer%readdata(j, R)
                t = R(1)
                ideb = 2
                ifin = 1+3*nplan
                Phpla = reshape(R(ideb:ifin), shape(Phpla))
                ideb = 2+6*nplan
                ifin = 1+3*npart+6*nplan
                Phpart = reshape(R(ideb:ifin), shape(Phpart))
                do k=1, npart
                 if (buffer%m_mce(k)%m_stop.eq.0) then
                    do p=1, nplan
                         XPP = Phpart(:,k)-Phpla(:,p)
                         d=dot_product(XPP,XPP)
                         if (d.lt.dmin2) then
                          ! valeur trop petite de la distance  a la planete p 
                          ! => on notifie le predecesseur
                          call arret%set_error(-9)
                          call arret%set_time(t)
                          call arret%set_value(sqrt(d))
                          call arret%set_body(p)
                          buffer%m_mce(k) = arret
                         endif
                     enddo
                 endif
                enddo
               enddo
            
            end subroutine ctrl_distpartpla_processbuffer

!***********************************************************************
!> @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
            subroutine ctrl_distpartpla_onbufferfull(this, userdata,           &
     &              poswrite)
             implicit none
             class(t_buffer), intent(inout) :: this !< buffer de sortie
             class(*), intent(inout) :: userdata    !< de type t_ctrl_distpartpla
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)

             select type(userdata)
              class is(t_ctrl_distpartpla)
               call ctrl_distpartpla_processbuffer(userdata, this,          &
     &               poswrite)

              class default
               stop 'ctrl_distpartpla_onbufferfull : bad class'
             end select 
             
             call this%empty()
            
            end subroutine ctrl_distpartpla_onbufferfull
 
      end module mod_ctrl_distpartpla

!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ctrl_diststar.f 
!!  \brief Controle de la distance des corps a l'etoile 
!!
!!    --------------   ctrl_diststar    ---------------------\n
!!    | buffer src |  ----------------> | arret+ buffer dst |\n
!!    --------------                    ---------------------\n
!!
!!   on envoie les donnees uniquement en cas d'erreur au buffer de sortie.
!!   le buffer de sortie est optionnel. 
!!
! history : creation 09/04/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le controle de la distance des corps a l'etoile
!***********************************************************************
      module mod_ctrl_diststar
       use mod_ctrl_base
       

!***********************************************************************
!> @class t_ctrl_diststar
!! Controle de la distance des planetes a l'etoile.\n
!! Notifie le buffer predecesseur si la planete est trop proche ou trop loin.\n
!! Cela provoque une erreur de code -6 ou -7.\n
!! Le buffer d'entree doit contenir les positions heliocentriques  + le temps.\n
!! les vitesses n'ont pas d'importance, car elles sont ignorees.\n
!!  
!***********************************************************************
      type, extends(t_ctrl_base) :: t_ctrl_diststar
          private
          
         real(TREAL) ::   m_distmin !< distance minimale a l'etoile. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
         real(TREAL) ::   m_distmax !< distance maximale a l'etoile. si m_distmax =-1, pas de controle, si m_distmax>=0, genere une erreur si distance>m_distmax
         integer     ::   m_plan_nb !< nombre de planetes a verifier (identique au nombre de planetes dans le buffer)

      contains
          procedure :: set_minmax => ctrl_diststar_set_minmax ! definit la distance minimale et maximale a l'etoile
          procedure :: processinputbuffer =>                            &
     &           ctrl_diststar_processbuffer ! procedure appellee lorsque le buffer d'entree est plein
          
      end type t_ctrl_diststar  

      contains
           

!***********************************************************************
!> @brief definit la distance minimale ou maximale des planetes a l'etoile
!***********************************************************************
            subroutine ctrl_diststar_set_minmax(this,nbplan,mins,maxs)
             implicit none
             class(t_ctrl_diststar), intent(inout):: this  !< dummy argument
             integer, intent(in) :: nbplan  !<nombre de planetes a verifier (identique au nombre de planetes dans le buffer)
             real(TREAL), intent(in) :: mins  !< distance minimale a l'etoile. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
             real(TREAL), intent(in) :: maxs  !< distance maximale a l'etoile. si m_distmax =-1, pas de controle, si m_distmax>=0, genere une erreur si distance>m_distmax
             procedure(t_procbufferfull), pointer :: pfcnfull
             
             this%m_plan_nb = nbplan
             this%m_distmin = mins
             this%m_distmax = maxs
             
             call this%init(1+6*nbplan, "CTRL diststar min/max")
             
            end  subroutine ctrl_diststar_set_minmax


!***********************************************************************
!> @brief verifie les distances dans le buffer, retourne
!!  l'indice de la derniere erreur
!***********************************************************************
           function ctrl_diststar_processbuffer(this,buffer,poswrite)
             implicit none
             class(t_ctrl_diststar), intent(inout) :: this    !< de type t_ctrl_diststar
             class(t_buffer), intent(inout) :: buffer !< buffer de sortie
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: ctrl_diststar_processbuffer
             integer(8) :: j
             integer k
             real(TREAL), dimension(1+6*this%m_plan_nb) :: R
             real(TREAL), dimension(3,this%m_plan_nb) :: Ph
             real(TREAL)  t, dmin2, dmax2, d
             
             integer nplan
             type(t_arret) :: arret
             
               ctrl_diststar_processbuffer = 0
               nplan = this%m_plan_nb
               dmin2 = this%m_distmin**2
               dmax2 = this%m_distmax**2
               do j=1, poswrite
                call buffer%readdata(j, R)
                t = R(1)
                Ph = reshape(R(2:1+3*nplan), shape(Ph))
                do k=1, nplan
                 d=dot_product(Ph(:,k),Ph(:,k))
                 if (d.lt.dmin2) then
                  ! valeur trop petite de la distance  => on notifie le predecesseur
                  call arret%set_error(-6)
                  call arret%set_time(t)
                  call arret%set_value(sqrt(d))
                  call arret%set_body(k)
                  call buffer%notify_successor_error(arret)
                  ctrl_diststar_processbuffer = j
                 else if (d.gt.dmax2) then
                  ! valeur trop grande de la distance  => on notifie le predecesseur
                  call arret%set_error(-7)
                  call arret%set_time(t)
                  call arret%set_value(sqrt(d))
                  call arret%set_body(k)
                  call buffer%notify_successor_error(arret)
                  ctrl_diststar_processbuffer = j
                 endif
                enddo
               enddo
            
            end function ctrl_diststar_processbuffer

 

      end module mod_ctrl_diststar

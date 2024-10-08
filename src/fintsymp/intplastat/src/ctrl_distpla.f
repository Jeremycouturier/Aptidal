!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ctrl_distpla.f 
!!  \brief Controle de la distance des corps entre eux (hors etoile)
!!
!!    --------------   ctrl_distpla     --------------\n
!!    | buffer src |  ----------------> | buffer dst |\n
!!    --------------                    --------------\n
!!
!!   on envoie les donnees uniquement en cas d'erreur au buffer de sortie.
!!   le buffer de sortie est optionnel. 
!!
! history : creation 24/01/2018
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le controle de la distance des corps entre eux (hors etoile)
!***********************************************************************
      module mod_ctrl_distpla
       use mod_ctrl_base
       
!***********************************************************************
!> @class t_ctrl_distpla
!! Controle de la distance des planetes entre elles.\n
!! Chaque planete controle si une autre planete n'est pas proche trop d'elle.\n
!! Notifie le buffer predecesseur si une autre planete est trop proche.\n
!! Cela provoque une erreur de code -8.\n
!! Le buffer d'entree doit contenir les positions heliocentriques  + le temps.\n
!! les vitesses n'ont pas d'importance, car elles sont ignorees.\n
!! le buffer de sortie est optionnel. \n
!!  
!***********************************************************************
      type, extends(t_ctrl_base) :: t_ctrl_distpla
          private
          
         real(TREAL), dimension(:), allocatable ::   m_distmin !< tableau de m_plan_nb distance minimale. si m_distmin(i) =-1, pas de controle pour la planete i, si m_distmin(i)>=0, genere une erreur si distance d'une autre planete<m_distmin
         integer     ::   m_plan_nb !< nombre de planetes a verifier (identique au nombre de planetes dans le buffer)

      contains
          procedure :: set_min => ctrl_distpla_set_min ! definit la distance minimale a chaque planete
          procedure :: processinputbuffer =>                            &
     &           ctrl_distpla_processbuffer ! procedure appellee lorsque le buffer d'entree est plein
          
      end type t_ctrl_distpla  

      contains
           

!***********************************************************************
!> @brief definit la distance minimale des planetes entre elles
!***********************************************************************
            subroutine ctrl_distpla_set_min(this,nbplan,mins)
             implicit none
             class(t_ctrl_distpla), intent(inout):: this  !< dummy argument
             integer, intent(in) :: nbplan  !<nombre de planetes a verifier (identique au nombre de planetes dans le buffer)
             real(TREAL), dimension(:), intent(in) :: mins  !< tableau de "nbplan" distances minimales a l'etoile. si m_distmin(k) =-1, pas de controle pour la planete i, si m_distmin(i)>=0, genere une erreur si distance<m_distmin(i)
             
             this%m_plan_nb = nbplan
             allocate(this%m_distmin(1:nbplan))
             this%m_distmin(1:nbplan) = mins(1:nbplan)
             
             call this%init(1+6*nbplan, "CTRL distpla min")
             
            end  subroutine ctrl_distpla_set_min

!***********************************************************************
!> @brief verifie les distances dans le buffer
!***********************************************************************
           function ctrl_distpla_processbuffer(this,buffer,poswrite)
             implicit none
             class(t_ctrl_distpla), intent(inout) :: this    !< de type t_ctrl_distpla
             class(t_buffer), intent(inout) :: buffer !< buffer de sortie
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: ctrl_distpla_processbuffer
             integer(8) :: j
             integer k, p
             real(TREAL), dimension(1+6*this%m_plan_nb) :: R
             real(TREAL), dimension(3,this%m_plan_nb) :: Ph
             real(TREAL)  t,  d, dmin2
             
             
             integer nplan
             type(t_arret) :: arret
             
               nplan = this%m_plan_nb
               ctrl_distpla_processbuffer = 0
               
               do j=1, poswrite
                call buffer%readdata(j, R)
                t = R(1)
                Ph = reshape(R(2:1+3*nplan), shape(Ph))
                do k=1, nplan
                  dmin2 = this%m_distmin(k)**2
                  do p=1, nplan
                  if (p.ne.k) then 
                     d=dot_product(Ph(:,k)-Ph(:,p),Ph(:,k)-Ph(:,p))
                     ! write(*,*) sqrt(dmin2), sqrt(d) , k, p
                     if (d.lt.dmin2) then
                      ! write(*,*) 'arret' dmin2, d , k, p
                      ! valeur trop petite de la distance  => on notifie le predecesseur
                      call arret%set_error(-9)
                      call arret%set_time(t)
                      call arret%set_value(p*1.D0)
                      call arret%set_body(k)
                      call buffer%notify_successor_error(arret)
                      ctrl_distpla_processbuffer = j
                     endif
                  endif
                 enddo
                enddo
               enddo
            
            end function ctrl_distpla_processbuffer

      end module mod_ctrl_distpla

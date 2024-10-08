!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ode.f 
!!  \brief type de base dont doivent derive les systemes integree 
!!   par un integrateur de type dopri ou odex
!!   les types derivant de celui-ci doivent surcharge la fonction fcn
!!
!
! history : creation 04/12/2015
!***********************************************************************
!***********************************************************************
! module pour une integration de type ODE
!***********************************************************************
      module mod_ode
       
!***********************************************************************
!> @class t_schemaode
!! classe de base decrivant 1 schema d'integration de type ODE
!!
!***********************************************************************
      type, abstract :: t_schemaode
       real(TREAL) :: m_eps      !< precision interne de l'integrateur
       real(TREAL) :: m_hmax     !< pas maximal de l'integrateur
       real(TREAL) :: m_h        !< pas de l'integrateur
       
       contains         
          procedure :: set_prec_step => schemaode_set_prec_step ! fixe la precision et les pas internes pour l'integrateur
          procedure(procdepart), deferred :: start  !< fonction appellee lors du demarrage (doit appeller la fonction t_sysode%fill pour remplir le vecteur d'etat)
          procedure(procrun), deferred :: run  !< fonction appellee lors d'un pas normal
          procedure :: check_error => schemaode_check_error  ! controle si une erreur s'est produite
      end type t_schemaode  

!***********************************************************************
!> @class t_sysode
!! classe de base decrivant 1 systeme integre avec un integrateur d'ODE de type dopri ou odex
!!  
!***********************************************************************
      type, abstract :: t_sysode
       contains
          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          procedure(t_odeallocatey), deferred :: allocate_y         
          procedure(t_odefilltoy), deferred :: fill_to_y         
          procedure(t_odefillfromy), deferred :: fill_from_y         
          procedure(t_odefcn), deferred ::  fcn        
      end type t_sysode  
      
           abstract interface
!***********************************************************************
! type de fonction pour allouer le vecteur d'etat du second membre
!***********************************************************************
            subroutine t_odeallocatey(this, y)
             import t_sysode
             class(t_sysode), intent(inout) :: this
             real(TREAL), dimension(:), allocatable, intent(inout) :: y   !< vecteur d'etat au temps t
            end subroutine t_odeallocatey

!***********************************************************************
! type de fonction pour remplir le vecteur d'etat du second membre
!***********************************************************************
            subroutine t_odefilltoy(this, y, t)
             import t_sysode
             class(t_sysode), intent(inout) :: this
             real(TREAL), intent(in) :: t                 !< temps ou est appelle t_odefill
             real(TREAL), dimension(:), intent(inout) :: y   !< vecteur d'etat au temps t
            end subroutine t_odefilltoy
            
!***********************************************************************
! type de fonction pour remplir this a partir du vecteur d'etat du second membre
!***********************************************************************
            subroutine t_odefillfromy(this, y)
             import t_sysode
             class(t_sysode), intent(inout) :: this
             real(TREAL), dimension(:), intent(in) :: y   !< vecteur d'etat au temps t
            end subroutine t_odefillfromy
            
!***********************************************************************
! type de fonction pour calculer le second membre de l'ODE
!***********************************************************************
            subroutine t_odefcn(this, n, t, y, dy)
             import t_sysode
             class(t_sysode), intent(inout) :: this
             integer, intent(in) :: n                     !< nombre decomposante de y
             real(TREAL), intent(in) :: t                 !< temps ou est appelle t_odefcn
             real(TREAL), dimension(:), intent(in) :: y   !< vecteur d'etat au temps t
             real(TREAL), dimension(:), intent(out) :: dy !< derivee du vecteur d'etat calcule au temps t
            end subroutine t_odefcn
            
!***********************************************************************
! type de fonction pour executer un pas normal d'integration avec l'integrateur
!***********************************************************************
            subroutine procrun(this, system, tstart, tend, dt)
            import t_schemaode
            import t_sysode
            class(t_schemaode), intent(inout) :: this
            class(t_sysode), intent(inout) :: system ! systeme qui est integre (doit fournir: fill et fcn)
            real(TREAL), intent(in) :: dt ! pas d'integration
            real(TREAL), intent(in) :: tstart ! temps initial
            real(TREAL), intent(in) :: tend ! temps final = tstart+dt
            end subroutine procrun

!***********************************************************************
! type de fonction pour demarrer une integration  au temps tstart  avec l'integrateur
!***********************************************************************
            subroutine procdepart(this, system, tstart, dt)
            import t_schemaode
            import t_sysode
            class(t_schemaode), intent(inout) :: this
            class(t_sysode), intent(inout) :: system ! systeme qui est integre (doit fournir: fill et  fcn)
            real(TREAL), intent(in) :: dt ! pas d'integration
            real(TREAL), intent(in) :: tstart ! temps initial
            end subroutine procdepart

          end interface
          
          
       contains   
          
!***********************************************************************
!> @brief fixe la precision et les pas internes
!***********************************************************************
      subroutine schemaode_set_prec_step(this, eps, hmax, h)
       implicit none
       class(t_schemaode), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: eps      !< precision interne de l'integrateur
       real(TREAL), intent(in) :: hmax     !< pas maximal de l'integrateur
       real(TREAL), intent(in) :: h        !< pas de l'integrateur
        
       this%m_eps = eps
       this%m_hmax = hmax
       this%m_h = h       
      end  subroutine schemaode_set_prec_step  

!***********************************************************************
!> @brief controle si une erreur s'est produite. 
!! retourne true en cas d'erreur
!***********************************************************************
      function schemaode_check_error(this) result(ret)
        implicit none
        class(t_schemaode), intent(inout) :: this  !< dummy argument
        logical ret ! code de retour
        ret = .false.
      end function schemaode_check_error

      end  module mod_ode
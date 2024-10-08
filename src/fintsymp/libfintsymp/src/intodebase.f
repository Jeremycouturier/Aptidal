!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file intsympbase.f 
!!  \brief integrateur d'un systeme derivant de t_syssymp   
!!         pouvant etre integre par un schema ODE
!!
!!     cette classe t_intodebase fournit les routines elementaires communes a l'ensemble 
!!      des systemes possibles
!
! history : creation 03/12/2015
!***********************************************************************
#include "realc.f"

!***********************************************************************
!> module pour l'integration ODE 
!***********************************************************************
      module mod_intodebase
      use mod_intsympbase
      use mod_ode
      use mod_sysplaode
      

!***********************************************************************
!> @class t_intodebase
!! classe de base decrivant l'integrateur ODE d'un systeme derivant de t_sysode
!!  
!!  
!***********************************************************************
      type, extends(t_intsympbase)  :: t_intodebase 
                  
          class(t_schemaode), pointer, private :: m_schemeode !< schema d'integration
          class (t_sysplaode), pointer     :: m_sysode  !< systeme planetaire integre

       contains
          private
                    
          procedure, public:: set_sysode => odeb_set_sys ! fixe le systeme a integrer
          procedure,public:: set_integ_scheme => odeb_set_integ_scheme ! fixe le type d'integrateur
          procedure,public:: set_integ_prec_step => odeb_set_prec_step ! fixe la precision et les pas internes pour l'integrateur
          procedure, public:: integ_run => odeb_integ_run  ! integre pendant niter iterations
                   
      end type t_intodebase  

     
      contains


!***********************************************************************
!> @brief fixe le systeme planetaire a integrer
!***********************************************************************
      subroutine odeb_set_sys(this, sysode)
       implicit none
       class(t_intodebase), intent(inout):: this  !< dummy argument
       class(t_sysplaode), target, intent(in) :: sysode  !< systeme planetaire
       
       this%m_sysode => sysode
      end  subroutine odeb_set_sys  
     

!***********************************************************************
!> @brief fixe le schema d'integration (e.g. 'DOPRI' ou 'ODEX') 
!***********************************************************************
      subroutine odeb_set_integ_scheme(this, schema)
       use mod_dopri
       use mod_odex1
       use mod_odex2
       implicit none
       class(t_intodebase), intent(inout):: this  !< dummy argument
       character(len=*), intent(in) :: schema      !< schema d'integration (e.g., 'DOPRI' ou 'ODEX' )
#if TREAL==8
       type(t_dopri),pointer  :: sysdopri
       type(t_odex1),pointer  :: sysodex1
       type(t_odex2),pointer  :: sysodex2
        
      ! cas de l'integrateur DOPRI
      if (schema.eq.'DOPRI') then
        allocate(sysdopri)
        this%m_schemeode => sysdopri
      ! cas de l'integrateur ODEX1
      else if (schema.eq.'ODEX1') then
        allocate(sysodex1)
        this%m_schemeode => sysodex1
      ! cas de l'integrateur ODEX2
      else if (schema.eq.'ODEX2') then
        allocate(sysodex2)
        this%m_schemeode => sysodex2
      else
#endif
        write(*,*) 'integrateur=', schema
        write(*,*) 'integrateur non supporte par odeb_set_integ_scheme'
        stop
#if TREAL==8
      endif  
#endif              
      end  subroutine odeb_set_integ_scheme  
      
!***********************************************************************
!> @brief fixe la precision et les pas internes
!***********************************************************************
      subroutine odeb_set_prec_step(this, eps, hmax, h)
       implicit none
       class(t_intodebase), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: eps      !< precision interne de l'integrateur
       real(TREAL), intent(in) :: hmax     !< pas maximal de l'integrateur
       real(TREAL), intent(in) :: h        !< pas de l'integrateur
        
       call this%m_schemeode%set_prec_step(eps, hmax, h)
       
      end  subroutine odeb_set_prec_step  
      
      
!***********************************************************************
!> @brief controle des erreurs pendant l'integration
!! @return retourne 0  si aucune action a faire \n
!! retourne 1 si on doit s'arreter d'integrer (mais poursuite de l'execution)
!***********************************************************************
      function odeb_integ_check_error(this, t) result(ret)
       implicit none
       class(t_intodebase), intent(inout):: this  !< dummy argument
       real(TREAL), intent (in) :: t !< temps actuel
       integer ret ! code de retour
       ret = 0
       
         if (this%m_sysode%check_error().or.                              &
     &    this%m_schemeode%check_error()) then
           select case (this%m_errorbehavior)
           case  (1)  
             write (*,*) " arret de l'integration ODE a t=", t
             ret = this%m_sysode%set_error_time(t)
             
           case  default  
             write (*,*) " arret de l'integration ODE a  t=", t
             write (*,*) " arret du programme !!"
             stop
           end select  
         endif
      end function  odeb_integ_check_error  
         

!***********************************************************************
!> @brief execute niter iterations du schema d'integration
!***********************************************************************
      subroutine odeb_integ_run(this, niter)
       implicit none
       class(t_intodebase), intent(inout):: this   !< dummy argument
       integer(8), intent(in) :: niter             !< nombre d'iterations
       integer(8) i, i0
       integer(8) minmodulo
       integer nstepscall,j
       real(TREAL) :: dt, t,ts
       integer(8),dimension(:), allocatable :: arstepscall
       logical isoutputstep
       integer code_err
       minmodulo = 0
       i = 1
       dt = this%m_dt
       
       call this%m_sysode%get_outputsteps(nstepscall,arstepscall)
       
       ! sortie dans le cas du demarrage standard
       i0 = 0
       call this%m_sysode%write_output(i0, this%m_t0)            
    
       
       ts = this%m_t0+this%m_it*dt
       call this%m_schemeode%start(this%m_sysode,ts,dt)

       do while (i.le.niter)
         ts = this%m_t0 + this%m_it*dt
         this%m_it = this%m_it+1
         t = this%m_t0 + this%m_it*dt
         call this%m_schemeode%run(this%m_sysode,ts,t,dt)
         ! controle de l'erreur
         code_err = odeb_integ_check_error(this, t)
         if (code_err.eq.1) then
          exit ! sortie de boucle 
         endif
         
         ! verifie si on est un pas de sortie
         isoutputstep = .false.
         do j=1, nstepscall
          if (mod(i,arstepscall(j)).eq.0) then
           isoutputstep = .true.
          endif
         enddo 
         ! sortie des donnees
         if (isoutputstep.eqv..true.) then
            !-- utilisation du systeme de sortie
            call this%m_sysode%write_output(this%m_it, t)     
            
            ! controle de l'erreur dans les sorties 
            code_err = odeb_integ_check_error(this,t)
            if (code_err.eq.1) then
             exit ! sortie de boucle 
            endif
         end if
         i = i + 1
       end do
       
       !-- flush des sorties 
       call this%m_sysode%flush_output()
       
       if (allocated(arstepscall)) then
          deallocate(arstepscall)
       endif
      end  subroutine odeb_integ_run  


      end module mod_intodebase

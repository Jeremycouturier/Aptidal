!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file sysplanewth.f 
!!  \brief systeme avec N planetes
!!         pouvant etre integre par un schema symplectique.
!!         les coordonnees des corps sont exprimees en heliocentriques (position/vitesse).
!!         seules les interactions newtonniennes sont prises en compte.
!!
!!      Le schema est adaptatif et suit la methode du Hamiltonien logarithmique (Preto and Tremaine 1999)
!!      On herite de la classe t_sysplanewtH
!!
!
! history : creation 20/07/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec N planetes 
! exprime en variables heliocentriques
! seules les interactions newtonniennes sont prises en compte
!***********************************************************************
      module mod_sysplanewtH_adapt
       use mod_sysplanewtH
        use kahansum
      
       integer, parameter :: NDATAH_SCH = 5 ! nombre de donnees dans le fichier de sortie (cf. add_adapt_out_step  et sysplanewtJ_adapt_write_output)
                                           ! 4 = rtime + temps(s) + E0 + E1 

!***********************************************************************
!> @class t_sysplanewtH
!! classe decrivant 1 systeme planetaire avec N planetes
!! exprimee en variables heliocentriques.
!! seules les interactions newtonniennes sont prises en compte.
!!  
!!  
!!  
!***********************************************************************
      type, extends(t_sysplanewtH) :: t_sysplanewtH_adapt 

          private 
          type(compensatedsum) :: rtime !< real time that is integrated by the system at the same time the Keplerian steps are done.
          type(compensatedsum) :: E0    !< Inital energy of the system, it is kept constant during the time of integration
          real(TREAL) :: E1     !< Threshold used in order to trigger the renormalization

       contains
          private
          
       procedure, public::set_plan_mass=>sysplanewtH_adapt_set_plan_mass ! fixe la masse des planetes 
       procedure, public :: set_plan_coord => sysplanewtH_set_plan_coord ! fixe les coordonnees des planetes (meme que celle integree)
          
          procedure, public :: pasA =>sysplanewtH_adapt_pasA
          procedure, public :: pasB =>sysplanewtH_adapt_pasB
          procedure, public :: pasC =>sysplanewtH_adapt_pasC
          
          procedure, public :: energie =>sysplanewtH_adapt_gamma_energie
          procedure, public :: energie_with_err   =>                    &
     &                         sysplanewtH_adapt_gamma_energie_with_err

          procedure, public :: create_output  =>                        &
     &             sysplanewtH_adapt_create_output! creation d'une copie pour effectuer les pas de sortie
          procedure, public :: write_output   =>                               &
     &             sysplanewtH_adapt_write_output
          procedure, public:: copy_output     =>                               &
     &             sysplanewtH_adapt_copy_output ! copie de this dans la copie de sortie

          procedure, public :: add_adapt_out_step                                &
     &             => sysplanewtH_adapt_add_adapt_out_step   ! fixe le pas de sortie, la taille du buffer des donnees suplementaires du schema
          
          procedure, public :: get_adaptativetime                                &
     &             => sysplanewtH_adapt_get_adaptativetime ! retourne le temps dans le cas d'un schemea adaptatif

          procedure, public :: dumpbin                                            &
     &             => sysplanewtH_adapt_dumpbin ! dump binaire vers un fichier
          procedure, public :: restorebin                                         &
     &             => sysplanewtH_adapt_restorebin ! restauration binaire depuis un fichier

      end type t_sysplanewtH_adapt
      
      contains
         
!***********************************************************************
!> @brief fixe la masse des planetes 
!!  la fonction set_star_mass doit etre appellee auparavant
!***********************************************************************
      subroutine sysplanewtH_adapt_set_plan_mass(this, m)
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:), intent(in) :: m  !< masse des planetes
       real(TREAL) :: m0,sum_mpl,E1
       integer :: i,j
       
       this%sys_pas_fixe = .false. ! indique que le systeme est a pas adapatif

       call sysplanewtH_set_plan_mass(this, m)
       call this%set_graphdotname("SYSPLA Helio adaptatif")
       sum_mpl = REALC(0.E0)
       E1 = REALC(0.E0)
       m0 = this%m_star_mass

       sum_mpl = 0.
       E1 = 0.
       do i=1,this%m_plan_nb

          sum_mpl = sum_mpl+m(i)
          do j=i+1,this%m_plan_nb
             E1 = E1 + m(i)*m(j)
          end do
       end do
       
       this%E1 =  2*E1/(sum_mpl*(m0+sum_mpl))
       
      end  subroutine sysplanewtH_adapt_set_plan_mass  

!***********************************************************************
!> @brief fixe les coordonnees des planetes. 
!! Il n'y a aucune conversion des coordonnees.
!! Ces coordoonnees doivent etre compatibles avec l'integrateur
!! calcule l'energie E0
!***********************************************************************
      subroutine sysplanewtH_set_plan_coord(this, X, XP)
       use kahansum
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:,:), intent(in) :: X  !< position des planetes
       real(TREAL), dimension(:,:), intent(in) :: XP  !< moment des planetes
       type(compensatedsum) :: E0
       
       call syspla_set_plan_coord(this, X, XP)
       
       call  sysplanewtH_energie_with_err(this, E0)
       this%E0 = E0
       this%E1 = -this%E1*E0%v
       !Setting the initial energy to zero
       this%m_Hinit%v = 0.
       this%m_Hinit%err = 0.
      end  subroutine sysplanewtH_set_plan_coord  


!***********************************************************************
!> @brief creation d'une copie pour effectuer les pas de sortie
!! en general, duplique simplement this
!***********************************************************************
      subroutine sysplanewtH_adapt_create_output(this, sysout)
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       class(t_syssymp), allocatable, intent(out) :: sysout  !< systeme cree pour la sortie
       
       call syspla_create_output(this, sysout)
             
      end  subroutine sysplanewtH_adapt_create_output  


!***********************************************************************
!> @brief ecriture de la sortie (appelee par l'integrateur)
!!   a l'iteration donne et au temps donne
!***********************************************************************
      subroutine sysplanewtH_adapt_write_output(this, it, t)
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: it !< iteration actuelle
       real(TREAL), intent(in) :: t !< temps correspondant a it
       real(TREAL) :: rtime
       type(compensatedsum) :: H,HmE0
       real(TREAL), dimension(1:NDATAH_SCH) :: schema_dataout
       integer j
       
       rtime = this%rtime%v
       call  sysplanewtH_energie_with_err(this, H)
       HmE0 = H-this%E0
       if (allocated(this%m_stepout)) then
       do j=1, size(this%m_stepout)
        if (this%isvalidoutput(it,j)) then 
               
         ! cas d'un buffer du schema adaptatif
         if (this%m_kind(j).eq.3) then
          schema_dataout(1) = rtime
          schema_dataout(2) = t
          schema_dataout(3) = this%E0%v
          schema_dataout(4) = this%E1
          schema_dataout(5) = HmE0%v
          ! attention, il faut remplir autant d'element que NDATAH_SCH
          call this%m_buffer(j)%writedata(schema_dataout)
         endif
         
        endif
       enddo
       endif

       ! appel pour les autres buffers
       call syspla_write_output(this, it, rtime)

       end  subroutine sysplanewtH_adapt_write_output


!***********************************************************************
!> @brief copie de this dans la copie de sortie (sysout)
!! en general, sysout = this 
!***********************************************************************
      subroutine sysplanewtH_adapt_copy_output(this, sysout)
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       class(t_syssymp), intent(inout) :: sysout  !< systeme cree pour la sortie

       select type(sysout)
       class is(t_sysplanewtH_adapt)
         
         call syspla_copy_output(this, sysout)

         sysout%rtime = this%rtime
         sysout%E0 = this%E0
         sysout%E1 = this%E1
         
        class default
         stop 'Erreur: this ne derive pas de syspla_copy_output' 
      end select
      end  subroutine sysplanewtH_adapt_copy_output

!***********************************************************************
!> @brief fixe le pas de sortie, la taille du buffer des informations du schema adaptatif
!! cela cree un nouveau buffer
!! lorsque le buffer est plein, le consommateur buffercs est appelle
!***********************************************************************
      subroutine sysplanewtH_adapt_add_adapt_out_step(this, stepout,                &
     &                  buffersize, buffercs)
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les inetgrales premieres
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
               
       call this%add_out_step(stepout, buffersize,                            &
     &             buffercs, NDATAH_SCH , 3) 
       
      end  subroutine sysplanewtH_adapt_add_adapt_out_step  
  
!***********************************************************************
!> @brief execute le pas A de l'integrateur
!***********************************************************************
      subroutine sysplanewtH_adapt_pasA(this, mdt)
       use kahansum
       use mod_kepsaut
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= tau
       
       real(TREAL) :: tau
       type(compensatedsum) :: H0, E0, diffH0mE0
       real(TREAL) ::  f_renorm
       
       E0 = this%E0
       
       call sysplanewtH_keplerian_energy(this, H0)
       diffH0mE0 = H0-E0    
       f_renorm = (diffH0mE0%v/this%E1)
       !if (f_renorm.lt.0.) then
       !  write(*,*) "sysplanewtJ_pasA", f_renorm
       !end if
         
       tau = (mdt)/sqrt(1.+(f_renorm)**2)
      
       this%rtime = this%rtime + tau
       
       call sysplanewtH_pasA(this, tau)

      end  subroutine sysplanewtH_adapt_pasA  


!***********************************************************************
!> @brief execute le pas B de l'integrateur = B1(tau)B2(tau)B1(tau) 
!***********************************************************************
      subroutine sysplanewtH_adapt_pasB(this, mdt)
       use kahansum
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       type(compensatedsum) :: H1
       real(TREAL) ::  f_renorm, tau

       call sysplanewtH_perturbation_energy(this, H1)
       f_renorm = ((-H1%v)/this%E1)
      
       tau = (mdt)/sqrt(1.+(f_renorm)**2)
       
       call sysplanewtH_pasB(this, tau)

      end  subroutine sysplanewtH_adapt_pasB  


!***********************************************************************
!> @brief execute le pas C de l'integrateur
!***********************************************************************
      subroutine sysplanewtH_adapt_pasC(this, mdt)
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       
       write(*,*) "sysplanewtH_adapt_pasC", mdt
       write(*,*) "pas de correcteur en heliocentrique"
       stop
       
      end  subroutine sysplanewtH_adapt_pasC  


      
!***********************************************************************
!> @brief calcule l'energie de this
!***********************************************************************
      subroutine sysplanewtH_adapt_gamma_energie_with_err(this, H)
      implicit none
      
       class(t_sysplanewtH_adapt), intent(inout):: this !< dummy argument
       type(compensatedsum), intent(out) :: H!< energie calculee
       real(TREAL) :: E1
       real(TREAL) :: fprime
       type(compensatedsum) :: H0, H1, E0 , DeltaE
       
       call sysplanewtH_keplerian_energy(this,H0)
       call sysplanewtH_perturbation_energy(this, H1)

       E0 = this%E0
       E1 = this%E1
       DeltaE = H0+H1-E0
       fprime = 1./sqrt(1.+((H1%v/E1)**2)) 
       H = DeltaE + H1%v/2*(fprime*DeltaE%v/E1)**2
       H = fprime*H
       !Too small so we expand to the second order
      ! H = log((diffH0mE0+(diffH0mE0**2+1)**.5)/                        &
  !   &             (-H1+(H1**2+1)**.5))
       !write(*,*) "sysplanewtJ_energie ",H,log((H0-E0)/(-H1)),(H0-E0)+H1
       ! write(*,*) "sysplanewtJ_energie offset ",2*abs(H1)
      end  subroutine  sysplanewtH_adapt_gamma_energie_with_err


   !***********************************************************************
!> @brief calcule l'energie de this
!***********************************************************************
      subroutine sysplanewtH_adapt_gamma_energie(this, H)
      implicit none
      
       class(t_sysplanewtH_adapt), intent(inout):: this !< dummy argument
       real(TREAL), intent(out) :: H !< energie calculee
       type(compensatedsum) :: Gamma_comp !< energie calculee

        call sysplanewtH_adapt_gamma_energie_with_err(this,Gamma_comp)
        H = Gamma_comp%v
        
      end  subroutine  sysplanewtH_adapt_gamma_energie   
  
!***********************************************************************
!> @brief retourne le temps du systeme dans le cas d'un schemea adaptatif
!***********************************************************************
      function sysplanewtH_adapt_get_adaptativetime(this) result(t)
       implicit none
       class(t_sysplanewtH_adapt), intent(in):: this  !< dummy argument
       
       real(TREAL) ::t ! temps du systeme
       t  = this%rtime%v
         
      end function  sysplanewtH_adapt_get_adaptativetime
        
!***********************************************************************
!> @brief dump binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (in) objet a sauvegarder
!***********************************************************************
      subroutine sysplanewtH_adapt_dumpbin(this, file) 
       use mod_io_dump
       implicit none
       class(t_sysplanewtH_adapt), intent(in):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       call sysplanewtH_dumpbin(this, file)

       write(file%m_nf) this%rtime 
       write(file%m_nf) this%E0 
       write(file%m_nf) this%E1 
       
      end subroutine  sysplanewtH_adapt_dumpbin  

!***********************************************************************
!> @brief restauration binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (inout) objet a sauvegarder
!***********************************************************************
      subroutine sysplanewtH_adapt_restorebin(this, file) 
       use mod_io_dump
       implicit none
       class(t_sysplanewtH_adapt), intent(inout):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       call sysplanewtH_restorebin(this, file)

       read(file%m_nf) this%rtime 
       read(file%m_nf) this%E0 
       read(file%m_nf) this%E1 
 
       end subroutine  sysplanewtH_adapt_restorebin  

      end module mod_sysplanewtH_adapt

!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file intgodex2016.f 
!!  \brief  routine d'integration : lance l'integration 
!!   et ecrit dans les fichiers de sortie
!!
! history : creation 11/01/2016
!***********************************************************************
#include "realc.f"
!-----------------------------------------------------------
! En entree: tinit: temps initial
!            dt : intervalle de temps pour une iteration
!            n_iter : nombre d'iterations a effectuer
!            n_out : nombre d'iterations apres lequel on appelle ecrit
!-----------------------------------------------------------
      module mod_intg_ode
      contains
      


!***********************************************************************
!> fonction appellee par le distributeur de condition initiale pour la condition ligneci
!***********************************************************************
      subroutine cidist_integfunc(ligneci, filepar)
        use mod_io_int
        use mod_io_txt_ncomp
        use mod_io_txt_car
        use mod_sysplaodehelio
        use mod_converter
        use mod_ci_planewt
        use mod_filepar
        use mod_gmsun
        implicit none
        
        class(t_ci_base), intent(inout) :: ligneci !< condition initiale recue
        class(*), intent(inout) :: filepar !< donnee utilisateur : fichier de parametres
        
        real(TREAL), dimension(:), allocatable :: mpl ! masse normalise (masse de l'etoile a 1)
        real(TREAL), dimension(:,:), allocatable :: CI
        real(TREAL),allocatable,dimension(:,:) :: Xplan,XPplan
        real(TREAL),allocatable,dimension(:,:) :: Ph,Vh, Vb
        real(TREAL),allocatable,dimension(:,:) :: Phinv, Vbinv, Vhinv
        type(t_io_int) ::io_int
        type(t_io_txt_car) ::io_cart
        type(t_io_txt_ncomp) ::io_ell
        type(t_sysplaodehelio),allocatable, target  :: sysH ! systeme en helio
        class(t_sysplaodehelio), pointer  :: psys ! systeme en jacobi ou helio
        character(100)  :: int_type ! type de l'integrateur
        integer :: nplan ! nombre de planetes
        integer k
        integer ci_type ! type des coordonnnees des conditions initiales
        real(TREAL) :: cG
        
! --- initialisation des constantes      
        ! constante GM du soleil
        select type(filepar)
         class is(t_extfilepar)
            cG = calcGM_sun(filepar%ref_gmsun)            
         class default
          stop 'filepar in cidist_integfunc : bad class'
        end select 

        
        
        write(*,*) "start cidist_integfunc", trim(ligneci%m_id) 
        call ligneci%debug()

        select type(ligneci)
         class is(t_ci_planewt)
          nplan = ligneci%m_plan_nb
          ci_type = ligneci%m_ci_type
          allocate(mpl(0:nplan))
          allocate(CI(1:6,1:nplan))
          mpl(0:nplan) = ligneci%m_plan_mass(0:nplan)
! Masses normalisees a 1
          cG=cG*mpl(0)
          write(*,*) 'cG*m0=', cG
          do k=1,nplan
            mpl(k)=mpl(k)/mpl(0)
          enddo
          mpl(0)=1.D0
          ! ici le 1:nplan au lieu de 0:nplan car ci_type: heliocentrique
          CI(1:6,1:nplan) = ligneci%m_plan_coord_ci(1:6,1:nplan)
         
         class default
          stop 'ligneci in cidist_integfunc : bad class'
        end select 

        select type(filepar)
         class is(t_extfilepar)

      int_type = trim(filepar%int_type)
      
 ! --- create the filename of output files  
      call filepar%create_filename(ligneci%m_id)

      if (filepar%if_int.eq.1) then
         write(*,*)  ' nf_int ', trim(filepar%nf_int)
         call io_int%set_name(trim(filepar%nf_int),60)
         call io_int%fcreate()
      end if
      if (filepar%if_ell.eq.1) then
         write(*,*)  ' nf_ell ', trim(filepar%nf_ell)
         io_ell = t_io_txt_ncomp("I/O ell")
         call io_ell%set_comp_nb(1+6*nplan)
         call io_ell%set_name(trim(filepar%nf_ell),61)
         call io_ell%fcreate()
      end if
      if (filepar%if_car.eq.1) then
         write(*,*) ' nf_car ', trim(filepar%nf_car)
         io_cart = t_io_txt_car("I/O cart")
         call io_cart%set_comp_nb(1+6*nplan)
         call io_cart%set_name(trim(filepar%nf_car),62)
         call io_cart%set_plan_nb(nplan)
         call io_cart%fcreate()
      end if


! --- Creation des tableaux
      allocate(Xplan(3,nplan))
      allocate(XPplan(3,nplan)) 
      allocate(Ph(3,nplan))
      allocate(Vh(3,nplan)) 
      allocate(Vb(3,nplan)) 
      allocate(Phinv(3,nplan))
      allocate(Vhinv(3,nplan)) 
      allocate(Vbinv(3,nplan)) 
        write(*,*)'Heliocentric Coord'
        filepar%int_var = 1               ! variables heliocentriques canoniques
        allocate(sysH)
        psys=>sysH

      
      ! fixe le systeme a integrer
      call psys%set_plan_nb(nplan)
      call psys%set_star_mass(cG, mpl(0))
      call psys%set_plan_mass(mpl(1:))

      !
      ! changement de coordonnees depuis les CI lues
      !
      write(*,*) 'CI='
      write(*,*) CI
      

       call CI2PhVh(nplan,mpl,cG, ci_type,CI,Ph,Vh)
       if (filepar%if_invar.eq.1) then 
        ! passage vers invariant 
        call coord_vh2vb(nplan,mpl,Vh,Vb)
        call PhVb2inv(nplan,mpl,Ph,Vb, Phinv, Vbinv)
        call coord_vb2vh(nplan,mpl,Vbinv,Vhinv)
        Xplan = Phinv
        XPplan = Vhinv
       else
        ! reste dans le repere des ci
        Xplan = Ph
        XPplan = Vh
       endif


      write(*,*) ' X , V dans cidist_integfunc'
      write(*,*)  Xplan,XPplan

      ! -- Conditions initiales des planetes
         write(*,*) CI
         
         call psys%set_plan_coord(Xplan,XPplan)


         call intg_symp(filepar,psys, cG, io_int, io_cart, io_ell)

         write(*,*) 'fin d integration'

        write(*,*) "end cidist_integfunc"
        
       deallocate(sysH)
       deallocate(Xplan)
       deallocate(XPplan) 
       deallocate(Ph)
       deallocate(Vh) 
       deallocate(Vb) 
       deallocate(Phinv)
       deallocate(Vhinv) 
       deallocate(Vbinv) 
       deallocate(mpl)
       deallocate(CI)
       
         class default
          stop 'filepar in cidist_integfunc : bad class'
        end select 

       end subroutine cidist_integfunc
      
      
!------------------------------------------------------------
!> routine d'integration
!------------------------------------------------------------
      subroutine intg_symp(filepar, sys, cG, io_int,io_cart, io_ell)
      use mod_sysplaodehelio
      use mod_intodebase
      use mod_io_int
      use mod_converter
      use mod_io_txt_ncomp
      use mod_io_txt_car
      use mod_filepar
      
      implicit none
      type (t_extfilepar),target, intent(inout) :: filepar !< fichier de parametres + variables derivees
      type (t_io_int),intent(inout) ::io_int  !< sortie des integrales premieres
      type (t_io_txt_car), intent(inout) :: io_cart  !< sortie des elements cartesiens
      type (t_io_txt_ncomp), intent(inout) :: io_ell !< sortie des elements elliptiques
      class(t_sysplaodehelio), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
     
      integer(8) :: n_iter
      integer(8) :: n_out
      real(TREAL) :: tinit,dt
      
      type (t_plan_PhVh2ell) :: PhVh2ell

      type(t_intodebase) :: integrator
      integer(8) :: i
      character(100)  :: int_type ! type de l'integrateur
      integer graphnode ! pour la generation du graphe
     
      n_iter = filepar%n_iter
      n_out = filepar%n_out
      tinit = filepar%tinit
      dt = filepar%dt
      int_type = filepar%int_type
      
      i = 0
      
      write(*,*) 'int_type = ', int_type

      call sys%clear_out_step()
      
      ! sortie integrales premieres
      if (filepar%if_int.eq.1) then
        call sys%add_int_out_step(n_out,1_8,io_int%getconsumer())  
      endif
      
      !--------------------------------
      ! cas integration en helio
      !--------------------------------

      ! sortie cartesienne
      if (filepar%if_car.eq.1) then
       call sys%add_plan_out_step(n_out,1_8, io_cart%getconsumer())
      endif
     
      ! sortie elements elliptiques
      if (filepar%if_ell.eq.1) then
       call PhVh2ell%set_kind(filepar%out_ell)      
       call PhVh2ell%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVh2ell%set_output(1_8, io_ell%getconsumer())

       call sys%add_plan_out_step(n_out,1_8,PhVh2ell%getconsumer())
      endif
      

      ! lance l'integration
      call integrator%set_error_behavior(1)
      call integrator%set_sysode(sys)
      call integrator%set_integ_scheme(int_type) 
      call integrator%set_integ_time(tinit, dt, i)
      call integrator%set_integ_prec_step(filepar%dopri_eps,              &
     &        filepar%dopri_hmax, filepar%dopri_h)
      call integrator%integ_run(n_iter)
      
      if (sys%check_error()) then
       write(*,*) 'une erreur s est produite'
       write(*,*) 'code :', sys%get_error_code()
       write(*,*) 'time :', sys%get_error_time()
       write(*,*) 'body :', sys%get_error_body()
      endif

      ! genere le graphe d'appel 
      write(6,*) "digraph g {"
      graphnode = 0
      call sys%graphdot(0, graphnode, 6)
      write(6,*) "}"

      end subroutine
      end 

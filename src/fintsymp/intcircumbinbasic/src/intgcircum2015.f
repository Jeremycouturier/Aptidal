!!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file intgcircum2015.f 
!!  \brief  routine d'integration : lance l'integration 
!!   et ecrit dans les fichiers de sortie
!!
! history : creation 19/11/2015
!***********************************************************************
#include "realc.f"
!-----------------------------------------------------------
! En entree: tinit: temps initial
!            dt : intervalle de temps pour une iteration
!            n_iter : nombre d'iterations a effectuer
!            n_out : nombre d'iterations apres lequel on appelle ecrit
!-----------------------------------------------------------
      module mod_intg_symp
      contains
      
!-----------------------------------------------
!!> calcul de la constante de gauss en UA/an et masse solaire 
!----------------------------------------------------
      function calcG() result(cG)
      use constantes
      implicit none
      real(TREAL) :: cG !< constante de Gauss  en AU, an et masse solaire 
      cG = kgauss**2*(REALC(365.25e0))**2
      write (*,*) ' cG =', cG
      end function calcG

      

!***********************************************************************
!> fonction appellee par le distributeur de condition initiale pour la condition ligneci
!***********************************************************************
      subroutine cidist_integfunc(ligneci, filepar)
        use mod_io_int
        use mod_io_txt_ncomp
        use mod_io_txt_car
        use mod_syscircumnewt
        use mod_convertercircumbin
        use mod_ci_planewt
        use mod_ci_planewt
        use mod_filepar
        implicit none
        
        class(t_ci_base), intent(inout) :: ligneci !< condition initiale recue
        class(*), intent(inout) :: filepar !< donnee utilisateur : fichier de parametres
        
        real(TREAL), dimension(:), allocatable :: mpl ! masse normalise (somme des masses des etoiles a 1)
        real(TREAL), dimension(:,:), allocatable :: CI
        real(TREAL),allocatable,dimension(:,:) :: Xplan,XPplan
        real(TREAL),allocatable,dimension(:,:) :: V,Vhat
        real(TREAL),allocatable,dimension(:,:) :: Vinv, Vhatinv
        type(t_io_int) ::io_int
        type(t_io_txt_car) ::io_cart
        type(t_io_txt_ncomp) ::io_ell
        type(t_syscircumnewt),allocatable, target  :: sysC ! systeme en circumbinaire
        class(t_syscircumnewt), pointer  :: psys ! systeme en circumbinaire ou ...
        character(100)  :: int_type ! type de l'integrateur
        integer :: nplan ! nombre de planetes + 1 (car 2nde etoile a l'indice 1)
        integer k
        integer ci_type ! type des coordonnnees des conditions initiales
        real(TREAL) :: cG
        real(TREAL) :: mstar ! somme de la masse des 2 etoiles
        real(TREAL)  :: factmassorig !< facteur de masse pour se ramener dans les masses originales
        
! --- initialisation des constantes      
        cG = calcG()
        
        write(*,*) "start cidist_integfunc", trim(ligneci%m_id) 
        call ligneci%debug()

        select type(ligneci)
         class is(t_ci_planewt)
          nplan = ligneci%m_plan_nb ! attention : ici inclus la 1ere etoile
          ci_type = ligneci%m_ci_type
          allocate(mpl(0:nplan))
          allocate(CI(1:6,1:nplan))
          mpl(0:nplan) = ligneci%m_plan_mass(0:nplan)
! Masses normalisees a 1 = (somme des masses des 2 etoiles)
          mstar = mpl(0)+mpl(1)
          factmassorig = mstar 
          cG=cG*mstar
          do k=0,nplan
            mpl(k)=mpl(k)/mstar
          enddo
          ! ici le 1:nplan au lieu de 0:nplan car pas de coordonnees pour la 1ere etoile
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
         io_cart = t_io_txt_car("I/O cart circum")
         call io_cart%set_comp_nb(1+6*nplan)
         call io_cart%set_name(trim(filepar%nf_car),62)
         call io_cart%set_plan_nb(nplan)
         call io_cart%fcreate()
      end if


! --- Creation des tableaux
      allocate(Xplan(3,nplan))
      allocate(XPplan(3,nplan)) 
      allocate(V(3,nplan))
      allocate(Vhat(3,nplan)) 
      allocate(Vinv(3,nplan))
      allocate(Vhatinv(3,nplan)) 
      if (int_type(4:4).EQ.'H') then
        write(*,*)'coord. circumbinaire'
        filepar%int_var = 2               ! variables circum-binaires
        allocate(sysC)
        psys=>sysC
      else 
        write(*,*)'coord. integration non geres (mauvais integrateur)'
        stop
      endif

      
      ! fixe le systeme a integrer
      call psys%set_plan_nb(nplan-1) ! nombre de planetes 
       write(*,*) "mpl=", mpl
      call psys%set_mass(cG, mpl)

      !
      ! changement de coordonnees depuis les CI lues
      !
      write(*,*) 'CI='
      write(*,*) CI
      

      if (filepar%int_var.ne.2) then
       write(*,*) " cas non implementee : coord integree inconnu"
       stop
      
      else
      !
      !! integration circumbinaire : 2eme etoile  reperee par rapport a la premiere etoile
      !! les planetes sont reperees par rapport au barycentre des 2 etoiles
      !
       call circum_CI2VVhat(nplan,mpl,cG, ci_type,CI,V,Vhat)
       if (filepar%if_invar.eq.1) then 
        ! passage vers invariant 
        call circum_VVhat2inv(nplan,mpl,V,Vhat, Vinv, Vhatinv)
        Xplan = Vinv
        XPplan = Vhatinv
       else
        ! reste dans le repere des ci
        Xplan = V
        XPplan = Vhat
       endif
      endif


      write(*,*) ' V, Vhat dans cidist_integfunc'
      write(*,*)  Xplan,XPplan

      ! -- Conditions initiales des planetes et 2nde etoile
         write(*,*) CI
         
         call psys%set_plan_coord(Xplan,XPplan)


         call intg_symp(filepar,psys,cG,factmassorig,io_int,io_cart,      &
     &             io_ell)

         write(*,*) 'fin d integration'

        write(*,*) "end cidist_integfunc"
        
       deallocate(Xplan)
       deallocate(XPplan) 
       deallocate(V)
       deallocate(Vhat) 
       deallocate(Vinv)
       deallocate(Vhatinv) 
       deallocate(mpl)
       deallocate(CI)
       
         class default
          stop 'filepar in cidist_integfunc : bad class'
        end select 

       end subroutine cidist_integfunc
      
!------------------------------------------------------------
!> routine d'integration
!------------------------------------------------------------
      subroutine intg_symp(filepar, sys, cG, factmassorig, io_int,        &
     &                        io_cart, io_ell)
      use mod_syscircumnewt
      use mod_intsympbase
      use mod_schema
      use mod_io_int
      use mod_convertercircumbin
      use mod_io_txt_ncomp
      use mod_io_txt_car
      use mod_filepar
      
      implicit none
      type (t_extfilepar),target, intent(inout) :: filepar !< fichier de parametres + variables derivees
      type (t_io_int),intent(inout) ::io_int  !< sortie des integrales premieres
      type (t_io_txt_car), intent(inout) :: io_cart  !< sortie des elements cartesiens
      type (t_io_txt_ncomp), intent(inout) :: io_ell !< sortie des elements elliptiques
      class(t_syscircumnewt), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
      real(TREAL), intent(in) :: factmassorig !< facteur de masse pour se ramener dans les masses originales (notamment les integrales premieres)
      
     
      integer(8) :: n_iter
      integer(8) :: n_out
      real(TREAL) :: tinit,dt
      
      type (t_circum_VVhat2ell) :: VVhat2ell
      type (t_circum_VVhat2PhVh) :: VVhat2PhVh

      type(t_intsympbase) :: integrator
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
        ! on se ramene dans le systeme de masse originale
        call sys%set_factor_int(factmassorig, factmassorig)
      endif
      
      !! integration circumbinaire : 2eme etoile  reperee par rapport a la premiere etoile
      !! les planetes sont reperees par rapport au barycentre des 2 etoiles

      ! sortie cartesienne en coordonnes circum-binaire (V,Vhat)
      if (filepar%if_car.eq.1) then
       call VVhat2PhVh%set_mass(sys%m_star_mass,sys%m_plan_mass)      
       call VVhat2PhVh%set_output(1_8, io_cart%getconsumer())
       call sys%add_plan_out_step(n_out,1_8, VVhat2PhVh%getconsumer())
      endif
     
      ! sortie elements elliptiques
      if (filepar%if_ell.eq.1) then
       call VVhat2ell%set_kind(filepar%out_ell)      
       call VVhat2ell%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call VVhat2ell%set_output(1_8, io_ell%getconsumer())

       call sys%add_plan_out_step(n_out,1_8,VVhat2ell%getconsumer())
      endif
      

      ! lance l'integration
      call integrator%set_error_behavior(1)
      call integrator%set_syssymp(sys)
      call integrator%set_integ_scheme(int_type) 
      call integrator%set_integ_time(tinit, dt, i)
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

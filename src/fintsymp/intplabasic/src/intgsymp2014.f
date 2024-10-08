!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file filepar.f 
!!  \brief  routine d'integration : lance l'integration 
!!   et ecrit dans les fichiers de sortie
!!
! history : creation 17/03/2014
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


!***********************************************************************
!> fonction appellee par le distributeur de condition initiale pour la condition ligneci
!***********************************************************************
      subroutine cidist_integfunc(ligneci, filepar)
        use mod_io_int
        use mod_io_txt_ncomp
        use mod_io_txt_car
        use mod_sysplanewth
        use mod_sysplanewth_adapt
        use mod_sysplanewtj
        use mod_sysplanewtj_adapt
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
        real(TREAL),allocatable,dimension(:,:) :: Ph,Vh, Vb, Pj,Vj
        real(TREAL),allocatable,dimension(:,:) :: Phinv, Vbinv, Vhinv
        type(t_io_int) ::io_int
        type(t_io_txt_car) ::io_cart
        type(t_io_txt_ncomp) ::io_ell
        type(t_sysplanewtH),allocatable, target  :: sysH ! systeme en helio
        type(t_sysplanewtJ),allocatable, target  :: sysJ ! systeme en jacobi
        type(t_sysplanewtJ_adapt),allocatable, target  :: sysJ_adapt ! systeme en jacobi et pas adaptatif
        type(t_sysplanewtH_adapt),allocatable, target  :: sysH_adapt ! systeme en heliocentrique et pas adaptatif
        class(t_syspla), pointer  :: psys ! systeme en jacobi ou helio
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


! --- We always use Jacobi coordinates !!
! --- Creation des tableaux
      allocate(Xplan(3,nplan))
      allocate(XPplan(3,nplan)) 
      allocate(Ph(3,nplan))
      allocate(Vh(3,nplan)) 
      allocate(Vb(3,nplan)) 
      allocate(Pj(3,nplan))
      allocate(Vj(3,nplan)) 
      allocate(Phinv(3,nplan))
      allocate(Vhinv(3,nplan)) 
      allocate(Vbinv(3,nplan)) 
      if (int_type(4:4).NE.'H') then
        filepar%int_var = 0 ! variable de jacobi             
        if (filepar%type_pas.eq.0) then 
          write(*,*)'Jacobi Coord'
          allocate(sysJ)
          psys=>sysJ  
        else
          write(*,*)'Adaptative Jacobi Coord'
          allocate(sysJ_adapt)
          psys=>sysJ_adapt  
        endif  
      else if(int_type(5:5).NE.'W') then
          filepar%int_var = 1               ! variables heliocentriques canoniques
        if (filepar%type_pas.eq.0) then 
          write(*,*)'Heliocentric Coord'
          allocate(sysH)
          psys=>sysH
        else
          write(*,*)'Adaptative Heliocentric Coord'
          allocate(sysH_adapt)
          psys=>sysH_adapt
        endif   
      else 
        write(*,*)'Heliocentric Coord (Wisdom)'
        filepar%int_var = 2               ! variables heliocentriques canoniques
        if (filepar%type_pas.ne.0) then 
          write(*,*) 'Adapatitive Heliocentric Coord (Wisdom) not impl'
          stop 1
        endif
      endif

      
      ! fixe le systeme a integrer
      call psys%set_plan_nb(nplan)
      call psys%set_star_mass(cG, mpl(0))
      call psys%set_plan_mass(mpl(1:))

      !
      ! changement de coordonnees depuis les CI lues
      !
      write(*,*) 'CI='
      write(*,*) CI
      

      if (filepar%int_var.eq.0) then
      !
      ! integration en jacobi
      !
       call CI2PhVh(nplan,mpl,cG, ci_type,CI,Ph,Vh)
       if (filepar%if_invar.eq.1) then 
        ! passage vers invariant
        call coord_vh2vb(nplan,mpl,Vh,Vb)
        call PhVb2inv(nplan,mpl,Ph,Vb, Phinv, Vbinv)
        call coord_vb2vh(nplan,mpl,Vbinv,Vhinv)
        call PhVh2PjVj(nplan,mpl,Phinv,Vhinv, Pj, Vj)
       else
        ! reste dans le repere des ci
        call PhVh2PjVj(nplan,mpl,Ph,Vh, Pj, Vj)
       endif
       Xplan = Pj
       XPplan = Vj
      
      else
      !
      ! integration en helio
      !
       call CI2PhVh(nplan,mpl,cG, ci_type,CI,Ph,Vh)
       call coord_vh2vb(nplan,mpl,Vh,Vb)
       if (filepar%if_invar.eq.1) then 
        ! passage vers invariant 
        call PhVb2inv(nplan,mpl,Ph,Vb, Phinv, Vbinv)
        Xplan = Phinv
        XPplan = Vbinv
       else
        ! reste dans le repere des ci
        Xplan = Ph
        XPplan = Vb
       endif
      endif


      write(*,*) ' X , V dans cidist_integfunc'
      write(*,*)  Xplan,XPplan

      ! -- Conditions initiales des planetes
         write(*,*) CI
         
         call psys%set_plan_coord(Xplan,XPplan)


         call intg_symp(filepar,psys, cG, io_int, io_cart, io_ell)

         write(*,*) 'fin d integration'

        write(*,*) "end cidist_integfunc"
        
      if (allocated(sysJ)) then
         deallocate(sysJ)
      end if
      if (allocated(sysH)) then
         deallocate(sysH)
      end if
       deallocate(Xplan)
       deallocate(XPplan) 
       deallocate(Ph)
       deallocate(Vh) 
       deallocate(Vb) 
       deallocate(Pj)
       deallocate(Vj) 
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
      use mod_sysplanewth
      use mod_intsympbase
      use mod_sysplanewtj_adapt
      use mod_sysplanewth_adapt
      use mod_schema
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
      class(t_syspla), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
     
      integer(8) :: n_iter
      integer(8) :: n_out
      real(TREAL) :: tinit,dt
      
      type (t_plan_PhVb2PhVh) :: PhVb2PhVh
      type (t_plan_PjVj2PhVh) :: PjVj2PhVh
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb
      type (t_plan_PhVb2ell) :: PhVb2ell
      type (t_io_txt_ncomp) :: io_sch ! donnes supplemntaires du schema adaptatifs

      type(t_intsympbase) :: integrator
      
      integer(8) :: i
      logical :: issyshelio
      character(100)  :: int_type ! type de l'integrateur
      integer graphnode ! pour la generation du graphe
     
      n_iter = filepar%n_iter
      n_out = filepar%n_out
      tinit = filepar%tinit
      dt = filepar%dt
      int_type = filepar%int_type
      
      issyshelio = .false.
      i = 0
      
      if (trim(int_type(3:3)).eq.'H') then
       issyshelio = .true.
      endif 
      if (trim(int_type(4:4)).eq.'H') then
       issyshelio = .true.
      endif 
      
      write(*,*) 'int_type = ', int_type
      write(*,*) 'issyshelio = ', issyshelio

      call sys%clear_out_step()
      
      call integrator%set_type_pas(filepar%type_pas)
      
      ! sortie integrales premieres
      if (filepar%if_int.eq.1) then
        call sys%add_int_out_step(n_out,1_8,io_int%getconsumer())  
      endif
      
      ! sortie schema adaptatif
      if (filepar%type_pas.ne.0) then
        write(*,*)  ' nf_sch ', trim(filepar%nf_sch)
        if (NDATAJ_SCH.ne.NDATAH_SCH) then
          write(*,*) 'NDATAJ_SCH doit etre egal a NDATAH_SCH'
          stop 1
        endif
        io_sch = t_io_txt_ncomp("I/O sch adapatif")
        call io_sch%set_comp_nb(NDATAJ_SCH)
        call io_sch%set_name(trim(filepar%nf_sch),63)
        call io_sch%fcreate()
        select type(sys)
         class is(t_sysplanewtJ_adapt)
           call sys%add_adapt_out_step(n_out,1_8,io_sch%getconsumer()) 
         class is(t_sysplanewtH_adapt)
           call sys%add_adapt_out_step(n_out,1_8,io_sch%getconsumer()) 
        end select    
      endif


      if (issyshelio) then
      !--------------------------------
      ! cas integration en helio
      !--------------------------------

      ! sortie cartesienne
      if (filepar%if_car.eq.1) then
       call PhVb2PhVh%set_mass(sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2PhVh%set_output(1_8,io_cart%getconsumer())
       call sys%add_plan_out_step(n_out,1_8, PhVb2PhVh%getconsumer())
      endif
     
      ! sortie elements elliptiques
      if (filepar%if_ell.eq.1) then
       call PhVb2ell%set_kind(filepar%out_ell)      
       call PhVb2ell%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell%set_output(1_8, io_ell%getconsumer())

       call sys%add_plan_out_step(n_out,1_8,PhVb2ell%getconsumer())
      endif
      
      else 
      !--------------------------------
      ! cas integration en jacobi
      !--------------------------------

      ! sortie cartesienne 
      if (filepar%if_car.eq.1) then
        call PjVj2PhVh%set_mass(sys%m_star_mass,sys%m_plan_mass)    
        call PjVj2PhVh%set_output(1_8, io_cart%getconsumer())
        call sys%add_plan_out_step(n_out,1_8, PjVj2PhVh%getconsumer())
      endif
             
      ! sortie elliptique 
      if (filepar%if_ell.eq.1) then
       call PhVb2ell%set_kind(filepar%out_ell)      
       call PhVb2ell%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell%set_output(1_8,io_ell%getconsumer())

        call PjVj2PhVb%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb%set_output(1_8,PhVb2ell%getconsumer())
        call sys%add_plan_out_step(n_out,1_8,PjVj2PhVb%getconsumer())
      endif

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

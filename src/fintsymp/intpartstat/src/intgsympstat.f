!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file intgsympstat.f 
!!  \brief  routine d'integration : construit le graphe de la simulation 
!! et lance l'integration 
!!
! history : creation 15/09/2014
!***********************************************************************
#include "realc.f"
      module mod_intg_symp
      
        use mod_filepar
        use mod_ci_planewt

       private init_minmax_diffalp
       private intg_symp
       private init_naf_alkhqp
       private init_naf_diffalp

!-----------------------------------------------------------------------
!!> type aggregateur pour le passage d'argument a  cidist_integfunc
!-----------------------------------------------------------------------
        type t_cidist_integfunc_arg
            type(t_extfilepar) :: m_filepar !< fichier de parametres
            class(t_ci_base), pointer :: m_cipla   !< condition du systeme planetaire
            logical :: m_boutpla  !< = true => sortir les donnees des planetes
        end type t_cidist_integfunc_arg


      contains
      
!***********************************************************************
!> fonction appellee par le distributeur de condition initiale pour les conditions tabci
!***********************************************************************
      subroutine cidist_integfunc(nbci, tabci, multipledata)
        use mod_io_int
        use mod_io_txt_ncomp
        use mod_io_txt_car
        use mod_forplatab
        use mod_converterpart 
        use mod_ci_partnewt
        use mod_ci_planewt
        use mod_filepar
        use mod_coordpart
        use mod_syspartnewtH
        use mod_io_txt_carpart
        use mod_io_txt_ellpart
        use mod_gmsun
        implicit none
        
        integer, intent(in) :: nbci !< nombre d'elements valides dans ci
        type(t_ci_base_ptr), dimension(:), intent(inout) ::   tabci !< tableau de conditions initiales
        class(*), intent(inout) :: multipledata !< donnee utilisateur : fichier de parametres + autres
        
        real(TREAL), dimension(:), allocatable :: mpl ! masse normalise (masse de l'etoile a 1)
        real(TREAL), dimension(:,:), allocatable :: COORDpla, COORDpart
        real(TREAL),allocatable,dimension(:,:) :: Xplan,XPplan
        real(TREAL),allocatable,dimension(:,:) :: Ph,Vh, Vb
        real(TREAL),allocatable,dimension(:,:) :: Phinv, Vbinv, Vhinv
        real(TREAL),allocatable,dimension(:,:) :: Xpart,XPpart
        real(TREAL),allocatable,dimension(:,:) :: Phpart,Vhpart, Vbpart
        real(TREAL),allocatable,dimension(:,:) :: Phinvpart,Vhinvpart
        real(TREAL),allocatable,dimension(:,:) :: Vbinvpart
        class(t_sysplaextH), allocatable  :: sysH !<  systeme en helio
        class(t_forcing), pointer  :: forcedriver!<  systeme force
        type(t_syspartnewtH) :: syspartH ! particules en heliocentrique
        character(100)  :: int_type ! type de l'integrateur
        integer :: nplan ! nombre de planetes
        integer k
        integer cipla_type ! type des coordonnnees des conditions initiales des planetes
        integer cipart_type ! type des coordonnnees des conditions initiales des particules
        real(TREAL) :: cG
        type(t_extfilepar) :: filepar ! fichier de parametres
        class(t_ci_base), pointer :: cipla ! conditions initiales des planetes
        integer npart ! nombre de particules
        class(t_ci_base), pointer ::   ptrcitype 
                
        select type(multipledata)
         class is(t_cidist_integfunc_arg)

         filepar = multipledata%m_filepar
         cipla => multipledata%m_cipla

!--- desactive si necessaire les sorties des planetes
        if (multipledata%m_boutpla.eqv..false.) then
            filepar%if_int = 0
            filepar%if_ell_pla = 0
            filepar%if_car_pla = 0
            filepar%naf_alkhqp_pla%m_compute = 0 
        endif


! --- initialisation des constantes      
        npart = nbci
        cG = calcGM_sun(filepar%ref_gmsun)
        
        write(*,*) "start cidist_integfunc"
        do k=1, npart
          write(*,*) 'nom ci : ', trim(tabci(k)%m_ptr%m_id)
          call tabci(k)%m_ptr%debug()
        enddo
!
! conditions initiales des planetes
!
        select type(cipla)
         class is(t_ci_planewt)
          nplan = cipla%m_plan_nb
          cipla_type = cipla%m_ci_type
          allocate(mpl(0:nplan))
          allocate(COORDpla(1:6,1:nplan))
          mpl(0:nplan) = cipla%m_plan_mass(0:nplan)
! Masses normalisees a 1
          cG=cG*mpl(0)
          do k=1,nplan
            mpl(k)=mpl(k)/mpl(0)
          enddo
          mpl(0)=1.D0
          ! ici le 1:nplan au lieu de 0:nplan car ci_type: heliocentrique
          COORDpla(1:6,1:nplan) = cipla%m_plan_coord_ci(1:6,1:nplan)
         
         class default
          stop 'cipla in cidist_integfunc : bad class'
        end select 


!
! conditions initiales des particules
!
         allocate(COORDpart(1:6,1:npart))
        do k=1, npart
         ptrcitype => tabci(k)%m_ptr
         select type(ptrcitype)
           class is(t_ci_partnewt)
             cipart_type = ptrcitype%m_ci_type
             COORDpart(1:6,k) = ptrcitype%m_plan_coord_ci(1:6)
          
           class default
             stop 'ptrcitype in cidist_integfunc : bad class'
         end select 
        enddo

      int_type = trim(filepar%int_type)
      


! --- Creation des tableaux
      allocate(Xplan(3,nplan))
      allocate(XPplan(3,nplan)) 
      allocate(Ph(3,nplan))
      allocate(Vh(3,nplan)) 
      allocate(Vb(3,nplan)) 
      allocate(Phinv(3,nplan))
      allocate(Vhinv(3,nplan)) 
      allocate(Vbinv(3,nplan)) 
      
      allocate(Xpart(3,npart))
      allocate(XPpart(3,npart)) 
      allocate(Phpart(3,npart))
      allocate(Vhpart(3,npart)) 
      allocate(Vbpart(3,npart)) 
      allocate(Phinvpart(3,npart))
      allocate(Vhinvpart(3,npart)) 
      allocate(Vbinvpart(3,npart)) 
      if (int_type(4:4).NE.'H') then
        write(*,*)'Jacobi Coord'
        write(*,*) "integration en Jacobi impossible"
        stop                
      else if(int_type(5:5).NE.'W') then
        write(*,*)'Heliocentric Coord'
        filepar%int_var = 1               ! variables heliocentriques canoniques
      else 
        write(*,*)'Heliocentric Coord (Wisdom)'
        write(*,*) "integration en HW impossible"
        stop                
      endif

      
!--- initialisation du modele planetaire : integree symplectiquement ou solution forcee
       select case (filepar%if_orb_pla)
       ! systeme planetaire en helio integree symplectiquement
       case (0)  
        allocate(t_sysplaextH::sysH)
        
       ! systeme planetaire en helio fournie par une solution tabulee
       case default
         if (filepar%if_invar.eq.1) then 
         ! passage vers invariant non supporte
         write(*,*) 'if_orb_pla!=0 alors if_invar doit etre 0'
         stop
         endif
         select case (filepar%if_orb_pla)
           case (1)  
             allocate(t_forcingplatab::forcedriver)
             select type(forcedriver)
               class is(t_forcingplatab)
                 call forcedriver%load(filepar%orb_pla%m_coord,             &
     &               filepar%orb_pla%m_nf, nplan) 
                 cipla_type = filepar%orb_pla%m_coord
                 call forcedriver%compute(0.D0, Xplan, XPplan)
                 COORDpla(1:3,1:nplan) = Xplan(1:3,1:nplan)
                 COORDpla(4:6,1:nplan) = XPplan(1:3,1:nplan)
             end select
             
           case default
             write(*,*) 'if_orb_pla=', filepar%if_orb_pla
             stop 'if_orb_pla : cas inconnu'
         end select

         ! allocation du systeme force
         allocate(t_forplaextH::sysH)
         select type(sysH)
           class is(t_forplaextH)
             call sysH%set_forcing(forcedriver)
         end select
       end select

      !
      ! changement de coordonnees depuis les CI lues
      !
      write(*,*) 'COORDpla='
      write(*,*) COORDpla
      write(*,*) 'COORDpart='
      write(*,*) COORDpart
      

      !
      ! integration en helio
      !
       call CI2PhVh_part(nplan,mpl,cG, cipla_type,COORDpla,Ph,Vh,        &
     &                   npart, COORDpart, Phpart, Vhpart)
       call coord_vh2vb(nplan,mpl,Vh,Vb)
       call coord_vh2vb_part(nplan,mpl,Vh, npart, Vhpart, Vbpart)
       if (filepar%if_invar.eq.1) then 
        ! passage vers invariant 
        call PhVb2inv_part(nplan,mpl,Ph,Vb, Phinv, Vbinv,                &
     &    npart, Phpart,Vbpart, Phinvpart, Vbinvpart)
        Xplan = Phinv
        XPplan = Vbinv
        Xpart = Phinvpart
        XPpart = Vbinvpart
       else
        ! reste dans le repere des ci
        Xplan = Ph
        XPplan = Vb
        Xpart = Phpart
        XPpart = Vbpart
       endif


      write(*,*) ' X , V dans cidist_integfunc'
      write(*,*)  Xplan,XPplan
      write(*,*) ' Xpart , Vpart dans cidist_integfunc'
      write(*,*)  Xpart,XPpart

      ! fixe le systeme a integrer
      call sysH%set_plan_nb(nplan)
      call sysH%set_star_mass(cG, mpl(0))
      call sysH%set_plan_mass(mpl(1:))
         
      call sysH%set_plan_coord(Xplan,XPplan)
      call sysH%set_extension_nb(1)
      call syspartH%set_nb(nplan, npart)
      call syspartH%set_star_mass(cG, mpl(0))
      call syspartH%set_coord(Xpart, XPpart)
      call sysH%add_extension(1,syspartH)



       call intg_symp(filepar,sysH, cG, cipla, tabci)

       write(*,*) 'fin d integration'

       write(*,*) "end cidist_integfunc"
        
       deallocate(Xplan)
       deallocate(XPplan) 
       deallocate(Ph)
       deallocate(Vh) 
       deallocate(Vb) 
       deallocate(Phinv)
       deallocate(Vhinv) 
       deallocate(Vbinv) 
       deallocate(mpl)
       deallocate(COORDpla)
       deallocate(COORDpart)
       deallocate(Xpart)
       deallocate(XPpart) 
       deallocate(Phpart)
       deallocate(Vhpart) 
       deallocate(Vbpart) 
       deallocate(Phinvpart)
       deallocate(Vhinvpart) 
       deallocate(Vbinvpart) 
      
       if (associated(forcedriver)) then
         deallocate(forcedriver)
       endif
       deallocate(sysH)

       ! uniquement la premiere fois, on effectue les sorties des planetes
       multipledata%m_boutpla = .false.
 
         class default
          stop 'multipledata in cidist_integfunc : bad class'
        end select 

      end subroutine cidist_integfunc
      
      

!------------------------------------------------------------
!> routine d'integration
!------------------------------------------------------------
      subroutine intg_symp(filepar, sys, cG, cipla, tabci)
#if USE_MPI
      use mpi
#endif
      use mod_sysplaexth
      use mod_intsympbase
      use mod_schema
      use mod_io_int
      use mod_converter
      use mod_io_txt_ncomp
      use mod_io_txt_car
      use mod_filepar
      use mod_minmax
      use mod_minmax_mda
      use mod_naf_base
      use mod_naf_diffang
      use mod_ctrl_distpartstar
      use mod_ctrl_distpartpla
      use mod_coordpart
      use mod_syspartnewtH
      use mod_converterpart 

      implicit none
      type (t_extfilepar),target, intent(inout) :: filepar !< fichier de parametres + variables derivees
      class(t_sysplaexth), intent(inout) :: sys  !< systeme planetaire a integrer
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
      type(t_ci_base_ptr), dimension(:), intent(inout) ::   tabci !< tableau de conditions initiales
      class(t_ci_base), intent(inout) :: cipla ! conditions initiales des planetes
     
      integer(8) :: n_iter
      integer(8) :: n_out
      real(TREAL) :: tinit,dt
      integer :: nplan, k, ntot
      integer npart
      
      ! pour les planetes
      type (t_plan_PhVb2ell) :: PhVb2ell_nafk_pla ! pour naf en alkhqp
      integer(8) :: len_nafk_pla ! longueur des tranches de naf
      type (t_naf_base) :: naf_alkhqp_pla ! naf en alkhqp

      ! pour les particules
      class(t_sysextension), pointer :: psysext 
      class(t_syspartnewtH), pointer :: psyspartH =>null()

#if 0
      type (t_plan_PjVj2PhVh) :: PjVj2PhVh
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb, PjVj2PhVb_2
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_mda ! pour les min/max des differences d'angle
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_nafk ! pour naf en alkhqp
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_nafd ! pour naf des differences d'angle
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_ds ! pour controle des distances
#endif

      type (t_plan_PhVb2PhVh) :: PhVb2PhVh
      type (t_plan_PhVb2ell) :: PhVb2ell
      type (t_part_PhVb2ell) :: PhVb2ell_maei
      type (t_part_PhVb2ell) :: PhVb2ell_mda ! pour les min/max des differences d'angle
      type (t_part_PhVb2ell) :: PhVb2ell_nafk ! pour naf en alkhqp
      type (t_part_PhVb2ell) :: PhVb2ell_nafd ! pour naf des differences d'angle
      type (t_part_PhVb2PhVh) :: PhVb2PhVhpart
      type (t_part_PhVb2ell) :: PhVb2ellpart
   
      type (t_minmax) :: minmax_aei
      integer(8) :: len_minmax_aei
      integer, dimension(:), allocatable :: aeicomp
      
      type(t_minmax_mda) :: minmax_diffalp ! min/moy/max en l1-l2,pibar1-pibar2, a1-a2


      type (t_naf_base) :: naf_alkhqp ! naf en alkhqp
      type (t_naf_diffang) :: naf_diffalp ! naf en differences d'angle
      integer(8) :: len_nafk, len_nafd ! longueur des tranches de naf

      type (t_io_int) ::io_int  ! sortie des integrales premieres
      type (t_io_txt_car) :: io_cart_pla  ! sortie des elements cartesiens des planetes
      type (t_io_txt_ncomp) :: io_ell_pla ! sortie des elements elliptiques des planetes
      
      type (t_ctrl_distpartstar) :: ctrl_diststar ! controle des distances a l'etoile
      
      type (t_ctrl_distpartpla) :: ctrl_distpla ! controle des distances aux planetes

      type(t_arret) :: arret ! code de l'arret
      type(t_intsympbase) :: integrator
      integer(8) :: i
      logical :: issyshelio
      character(100)  :: int_type ! type de l'integrateur
      character(256)  :: nameid
      integer graphnode ! pour la generation du graphe
      integer writegraph ! pour la generation du graphe
      integer rankmpi ! rang mpi
      integer code ! code mpi
      
      integer(8) :: len_malp
      

      rankmpi = 0
      code = 0
#if USE_MPI
      call MPI_COMM_RANK ( MPI_COMM_WORLD,rankmpi,code)
#endif

      n_iter = filepar%n_iter
      n_out = filepar%n_out
      tinit = filepar%tinit
      dt = filepar%dt
      int_type = filepar%int_type
      nplan = sys%m_plan_nb
      
      issyshelio = .false.
      i = 0
      
      if (filepar%int_var.eq.1) then
       issyshelio = .true.
      endif 
      
      write(*,*) 'int_type = ', int_type
      write(*,*) 'issyshelio = ', issyshelio

      call sys%clear_out_step()
      
      call sys%get_extension(1, psysext) ! recupere les particules
      select type(psysext)
        class is(t_syspartnewtH)
           psyspartH=>psysext
           npart = psyspartH%get_nb()
        class default
               stop 'intg_symp : bad class t_syspartnewtH for psysext'
      end select 

      ntot=npart+nplan

      ! sortie elliptique des planetes
      if (filepar%if_ell_pla.eq.1) then
         write(*,*)  ' nf_ell ', trim(filepar%nf_ell_pla)
         io_ell_pla = t_io_txt_ncomp("I/O ell")
         call io_ell_pla%set_comp_nb(1+6*nplan)
         call io_ell_pla%set_name(trim(filepar%nf_ell_pla),61)
         call io_ell_pla%fcreate()
      end if
      
      ! sortie cartesienne des planetes
      if (filepar%if_car_pla.eq.1) then
         write(*,*) ' nf_car_pla ', trim(filepar%nf_car_pla)
         io_cart_pla = t_io_txt_car("I/O cart")
         call io_cart_pla%set_comp_nb(1+6*nplan)
         call io_cart_pla%set_name(trim(filepar%nf_car_pla),62)
         call io_cart_pla%set_plan_nb(nplan)
         call io_cart_pla%fcreate()
      end if

      ! sortie integrales premieres
      if (filepar%if_int.eq.1) then
        write(*,*)  ' nf_int ', trim(filepar%nf_int)
        call io_int%set_name(trim(filepar%nf_int),60)
        call io_int%fcreate()
        call sys%add_int_out_step(n_out,1_8,io_int%getconsumer())  
      endif
      
      ! sortie elliptique des particules
      if (filepar%if_ell_part.eq.1) then
         call filepar%io_ell%set_col_nb(1+6, npart, 1+6*nplan+6*npart)
         call filepar%io_ell%set_part_nb(nplan, npart)
         do k=1, npart
           call filepar%io_ell%set_pref2dline(k,tabci(k)%m_ptr%m_id) 
         enddo
      end if
      
      ! sortie cartesienne des particules
      if (filepar%if_car_part.eq.1) then
         call filepar%io_cart%set_col_nb(1+6, npart, 1+6*nplan+6*npart)
         call filepar%io_cart%set_part_nb(nplan, npart)
         do k=1, npart
           call filepar%io_cart%set_pref2dline(k,tabci(k)%m_ptr%m_id) 
         enddo
      end if


      ! sortie min/moy/max a,e,i
      if (filepar%minmax_aei%m_compute.eq.1) then
       call filepar%io_minmax_aei%set_col_nb(10,npart,1+9*npart)
       do k=1, npart
           call filepar%io_minmax_aei%set_pref2dline(k,                 &
     &                  tabci(k)%m_ptr%m_id) 
       enddo
       len_minmax_aei = filepar%minmax_aei%m_stepout                    &
     &                  /filepar%minmax_aei%m_stepcalc
       call minmax_aei%set_stepout(len_minmax_aei)
       allocate(aeicomp(1:3*npart))
       do k=0, npart-1
        aeicomp(3*k+1) = 6*nplan+6*k+1 ! pour a
        aeicomp(3*k+2) = 6*nplan+6*k+2 ! pour e
        aeicomp(3*k+3) = 6*nplan+6*k+3 ! pour I
       enddo
       call minmax_aei%set_component(aeicomp,6*ntot)
       deallocate(aeicomp)
       call minmax_aei%set_output(1_8,                                   &
     &       filepar%io_minmax_aei%getconsumer())
       call minmax_aei%set_graphdotname("min/max a,e,i")
      
       call PhVb2ell_maei%set_kind(filepar%minmax_aei%m_elltype)      
       call PhVb2ell_maei%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_maei%set_nb(npart)      
       call PhVb2ell_maei%set_output(1_8, minmax_aei%getconsumer())
      endif

      ! sortie min/moy/max difference d'angle
      if (filepar%minmax_diffalp%m_compute.eq.1) then
        len_malp = filepar%minmax_diffalp%m_stepout                     &
     &                  /filepar%minmax_diffalp%m_stepcalc
        call init_minmax_diffalp(filepar%minmax_diffalp,minmax_diffalp, &
     &         filepar%io_minmax_diffalp, filepar%minmax_diffalp_pla,   &
     &         sys, psyspartH, cG, tabci, PhVb2ell_mda)
      endif

      ! sortie naf en alkhqp
      if (filepar%naf_alkhqp%m_compute.ne.0) then
        len_nafk = filepar%naf_alkhqp%m_stepout                         &
     &                  /filepar%naf_alkhqp%m_stepcalc
        call init_naf_alkhqp(filepar%naf_alkhqp,naf_alkhqp,             &
     &         filepar%io_naf_alkhqp,                                   &
     &         sys, psyspartH, cG, tabci, PhVb2ell_nafk)
      endif

      ! sortie naf en difference d'angle
      if (filepar%naf_diffalp%m_compute.eq.1) then
        len_nafd = filepar%naf_diffalp%m_stepout                        &
     &                  /filepar%naf_diffalp%m_stepcalc
        call init_naf_diffalp(filepar%naf_diffalp,naf_diffalp,          &
     &         filepar%io_naf_diffalp, filepar%naf_diffalp_pla,         &
     &         sys, psyspartH, cG, tabci, PhVb2ell_nafd)
      endif

      ! controle des distances a l'etoile
      if (filepar%ctrl_diststar%m_compute.eq.1) then
        call ctrl_diststar%set_minmax(nplan, npart,                     &
     &         filepar%ctrl_diststar%m_distmin,                         &
     &         filepar%ctrl_diststar%m_distmax)
      endif

      ! sortie naf en alkhqp aux planetes
      if (filepar%naf_alkhqp_pla%m_compute.ne.0) then
        len_nafk_pla = filepar%naf_alkhqp_pla%m_stepout                 &
     &                  /filepar%naf_alkhqp_pla%m_stepcalc
        call init_naf_alkhqp_pla(filepar%naf_alkhqp_pla,naf_alkhqp_pla, &
     &         filepar%io_naf_alkhqp_pla,                               &
     &         sys, cG, cipla%m_id, PhVb2ell_nafk_pla)
      endif

      ! controle des distances aux planetes
      if (filepar%ctrl_distpla%m_compute.eq.1) then
        call ctrl_distpla%set_minmax(nplan, npart,                      &
     &         filepar%ctrl_distpla%m_distmin)
      endif

      !--------------------------------
      !--------------------------------
      ! code dependant de l'integrateur
      !--------------------------------
      !--------------------------------
      if (issyshelio) then
      !--------------------------------
      ! cas integration en helio
      !--------------------------------

      ! sortie cartesienne des planetes
      if (filepar%if_car_pla.eq.1) then
       call PhVb2PhVh%set_mass(sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2PhVh%set_output(1_8,io_cart_pla%getconsumer())
       call sys%add_plan_out_step(n_out,1_8, PhVb2PhVh%getconsumer())
      endif
     
      ! sortie elements elliptiques des planetes
      if (filepar%if_ell_pla.eq.1) then
       call PhVb2ell%set_kind(filepar%out_ell_pla)      
       call PhVb2ell%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell%set_output(1_8, io_ell_pla%getconsumer())
       call sys%add_plan_out_step(n_out,1_8,PhVb2ell%getconsumer())
      endif
      
      ! sortie cartesienne des particules
      if (filepar%if_car_part.eq.1) then
       call PhVb2PhVhpart%set_mass(sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2PhVhpart%set_nb(npart)      
       call PhVb2PhVhpart%set_output(1_8,filepar%io_cart%getconsumer())
       call psyspartH%add_part_out_step(n_out,1_8,                         &
     &           PhVb2PhVhpart%getconsumer())
      endif
     
      ! sortie elements elliptiques des particules
      if (filepar%if_ell_part.eq.1) then
       call PhVb2ellpart%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ellpart%set_nb(npart)      
       call PhVb2ellpart%set_output(1_8, filepar%io_ell%getconsumer())
       call PhVb2ellpart%set_kind(filepar%out_ell_part)      
       call psyspartH%add_part_out_step(n_out,1_8,                         &
     &           PhVb2ellpart%getconsumer())
      endif
      

      ! sortie min/max a,e,i
      if (filepar%minmax_aei%m_compute.eq.1) then
       call psyspartH%add_part_out_step(filepar%minmax_aei%m_stepcalc,  &
     &          1_8, PhVb2ell_maei%getconsumer())
      endif
      
      ! sortie min/moy/max difference d'angle
      if (filepar%minmax_diffalp%m_compute.eq.1) then
       call psyspartH%add_part_out_step(                                &
     &          filepar%minmax_diffalp%m_stepcalc,                      &
     &          1_8, PhVb2ell_mda%getconsumer())
      endif

      ! sortie naf alkhqp
      if (filepar%naf_alkhqp%m_compute.ne.0) then
       call psyspartH%add_part_out_step(filepar%naf_alkhqp%m_stepcalc,   &
     &          len_nafk, PhVb2ell_nafk%getconsumer())
      endif
      
      ! sortie naf difference d'angle
      if (filepar%naf_diffalp%m_compute.eq.1) then
       call psyspartH%add_part_out_step(filepar%naf_diffalp%m_stepcalc,   &
     &          len_nafd, PhVb2ell_nafd%getconsumer())
      endif

      ! controle des distances avec l'etoile
      if (filepar%ctrl_diststar%m_compute.eq.1) then
       call psyspartH%add_part_out_step(                                   &
     &          filepar%ctrl_diststar%m_stepcalc, 1_8,                     &
     &          ctrl_diststar%getconsumer())
      endif

      ! sortie naf alkhqp des planetes
      if (filepar%naf_alkhqp_pla%m_compute.ne.0) then
       call sys%add_plan_out_step(filepar%naf_alkhqp_pla%m_stepcalc,         &
     &          len_nafk_pla, PhVb2ell_nafk_pla%getconsumer())
      endif

      ! controle des distances aux planetes
      if (filepar%ctrl_distpla%m_compute.eq.1) then
       call psyspartH%add_part_out_step(                                 &
     &          filepar%ctrl_distpla%m_stepcalc, 1_8,                    &
     &          ctrl_distpla%getconsumer())
      endif

      else 
      !--------------------------------
      ! cas integration en jacobi
      !--------------------------------
      write(*,*) "cas integration en jacobi non implemente"
      stop
      
#if 0
      ! sortie cartesienne 
      if (filepar%if_car_pla.eq.1) then
        call PjVj2PhVh%set_mass(sys%m_star_mass,sys%m_plan_mass)    
        call PjVj2PhVh%set_output(1, io_cart_pla%getconsumer())
        call sys%add_plan_out_step(n_out,1, PjVj2PhVh%getconsumer())
      endif
             
      ! sortie elliptique 
      if (filepar%if_ell_pla.eq.1) then
       call PhVb2ell%set_kind(filepar%out_ell_pla)      
       call PhVb2ell%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell%set_output(1,io_ell_pla%getconsumer())

        call PjVj2PhVb%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb%set_output(1,PhVb2ell%getconsumer())
        call sys%add_plan_out_step(n_out,1,PjVj2PhVb%getconsumer())
      endif

      ! sortie min/max a,e,i
      if (filepar%minmax_aei%m_compute.eq.1) then
        call PjVj2PhVb_2%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_2%set_output(len_malp,PhVb2ell_2%getconsumer())
        call sys%add_plan_out_step(filepar%minmax_aei%m_stepcalc,       &
     &          len_minmax_aei, PjVj2PhVb_2%getconsumer())
      endif

      ! sortie min/moy/max difference d'angle
      if (filepar%minmax_diffalp%m_compute.eq.1) then
        call PjVj2PhVb_mda%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_mda%set_output(1,PhVb2ell_mda%getconsumer())
        call sys%add_plan_out_step(filepar%minmax_diffalp%m_stepcalc,   &
     &          len_malp, PjVj2PhVb_mda%getconsumer())
      endif

      ! sortie naf alkhqp
      if (filepar%naf_alkhqp%m_compute.ne.0) then
        call PjVj2PhVb_nafk%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_nafk%set_output(len_nafk,                        &
     &          PhVb2ell_nafk%getconsumer())
        call sys%add_plan_out_step(filepar%naf_alkhqp%m_stepcalc,       &
     &          len_nafk, PjVj2PhVb_nafk%getconsumer())
      endif
      
      ! sortie naf difference d'angle
      if (filepar%naf_diffalp%m_compute.eq.1) then
        call PjVj2PhVb_nafd%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_nafd%set_output(len_nafd,                        &
     &          PhVb2ell_nafd%getconsumer())
        call sys%add_plan_out_step(filepar%naf_diffalp%m_stepcalc,      &
     &          len_nafd, PjVj2PhVb_nafd%getconsumer())
      endif

      ! controle des distances
      if (filepar%ctrl_diststar%m_compute.eq.1) then
        call PjVj2PhVb_ds%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_ds%set_output(1, ctrl_dist%getconsumer())
        call sys%add_plan_out_step(filepar%ctrl_diststar%m_stepcalc,    &
     &          1, PjVj2PhVb_ds%getconsumer())
      endif

#endif
      endif
      !--------------------------------
      !--------------------------------
      ! fin du code dependant de l'integrateur
      !--------------------------------
      !--------------------------------

      ! lance l'integration
      call integrator%set_error_behavior(1)
      call integrator%set_syssymp(sys)
      call integrator%set_integ_scheme(int_type) 
      call integrator%set_integ_time(tinit, dt, i)
      call integrator%integ_run(n_iter)
      
      ! verifie le code controle global de l'integration
      if (sys%check_error()) then
       arret=sys%get_error()
       nameid = "SYSTEM"//tabci(1)%m_ptr%m_id
       if (arret%m_cause.ne.-8) then
       ! cas ou ce n'est pas lie a l'absence de particules
           write(*,*) 'une erreur s est produite'
           write(*,*) 'code :', arret%m_cause
           write(*,*) 'time :', arret%m_time
           write(*,*) 'body :', arret%m_body
           write(*,*) 'value :', arret%m_value
           write(filepar%io_control_pla,1004)  trim(nameid),                  &
     &  arret%m_cause, tinit, arret%m_time, arret%m_body, arret%m_value
       else
           write(filepar%io_control_pla,1004)  trim(nameid),                  &
     &   0, tinit, integrator%get_integ_current_time(), 0, REALC(0.)
       endif
      else
           nameid = "SYSTEM"//tabci(1)%m_ptr%m_id
           write(filepar%io_control_pla,1004)  trim(nameid),                  &
     &   0, tinit, integrator%get_integ_current_time(), 0, REALC(0.)
      endif
      
      ! ecriture du controle de l'integration
      do k=1, npart
       arret = psyspartH%m_part_mce(k)
       if (arret%m_stop.ne.0) then
       write(filepar%io_control,1004) trim(tabci(k)%m_ptr%m_id),           &
     &  arret%m_cause, tinit, arret%m_time, arret%m_body,arret%m_value
      else
       write(filepar%io_control,1004) trim(tabci(k)%m_ptr%m_id),           &
     &   0, tinit, integrator%get_integ_current_time(), 0, REALC(0.)
      endif 
      enddo

1004  format(A,1X,I2,1X,FMT_TREAL,1X,FMT_TREAL,I4,1X,FMT_TREAL)  

      ! genere le graphe d'appel 
#if USE_MPI
      if (rankmpi.eq.1) then 
#else
      if (rankmpi.eq.0) then 
#endif
       writegraph = 0
       if (writegraph.eq.1) then
        write(6,*) "digraph g {"
        graphnode = 0
        call sys%graphdot(0, graphnode, 6)
        write(6,*) "}"
       endif
      endif
      
      end subroutine

!------------------------------------------------------------
!> initialise l'analyse en min/max/moy des differences d'angle des particules
!!  precise quelles composantes sont utilisees
!------------------------------------------------------------
      subroutine init_minmax_diffalp(par_minmax_diffalp,                    &
     &              comp_minmax_diffalp,  io_minmax_diffalp,                &
     &              minmax_diffalp_pla, sys, psyspartH, cG, tabci,          &
     &              PhVb2ell_mda)
        use mod_filepar
        use mod_syspla
        use mod_converter
        use mod_minmax_mda
        use mod_converterpart 
        use mod_sysplaexth
        use mod_syspartnewtH
      implicit none
      type(t_minmax_par), intent(inout) :: par_minmax_diffalp !< donnee du fichier de parametres pour min/moy/max en l1-l2,pibar1-pibar2, a1-a2
      type(t_minmax_mda), intent(inout) :: comp_minmax_diffalp !< module de calcul des min/moy/max en l1-l2,pibar1-pibar2, a1-a2
      type (t_io_txt_mligncomp), intent(inout) :: io_minmax_diffalp  !< sortie des min/max des differences a_p1-a_p2,....
      integer,dimension(1:1), intent(in) :: minmax_diffalp_pla !< numero des 2 planetes utilisees pour minmax_diffalp
      class(t_sysplaexth), intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_part_PhVb2ell), intent(inout) :: PhVb2ell_mda !< buffer des elements elliptiques en entree
      type(t_ci_base_ptr), dimension(:), intent(inout) ::   tabci !< tableau de conditions initiales
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
      class(t_syspartnewtH), intent(inout) :: psyspartH !< systemes des particules

      integer(8) len_minmax
      integer, dimension(1:3*psyspartH%get_nb()) :: comp_p1, comp_p2
      logical, dimension(1:3*psyspartH%get_nb()) :: angle
      integer p1, nplan, npart, k, ntot

       nplan = sys%m_plan_nb
       npart = psyspartH%get_nb()
       ntot = nplan+npart

       call io_minmax_diffalp%set_col_nb(1+15,npart,1+15*npart)
       do k=1, npart
           call io_minmax_diffalp%set_pref2dline(k,tabci(k)%m_ptr%m_id) 
       enddo
       len_minmax = par_minmax_diffalp%m_stepout                                  &
     &                  /par_minmax_diffalp%m_stepcalc
       call comp_minmax_diffalp%set_stepout(len_minmax)

       p1 = minmax_diffalp_pla(1)
       do k=1, npart
            angle(3*(k-1)+1) = .false. ! pour a 
            angle(3*(k-1)+2) = .true. ! pour la
            angle(3*(k-1)+3) = .true. ! pour pi
            ! pour la  planete
            comp_p2(3*(k-1)+1) = (p1-1)*6+1 ! pour a
            comp_p2(3*(k-1)+2) = (p1-1)*6+4 ! pour la
            comp_p2(3*(k-1)+3) = (p1-1)*6+5 ! pour pi
            ! pour la particule
            comp_p1(3*(k-1)+1) = 6*nplan+(k-1)*6+1 ! pour a
            comp_p1(3*(k-1)+2) = 6*nplan+(k-1)*6+4 ! pour la
            comp_p1(3*(k-1)+3) = 6*nplan+(k-1)*6+5 ! pour pi
       enddo

       call comp_minmax_diffalp%set_component(comp_p1,comp_p2,angle,             &
     &        6*ntot)
       call comp_minmax_diffalp%set_output(1_8,                                  &
     &       io_minmax_diffalp%getconsumer())
       call comp_minmax_diffalp%set_graphdotname(                                &
     &     "min/max a_partk-a_plan,.")
      
       call PhVb2ell_mda%set_kind(par_minmax_diffalp%m_elltype)      
       call PhVb2ell_mda%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_mda%set_nb(npart)      
       call PhVb2ell_mda%set_output(1_8,                                 &
     &     comp_minmax_diffalp%getconsumer())
      end subroutine init_minmax_diffalp

!------------------------------------------------------------
!> initialise l'analyse en frequence en alkhqp des particules
!!  precise quelles composantes sont utilisees
!------------------------------------------------------------
      subroutine init_naf_alkhqp(par_naf_alkhqp, comp_naf_alkhqp,        &
     &               io_naf_alkhqp, sys, psyspartH, cG, tabci,           &
     &              PhVb2ell_nafk)
        use mod_filepar
        use mod_converter
        use mod_naf_base
        use mod_converterpart 
        use mod_sysplaexth
        use mod_syspartnewtH
      implicit none
      type(t_naf_par), intent(inout) :: par_naf_alkhqp !< donnee du fichier de parametres pour l'analyse en frequence en a,l,k,h,q,p
      type(t_naf_base), intent(inout) :: comp_naf_alkhqp !< module de calcul de l'analyse en frequence en a,l,k,h,q,p
      type (t_io_txt_mligncomp), intent(inout) :: io_naf_alkhqp  !< sortie de l'analyse en frequence en a,l,k,h,q,p.
      class(t_sysplaexth), intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_part_PhVb2ell), intent(inout) :: PhVb2ell_nafk !< buffer des elements elliptiques en entree
      type(t_ci_base_ptr), dimension(:), intent(inout) ::   tabci !< tableau de conditions initiales
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
      class(t_syspartnewtH), intent(inout) :: psyspartH !< systemes des particules

      integer(8) len_naf
      integer nplan, ncompout, npart, ntot
      integer ncompinput ! nombre de composante par planete
      integer ncompouttotal ! nombre de composante total en sortie
      integer elltypenaf ! type de ell pour naf
      integer j,k
      integer, dimension(:), allocatable :: compinput, comp_err
      
       nplan = sys%m_plan_nb
       npart = psyspartH%get_nb()
       ntot = nplan+npart
       
       ncompout = npart*par_naf_alkhqp%m_naf_nterm
       select case(par_naf_alkhqp%m_compute)
        case (1) 
          ncompinput = 3  ! cas standard a*exp(i*l), k+i*h et q+i*p
        case (2) 
          ncompinput = 2  ! cas standard a*exp(i*l), k+i*h 
        case default 
         stop "init_naf_alkhqp  erreur m_compute"
       end select 
       ncompouttotal = ncompinput*ncompout
       
       call io_naf_alkhqp%set_col_nb(1+3*ncompouttotal/npart,               &
     &               npart, 1+3*ncompouttotal)
        do k=1, npart
           call io_naf_alkhqp%set_pref2dline(k,tabci(k)%m_ptr%m_id) 
       enddo
       
       len_naf = par_naf_alkhqp%m_stepout/par_naf_alkhqp%m_stepcalc
       call comp_naf_alkhqp%set_stepout(len_naf)
       
       ! pour chaque particule 
       allocate(compinput(1:npart*ncompinput))
       allocate(comp_err(1:npart*ncompinput))
       do j=1, npart
        do k=1, ncompinput
         compinput((j-1)*ncompinput+k) = 6*nplan+(j-1)*6+2*(k-1)+1
         comp_err((j-1)*ncompinput+k) = j
        enddo
       enddo
       write(*,*) 'compinput naf=',compinput
       
       call comp_naf_alkhqp%set_component(compinput, 6*ntot)
       call comp_naf_alkhqp%set_multictrlerr(comp_err)
       deallocate(compinput)
       deallocate(comp_err)
       call comp_naf_alkhqp%set_param_naf(par_naf_alkhqp%m_naf_nterm,     &
     &   par_naf_alkhqp%m_naf_iw, par_naf_alkhqp%m_naf_isec,              &
     &   par_naf_alkhqp%m_naf_dtour, par_naf_alkhqp%m_naf_tol)
       call comp_naf_alkhqp%set_output(1_8,io_naf_alkhqp%getconsumer())
       call comp_naf_alkhqp%set_graphdotname("NAF alkhqp")
      
       ! conversion des elements d'entree dans le bon type
       elltypenaf = 0
       select case(par_naf_alkhqp%m_elltype)
        case (3) 
         elltypenaf = 7
        case (4) 
         elltypenaf = 8
        case default 
         stop "init_naf_alkhqp  erreur elltype"
       end select 
       call PhVb2ell_nafk%set_kind(elltypenaf)      
       call PhVb2ell_nafk%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_nafk%set_nb(npart)      
       call PhVb2ell_nafk%set_output(len_naf,                               &
     &       comp_naf_alkhqp%getconsumer())

      end subroutine init_naf_alkhqp

!------------------------------------------------------------
!> initialise l'analyse en frequence en difference d'angles des particules
!!  precise quelles composantes sont utilisees
!------------------------------------------------------------
      subroutine init_naf_diffalp(par_naf_diffalp,                               &
     &              comp_naf_diffalp,  io_naf_diffalp,                           &
     &              naf_diffalp_pla, sys, psyspartH, cG, tabci,                  &
     &              PhVb2ell_nafd)
        use mod_filepar
        use mod_converter
        use mod_naf_diffang
        use mod_converterpart 
        use mod_sysplaexth
        use mod_syspartnewtH
      implicit none
      type(t_naf_par), intent(inout) :: par_naf_diffalp !< donnee du fichier de parametres pour l'analyse en frequence en a,l,k,h,q,p
      type(t_naf_diffang), intent(inout) :: comp_naf_diffalp !< module de calcul de l'analyse en frequence en a,l,k,h,q,p
      type (t_io_txt_mligncomp), intent(inout) :: io_naf_diffalp  !< sortie de l'analyse en frequence en a,l,k,h,q,p.
      integer,dimension(1:1), intent(in) :: naf_diffalp_pla !< numero des 2 planetes utilisees pour naf_diffalp
      class(t_sysplaexth), intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_part_PhVb2ell), intent(inout) :: PhVb2ell_nafd !< buffer des elements elliptiques en entree
      type(t_ci_base_ptr), dimension(:), intent(inout) ::   tabci !< tableau de conditions initiales
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
      class(t_syspartnewtH), intent(inout) :: psyspartH !< systemes des particules

      integer(8) len_naf
      integer, dimension(1:2*psyspartH%get_nb()) :: comp_p1, comp_p2
      integer, dimension(1:2*psyspartH%get_nb()) :: comp_err
      integer p1, k
      integer nplan, ncompout, npart, ntot
      integer ncompinput ! nombre de composante par planete
      integer ncompouttotal ! nombre de composante total en sortie

       nplan = sys%m_plan_nb
       npart = psyspartH%get_nb()
       ntot = nplan+npart

       ncompout = npart*par_naf_diffalp%m_naf_nterm
       ncompinput = 2  ! cas standard exp(i*(l1-l2)), exp(i*(pi1-pi2))
       ncompouttotal = ncompinput*ncompout
       
       call io_naf_diffalp%set_col_nb(1+3*ncompouttotal/npart,            &
     &               npart, 1+3*ncompouttotal)
       do k=1, npart
           call io_naf_diffalp%set_pref2dline(k,tabci(k)%m_ptr%m_id) 
       enddo
       len_naf = par_naf_diffalp%m_stepout/par_naf_diffalp%m_stepcalc
       call comp_naf_diffalp%set_stepout(len_naf)
       
       p1 = naf_diffalp_pla(1)
       do k=1, npart
            ! pour la  planete
            comp_p2(2*(k-1)+1) = (p1-1)*6+4 ! pour la
            comp_p2(2*(k-1)+2) = (p1-1)*6+5 ! pour pi
            ! pour la particule
            comp_p1(2*(k-1)+1) = 6*nplan+(k-1)*6+4 ! pour la
            comp_p1(2*(k-1)+2) = 6*nplan+(k-1)*6+5 ! pour pi
            comp_err(2*(k-1)+1) = k
            comp_err(2*(k-1)+2) = k
       enddo

       call comp_naf_diffalp%set_component12(comp_p1, comp_p2, 6*ntot)
       call comp_naf_diffalp%set_multictrlerr(comp_err)
       call comp_naf_diffalp%set_param_naf(par_naf_diffalp%m_naf_nterm,     &
     &   par_naf_diffalp%m_naf_iw, par_naf_diffalp%m_naf_isec,              &
     &   par_naf_diffalp%m_naf_dtour, par_naf_diffalp%m_naf_tol)
       call comp_naf_diffalp%set_output(1_8,                                &
     &   io_naf_diffalp%getconsumer())
       call comp_naf_diffalp%set_graphdotname(                              &
     &    "NAF exp(i*(l_partk-l_plan)),...")
      
       call PhVb2ell_nafd%set_kind(par_naf_diffalp%m_elltype)      
       call PhVb2ell_nafd%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_nafd%set_nb(npart)      
       call PhVb2ell_nafd%set_output(len_naf,                                 &
     &       comp_naf_diffalp%getconsumer())
      
      end subroutine init_naf_diffalp
      
!------------------------------------------------------------
!> initialise l'analyse en frequence en alkhqp pour les planetes
!!  precise quelles composantes sont utilisees
!------------------------------------------------------------
      subroutine init_naf_alkhqp_pla(par_naf_alkhqp, comp_naf_alkhqp,        &
     &               io_naf_alkhqp, sys, cG, id_ci, PhVb2ell_nafk)
        use mod_filepar
        use mod_syspla
        use mod_converter
        use mod_naf_base
      implicit none
      type(t_naf_par), intent(inout) :: par_naf_alkhqp !< donnee du fichier de parametres pour l'analyse en frequence en a,l,k,h,q,p
      type(t_naf_base), intent(inout) :: comp_naf_alkhqp !< module de calcul de l'analyse en frequence en a,l,k,h,q,p
      type (t_io_txt_ncomp), intent(inout) :: io_naf_alkhqp  !< sortie de l'analyse en frequence en a,l,k,h,q,p.
      class(t_syspla), intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_plan_PhVb2ell), intent(inout) :: PhVb2ell_nafk !< buffer des elements elliptiques en entree
      character(len=*), intent(in) :: id_ci  !< identifiant de la condition initiale
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an

      integer(8) len_naf
      integer nplan, ncompout
      integer ncompinput ! nombre de composante par planete
      integer ncompouttotal ! nombre de composante total en sortie
      integer elltypenaf ! type de ell pour naf
      integer j,k
      integer, dimension(:), allocatable :: compinput
      
       nplan = sys%m_plan_nb
       ncompout = nplan*par_naf_alkhqp%m_naf_nterm
       select case(par_naf_alkhqp%m_compute)
        case (1) 
          ncompinput = 3  ! cas standard a*exp(i*l), k+i*h et q+i*p
        case (2) 
          ncompinput = 2  ! cas standard a*exp(i*l), k+i*h 
        case default 
         stop "init_naf_alkhqp  erreur m_compute"
       end select 
       ncompouttotal = ncompinput*ncompout
       
       call io_naf_alkhqp%set_comp_nb(1+ ncompouttotal*3)
       call io_naf_alkhqp%set_prefixline(id_ci)
       len_naf = par_naf_alkhqp%m_stepout/par_naf_alkhqp%m_stepcalc
       call comp_naf_alkhqp%set_stepout(len_naf)
       
       ! pour chaque planete 
       allocate(compinput(1:nplan*ncompinput))
       do j=1, nplan
        do k=1, ncompinput
         compinput((j-1)*ncompinput+k) = (j-1)*6+2*(k-1)+1
        enddo
       enddo
       write(*,*) 'compinput=',compinput
       
       call comp_naf_alkhqp%set_component(compinput, 6*nplan)
       deallocate(compinput)
       call comp_naf_alkhqp%set_param_naf(par_naf_alkhqp%m_naf_nterm,     &
     &   par_naf_alkhqp%m_naf_iw, par_naf_alkhqp%m_naf_isec,              &
     &   par_naf_alkhqp%m_naf_dtour, par_naf_alkhqp%m_naf_tol)
       call comp_naf_alkhqp%set_output(1_8,io_naf_alkhqp%getconsumer())
       call comp_naf_alkhqp%set_graphdotname("NAF alkhqp pla")
      
       ! conversion des elements d'entree dans le bon type
       elltypenaf = 0
       select case(par_naf_alkhqp%m_elltype)
        case (3) 
         elltypenaf = 7
        case (4) 
         elltypenaf = 8
        case default 
         stop "init_naf_alkhqp_pla  erreur elltype"
       end select 
       call PhVb2ell_nafk%set_kind(elltypenaf)      
       call PhVb2ell_nafk%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_nafk%set_output(len_naf,                               &
     &       comp_naf_alkhqp%getconsumer())
     
      end subroutine init_naf_alkhqp_pla

      end 


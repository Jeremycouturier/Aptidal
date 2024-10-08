!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file intgsympstat.f 
!!  \brief  routine d'integration : construit le graphe de la simulation 
!! et lance l'integration 
!!
! history : creation 17/03/2014
!***********************************************************************
#include "realc.f"
      module mod_intg_symp
      
       private init_minmax_diffalp
       private init_naf_alkhqp
       private init_naf_diffalp
       private init_naf_difflpm
      contains


!------------------------------------------------------------
!> routine d'integration
!------------------------------------------------------------
      subroutine intg_symp(filepar, sys, cG, ctrlpladmin, id_ci, i, si)
#if USE_MPI
      use mpi
#endif
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
      use mod_minmax
      use mod_minmax_mda
      use mod_minmax_mdc
      use mod_minmax_dae2
      use mod_naf_base
      use mod_naf_diffang
      use mod_naf_difflpm
      use mod_ctrl_ener
      use mod_ctrl_diststar
      use mod_ctrl_distpla
      use mod_io_dump

      implicit none
      type (t_extfilepar),target, intent(inout) :: filepar !< fichier de parametres + variables derivees
      class(t_syspla), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
      character(len=*), intent(in) :: id_ci  !< identifiant de la condition initiale
      real(TREAL), dimension(:), pointer, intent(in) :: ctrlpladmin ! distance minimale entre planetes
      integer(8), intent(in) :: i !< iteration
      integer(8), intent(in) :: si !< iteration fixe
     
      integer(8) :: n_iter
      integer(8) :: n_out
      real(TREAL) :: tinit,dt
      integer :: nplan, k
      
      type (t_plan_PhVb2PhVh) :: PhVb2PhVh
      type (t_plan_PjVj2PhVh) :: PjVj2PhVh
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb, PjVj2PhVb_2
      type (t_plan_PhVb2ell) :: PhVb2ell, PhVb2ell_2
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_mda ! pour les min/max des differences d'angle entre [-pi,pi] et [0,2*pi]
      type (t_plan_PhVb2ell) :: PhVb2ell_mda ! pour les min/max des differences d'angle  [-pi,pi] et [0,2*pi]
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_mdc ! pour les min/max des differences d'angle continu
      type (t_plan_PhVb2ell) :: PhVb2ell_mdc ! pour les min/max des differences d'angle  continu
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_ae2 ! pour les min/max des differences de a**2,e**2,....
      type (t_plan_PhVb2ell) :: PhVb2ell_ae2 ! pour les min/max des differences de a**2,e**2,....
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_nafk ! pour naf en alkhqp
      type (t_plan_PhVb2ell) :: PhVb2ell_nafk ! pour naf en alkhqp
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_nafd ! pour naf des differences d'angle en exp(i*....)
      type (t_plan_PhVb2ell) :: PhVb2ell_nafd ! pour naf des differences d'angle en exp(i*....)
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_nafm ! pour naf des differences d'angle en (...) mod 2pi ou +-pi
      type (t_plan_PhVb2ell) :: PhVb2ell_nafm ! pour naf des differences d'angle en (...) mod 2pi ou +-pi
    
      type (t_minmax) :: minmax_aei
      integer(8) :: len_minmax_aei
      integer, dimension(:), allocatable :: aeicomp
      
      type(t_minmax_mda) :: minmax_diffalp ! min/moy/max en l1-l2,pibar1-pibar2, a1-a2 (angle entre [-pi,pi] et [0,2*pi])
      type(t_minmax_mdc) :: minmax_diffalc ! min/moy/max en l1-l2,pibar1-pibar2, a1-a2 (angle continu)
      type(t_minmax_dae2) :: minmax_diffae2 ! min/moy/max en (a1-a2)**2,(e1-e2)**2, (a1-a1)**2+(e1-e2)**2


      type (t_naf_base) :: naf_alkhqp ! naf en alkhqp
      type (t_naf_diffang) :: naf_diffalp ! naf en differences d'angle exp i*....
      type (t_naf_difflpm) :: naf_difflpm ! naf en differences d'angle ... mod 2*pi
      integer(8) :: len_nafk, len_nafd, len_nafm ! longueur des tranches de naf

      type (t_io_int) ::io_int  ! sortie des integrales premieres
      type (t_io_txt_car) :: io_cart  ! sortie des elements cartesiens
      type (t_io_txt_ncomp) :: io_ell ! sortie des elements elliptiques
      type (t_io_dump) :: io_dump ! dump de l'integration pour un redemarrage
      
      type (t_ctrl_diststar) :: ctrl_diststar ! controle des distances avec l'etoile
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_ds ! pour controle des distances a l'etoile
      type (t_plan_PhVb2PhVh) :: PhVb2PhVh_ds ! pour dump du controle des distances a l'etoile
      type (t_plan_PhVb2ell) :: PhVb2ell_ds   ! pour dump du controle des distances a l'etoile

      type (t_ctrl_distpla) :: ctrl_distpla ! controle des distances entre planetes
      type (t_plan_PjVj2PhVb) :: PjVj2PhVb_dp ! pour controle des distances entre planetes
      type (t_plan_PhVb2PhVh) :: PhVb2PhVh_dp ! pour dump du controle des distances entre planetes
      type (t_plan_PhVb2ell) :: PhVb2ell_dp   ! pour dump du controle des distances entre planetes

      type (t_io_txt_ncomp) :: io_sch ! donnes supplmentaires du schema adaptatifs
      
      type (t_ctrl_ener) :: ctrl_ener ! controle de l'energie

      type(t_arret) :: arret ! code de l'arret
      type(t_intsympbase) :: integrator
      logical :: issyshelio
      character(100)  :: int_type ! type de l'integrateur
      integer graphnode ! pour la generation du graphe
      integer writegraph ! pour la generation du graphe
      integer rankmpi ! rang mpi
      integer code ! code mpi
      
      integer(8) :: len_malp
      integer(8) :: len_malc
      integer(8) :: len_mae2
      integer(8) :: len_minmax_aei_buf
      integer(8) :: bufferlimit ! taille des buffers pour les min/max
      
      bufferlimit = 1000
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
      
      if (filepar%int_var.eq.1) then
       issyshelio = .true.
      endif 
      
      write(*,*) 'int_type = ', int_type
      write(*,*) 'issyshelio = ', issyshelio

      call sys%clear_out_step()
      write (*,*) 'type_pas = ',filepar%type_pas
      call integrator%set_type_pas(filepar%type_pas)
      call integrator%set_id(id_ci)

      
      ! sortie elliptique
      if (filepar%if_ell.eq.1) then
         write(*,*)  ' nf_ell ', trim(filepar%nf_ell)
         io_ell = t_io_txt_ncomp("I/O ell")
         call io_ell%set_comp_nb(1+6*nplan)
         call io_ell%set_name(trim(filepar%nf_ell),61)
         call io_ell%fcreate()
      end if
      
      ! sortie cartesienne
      if (filepar%if_car.eq.1) then
         write(*,*) ' nf_car ', trim(filepar%nf_car)
         io_cart = t_io_txt_car("I/O cart")
         call io_cart%set_comp_nb(1+6*nplan)
         call io_cart%set_name(trim(filepar%nf_car),62)
         call io_cart%set_plan_nb(nplan)
         call io_cart%fcreate()
      end if

      ! sortie integrales premieres
      if (filepar%if_int.eq.1) then
        write(*,*)  ' nf_int ', trim(filepar%nf_int)
        call io_int%set_name(trim(filepar%nf_int),60)
        call io_int%fcreate()
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
        call io_sch%set_name(trim(filepar%nf_sch),57)
        call io_sch%fcreate()
        select type(sys)
         class is(t_sysplanewtJ_adapt)
           call sys%add_adapt_out_step(n_out,1_8,io_sch%getconsumer()) 
         class is(t_sysplanewtH_adapt)
           call sys%add_adapt_out_step(n_out,1_8,io_sch%getconsumer()) 
        end select    
      endif
      
      ! sortie min/moy/max a,e,i
      if (filepar%minmax_aei%m_compute.eq.1) then
       call filepar%io_minmax_aei%set_comp_nb(1+9*nplan)
       call filepar%io_minmax_aei%set_prefixline(id_ci)
       len_minmax_aei = filepar%minmax_aei%m_stepout                    &
     &                  /filepar%minmax_aei%m_stepcalc
       call minmax_aei%set_stepout(len_minmax_aei)
       len_minmax_aei_buf = len_minmax_aei
       if (len_minmax_aei_buf.ge.bufferlimit) then
         len_minmax_aei_buf = bufferlimit
       endif
       allocate(aeicomp(1:3*nplan))
       do k=0, nplan-1
        aeicomp(3*k+1) = 6*k+1 ! pour a
        aeicomp(3*k+2) = 6*k+2 ! pour e
        aeicomp(3*k+3) = 6*k+3 ! pour I
       enddo
       call minmax_aei%set_component(aeicomp,6*nplan)
       deallocate(aeicomp)
       call minmax_aei%set_output(1_8,                                    &
     &       filepar%io_minmax_aei%getconsumer())
       call minmax_aei%set_graphdotname("min/max a,e,i")
      
       call PhVb2ell_2%set_kind(filepar%minmax_aei%m_elltype)      
       call PhVb2ell_2%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_2%set_output(1_8, minmax_aei%getconsumer())
      endif

      ! sortie min/moy/max difference d'angle entre [-pi,pi] et [0, 2*pi]
      if (filepar%minmax_diffalp%m_compute.eq.1) then
        len_malp = filepar%minmax_diffalp%m_stepout                     &
     &                  /filepar%minmax_diffalp%m_stepcalc
        if (len_malp.ge.bufferlimit) then
         len_malp = bufferlimit
        endif
        call init_minmax_diffalp(filepar%minmax_diffalp,minmax_diffalp, &
     &         filepar%io_minmax_diffalp, filepar%minmax_diffalp_pla,   &
     &         sys, cG, id_ci, PhVb2ell_mda)
      endif

      ! sortie min/moy/max difference d'angle continu
      if (filepar%minmax_diffalc%m_compute.eq.1) then
        len_malc = filepar%minmax_diffalc%m_stepout                     &
     &                  /filepar%minmax_diffalc%m_stepcalc
        if (len_malc.ge.bufferlimit) then
         len_malc = bufferlimit
        endif
        call init_minmax_diffalc(filepar%minmax_diffalc,minmax_diffalc, &
     &         filepar%io_minmax_diffalc, filepar%minmax_diffalc_pla,   &
     &         sys, cG, id_ci, PhVb2ell_mdc)
      endif

      ! sortie min/moy/max difference en a**2,e**2,a**2+e**2
      if (filepar%minmax_diffae2%m_compute.eq.1) then
        len_mae2 = filepar%minmax_diffae2%m_stepout                     &
     &                  /filepar%minmax_diffae2%m_stepcalc
        if (len_mae2.ge.bufferlimit) then
         len_mae2 = bufferlimit
        endif
       call init_minmax_diffae2(filepar%minmax_diffae2,minmax_diffae2,   &
     &         filepar%io_minmax_diffae2, filepar%minmax_diffae2_pla,   &
     &         sys, cG, id_ci, PhVb2ell_ae2)
      endif

      ! sortie naf en alkhqp
      if (filepar%naf_alkhqp%m_compute.ne.0) then
        len_nafk = filepar%naf_alkhqp%m_stepout                         &
     &                  /filepar%naf_alkhqp%m_stepcalc
        call init_naf_alkhqp(filepar%naf_alkhqp,naf_alkhqp,             &
     &         filepar%io_naf_alkhqp,                                   &
     &         sys, cG, id_ci, PhVb2ell_nafk)
      endif

      ! sortie naf en difference d'angle : exp(i*....)
      if (filepar%naf_diffalp%m_compute.eq.1) then
        len_nafd = filepar%naf_diffalp%m_stepout                        &
     &                  /filepar%naf_diffalp%m_stepcalc
        call init_naf_diffalp(filepar%naf_diffalp,naf_diffalp,          &
     &         filepar%io_naf_diffalp, filepar%naf_diffalp_pla,         &
     &         sys, cG, id_ci, PhVb2ell_nafd)
      endif

      ! sortie naf en difference d'angle : (...+i*0) mod2pi ou +-pi
      if (filepar%naf_difflpm%m_compute.eq.1) then
        len_nafm = filepar%naf_difflpm%m_stepout                        &
     &                  /filepar%naf_difflpm%m_stepcalc
        call init_naf_difflpm(filepar%naf_difflpm,naf_difflpm,          &
     &         filepar%io_naf_difflpm, filepar%naf_difflpm_pla,         &
     &         sys, cG, id_ci, PhVb2ell_nafm)
      endif


      ! controle des distances avec l'etoile
      if (filepar%ctrl_diststar%m_compute.ge.1) then
        call ctrl_diststar%set_minmax(nplan,                            &
     &         filepar%ctrl_diststar%m_distmin,                         &
     &         filepar%ctrl_diststar%m_distmax)
     
      ! dump du controle de distance a l'etoile
       select case (filepar%ctrl_diststar%m_compute)  
       case (0,1)
       case (2) ! sortie cartesien
       call filepar%io_ctrlstar_dump%set_comp_nb(1+6*nplan)
       call filepar%io_ctrlstar_dump%set_prefixline(id_ci)
       call PhVb2PhVh_ds%set_mass(sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2PhVh_ds%set_output(1_8,                                 &
     &          filepar%io_ctrlstar_dump%getconsumer())
       call ctrl_diststar%set_output(PhVb2PhVh_ds%getconsumer())

       case (3) ! sortie elliptique
       call filepar%io_ctrlstar_dump%set_comp_nb(1+6*nplan)
       call filepar%io_ctrlstar_dump%set_prefixline(id_ci)
       call PhVb2ell_ds%set_kind(filepar%out_ell)      
       call PhVb2ell_ds%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_ds%set_output(1_8,                                  &
     &           filepar%io_ctrlstar_dump%getconsumer())
       call ctrl_diststar%set_output(PhVb2ell_ds%getconsumer())

        case default
          stop "intg_symp : ctrl_diststar_compute inconnu"
       end select

      endif
      
      ! controle des distances entre planetes
      if (filepar%ctrl_distpla%m_compute.ne.0) then
        call ctrl_distpla%set_min(nplan,ctrlpladmin)
     
      ! dump du controle de distance entre planetes
       select case (filepar%ctrl_distpla%m_compute)  
       case (0,1)
       case (2) ! sortie cartesien
       call filepar%io_ctrlpla_dump%set_comp_nb(1+6*nplan)
       call filepar%io_ctrlpla_dump%set_prefixline(id_ci)
       call PhVb2PhVh_dp%set_mass(sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2PhVh_dp%set_output(1_8,                                 &
     &          filepar%io_ctrlpla_dump%getconsumer())
       call ctrl_distpla%set_output(PhVb2PhVh_dp%getconsumer())

       case (3) ! sortie elliptique
       call filepar%io_ctrlpla_dump%set_comp_nb(1+6*nplan)
       call filepar%io_ctrlpla_dump%set_prefixline(id_ci)
       call PhVb2ell_dp%set_kind(filepar%out_ell)      
       call PhVb2ell_dp%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_dp%set_output(1_8,                                  &
     &           filepar%io_ctrlpla_dump%getconsumer())
       call ctrl_distpla%set_output(PhVb2ell_dp%getconsumer())

        case default
          stop "intg_symp : ctrl_distpla_compute inconnu"
       end select

      endif

      ! controle de l'energie
      if (filepar%ctrl_energie%m_compute.eq.1) then
       call ctrl_ener%set_relenermax(filepar%ctrl_energie%m_relenermax)
       call sys%add_int_out_step(filepar%ctrl_energie%m_stepcalc,1_8,   &
     &         ctrl_ener%getconsumer())  
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
      
      ! sortie min/max a,e,i
      if (filepar%minmax_aei%m_compute.eq.1) then
       call sys%add_plan_out_step(filepar%minmax_aei%m_stepcalc,         &
     &          len_minmax_aei_buf, PhVb2ell_2%getconsumer())
      endif
      
      ! sortie min/moy/max difference d'angle entre [-pi,pi] et [0,2*pi]
      if (filepar%minmax_diffalp%m_compute.eq.1) then
       call sys%add_plan_out_step(filepar%minmax_diffalp%m_stepcalc,     &
     &          len_malp, PhVb2ell_mda%getconsumer())
      endif

      ! sortie min/moy/max difference d'angle continu
      if (filepar%minmax_diffalc%m_compute.eq.1) then
       call sys%add_plan_out_step(filepar%minmax_diffalc%m_stepcalc,     &
     &          len_malc, PhVb2ell_mdc%getconsumer())
      endif

      ! sortie min/moy/max difference de a**2,e**2, a**2+e**2
      if (filepar%minmax_diffae2%m_compute.eq.1) then
       call sys%add_plan_out_step(filepar%minmax_diffae2%m_stepcalc,     &
     &          len_mae2, PhVb2ell_ae2%getconsumer())
      endif

      ! sortie naf alkhqp
      if (filepar%naf_alkhqp%m_compute.ne.0) then
       call sys%add_plan_out_step(filepar%naf_alkhqp%m_stepcalc,         &
     &          len_nafk, PhVb2ell_nafk%getconsumer())
      endif
      
      ! sortie naf difference d'angle en exp(i*....)
      if (filepar%naf_diffalp%m_compute.eq.1) then
       call sys%add_plan_out_step(filepar%naf_diffalp%m_stepcalc,        &
     &          len_nafd, PhVb2ell_nafd%getconsumer())
      endif

      ! sortie naf difference d'angle en (...+i*0) mod 2pi ou +-pi
      if (filepar%naf_difflpm%m_compute.eq.1) then
       call sys%add_plan_out_step(filepar%naf_difflpm%m_stepcalc,        &
     &          len_nafm, PhVb2ell_nafm%getconsumer())
      endif

      ! controle des distances avec l'etoile
      if (filepar%ctrl_diststar%m_compute.ge.1) then
       call sys%add_plan_out_step(filepar%ctrl_diststar%m_stepcalc,      &
     &          1_8, ctrl_diststar%getconsumer())
      endif

      ! controle des distances entre planetes
      if (filepar%ctrl_distpla%m_compute.ne.0) then
       call sys%add_plan_out_stepmid(filepar%ctrl_distpla%m_stepcalc,       &
     &          1_8, ctrl_distpla%getconsumer())
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

      ! sortie min/max a,e,i
      if (filepar%minmax_aei%m_compute.eq.1) then
        call PjVj2PhVb_2%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_2%set_output(len_malp,PhVb2ell_2%getconsumer())
        call sys%add_plan_out_step(filepar%minmax_aei%m_stepcalc,       &
     &          len_minmax_aei_buf, PjVj2PhVb_2%getconsumer())
      endif

      ! sortie min/moy/max difference d'angle entre [-pi,pi] et [0, 2*pi]
      if (filepar%minmax_diffalp%m_compute.eq.1) then
        call PjVj2PhVb_mda%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_mda%set_output(1_8,PhVb2ell_mda%getconsumer())
        call sys%add_plan_out_step(filepar%minmax_diffalp%m_stepcalc,   &
     &          len_malp, PjVj2PhVb_mda%getconsumer())
      endif

      ! sortie min/moy/max difference d'angle continu
      if (filepar%minmax_diffalc%m_compute.eq.1) then
        call PjVj2PhVb_mdc%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_mdc%set_output(1_8,PhVb2ell_mdc%getconsumer())
        call sys%add_plan_out_step(filepar%minmax_diffalc%m_stepcalc,   &
     &          len_malc, PjVj2PhVb_mdc%getconsumer())
      endif

      ! sortie min/moy/max difference de a**2,e**2, a**2+e**2
      if (filepar%minmax_diffae2%m_compute.eq.1) then
        call PjVj2PhVb_ae2%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_ae2%set_output(1_8,PhVb2ell_ae2%getconsumer())
        call sys%add_plan_out_step(filepar%minmax_diffae2%m_stepcalc,   &
     &          len_mae2, PjVj2PhVb_ae2%getconsumer())
      endif

      ! sortie naf alkhqp
      if (filepar%naf_alkhqp%m_compute.ne.0) then
        call PjVj2PhVb_nafk%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_nafk%set_output(len_nafk,                        &
     &          PhVb2ell_nafk%getconsumer())
        call sys%add_plan_out_step(filepar%naf_alkhqp%m_stepcalc,       &
     &          len_nafk, PjVj2PhVb_nafk%getconsumer())
      endif
      
      ! sortie naf difference d'angle en exp(i*(....))
      if (filepar%naf_diffalp%m_compute.eq.1) then
        call PjVj2PhVb_nafd%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_nafd%set_output(len_nafd,                        &
     &          PhVb2ell_nafd%getconsumer())
        call sys%add_plan_out_step(filepar%naf_diffalp%m_stepcalc,      &
     &          len_nafd, PjVj2PhVb_nafd%getconsumer())
      endif

      ! sortie naf difference d'angle en (...)+i*0 mod 2pi ou +-pi
      if (filepar%naf_difflpm%m_compute.eq.1) then
        call PjVj2PhVb_nafm%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_nafm%set_output(len_nafm,                        &
     &          PhVb2ell_nafm%getconsumer())
        call sys%add_plan_out_step(filepar%naf_difflpm%m_stepcalc,      &
     &          len_nafm, PjVj2PhVb_nafm%getconsumer())
      endif


      ! controle des distances avec l'etoile
      if (filepar%ctrl_diststar%m_compute.ge.1) then
        call PjVj2PhVb_ds%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_ds%set_output(1_8, ctrl_diststar%getconsumer())
        call sys%add_plan_out_step(filepar%ctrl_diststar%m_stepcalc,    &
     &          1_8, PjVj2PhVb_ds%getconsumer())
      endif

      ! controle des distances entre planetes
      if (filepar%ctrl_distpla%m_compute.ne.0) then
        call PjVj2PhVb_dp%set_mass(sys%m_star_mass,sys%m_plan_mass)      
        call PjVj2PhVb_dp%set_output(1_8, ctrl_distpla%getconsumer())
        call sys%add_plan_out_stepmid(filepar%ctrl_distpla%m_stepcalc,     &
     &          1_8, PjVj2PhVb_dp%getconsumer())
      endif

      endif
      !--------------------------------
      !--------------------------------
      ! fin du code dependant de l'integrateur
      !--------------------------------
      !--------------------------------

      ! controle du dump de l'integration
      if (filepar%if_dump.ne.0) then
         write(*,*)  ' nf_dump ', trim(filepar%nf_dump)
         call io_dump%set_name(trim(filepar%nf_dump),58)
         call io_dump%set_output(filepar%n_dump)
         call io_dump%fcreate()
         call integrator%set_dump(io_dump)
      end if


      ! lance l'integration
      call integrator%set_error_behavior(1)
      call integrator%set_syssymp(sys)
      call integrator%set_integ_scheme(int_type) 
      call integrator%set_integ_time(tinit, dt, i)
      call integrator%set_integ_time_si(si)
      call integrator%integ_run(n_iter)
      
      ! ecriture du controle de l'integration
      if (sys%check_error()) then
       arret=sys%get_error()
       write(*,*) 'une erreur s est produite'
       write(*,*) 'code :', arret%m_cause
       write(*,*) 'time :', arret%m_time
       write(*,*) 'body :', arret%m_body
       write(*,*) 'value :', arret%m_value
       write(filepar%io_control,1004) trim(id_ci),arret%m_cause,           &
     &          tinit, arret%m_time, arret%m_body, arret%m_value
      else
       write(filepar%io_control,1004) trim(id_ci),  0, tinit,              &
     &          integrator%get_integ_current_time(), 0, 0.D0
      endif 

1004  format(A,1X,I2,1X,FMT_TREAL,1X,FMT_TREAL,I4,1X,FMT_TREAL)  

      ! genere le graphe d'appel 
#if USE_MPI
      if (rankmpi.eq.1) then 
#else
      if (rankmpi.eq.0) then 
#endif
       writegraph = 1
       if (writegraph.eq.1) then
        write(6,*) "digraph g {"
        graphnode = 0
        call sys%graphdot(0, graphnode, 6)
        write(6,*) "}"
       endif
      endif
      
      end subroutine

!----------------------------------------------------------------------------------------
!> initialise l'analyse en min/max/moy des differences d'angle entre [-pi,pi] et [0,2*pi]
!!  precise quelles composantes sont utilisees
!----------------------------------------------------------------------------------------
      subroutine init_minmax_diffalp(par_minmax_diffalp,                    &
     &              comp_minmax_diffalp,  io_minmax_diffalp,                &
     &              minmax_diffalp_pla, sys, cG, id_ci,  PhVb2ell_mda)
        use mod_filepar
        use mod_syspla
        use mod_converter
        use mod_minmax_mda
      implicit none
      type(t_minmax_par), intent(inout) :: par_minmax_diffalp !< donnee du fichier de parametres pour min/moy/max en l1-l2,pibar1-pibar2, a1-a2
      type(t_minmax_mda), intent(inout) :: comp_minmax_diffalp !< module de calcul des min/moy/max en l1-l2,pibar1-pibar2, a1-a2
      type (t_io_txt_ncomp), intent(inout) :: io_minmax_diffalp  !< sortie des min/max des differences a_p1-a_p2,....
      integer,dimension(1:2), intent(in) :: minmax_diffalp_pla !< numero des 2 planetes utilisees pour minmax_diffalp
      class(t_syspla), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_plan_PhVb2ell), intent(inout) :: PhVb2ell_mda !< buffer des elements elliptiques en entree
      character(len=*), intent(in) :: id_ci  !< identifiant de la condition initiale
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an

      integer(8) len_minmax
      integer, dimension(1:3) :: comp_p1, comp_p2
      logical, dimension(1:3) :: angle
      integer p1, p2, nplan

       nplan = sys%m_plan_nb

       call io_minmax_diffalp%set_comp_nb(1+15)
       call io_minmax_diffalp%set_prefixline(id_ci)
       len_minmax = par_minmax_diffalp%m_stepout                                  &
     &                  /par_minmax_diffalp%m_stepcalc
       call comp_minmax_diffalp%set_stepout(len_minmax)
       angle(1) = .false. ! pour a 
       angle(2:3) = .true. ! pour la et pi
       ! premiere planete
       p1 = minmax_diffalp_pla(1)
       comp_p1(1) = (p1-1)*6+1 ! pour a
       comp_p1(2) = (p1-1)*6+4 ! pour la
       comp_p1(3) = (p1-1)*6+5 ! pour pi
       ! deuxieme planete
       p2 = minmax_diffalp_pla(2)
       comp_p2(1) = (p2-1)*6+1 ! pour a
       comp_p2(2) = (p2-1)*6+4 ! pour la
       comp_p2(3) = (p2-1)*6+5 ! pour pi
       call comp_minmax_diffalp%set_component(comp_p1,comp_p2,angle,             &
     &        6*nplan)
       call comp_minmax_diffalp%set_output(1_8,                                  &
     &       io_minmax_diffalp%getconsumer())
       call comp_minmax_diffalp%set_graphdotname(                                &
     &   "min/max a_p1-a_p2,...")
      
       call PhVb2ell_mda%set_kind(par_minmax_diffalp%m_elltype)      
       call PhVb2ell_mda%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_mda%set_output(1_8,                                         &
     &       comp_minmax_diffalp%getconsumer())
      end subroutine init_minmax_diffalp

!----------------------------------------------------------------------------------------
!> initialise l'analyse en min/max/moy des differences d'angle continu
!!  precise quelles composantes sont utilisees
!----------------------------------------------------------------------------------------
      subroutine init_minmax_diffalc(par_minmax_diffalc,                    &
     &              comp_minmax_diffalc,  io_minmax_diffalc,                &
     &              minmax_diffalc_pla, sys, cG, id_ci,  PhVb2ell_mdc)
        use mod_filepar
        use mod_syspla
        use mod_converter
        use mod_minmax_mdc
      implicit none
      type(t_minmax_par), intent(inout) :: par_minmax_diffalc !< donnee du fichier de parametres pour min/moy/max en l1-l2,pibar1-pibar2, a1-a2
      type(t_minmax_mdc), intent(inout) :: comp_minmax_diffalc !< module de calcul des min/moy/max en l1-l2,pibar1-pibar2, a1-a2
      type (t_io_txt_ncomp), intent(inout) :: io_minmax_diffalc  !< sortie des min/max des differences a_p1-a_p2,....
      integer,dimension(1:2), intent(in) :: minmax_diffalc_pla !< numero des 2 planetes utilisees pour minmax_diffalc
      class(t_syspla), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_plan_PhVb2ell), intent(inout) :: PhVb2ell_mdc !< buffer des elements elliptiques en entree
      character(len=*), intent(in) :: id_ci  !< identifiant de la condition initiale
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an

      integer(8) len_minmax
      integer, dimension(1:3) :: comp_p1, comp_p2
      logical, dimension(1:3) :: angle
      integer p1, p2, nplan

       nplan = sys%m_plan_nb

       call io_minmax_diffalc%set_comp_nb(1+9)
       call io_minmax_diffalc%set_prefixline(id_ci)
       len_minmax = par_minmax_diffalc%m_stepout                                  &
     &                  /par_minmax_diffalc%m_stepcalc
       call comp_minmax_diffalc%set_stepout(len_minmax)
       angle(1) = .false. ! pour a 
       angle(2:3) = .true. ! pour la et pi
       ! premiere planete
       p1 = minmax_diffalc_pla(1)
       comp_p1(1) = (p1-1)*6+1 ! pour a
       comp_p1(2) = (p1-1)*6+4 ! pour la
       comp_p1(3) = (p1-1)*6+5 ! pour pi
       ! deuxieme planete
       p2 = minmax_diffalc_pla(2)
       comp_p2(1) = (p2-1)*6+1 ! pour a
       comp_p2(2) = (p2-1)*6+4 ! pour la
       comp_p2(3) = (p2-1)*6+5 ! pour pi
       call comp_minmax_diffalc%set_component(comp_p1,comp_p2,angle,             &
     &        6*nplan)
       call comp_minmax_diffalc%set_output(1_8,                                  &
     &       io_minmax_diffalc%getconsumer())
       call comp_minmax_diffalc%set_graphdotname(                                &
     &   "min/max a_p1-a_p2,...")
      
       call PhVb2ell_mdc%set_kind(par_minmax_diffalc%m_elltype)      
       call PhVb2ell_mdc%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_mdc%set_output(1_8,                                         &
     &       comp_minmax_diffalc%getconsumer())
      end subroutine init_minmax_diffalc

!------------------------------------------------------------
!> initialise l'analyse en min/max/moy des differences d'a**2,e**2
!!  precise quelles composantes sont utilisees
!------------------------------------------------------------
      subroutine init_minmax_diffae2(par_minmax_diffae2,                    &
     &              comp_minmax_diffae2,  io_minmax_diffae2,                &
     &              minmax_diffae2_pla, sys, cG, id_ci,  PhVb2ell_ae2)
        use mod_filepar
        use mod_syspla
        use mod_converter
        use mod_minmax_dae2
      implicit none
      type(t_minmax_par), intent(inout) :: par_minmax_diffae2 !< donnee du fichier de parametres pour min/moy/max en (a1-a2)**2, (e1-e2)**2, (a1-a2)**2+(e1-e2)**2
      type(t_minmax_dae2), intent(inout) :: comp_minmax_diffae2 !< module de calcul des min/moy/max en (a1-a2)**2, (e1-e2)**2, (a1-a2)**2+(e1-e2)**2
      type (t_io_txt_ncomp), intent(inout) :: io_minmax_diffae2  !< sortie des min/max des differences (a1-a2)**2, (e1-e2)**2, (a1-a2)**2+(e1-e2)**2
      integer,dimension(1:2), intent(in) :: minmax_diffae2_pla !< numero des 2 planetes utilisees pour minmax_diffalp
      class(t_syspla), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_plan_PhVb2ell), intent(inout) :: PhVb2ell_ae2 !< buffer des elements elliptiques en entree
      character(len=*), intent(in) :: id_ci  !< identifiant de la condition initiale
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an

      integer(8) len_minmax
      integer, dimension(1:2) :: comp_p1, comp_p2
      integer p1, p2, nplan

       nplan = sys%m_plan_nb

       call io_minmax_diffae2%set_comp_nb(1+9)
       call io_minmax_diffae2%set_prefixline(id_ci)
       len_minmax = par_minmax_diffae2%m_stepout                                  &
     &                  /par_minmax_diffae2%m_stepcalc
       call comp_minmax_diffae2%set_stepout(len_minmax)
       ! premiere planete
       p1 = minmax_diffae2_pla(1)
       comp_p1(1) = (p1-1)*6+1 ! pour a
       comp_p1(2) = (p1-1)*6+2 ! pour e
       ! deuxieme planete
       p2 = minmax_diffae2_pla(2)
       comp_p2(1) = (p2-1)*6+1 ! pour a
       comp_p2(2) = (p2-1)*6+2 ! pour e
       call comp_minmax_diffae2%set_component(comp_p1,comp_p2,6*nplan)
       call comp_minmax_diffae2%set_output(1_8,                                  &
     &       io_minmax_diffae2%getconsumer())
       call comp_minmax_diffae2%set_graphdotname(                                &
     &   "min/max (a_p1-a_p2)**2,...")
      
       call PhVb2ell_ae2%set_kind(par_minmax_diffae2%m_elltype)      
       call PhVb2ell_ae2%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_ae2%set_output(1_8,                                  &
     &       comp_minmax_diffae2%getconsumer())
      end subroutine init_minmax_diffae2

!------------------------------------------------------------
!> initialise l'analyse en frequence en alkhqp
!!  precise quelles composantes sont utilisees
!------------------------------------------------------------
      subroutine init_naf_alkhqp(par_naf_alkhqp, comp_naf_alkhqp,        &
     &               io_naf_alkhqp, sys, cG, id_ci, PhVb2ell_nafk)
        use mod_filepar
        use mod_syspla
        use mod_converter
        use mod_naf_base
      implicit none
      type(t_naf_par), intent(inout) :: par_naf_alkhqp !< donnee du fichier de parametres pour l'analyse en frequence en a,l,k,h,q,p
      type(t_naf_base), intent(inout) :: comp_naf_alkhqp !< module de calcul de l'analyse en frequence en a,l,k,h,q,p
      type (t_io_txt_ncomp), intent(inout) :: io_naf_alkhqp  !< sortie de l'analyse en frequence en a,l,k,h,q,p.
      class(t_syspla), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
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
       call PhVb2ell_nafk%set_output(len_naf,                               &
     &       comp_naf_alkhqp%getconsumer())
     
      end subroutine init_naf_alkhqp

!------------------------------------------------------------
!> initialise l'analyse en frequence en difference d'angles exp(i*(..-..))
!!  precise quelles composantes sont utilisees
!------------------------------------------------------------
      subroutine init_naf_diffalp(par_naf_diffalp,                               &
     &              comp_naf_diffalp,  io_naf_diffalp,                           &
     &              naf_diffalp_pla, sys, cG, id_ci,  PhVb2ell_nafd)
        use mod_filepar
        use mod_syspla
        use mod_converter
        use mod_naf_diffang
      implicit none
      type(t_naf_par), intent(inout) :: par_naf_diffalp !< donnee du fichier de parametres pour l'analyse en frequence en a,l,k,h,q,p
      type(t_naf_diffang), intent(inout) :: comp_naf_diffalp !< module de calcul de l'analyse en frequence en a,l,k,h,q,p
      type (t_io_txt_ncomp), intent(inout) :: io_naf_diffalp  !< sortie de l'analyse en frequence en a,l,k,h,q,p.
      integer,dimension(1:2), intent(in) :: naf_diffalp_pla !< numero des 2 planetes utilisees pour naf_diffalp
      class(t_syspla), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_plan_PhVb2ell), intent(inout) :: PhVb2ell_nafd !< buffer des elements elliptiques en entree
      character(len=*), intent(in) :: id_ci  !< identifiant de la condition initiale
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an

      integer(8) len_naf
      integer, dimension(1:2) :: comp_p1, comp_p2
      integer p1, p2
      integer nplan, ncompout
      integer ncompinput ! nombre de composante par planete
      integer ncompouttotal ! nombre de composante total en sortie

       nplan = sys%m_plan_nb

       ncompout = par_naf_diffalp%m_naf_nterm
       ncompinput = 2  ! cas standard exp(i*(l1-l2)), exp(i*(pi1-pi2))
       ncompouttotal = ncompinput*ncompout
       
       call io_naf_diffalp%set_comp_nb(1+ncompouttotal*3)
       call io_naf_diffalp%set_prefixline(id_ci)
       len_naf = par_naf_diffalp%m_stepout/par_naf_diffalp%m_stepcalc
       call comp_naf_diffalp%set_stepout(len_naf)
       
       ! premiere planete
       p1 = naf_diffalp_pla(1)
       comp_p1(1) = (p1-1)*6+4 ! pour la
       comp_p1(2) = (p1-1)*6+5 ! pour pi
       ! deuxieme planete
       p2 = naf_diffalp_pla(2)
       comp_p2(1) = (p2-1)*6+4 ! pour la
       comp_p2(2) = (p2-1)*6+5 ! pour pi
       
       call comp_naf_diffalp%set_component12(comp_p1, comp_p2, 6*nplan)
       call comp_naf_diffalp%set_param_naf(par_naf_diffalp%m_naf_nterm,     &
     &   par_naf_diffalp%m_naf_iw, par_naf_diffalp%m_naf_isec,              &
     &   par_naf_diffalp%m_naf_dtour, par_naf_diffalp%m_naf_tol)
       call comp_naf_diffalp%set_output(1_8,                                &
     &   io_naf_diffalp%getconsumer())
       call comp_naf_diffalp%set_graphdotname("NAF exp(i*(l1-l2)),...")
      
       call PhVb2ell_nafd%set_kind(par_naf_diffalp%m_elltype)      
       call PhVb2ell_nafd%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_nafd%set_output(len_naf,                                 &
     &       comp_naf_diffalp%getconsumer())
      
      end subroutine init_naf_diffalp
      
!------------------------------------------------------------
!> initialise l'analyse en frequence en difference d'angles exp(i*(..-..))
!!  precise quelles composantes sont utilisees
!------------------------------------------------------------
      subroutine init_naf_difflpm(par_naf_difflpm,                               &
     &              comp_naf_difflpm,  io_naf_difflpm,                           &
     &              naf_difflpm_pla, sys, cG, id_ci,  PhVb2ell_nafd)
        use mod_filepar
        use mod_syspla
        use mod_converter
        use mod_naf_difflpm
      implicit none
      type(t_naf_par), intent(inout) :: par_naf_difflpm !< donnee du fichier de parametres pour l'analyse en frequence en a,l,k,h,q,p
      type(t_naf_difflpm), intent(inout) :: comp_naf_difflpm !< module de calcul de l'analyse en frequence en a,l,k,h,q,p
      type (t_io_txt_ncomp), intent(inout) :: io_naf_difflpm  !< sortie de l'analyse en frequence en a,l,k,h,q,p.
      integer,dimension(1:2), intent(in) :: naf_difflpm_pla !< numero des 2 planetes utilisees pour naf_difflpm
      class(t_syspla), pointer, intent(inout) :: sys  !< systeme planetaire a integrer
      type (t_plan_PhVb2ell), intent(inout) :: PhVb2ell_nafd !< buffer des elements elliptiques en entree
      character(len=*), intent(in) :: id_ci  !< identifiant de la condition initiale
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an

      integer(8) len_naf
      integer, dimension(1:2) :: comp_p1, comp_p2
      integer p1, p2
      integer nplan, ncompout
      integer ncompinput ! nombre de composante par planete
      integer ncompouttotal ! nombre de composante total en sortie

       nplan = sys%m_plan_nb

       ncompout = par_naf_difflpm%m_naf_nterm
       ncompinput = 4  ! cas standard ((l1-l2)+i*0) mod [0,2pi] et [-PI,PI] , ((pi1-pi2)+i*0)  mod [0,2pi] et [-PI,PI]
       ncompouttotal = ncompinput*ncompout
       
       call io_naf_difflpm%set_comp_nb(1+ncompouttotal*3)
       call io_naf_difflpm%set_prefixline(id_ci)
       len_naf = par_naf_difflpm%m_stepout/par_naf_difflpm%m_stepcalc
       call comp_naf_difflpm%set_stepout(len_naf)
       
       ! premiere planete
       p1 = naf_difflpm_pla(1)
       comp_p1(1) = (p1-1)*6+4 ! pour la
       comp_p1(2) = (p1-1)*6+5 ! pour pi
       ! deuxieme planete
       p2 = naf_difflpm_pla(2)
       comp_p2(1) = (p2-1)*6+4 ! pour la
       comp_p2(2) = (p2-1)*6+5 ! pour pi
       
       call comp_naf_difflpm%set_component12(comp_p1, comp_p2, 6*nplan)
       call comp_naf_difflpm%set_param_naf(par_naf_difflpm%m_naf_nterm,     &
     &   par_naf_difflpm%m_naf_iw, par_naf_difflpm%m_naf_isec,              &
     &   par_naf_difflpm%m_naf_dtour, par_naf_difflpm%m_naf_tol)
       call comp_naf_difflpm%set_output(1_8,                                &
     &   io_naf_difflpm%getconsumer())
       call comp_naf_difflpm%set_graphdotname(                              &
     &   "NAF ((l1-l2) +i*0) mod 2PI ")
      
       call PhVb2ell_nafd%set_kind(par_naf_difflpm%m_elltype)      
       call PhVb2ell_nafd%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell_nafd%set_output(len_naf,                                 &
     &       comp_naf_difflpm%getconsumer())
      
      end subroutine init_naf_difflpm

      end module

      module mod_intg_symp_start
      contains
!***********************************************************************
!> fonction appellee par le distributeur de condition initiale pour la condition ligneci
!***********************************************************************
      subroutine cidist_integfunc(ligneci, filepar)
        use mod_io_int
        use mod_sysplanewth
        use mod_sysplanewth_adapt
        use mod_sysplanewtj
        use mod_sysplanewtj_adapt
        use mod_converter
        use mod_ci_planewt
        use mod_ci_pladistmin
        use mod_ci_multibase
        use mod_filepar
        use mod_gmsun
        use mod_intg_symp
        implicit none
        
        class(t_ci_base), intent(inout) :: ligneci !< condition initiale recue
        class(*), intent(inout) :: filepar !< donnee utilisateur : fichier de parametres
        
        class(t_ci_base), pointer :: pci1 !< condition initiale recue 1
        real(TREAL), dimension(:), allocatable :: mpl ! masse normalise (masse de l'etoile a 1)
        real(TREAL), dimension(:,:), allocatable :: CI
        real(TREAL), dimension(:), allocatable, target :: ctrlpladmin ! distance minimale entre planetes
        real(TREAL), dimension(:), pointer :: pctrlpladmin => NULL()! distance minimale entre planetes
        real(TREAL),allocatable,dimension(:,:) :: Xplan,XPplan
        real(TREAL),allocatable,dimension(:,:) :: Ph,Vh, Vb, Pj,Vj
        real(TREAL),allocatable,dimension(:,:) :: Phinv, Vbinv, Vhinv
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
        
        select type(filepar)
         class is(t_extfilepar)
        
        select type(ligneci)
         class is(t_ci_multibase)
         ! condition initiale du systeme
         pci1 => ligneci%m_arci(1)%m_ptr
         select type(pci1)
         class is(t_ci_planewt)
          nplan = pci1%m_plan_nb
          ci_type = pci1%m_ci_type
          allocate(mpl(0:nplan))
          allocate(CI(1:6,1:nplan))
          mpl(0:nplan) = pci1%m_plan_mass(0:nplan)
! Masses normalisees a 1
          cG=cG*mpl(0)
          do k=1,nplan
            mpl(k)=mpl(k)/mpl(0)
          enddo
          mpl(0)=1.D0
          ! ici le 1:nplan au lieu de 0:nplan car ci_type: heliocentrique
          CI(1:6,1:nplan) = pci1%m_plan_coord_ci(1:6,1:nplan)
         
         class default
          stop 'pci1 in cidist_integfunc : bad class 1'
         end select 

          ! recuperation du controle des distance minimales entre planetes
          if (filepar%ctrl_distpla%m_compute.ne.0) then
             nullify(pci1)
             pci1 => ligneci%m_arci(2)%m_ptr
             select type(pci1)
             class is(t_ci_pladistmin)
                allocate(ctrlpladmin(1:nplan))
                ctrlpladmin = pci1%m_plan_dmin
                pctrlpladmin => ctrlpladmin
             class default
               stop 't_ci_pladistmin in cidist_integfunc : bad class'
             end select 
          endif 
         

         class default
          stop 'ligneci in cidist_integfunc : bad class all'
        end select 


      int_type = trim(filepar%int_type)
      

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

! --- create the filename of output files  
      call filepar%create_filename_ci(ligneci%m_id)

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


      call intg_symp(filepar,psys,cG,pctrlpladmin,ligneci%m_id,0_8,0_8)

      write(*,*) 'fin d integration'

      write(*,*) "end cidist_integfunc"
        
      if (allocated(sysJ)) then
         deallocate(sysJ)
      end if
      if (allocated(sysH)) then
         deallocate(sysH)
      end if
      if (allocated(sysJ_adapt)) then
         deallocate(sysJ_adapt)
      end if
      if (allocated(sysH_adapt)) then
         deallocate(sysH_adapt)
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
       if (allocated(ctrlpladmin)) then
          deallocate(ctrlpladmin)
       end if
       
       
         class default
          stop 'filepar in cidist_integfunc : bad class'
        end select 

       end subroutine cidist_integfunc
      end module

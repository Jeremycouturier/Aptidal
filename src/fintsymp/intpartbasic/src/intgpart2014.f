!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file filepar.f 
!!  \brief  routine d'integration des particules : lance l'integration 
!!   et ecrit dans les fichiers de sortie
!!
! history : creation 16/06/2014
!***********************************************************************
#include "realc.f"
!-----------------------------------------------------------
! En entree: tinit: temps initial
!            dt : intervalle de temps pour une iteration
!            n_iter : nombre d'iterations a effectuer
!            n_out : nombre d'iterations apres lequel on appelle ecrit
!-----------------------------------------------------------
      module mod_intg_symp
        
        use mod_filepar
        use mod_ci_planewt
!-----------------------------------------------------------------------
!!> type aggregateur pour le passage d'argument a  cidist_integfunc
!-----------------------------------------------------------------------
        type t_cidist_integfunc_arg
            type(t_extfilepar) :: m_filepar !< fichier de parametres
            class(t_ci_base), pointer :: m_cipla   !< condition du systeme planetaire
        end type t_cidist_integfunc_arg
        
      contains


!***********************************************************************
!> fonction appellee par le distributeur de condition initiale pour les conditions tabci
!***********************************************************************
      subroutine cidist_integfunc(nbci, tabci, multipledata)
        use mod_io_int
        use mod_io_txt_ncomp
        use mod_io_txt_car
        use mod_sysplaexth
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
        type(t_io_int) ::io_int
        type(t_io_txt_car) ::iopla_cart
        type(t_io_txt_car_part) ::iopart_cart
        type(t_io_txt_ncomp) ::iopla_ell
        type(t_io_txt_ell_part) :: iopart_ell 
        type(t_sysplaextH)  :: sysH ! systeme en helio
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

         class default
          stop 'multipledata in cidist_integfunc : bad class'
        end select 

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
      
 ! --- create the filename of output files  of the particules
      call filepar%create_filename_part(tabci(1)%m_ptr%m_id)
      call filepar%create_filename_pla(cipla%m_id)

      if (filepar%if_int.eq.1) then
         write(*,*)  ' nf_int ', trim(filepar%nf_int)
         call io_int%set_name(trim(filepar%nf_int),60)
         call io_int%fcreate()
      end if
      if (filepar%if_ell.eq.1) then
         write(*,*)  ' nfpla_ell ', trim(filepar%nfpla_ell)
         iopla_ell = t_io_txt_ncomp("I/O PLA ell")
         call iopla_ell%set_comp_nb(1+6*nplan)
         call iopla_ell%set_name(trim(filepar%nfpla_ell),61)
         call iopla_ell%fcreate()
         write(*,*)  ' nfpart_ell ', trim(filepar%nfpart_ell)
         iopart_ell = t_io_txt_ell_part("I/O PART ell")
         call iopart_ell%set_col_nb(1+6, npart, 1+6*nplan+6*npart)
         call iopart_ell%set_name(trim(filepar%nfpart_ell),63)
         call iopart_ell%set_part_nb(nplan, npart)
         do k=1, npart
           call iopart_ell%set_pref2dline(k,tabci(k)%m_ptr%m_id) 
         enddo
         call iopart_ell%fcreate()
      end if
      if (filepar%if_car.eq.1) then
         write(*,*) ' nfpla_car ', trim(filepar%nfpla_car)
         iopla_cart = t_io_txt_car("I/O PLA cart")
         call iopla_cart%set_comp_nb(1+6*nplan)
         call iopla_cart%set_name(trim(filepar%nfpla_car),62)
         call iopla_cart%set_plan_nb(nplan)
         call iopla_cart%fcreate()
         write(*,*) ' nfpart_car ', trim(filepar%nfpart_car)
         iopart_cart = t_io_txt_car_part("I/O PART cart")
         call iopart_cart%set_col_nb(1+6, npart, 1+6*nplan+6*npart)
         call iopart_cart%set_name(trim(filepar%nfpart_car),64)
         call iopart_cart%set_part_nb(nplan, npart)
         do k=1, npart
           call iopart_cart%set_pref2dline(k,tabci(k)%m_ptr%m_id) 
         enddo
         call iopart_cart%fcreate()
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



       call intg_symp(filepar,sysH, cG, io_int,iopla_cart,iopla_ell,       &
     &           iopart_cart, iopart_ell)

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
      
       ! fermeture des fichiers pour les planetes
      if (filepar%if_int.eq.1) then
         call io_int%fclose()
      endif
      if (filepar%if_ell.eq.1) then
         call iopla_ell%fclose()
      end if
      if (filepar%if_car.eq.1) then
         call iopla_cart%fclose()
      end if

      end subroutine cidist_integfunc
      
      
!------------------------------------------------------------
!> routine d'integration
!------------------------------------------------------------
      subroutine intg_symp(filepar,sys,cG,io_int,iopla_cart,iopla_ell,    &
     &           iopart_cart, iopart_ell)
      use mod_sysplaexth
      use mod_intsympbase
      use mod_schema
      use mod_io_int
      use mod_converter
      use mod_io_txt_ncomp
      use mod_io_txt_car
      use mod_syspartnewth
      use mod_filepar
      use mod_converterpart
      use mod_io_txt_carpart
        use mod_io_txt_ellpart
      
      implicit none
      type (t_extfilepar),target, intent(inout) :: filepar !< fichier de parametres + variables derivees
      type (t_io_int),intent(inout) ::io_int  !< sortie des integrales premieres
      type (t_io_txt_car), intent(inout) :: iopla_cart  !< sortie des elements cartesiens des planetes
      type (t_io_txt_ncomp), intent(inout) :: iopla_ell !< sortie des elements elliptiques des planetes
      type (t_io_txt_car_part), intent(inout) :: iopart_cart  !< sortie des elements cartesiens des particules
      type (t_io_txt_ell_part), intent(inout) :: iopart_ell !< sortie des elements elliptiques des particules
      class(t_sysplaexth), intent(inout) :: sys  !< systeme planetaire a integrer
      real(TREAL), intent(in) :: cG !< constante de gauss en UA/an
     
      integer(8) :: n_iter
      integer(8) :: n_out
      real(TREAL) :: tinit,dt
      
      ! pour les planetes
      type (t_plan_PhVb2PhVh) :: PhVb2PhVh 
      type (t_plan_PhVb2ell) :: PhVb2ell   
      
      ! pour les particules
      class(t_sysextension), pointer :: psysext 
      class(t_syspartnewtH), pointer :: psyspartH =>null()
      type (t_part_PhVb2PhVh) :: PhVb2PhVhpart
      type (t_part_PhVb2ell) :: PhVb2ellpart
      integer npart
      
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
      
      call sys%clear_out_step()
      
      call sys%get_extension(1, psysext) ! recupere les particules
      select type(psysext)
        class is(t_syspartnewtH)
           psyspartH=>psysext
           npart = psyspartH%get_nb()
        class default
               stop 'intg_symp : bad class t_syspartnewtH for psysext'
      end select 

      
      ! sortie integrales premieres
      if (filepar%if_int.eq.1) then
        call sys%add_int_out_step(n_out,1_8,io_int%getconsumer())  
      endif
      
      !--------------------------------
      ! cas integration en helio
      !--------------------------------

      ! sortie cartesienne
      if (filepar%if_car.eq.1) then
       call PhVb2PhVh%set_mass(sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2PhVh%set_output(1_8,iopla_cart%getconsumer())
       call sys%add_plan_out_step(n_out,1_8, PhVb2PhVh%getconsumer())
       
       call PhVb2PhVhpart%set_mass(sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2PhVhpart%set_nb(npart)      
       call PhVb2PhVhpart%set_output(1_8,iopart_cart%getconsumer())
       call psyspartH%add_part_out_step(n_out,1_8,                         &
     &           PhVb2PhVhpart%getconsumer())
      endif
     
      ! sortie elements elliptiques
      if (filepar%if_ell.eq.1) then
       call PhVb2ell%set_kind(filepar%out_ell)      
       call PhVb2ell%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ell%set_output(1_8, iopla_ell%getconsumer())

       call sys%add_plan_out_step(n_out,1_8,PhVb2ell%getconsumer())

       call PhVb2ellpart%set_mass(cG,sys%m_star_mass,sys%m_plan_mass)      
       call PhVb2ellpart%set_nb(npart)      
       call PhVb2ellpart%set_output(1_8, iopart_ell%getconsumer())
       call PhVb2ellpart%set_kind(filepar%out_ell)      
       call psyspartH%add_part_out_step(n_out,1_8,                         &
     &           PhVb2ellpart%getconsumer())
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

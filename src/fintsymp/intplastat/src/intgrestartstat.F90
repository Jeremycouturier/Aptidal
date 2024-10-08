!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file intgrestartstat.F90
!!  \brief  redemarre une integration a partir d'une date fixee
!! et lance l'integration 
!!
! history : creation 23/11/2018
!***********************************************************************
#include "realc.f"
module mod_intg_symp_restart

contains

!***********************************************************************
!> fonction appellee pour le redemarrage
!***********************************************************************
subroutine cidist_integfunc_restart(filepar, readerpladmin)
  use mod_intg_symp
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
  use mod_ciread_pladistmin
  use mod_intsympbase
  use mod_io_dump
  implicit none
  
  class(*), target, intent(inout) :: filepar !< donnee utilisateur : fichier de parametres
  type(t_ciread_pladistmin), intent(inout) :: readerpladmin !< lecteur des distances minimales
  
  class(t_ci_base), pointer :: cipladistmin !< distances miniales lues
  real(TREAL), dimension(:), allocatable, target :: ctrlpladmin ! distance minimale entre planetes
  real(TREAL), dimension(:), pointer :: pctrlpladmin => NULL()! distance minimale entre planetes
  type(t_sysplanewtH),allocatable, target  :: sysH ! systeme en helio
  type(t_sysplanewtJ),allocatable, target  :: sysJ ! systeme en jacobi
  type(t_sysplanewtJ_adapt),allocatable, target  :: sysJ_adapt ! systeme en jacobi et pas adaptatif
  type(t_sysplanewtH_adapt),allocatable, target  :: sysH_adapt ! systeme en heliocentrique et pas adaptatif
  class(t_syspla), pointer  :: psys ! systeme en jacobi ou helio
  character(100)  :: int_type ! type de l'integrateur
  character(len=64)  :: idrestart ! numero de la condition a redemarrer
  integer :: nplan ! nombre de planetes
  real(TREAL) :: cG
  type(t_intsympbase) :: integrator
  logical cifoundpladmin
  type(t_io_dump) :: io_restart
  
  ! constante GM du soleil
select type(filepar)
class is(t_extfilepar)

cG = calcGM_sun(filepar%ref_gmsun)
    
int_type = trim(filepar%int_type)

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


! restaure le systeme integre
call integrator%set_syssymp(psys)
call io_restart%set_name(filepar%nf_initext, 58)
call io_restart%fopen()
call integrator%restorebin(io_restart, filepar%tinit)
call io_restart%fclose()
idrestart = integrator%m_id
nplan = psys%get_plan_nb()
write(*,*) 'restart with ', idrestart
write(*,*) 'number of planets ', nplan


! --- create the filename of output files  
call filepar%create_filename_ci(idrestart)

! controle des distances minimales
if (filepar%ctrl_distpla%m_compute.ne.0) then
    cifoundpladmin = .false.
    call readerpladmin%fopen();
    do while ((cifoundpladmin.eqv..false.) .and. (readerpladmin%freadline(cipladistmin)))
      cifoundpladmin = cipladistmin%m_id.eq.idrestart
      if (cifoundpladmin.eqv..true.) goto 500
    end do
    call cipladistmin%debug()

500 if (cifoundpladmin.eqv..false.) then
      write(*,*) 'distance minimale non trouve pour',idrestart 
      stop 1
    endif

    select type(cipladistmin)
    class is(t_ci_pladistmin)
      allocate(ctrlpladmin(1:nplan)) 
      ctrlpladmin = cipladistmin%m_plan_dmin
      pctrlpladmin => ctrlpladmin
    class default
       stop 't_ci_pladistmin in cidist_integfunc : bad class'
    end select 
endif 

filepar%tinit = integrator%m_t0
call intg_symp(filepar,psys, cG, pctrlpladmin, idrestart, integrator%m_it, integrator%m_si)

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

if (allocated(ctrlpladmin)) then
    deallocate(ctrlpladmin)
end if
 
 class default
 stop 'filepar in cidist_integfunc_restart : bad class'
end select 

end subroutine cidist_integfunc_restart

end module mod_intg_symp_restart
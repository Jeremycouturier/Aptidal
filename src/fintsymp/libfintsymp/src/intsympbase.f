!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file intsympbase.f 
!!  \brief integrateur d'un systeme derivant de t_syssymp  
!!         pouvant etre integre par un schema symplectique
!!
!!     cette classe t_intsympbase fournit les routines elementaires communes a l'ensemble 
!!      des systemes possibles
!!
!!
!
! history : creation 20/07/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
!> module pour l'integration symplectique de 1 systeme avec N planetes 
!! et M particules 
!***********************************************************************
      module mod_intsympbase
      use mod_syssympbase
      use mod_schema
      use mod_io_dump

!***********************************************************************
!> @class t_intsympbase
!! classe de base decrivant l'integrateur symplectique d'un systeme derivant de t_syssymp
!!
!!  
!***********************************************************************
      type :: t_intsympbase 
          ! les variables suivantes toutes privees ou reservees au type derivant de celui-ci 
          ! (ne doivent pas etre utilisees par les applications)
          
          character(len=64) :: m_id  !< identifiant 
          
          class (t_syssymp), pointer :: m_syspla !< systeme planetaire integre
          type (t_io_dump):: m_dump !< dump de l'ensemble de l'integration
          
          ! description de l'integration
          type(t_schema)    :: m_scheme !< schema d'integration
          real(TREAL)       :: m_t0     !< temps initial
          integer(8)        :: m_it     !< iteration courante
          real(TREAL)       :: m_dt     !< pas d'integration
          integer(8)        :: m_si     !< iteration courante
          
          ! gestion des erreurs
          !>comportement en en cas d'erreur
          !! * = 0 => arret immediat en cas d'erreur : sortie du programme
          !! * = 1 => arret immediat en cas d'erreur et retourne immediatement a la fonction appellante
          integer           :: m_errorbehavior  = 0
          
          !> gestion du type de pas
          !! * = 0 => integration a pas fixe
          !! * = 1 => integration a pas adaptatif avec sortie au temps fixe
          !! * = 2 => integration a pas adaptatif avec sortie au temps non fixe
          integer           :: m_type_pas  = 0
                  
       contains
          private
                    
          procedure, public:: set_syssymp ! fixe le systeme a integrer
          procedure, public:: set_id ! fixe le nom du systeme
          procedure, public:: set_integ_scheme ! fixe le schema de l'integrateur
          procedure, public:: set_integ_time ! fixe le temps initial, le pas et l'iteration courante d'integrateur
          procedure, public:: set_integ_time_si ! fixe l'iteration du pas fixe d'integrateur
          procedure, public:: get_integ_current_time ! recupere le temps actuel de l'integrateur
          procedure, public:: set_type_pas ! fixe le  type de pas
          procedure, public:: set_dump  => intsympbase_setdump ! fixe le dump
          procedure, public:: integ_run  ! integre pendant niter iterations
          procedure, NON_OVERRIDABLE, public:: set_error_behavior ! fixe le comportement en cas d'erreur
                   
          procedure, private:: dumpbin  => intsympbase_dumpbin ! dump binaire vers un fichier 
          procedure, public:: restorebin  => intsympbase_restorebin ! restaure binaire depuis un fichier

      end type t_intsympbase  

     
      contains
         
!***********************************************************************
!> @brief fixe le systeme planetaire a integrer
!***********************************************************************
      subroutine set_syssymp(this, syspla)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       class(t_syssymp), target, intent(in) :: syspla  !< systeme planetaire
       
       this%m_syspla => syspla
      end  subroutine set_syssymp  

!***********************************************************************
!> @brief fixe le nom du systeme planetaire 
!***********************************************************************
      subroutine set_id(this, id)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       character(len=64), intent(in) :: id  !< nom du systeme
       
       this%m_id = id
      end  subroutine set_id  


!***********************************************************************
!> @brief fixe le schema d'integration (e.g. 'ABA82' ou 'ABAH82') 
!***********************************************************************
      subroutine set_integ_scheme(this, schema)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       character(len=*), intent(in) :: schema      !< schema d'integration (e.g., 'ABA82' ou 'BAB82' )
       
       call create_schema(schema,this%m_scheme)
       
      end  subroutine set_integ_scheme  

!***********************************************************************
!> @brief fixe l'iteration courante pour le fixe fixe des integrateurs a pas adapatifs
!***********************************************************************
      subroutine set_integ_time_si(this, si)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       integer(8), intent(in)  :: si      !< iteration courante (en general, 0) 
       this%m_si = si
      end  subroutine set_integ_time_si  

!***********************************************************************
!> @brief fixe le temps initial (en genreal, 0.0) , le pas d'integration 
!! et l'iteration courante (en general 0)
!***********************************************************************
      subroutine set_integ_time(this, t0, dt, it)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: t0      !< temps initial (en general, 0.D0)
       real(TREAL), intent(in) :: dt      !< pas d'integration 
       integer(8), intent(in)  :: it      !< iteration courante (en general, 0) 
       
       this%m_t0 = t0
       this%m_dt = dt
       this%m_it = it
      end  subroutine set_integ_time  



!***********************************************************************
!> @brief recupere le temps actuel de l'integration
!***********************************************************************
      function get_integ_current_time(this) result(t)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       real(TREAL):: t     !< temps actuel de l'integration
       
       if (this%m_syspla%sys_pas_fixe) then
         t = this%m_t0 + this%m_it*this%m_dt
       else
         t = this%m_syspla%get_adaptativetime()
       endif
         
      end  function get_integ_current_time  


!***********************************************************************
!> @brief fixe le pas du dump
!***********************************************************************
      subroutine intsympbase_setdump(this, iodump)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       type(t_io_dump),  intent(in) ::  iodump  !< fichier de dump
      
       this%m_dump = iodump
         
      end  subroutine intsympbase_setdump  


!***********************************************************************
!> @brief controle des erreurs pendant l'integration
!! @return retourne 0  si aucune action a faire \n
!! retourne 1 si on doit s'arreter d'integrer (mais poursuite de l'execution)
!***********************************************************************
      function integ_check_error(this, t) result(ret)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       real(TREAL), intent (in) :: t !< temps actuel
       integer ret ! code de retour
       ret = 0
       
         if (this%m_syspla%check_error()) then
           select case (this%m_errorbehavior)
           case  (1)  
             write (*,*) " arret de l'integration a t=", t
             ret = this%m_syspla%set_error_time(t)
             
           case  default  
             write (*,*) " arret de l'integration a  t=", t
             write (*,*) " arret du programme !!"
             stop
           end select  
         endif
      end function  integ_check_error  
         
!***********************************************************************
!> @brief controle des erreurs dans les sucesseurs 
!! @return retourne 0  si aucune action a faire \n
!! retourne 1 si on doit s'arreter d'integrer (mais poursuite de l'execution)
!***********************************************************************
      function integ_check_sucessor_error(this,t,sysplaout) result(ret)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       real(TREAL), intent (in) :: t !< temps actuel
       class(t_syssymp), intent(inout) :: sysplaout
       integer ret ! code de retour
       ret = 0
       
         if (sysplaout%check_error()) then
           select case (this%m_errorbehavior)
           case  (1)  
             write (*,*) " arret de l'integration a t=", t
             call this%m_syspla%set_error(sysplaout%get_error())
             ret = sysplaout%set_error_time(t)
             
           case  default  
             write (*,*) " arret de l'integration a  t=", t
             write (*,*) " arret du programme !!"
             stop
           end select  
         endif
      end function  integ_check_sucessor_error  

!***********************************************************************
!> @brief execute niter iterations du schema d'integration
!***********************************************************************
      subroutine integ_run(this, niter)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: niter             !< nombre d'iterations

       ! verification de la coherence
       if ((this%m_type_pas.eq.0) .and.                                       &
     &      (this%m_syspla%sys_pas_fixe.eqv..false.)) then
          write(*,*) 'integrateur pas fixe avec systeme pas variable'
          write(*,*) 'bug interne'
          stop 1
       endif
       if ((this%m_type_pas.ne.0) .and.                                       &
     &     (this%m_syspla%sys_pas_fixe.eqv..true.)) then
          write(*,*) 'integrateur pas variable avec systeme pas fixe'
          write(*,*) 'bug interne'
          stop 1
       endif
       
       ! choix de l'integrateur       
       if (this%m_type_pas.eq.0 .or. this%m_type_pas.eq.2) then
         call integ_run_pas_fixe(this, niter)
       else
         call integ_run_pas_adaptatif(this, niter)
       end if

      end  subroutine integ_run 

!***********************************************************************
!> @brief execute niter iterations du schema d'integration pour un schema a pas fixe
!***********************************************************************
      subroutine integ_run_pas_fixe(this, niter)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: niter             !< nombre d'iterations
       integer(8) i, i0
       integer(8) minmodulo
       integer nstepscall,j, nstepscallmid
       real(TREAL) :: dt, t
       class(t_syssymp), allocatable :: sysplaout
       integer(8),dimension(:), allocatable :: arstepscall
       integer(8),dimension(:), allocatable :: arstepscallmid
       logical isoutputstep
       integer code_err
       
       write(*,*) 'niter=',niter, this%m_t0 , this%m_it, this%m_dt       
       minmodulo = 0
       i = 1
       dt = this%m_dt
       
       
       call this%m_syspla%get_outputsteps(nstepscallmid,arstepscallmid)
       ! creation d'une copie pour la sortie
       call this%m_syspla%create_output(sysplaout)
       call sysplaout%get_outputsteps(nstepscall,arstepscall)
       
       ! sortie dans le cas du demarrage standard
       call this%m_syspla%copy_output(sysplaout)   
       if (this%m_it.ne.0_8) then
          call this%m_scheme%m_procdepart(sysplaout,-dt)
       endif
       i0 = this%m_it
       call sysplaout%write_output(i0, this%m_t0+i0*dt)            
    
       if (this%m_it.eq.0_8) then
        call this%m_scheme%m_procdepart(this%m_syspla,dt)
       endif

       do while (i.le.niter)
         this%m_it = this%m_it+1
         t = this%m_t0 + this%m_it*dt
         call this%m_scheme%m_procrun(this%m_syspla,dt)

         ! dump binaire si besoin
         call this%dumpbin(this%m_it, this%m_dump) 
         
         ! controle de l'erreur
         code_err = integ_check_error(this, t)
         if (code_err.eq.1) then
          exit ! sortie de boucle 
         endif
         
         ! verifie si on est un pas de sortie symplectique
           isoutputstep = .false.
           i = this%m_it
           do j=1, nstepscallmid
            if (mod(i,arstepscallmid(j)).eq.0) then
             isoutputstep = .true.
            endif
           enddo
           if (isoutputstep) then
              call this%m_syspla%write_output(i, t)     
              ! controle de l'erreur dans les sorties 
              code_err = integ_check_sucessor_error(this, t,            &
     &                                              this%m_syspla)
              if (code_err.eq.1) then
               exit ! sortie de boucle 
              endif
           endif

         ! verifie si on est un pas de sortie normale
         isoutputstep = .false.
         do j=1, nstepscall
          if (mod(i,arstepscall(j)).eq.0) then
           isoutputstep = .true.
          endif
         enddo 
         ! sortie des donnees
         if (isoutputstep.eqv..true.) then
            !-- utilisation du systeme de sortie
            call this%m_syspla%copy_output(sysplaout)           
            call this%m_scheme%m_procdepart(sysplaout,-dt)
            call sysplaout%write_output(this%m_it, t)     
            call sysplaout%copy_output_feedback(this%m_syspla)     
            
            ! controle de l'erreur dans les sorties 
            code_err = integ_check_sucessor_error(this, t, sysplaout)
            if (code_err.eq.1) then
             exit ! sortie de boucle 
            endif
         end if
         i = i + 1
       end do
       
       !-- flush des sorties 
       call sysplaout%flush_output()
       
       !-- suppression du systeme de sortie
       call sysplaout%delete_output(sysplaout)
       
       if (allocated(arstepscall)) then
          deallocate(arstepscall)
       endif
       if (allocated(arstepscallmid)) then
          deallocate(arstepscallmid)
       endif
      end  subroutine integ_run_pas_fixe  


!***********************************************************************
!> @brief execute niter iterations du schema d'integration pour un schema a pas adaptatif
! gere le cas du pas adaptatif pour les sorties a temps fixes
!***********************************************************************
      subroutine integ_run_pas_adaptatif(this, niter)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: niter             !< nombre d'iterations
       integer(8) i, i0
       integer(8) prev_it 
       real(TREAL) :: st, ds, out_ds
       integer(8) minmodulo
       integer nstepscall,j, iterout, nstepscallmid
       real(TREAL) :: dt, t, tfinal, prev_t, exact_t, out_t
       real(TREAL) :: tlen
       class(t_syssymp), allocatable :: sysplaout
       integer(8),dimension(:), allocatable :: arstepscall
       integer(8),dimension(:), allocatable :: arstepscallmid
       logical isoutputstep
       integer code_err
       real(TREAL) eps

       
       eps = epsilon(REALC(1.e0))
       
       minmodulo = 0
       dt = this%m_dt
       ds = this%m_dt
       
       call this%m_syspla%get_outputsteps(nstepscallmid,arstepscallmid)
       
       ! creation d'une copie pour la sortie
       call this%m_syspla%create_output(sysplaout)
       call sysplaout%get_outputsteps(nstepscall,arstepscall)
       
       ! sortie dans le cas du demarrage standard
       call this%m_syspla%copy_output(sysplaout)           
       if (this%m_it.eq.0_8) then
          i0 = this%m_it
          call sysplaout%write_output(i0, this%get_integ_current_time())            
       else
       ! cas d'un redemerrage
              exact_t =  this%m_t0 + this%m_si*dt
              st = exact_t
              !-- utilisation du systeme de sortie
              call this%m_syspla%copy_output(sysplaout)
              tlen = 1 !sysplaout%get_adaptativetime() - prev_t 
              call this%m_scheme%m_procdepart(sysplaout,-ds)
              out_t = sysplaout%get_adaptativetime() 
              out_ds = ds    
              iterout = 0
              do while ((abs((out_t-exact_t)/exact_t).gt.(5.*eps)).and. &
     &               (iterout.le.5))
               out_ds = out_ds*(exact_t-out_t)/tlen
               call this%m_scheme%m_procdepart(sysplaout,out_ds)
               call this%m_scheme%m_procrun(sysplaout,out_ds)
               call this%m_scheme%m_procdepart(sysplaout,-out_ds)
               st = st + out_ds
               prev_t = out_t
               out_t = sysplaout%get_adaptativetime() 
               tlen = out_t-prev_t  
               iterout = iterout+1    
              enddo
              call sysplaout%write_output(this%m_it, st)     
       endif
    
       if (this%m_it.eq.0_8) then
        call this%m_scheme%m_procdepart(this%m_syspla,dt)
       endif
    
       prev_it = this%m_it
       t = this%m_t0
       tfinal = this%m_t0 + (this%m_si+niter)*dt
       
       do while (t.lt.tfinal)
         this%m_si = this%m_si + 1
         st = this%m_t0 + this%m_si*ds
         call this%m_scheme%m_procrun(this%m_syspla,ds)
         
         ! recuperation du temps du systeme
         prev_t = t
         t = this%get_integ_current_time()
         this%m_it = t/dt
         
         ! controle de l'erreur
         code_err = integ_check_error(this, t)
         if (code_err.eq.1) then
          exit ! sortie de boucle 
         endif


          ! verifie si on est un pas de sortie symplectique
           isoutputstep = .false.
           i = this%m_it
           do j=1, nstepscallmid
            if (mod(i,arstepscallmid(j)).eq.0) then
             isoutputstep = .true.
            endif
           enddo
           if (isoutputstep) then
              call this%m_syspla%write_output(i, t)     
              ! controle de l'erreur dans les sorties 
              out_t = this%m_syspla%get_adaptativetime()        
              code_err = integ_check_sucessor_error(this, out_t,           &
     &                                              this%m_syspla)
              if (code_err.eq.1) then
               exit ! sortie de boucle 
              endif
           endif


         ! verifie si on est un pas de sortie
         if (prev_it.ne.this%m_it) then
           ! dump binaire si besoin
           call this%dumpbin(this%m_it, this%m_dump) 

           prev_it = this%m_it
           isoutputstep = .false.
           i = this%m_it
           do j=1, nstepscall
            if (mod(i,arstepscall(j)).eq.0) then
             isoutputstep = .true.
            endif
           enddo
           ! sortie des donnees
           if (isoutputstep.eqv..true.) then
              exact_t =  this%m_t0 + this%m_it*dt
              !-- utilisation du systeme de sortie
              call this%m_syspla%copy_output(sysplaout)
              tlen = sysplaout%get_adaptativetime() - prev_t 
              call this%m_scheme%m_procdepart(sysplaout,-ds)
              out_t = sysplaout%get_adaptativetime() 
              out_ds = ds    
              iterout = 0
              do while ((abs((out_t-exact_t)/exact_t).gt.(5.*eps)).and. &
     &               (iterout.le.5))
               out_ds = out_ds*(exact_t-out_t)/tlen
               call this%m_scheme%m_procdepart(sysplaout,out_ds)
               call this%m_scheme%m_procrun(sysplaout,out_ds)
               call this%m_scheme%m_procdepart(sysplaout,-out_ds)
               st = st + out_ds
               prev_t = out_t
               out_t = sysplaout%get_adaptativetime() 
               tlen = out_t-prev_t  
               iterout = iterout+1    
              enddo
              call sysplaout%write_output(this%m_it, st)     
              call sysplaout%copy_output_feedback(this%m_syspla)     
           
              ! controle de l'erreur dans les sorties 
              out_t = sysplaout%get_adaptativetime()        
              code_err = integ_check_sucessor_error(this, out_t,           &
     &                                              sysplaout)
              if (code_err.eq.1) then
               exit ! sortie de boucle 
              endif
           end if
         endif
       end do
       
       !-- flush des sorties 
       call sysplaout%flush_output()
       
       !-- suppression du systeme de sortie
       call sysplaout%delete_output(sysplaout)
       
       if (allocated(arstepscall)) then
          deallocate(arstepscall)
       endif
       if (allocated(arstepscallmid)) then
          deallocate(arstepscallmid)
       endif
      end  subroutine integ_run_pas_adaptatif 

!***********************************************************************
!> @brief fixe le type de pas pour l'integration
!!
!! valeur des codes 
!! * = 0 => integration a pas fixe
!! * = 1 => integration a pas adaptatif avec sortie au temps fixe
!! * = 2 => integration a pas adaptatif avec sortie au temps non fixe
!***********************************************************************
      subroutine set_type_pas(this, code)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       integer, intent(in)  :: code      !< comportement
       
       this%m_type_pas = code
        
      end  subroutine set_type_pas  

!***********************************************************************
!> @brief fixe le comportement en cas d'erreur
!!
!! valeur des codes 
!! * = 0 => arret immediat en cas d'erreur : sortie du programme
!! * = 1 => arret immediat en cas d'erreur et retourne immediatement a la fonction appellante
!***********************************************************************
      subroutine set_error_behavior(this, code)
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       integer, intent(in)  :: code      !< comportement
       
       this%m_errorbehavior = code
        
      end  subroutine set_error_behavior  

!***********************************************************************
!> @brief dump binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (in) objet a sauvegarder
!***********************************************************************
      subroutine intsympbase_dumpbin(this, curstep, file) 
       use mod_syssympbase
       use mod_io_dump
       implicit none
       class(t_intsympbase), intent(in):: this  !< dummy argument
       integer(8), intent(in) :: curstep !< pas actuel
       type(t_io_dump), intent(in) :: file

        if (this%m_dump%isoutputstep(curstep)) then
          write(file%m_nf) this%m_t0+this%m_it*this%m_dt
          write(file%m_nf) this%m_id
          write(file%m_nf) this%m_t0
          write(file%m_nf) this%m_it
          write(file%m_nf) this%m_dt
          write(file%m_nf) this%m_si
          write(file%m_nf) this%m_errorbehavior
          write(file%m_nf) this%m_type_pas
        
          call this%m_syspla%dumpbin(file)
        endif 

      end subroutine  intsympbase_dumpbin  

!***********************************************************************
!> @brief restauration binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (inout) objet a sauvegarder
!! tload (in) temps ou il faut recherager les donnees
!***********************************************************************
      subroutine intsympbase_restorebin(this, file, tload) 
       use mod_syssympbase
       implicit none
       class(t_intsympbase), intent(inout):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file
       real(TREAL), intent(in) :: tload
       real(TREAL) :: tcur 

       tcur = -1D100
       this%m_dt = 0
       do while (abs(tcur-tload).gt.abs(2.*this%m_dt)) 
        read(file%m_nf,end=1000) tcur
        if (tload<0. .and. tcur.lt.tload) then
           tcur = tload
           goto 1000
        endif   
        if (tload>0. .and. tload.lt.tcur) then
           tcur = tload
           goto 1000
        endif   
        read(file%m_nf) this%m_id
        read(file%m_nf) this%m_t0
        read(file%m_nf) this%m_it
        read(file%m_nf) this%m_dt
        read(file%m_nf) this%m_si
        read(file%m_nf) this%m_errorbehavior
        read(file%m_nf) this%m_type_pas


        call this%m_syspla%restorebin(file)
       end do 

1000       if (abs(tcur-tload).gt.abs(this%m_dt)) then
        write(*,*) 'time not found', tload
        write(*,*) 'last found time ', tcur
        stop 1
       else
          write(*,*) 'read data for time=', tcur
       endif 



      end subroutine  intsympbase_restorebin  

      end module mod_intsympbase

!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file sysplaode.f 
!!  \brief systeme avec N planetes
!!         pouvant etre integre par un integrateur ODE
!!
!!     cette classe t_sysplaodeode fournit les routines elementaires communes a l'ensemble 
!!      des systemes possibles
!!
!!
!
! history : creation 07/12/2015
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec N planetes 
!***********************************************************************
      module mod_sysplaode
      use mod_buffer
      use mod_ode
      use mod_arret
      
!***********************************************************************
!> @class t_sysplaode
!! classe de base decrivant 1 systeme planetaire avec N planetes
!!
!!  
!***********************************************************************
      type, extends(t_sysode) , abstract :: t_sysplaode 
          ! ces champs ne doivent etre acceder en lecture que par les fonctions exterieures au module
          
          ! description du systeme planetaire
          
          integer :: m_plan_nb  !< nombre de planetes 
          real(TREAL), dimension(:), allocatable :: m_plan_mass !< masse des planetes
          real(TREAL), dimension(:), allocatable :: m_plan_mpsbeta !< (m_p+m_star)/m_star pour chaque  planete
          real(TREAL), dimension(:,:), allocatable :: m_plan_coordX !< coordonnees (position) des planetes : (3, m_plan_nb)
          real(TREAL), dimension(:,:), allocatable :: m_plan_coordXP !< coordonnees (vitesse) des planetes : (3, m_plan_nb)
          real(TREAL), dimension(:,:), allocatable :: m_plan_errX !< erreur sur les coordonnees (position) des planetes : (3, m_plan_nb)
          real(TREAL), dimension(:,:), allocatable :: m_plan_errXP !< erreur sur les  coordonnees (vitesse) des planetes : (3, m_plan_nb)
          real(TREAL), dimension(:), allocatable :: m_plan_dataout !< tableau temporaire pour copier les planetes+le temps dans le buffer
                  
          real(TREAL) :: m_star_mass !< masse de l'etoile
          real(TREAL) :: m_cG !< constante de Gauss AU**3/an**2
          real(TREAL) :: m_fact_H = 1 !< facteur multiplicatif pour l'energie (par defaut = 1.0)
          real(TREAL) :: m_fact_C = 1 !< facteur multiplicatif pour le moment cinetique (par defaut = 1.0)

          !private
          ! buffers de sortie
          integer(8), dimension(:), allocatable     :: m_stepout !< pas de sortie (planete et/ou particule)
          type(t_buffer), dimension(:), allocatable :: m_buffer !< buffer de sortie (planete et/ou particule)
          !> buffer de sortie
          !! = 0 : buffer pour les planetes
          !! = 1 : buffer pour les particules
          !! = 2 : buffer pour les integrales premieres
          integer(8), dimension(:), allocatable     :: m_kind 
          
          type(t_arret) :: m_plan_arret !< cause de l'arret des planetes

          character(len=30) :: m_dotname = "SYSPLAODE" !< nom du convertisseur pour le graphique dot

       contains
          
          
          procedure:: set_plan_nb => sysplaode_set_plan_nb ! fixe le nombre de planetes 
          procedure:: get_plan_nb => sysplaode_get_plan_nb ! retourne le nombre de planetes 
          procedure:: set_plan_mass =>sysplaode_set_plan_mass ! fixe la masse des planetes 
          procedure:: set_plan_coord => sysplaode_set_plan_coord ! fixe les coordonnees des planetes (meme que celle integree)
          
          procedure:: set_star_mass => sysplaode_set_star_mass ! fixe la masse de l'etoile
          
          procedure :: get_outputsteps => sysplaode_get_outputsteps ! fourni les pas de sorties  (appellee par l'integrateur)
          procedure :: write_output    => sysplaode_write_output ! ecriture de la sortie (appellee par l'integrateur)
          procedure :: flush_output    => sysplaode_flush_output ! forcer a transmettre les buffers (appellee par l'integrateur)
          procedure:: copy_output_feedback=>sysplaode_copy_output_feedb ! copie de la copie de sortie vers this des informations necessaires 
         
          procedure :: add_plan_out_step => sysplaode_add_plan_out_step ! fixe le pas de sortie, la taille du buffer des planetes
          procedure :: add_int_out_step => sysplaode_add_int_out_step   ! fixe le pas de sortie, la taille du buffer des integrales

          procedure, private:: sysplaode_add_out_step ! fixe le pas de sortie, la taille du buffer 

          procedure :: clear_out_step => sysplaode_clear_out_step ! supprime tous les buffers lies en sortie au systeme
         
          procedure :: set_factor_int => sysplaode_set_factor_int  ! fixe les facteurs multiplicatifs des integrales

         ! controle des erreurs
          procedure :: set_error_time => sysplaode_set_error_time ! fixe le temps pour le controle des erreurs
          procedure :: check_error => sysplaode_check_error ! retourne true si une erreur s'est produite
          procedure :: get_error_time  => sysplaode_get_error_time ! recupere le temps de l'erreur
          procedure :: get_error_code  => sysplaode_get_error_code ! recupere le code d'erreur
          procedure :: get_error_body  => sysplaode_get_error_body ! recupere le corps generant l'erreur
          procedure :: get_error  => sysplaode_get_error ! recupere l'erreur
          procedure :: set_error  => sysplaode_set_error ! fixe l'erreur
         
          procedure :: graphdot => sysplaode_graphdot ! affiche le graphique dot
          procedure :: set_graphdotname => sysplaode_set_graphdotname ! fixe le nom du graphique dot

          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          procedure(energie), deferred :: energie
          procedure(mom_cin), deferred :: mom_cin
          
      end type t_sysplaode  

           abstract interface

            subroutine energie(this, H)
             import t_sysplaode
             class(t_sysplaode), intent(inout) :: this
             real(TREAL), intent(out) :: H 
            end subroutine energie

            subroutine mom_cin(this, C)
             import t_sysplaode
             class(t_sysplaode), intent(inout) :: this
             real(TREAL), dimension(1:3), intent(out) :: C 
            end subroutine mom_cin

          end interface

     
      contains
         

!***********************************************************************
!> @brief fixe le nombre de planetes 
!***********************************************************************
      subroutine sysplaode_set_plan_nb(this, nplan)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       integer, intent(in) :: nplan                 !< nombre de planetes
       
       this%m_plan_nb = nplan
      end  subroutine sysplaode_set_plan_nb  
      
!***********************************************************************
!> @brief retourne le nombre de planetes 
!***********************************************************************
      function sysplaode_get_plan_nb(this) result(get_plan_nb)
       implicit none
       class(t_sysplaode), intent(in):: this  !< dummy argument
       integer :: get_plan_nb       
       
       get_plan_nb = this%m_plan_nb
      end  function sysplaode_get_plan_nb  

!***********************************************************************
!> @brief fixe la masse des planetes 
!!  la fonction set_star_mass doit etre appellee auparavant
!***********************************************************************
      subroutine sysplaode_set_plan_mass(this, m)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:), intent(in) :: m  !< masse des planetes
       real(TREAL) :: m0

       allocate(this%m_plan_mass(this%m_plan_nb))
       if( any( shape(this%m_plan_mass) /= shape(m) ) ) then
          stop "Erreur : set_plan_mass : dimensions differentes"
       endif
       this%m_plan_mass = m
       allocate(this%m_plan_mpsbeta(this%m_plan_nb))
       m0 = this%m_star_mass
       this%m_plan_mpsbeta = (this%m_plan_mass+m0)/m0
      end  subroutine sysplaode_set_plan_mass  

!***********************************************************************
!> @brief fixe les coordonnees des planetes. 
!! Il n'y a aucune conversion des coordonnees.
!! Ces coordoonnees doivent etre compatibles avec l'integrateur
!***********************************************************************
      subroutine sysplaode_set_plan_coord(this, X, XP)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:,:), intent(in) :: X  !< position des planetes
       real(TREAL), dimension(:,:), intent(in) :: XP  !< moment des planetes
       
       if (.not.allocated(this%m_plan_coordX)) then
        allocate(this%m_plan_coordX(3,this%m_plan_nb))
        if( any( shape(this%m_plan_coordX) /= shape(X) ) ) then
           stop "Erreur : set_plan_coord : X : dimensions differentes"
        endif
        allocate(this%m_plan_coordXP(3,this%m_plan_nb))
        if( any( shape(this%m_plan_coordXP) /= shape(XP) ) ) then
           stop "Erreur : set_plan_coord : XP : dimensions differentes"
        endif
        allocate(this%m_plan_errX(3,this%m_plan_nb))
        allocate(this%m_plan_errXP(3,this%m_plan_nb))
       endif
       this%m_plan_coordX = X
       this%m_plan_coordXP = XP
       ! tableau pour la sommation compensee
       this%m_plan_errX = REALC(0.E0)
       this%m_plan_errXP = REALC(0.E0)
      end  subroutine sysplaode_set_plan_coord  

!***********************************************************************
!> @brief fixe la masse de l'etoile 
!***********************************************************************
      subroutine sysplaode_set_star_mass(this, cG, m0)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: m0          !< masse de l'etoile
       real(TREAL), intent(in) :: cG          !< constante de Gauss AU**3/an**2
       
       this%m_star_mass = m0
       this%m_cG = cG
      end  subroutine sysplaode_set_star_mass  
      
!***********************************************************************
!> @brief met a jour le temps d'une erreur.
!! indique le comportement du programme qui l'appelle par la suite
!! @return retourne 0 si le programme doit se poursuivre \n
!! retourne 1 si on doit s'arreter d'integrer (mais poursuite de l'execution)
!***********************************************************************
      function sysplaode_set_error_time(this, t)  result(ret)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< systeme planetaire
       real(TREAL), intent(in):: t    !< temps a laquelle s'est produite l'erreur
       integer ret ! code de retour
       
       ret = 1 ! en l'absence de particule on s'arrete
       call this%m_plan_arret%set_time(t)
                    
      end function sysplaode_set_error_time  

!***********************************************************************
!> @brief controle des erreurs du systeme planetaire
!! @return retourne false  si aucune erreur \n
!! retourne true si une erreur s'est produite depuis le dernier appel a cette fonction
!***********************************************************************
      function sysplaode_check_error(this) result(ret)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       
       logical ret ! code de retour
       ret = .false.
       
         if (this%m_plan_arret%m_stop.ne.0) then
                ret = .true.
         endif
         
      end function  sysplaode_check_error  

!***********************************************************************
!> @brief recupere le temps de l'erreur dans les donnees
!***********************************************************************
            function sysplaode_get_error_time(this) result(t)
             class(t_sysplaode), intent(inout) :: this  !< dummy argument
             real(TREAL) :: t ! temps de l'erreur
             t = this%m_plan_arret%m_time 
            end function sysplaode_get_error_time

!***********************************************************************
!> @brief recupere le corps ayant genere l'erreur dans les donnees
!***********************************************************************
            function sysplaode_get_error_body(this) result(body)
             class(t_sysplaode), intent(inout) :: this  !< dummy argument
             integer :: body ! corps
             body = this%m_plan_arret%m_body 
            end function sysplaode_get_error_body

!***********************************************************************
!> @brief recupere le code de l'erreur
!***********************************************************************
            function sysplaode_get_error_code(this) result(code)
             class(t_sysplaode), intent(inout) :: this  !< dummy argument
             integer :: code ! code
             code = this%m_plan_arret%m_cause 
            end function sysplaode_get_error_code

!***********************************************************************
!> @brief recupere l'erreur du systeme
!***********************************************************************
      function sysplaode_get_error(this) result(arret)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       
       type(t_arret) :: arret ! code de retour
       arret = this%m_plan_arret
         
      end function  sysplaode_get_error  

!***********************************************************************
!> @brief fixe l'erreur du systeme
!***********************************************************************
      subroutine sysplaode_set_error(this, arret)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       type(t_arret), intent(in) :: arret !< erreur
       
       this%m_plan_arret = arret
       
      end subroutine  sysplaode_set_error  

 

!***********************************************************************
!> @brief copie de this dans dst des informations transmises 
!! par les buffers de sortie
!! en general, les erreurs
!***********************************************************************
      subroutine sysplaode_copy_output_feedb(this, dst)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       class(t_sysode), intent(inout) :: dst  !< systeme distant

      end  subroutine sysplaode_copy_output_feedb

!***********************************************************************
!> @brief retourne les pas de sortie ou l'integrateur doit appeller la fonction 
!!        write_output
!***********************************************************************
      subroutine sysplaode_get_outputsteps(this,nstepscall,arstepscall)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       integer, intent(out) :: nstepscall  !< taille de arstepscall
       integer(8), dimension(:), allocatable, intent(out) :: arstepscall !< tableau des pas de sortie

       nstepscall = size(this%m_stepout)
       allocate(arstepscall(1:nstepscall))
       arstepscall = this%m_stepout
       write(*,*) 'nstepscall',nstepscall, arstepscall

      end  subroutine sysplaode_get_outputsteps

!***********************************************************************
!> @brief ecriture de la sortie (appellee par l'integrateur)
!!   a l'iteration donne et au temps donne
!***********************************************************************
      subroutine sysplaode_write_output(this, it, t)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: it !< iteration actuelle
       real(TREAL), intent(in) :: t  !< temps correspondant a it
       integer :: j
       integer nplan
       real(TREAL), dimension(1:5) :: integrales
       real(TREAL) :: C(3), E
       nplan = this%get_plan_nb()

       do j=1, size(this%m_stepout)
        ! write(*,*) " sysplaode_write_output", j, it, this%m_stepout(j)
        if (mod(it,this%m_stepout(j)).eq.0) then
        
         ! cas d'un buffer de planetes
         if (this%m_kind(j).eq.0) then
          ! write(*,*) " sortie planete ",it,t
          this%m_plan_dataout(1) = t
          this%m_plan_dataout(2:1+3*nplan)=reshape(this%m_plan_coordX,      &
     &                 shape(this%m_plan_dataout(1:3*nplan)))
          this%m_plan_dataout(2+3*nplan:1+6*nplan)=reshape(                 &
     &                 this%m_plan_coordXP,                                 &
     &                 shape(this%m_plan_dataout(1:3*nplan)))
          call this%m_buffer(j)%writedata(this%m_plan_dataout)
         endif
         
         ! cas d'un buffer d' integrales premieres
         if (this%m_kind(j).eq.2) then
          ! write(*,*) " sortie int ",it,t
          call this%energie(E)
          call this%mom_cin(C)
          integrales(1) = t
          integrales(2) = E*this%m_fact_H ! energie
          integrales(3:5) = C*this%m_fact_C ! moment cinetique
          call this%m_buffer(j)%writedata(integrales)
         endif
         
         ! verifie si une erreur s'est produite dans le successeur
         if (this%m_buffer(j)%check_successor_error()) then
         call this%m_buffer(j)%get_successor_error(this%m_plan_arret)
         endif
        endif
       enddo
       
      end  subroutine sysplaode_write_output

!***********************************************************************
!> @brief force a transmettre le contenu des buffers (appellee par l'integrateur)
!***********************************************************************
      subroutine sysplaode_flush_output(this)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       integer :: j

       do j=1, size(this%m_stepout)
          call this%m_buffer(j)%flushdata()
       enddo
       
      end  subroutine sysplaode_flush_output

!***********************************************************************
!> @brief routine commune aux planetes et particules pour ajouter un 
!!  nouveau buffer au tableau de ceux existants
!***********************************************************************
      subroutine sysplaode_add_out_step(this, stepout, buffersize,           &
     &                             buffercs, ncomp, pkind)
       implicit none
       class(t_sysplaode), intent(inout):: this       !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les planetes
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
       integer, intent (in) :: pkind !< =0 => buffer de planete, = 1 =buffer de particule         
       integer, intent (in) ::  ncomp !< nombre de composantes par pas de sorties
       integer newind
       integer(8), dimension(:), allocatable       :: newstepout 
       type(t_buffer), dimension(:), allocatable :: newbuffer
       integer(8), dimension(:), allocatable       :: newkind
       
       newind = 1
       
       if (allocated(this%m_stepout)) then
           newind = size(this%m_stepout)+1
           allocate(newstepout(1:newind))
           allocate(newbuffer(1:newind))
           allocate(newkind(1:newind))
           newstepout(1:newind-1) = this%m_stepout
           newbuffer(1:newind-1) = this%m_buffer
           newkind(1:newind-1) = this%m_kind
           call move_alloc(TO=this%m_stepout, FROM=newstepout)
           call move_alloc(TO=this%m_buffer, FROM=newbuffer)
           call move_alloc(TO=this%m_kind, FROM=newkind)
       else
           allocate(this%m_stepout(1:1))
           allocate(this%m_buffer(1:1))
           allocate(this%m_kind(1:1))
       endif
      
       this%m_stepout(newind) = stepout
       call this%m_buffer(newind)%init(buffersize, ncomp, buffercs) 
       this%m_kind(newind) = pkind
       
      end  subroutine sysplaode_add_out_step  


!***********************************************************************
!> @brief routine pour supprimer tous les buffers lies au systeme en sortie
!***********************************************************************
      subroutine sysplaode_clear_out_step(this)
       implicit none
       class(t_sysplaode), intent(inout):: this       !< dummy argument
       
       if (allocated(this%m_stepout)) then
        deallocate(this%m_stepout)
        deallocate(this%m_buffer)
        deallocate(this%m_kind)
       endif
       
      end  subroutine sysplaode_clear_out_step  



!***********************************************************************
!> @brief fixe le pas de sortie, la taille du buffer des planetes
!! cela cree un nouveau buffer
!! lorsque le buffer est plein, le consommateur buffercs est appelle
!***********************************************************************
      subroutine sysplaode_add_plan_out_step(this, stepout, buffersize,           &
     &                             buffercs)
       implicit none
       class(t_sysplaode), intent(inout):: this       !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les planetes
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
       integer nplan
       
       nplan = this%get_plan_nb()
       call this%sysplaode_add_out_step(stepout, buffersize,                      &
     &                               buffercs, 1+nplan*6, 0)
       if (.not.allocated(this%m_plan_dataout)) then
           allocate(this%m_plan_dataout(1:1+nplan*6))
       endif    
       
      end  subroutine sysplaode_add_plan_out_step  


!***********************************************************************
!> @brief fixe le pas de sortie, la taille du buffer des integrales premieres
!! cela cree un nouveau buffer
!! lorsque le buffer est plein, le consommateur buffercs est appelle
!***********************************************************************
      subroutine sysplaode_add_int_out_step(this, stepout, buffersize,           &
     &                             buffercs)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les inetgrales premieres
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
               
       call this%sysplaode_add_out_step(stepout, buffersize,                      &
     &             buffercs, 5 , 2) ! 5 = temps+energie +moment cin.
       
      end  subroutine sysplaode_add_int_out_step  

!***********************************************************************
!> @brief fixe les facteurs multiplicatifs des integrales premieres
!! Ces facteurs multiplicatifs sont par defaut 1
!***********************************************************************
      subroutine sysplaode_set_factor_int(this, factH, factC)
       implicit none
       class(t_sysplaode), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: factH       !< facteur mutliplicatif pour l'energie
       real(TREAL), intent(in) :: factC       !< facteur mutliplicatif pour le moment cinetique
               
       this%m_fact_H = factH
       this%m_fact_C = factC
       
      end  subroutine sysplaode_set_factor_int  


!***********************************************************************
!> @brief  fonction pour generer le graph au format dot
!***********************************************************************
       subroutine sysplaode_graphdot(this, nodesrc,nodecounter,num_file)
        implicit none
        class(t_sysplaode), intent(inout) :: this  !< dummy argument
        integer, intent(in) :: nodesrc !< numero du noeud  source 
        integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
        integer, intent(in) :: num_file !<numero du fichier de sortie
        
        integer curnode, j

        curnode = nodecounter+1
        call buffer_dot_writelabel(num_file, curnode,                         &
     &         trim(this%m_dotname))
        call buffer_dot_writenode(num_file, nodesrc,curnode)
        nodecounter = curnode

        do j=1, size(this%m_buffer)
           call this%m_buffer(j)%graphdot(curnode,nodecounter,num_file)
        enddo
        
       end subroutine sysplaode_graphdot

!***********************************************************************
!> @brief specifie le nom pour le graphique DOT utilise par this
!***********************************************************************
       subroutine sysplaode_set_graphdotname (this, dotname)
        implicit none
        class(t_sysplaode), intent(inout) :: this  !< dummy argument
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        this%m_dotname = dotname
       end subroutine sysplaode_set_graphdotname

      end module mod_sysplaode

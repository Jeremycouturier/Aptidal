!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file syspla.f 
!!  \brief systeme avec N planetes
!!         pouvant etre integre par un schema symplectique
!!
!!     cette classe t_syspla fournit les routines elementaires communes a l'ensemble 
!!      des systemes possibles
!!
!!
!
! history : creation 20/07/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec N planetes 
!***********************************************************************
      module mod_syspla
      use mod_buffer
      use mod_syssympbase
      use mod_arret
      use kahansum
      
!***********************************************************************
!> @class t_syspla
!! classe de base decrivant 1 systeme planetaire avec N planetes
!!
!!  
!***********************************************************************
      type, extends(t_syssymp) , abstract :: t_syspla 
          ! ces champs ne doivent etre acceder en lecture que par les fonctions exterieures au module
          
          ! description du systeme planetaire
          
          integer :: m_plan_nb  !< nombre de planetes 
          real(TREAL), dimension(:), allocatable :: m_plan_mass !< masse des planetes
          real(TREAL), dimension(:), allocatable :: m_plan_mass0    !< masse etoile et masse planete 
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
          type(compensatedsum) :: m_Hinit      !< Energie initiale pour la sortie de delta Energie
          type(compensatedsum), dimension(1:3) :: m_Cinit      !< Energie initiale pour la sortie de delta moment cinetique
          !private
          ! buffers de sortie
          integer(8), dimension(:), allocatable     :: m_stepout !< pas de sortie (planete et/ou particule)
          type(t_buffer), dimension(:), allocatable :: m_buffer !< buffer de sortie (planete et/ou particule)
          !> buffer de sortie
          !! = 0 : buffer pour les planetes
          !! = 1 : buffer pour les particules
          !! = 2 : buffer pour les integrales premieres
          integer(8), dimension(:), allocatable     :: m_kind 
          !> sortie au pas de sortie ou au pas symplectique
          !! = true : sortie aux pas normaux
          !! = false  : sortie au pas symplectique
          logical, dimension(:), allocatable     :: m_stepkind 

          logical :: m_localstepkind = .false. !< indique quels m_stepkind seront selectionnes
          
          type(t_arret) :: m_plan_arret !< cause de l'arret des planetes

          character(len=30) :: m_dotname = "SYSPLA" !< nom du convertisseur pour le graphique dot

       contains
          
          
          procedure:: set_plan_nb ! fixe le nombre de planetes 
          procedure:: get_plan_nb ! retourne le nombre de planetes 
          procedure:: set_plan_mass =>syspla_set_plan_mass ! fixe la masse des planetes 
          procedure:: set_plan_coord => syspla_set_plan_coord ! fixe les coordonnees des planetes (meme que celle integree)
          
          procedure:: set_star_mass ! fixe la masse de l'etoile
          
          procedure:: create_output  => syspla_create_output! creation d'une copie pour effectuer les pas de sortie
          procedure:: copy_output    => syspla_copy_output ! copie de this dans la copie de sortie
          procedure, NOPASS:: delete_output  => syspla_delete_output ! destruction de l'objet cree par create_output
          
          procedure :: get_outputsteps => syspla_get_outputsteps ! fourni les pas de sorties  (appelee par l'integrateur)
          procedure :: write_output    => syspla_write_output ! ecriture de la sortie (appelee par l'integrateur)
          procedure :: flush_output    => syspla_flush_output ! forcer a transmettre les buffers (appelee par l'integrateur)
          procedure:: copy_output_feedback=>syspla_copy_output_feedback ! copie de la copie de sortie vers this des informations necessaires 
         
          procedure :: add_plan_out_step => syspla_add_plan_out_step ! fixe le pas de sortie, la taille du buffer des planetes (sortie normale)
          procedure:: add_plan_out_stepmid=>syspla_add_plan_out_stepmid ! fixe le pas de sortie, la taille du buffer des planetes (sortie symplectique)
          procedure :: add_int_out_step => syspla_add_int_out_step   ! fixe le pas de sortie, la taille du buffer des integrales

          procedure :: add_out_step  => syspla_add_out_step ! fixe le pas de sortie, la taille du buffer 

          procedure :: clear_out_step => syspla_clear_out_step ! supprime tous les buffers lies en sortie au systeme
          procedure :: isvalidoutput => syspla_isvalidoutput  ! retourne true si c'est un pas de sortie
         
          procedure :: set_factor_int => syspla_set_factor_int  ! fixe les facteurs multiplicatifs des integrales

         ! controle des erreurs
          procedure :: set_error_time => syspla_set_error_time ! fixe le temps pour le controle des erreurs
          procedure :: check_error => syspla_check_error ! retourne true si une erreur s'est produite
          procedure :: get_error_time  => syspla_get_error_time ! recupere le temps de l'erreur
          procedure :: get_error_code  => syspla_get_error_code ! recupere le code d'erreur
          procedure :: get_error_body  => syspla_get_error_body ! recupere le corps generant l'erreur
          procedure :: get_error  => syspla_get_error ! recupere l'erreur
          procedure :: set_error  => syspla_set_error ! fixe l'erreur
         
          procedure :: graphdot => syspla_graphdot ! affiche le graphique dot
          procedure :: set_graphdotname => syspla_set_graphdotname ! fixe le nom du graphique dot

          
          procedure :: energie_with_err => syspla_energie_with_err
          procedure :: mom_cin_with_err => syspla_momcin_with_err

       
          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
!          procedure(pasA), deferred :: pasA
!          procedure(pasA), deferred :: pasB
!          procedure(pasA), deferred :: pasC
          procedure(energie), deferred :: energie
          procedure(mom_cin), deferred :: mom_cin
          
          procedure :: dumpbin =>  syspla_dumpbin ! dump binaire vers un fichier
          procedure :: restorebin => syspla_restorebin ! restauration binaire depuis un fichier

      end type t_syspla  

           abstract interface

            subroutine energie(this, H)
             import t_syspla
             class(t_syspla), intent(inout) :: this
             real(TREAL), intent(out) :: H 
             end subroutine energie


            subroutine mom_cin(this, C)
             import t_syspla
             class(t_syspla), intent(inout) :: this
             real(TREAL), dimension(1:3), intent(out) :: C 
            end subroutine mom_cin

          end interface

     
      contains
         

!***********************************************************************
!> @brief fixe le nombre de planetes 
!***********************************************************************
      subroutine set_plan_nb(this, nplan)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       integer, intent(in) :: nplan                 !< nombre de planetes
       
       this%m_plan_nb = nplan
      end  subroutine set_plan_nb  
      
!***********************************************************************
!> @brief retourne le nombre de planetes 
!***********************************************************************
      function get_plan_nb(this)
       implicit none
       class(t_syspla), intent(in):: this  !< dummy argument
       integer :: get_plan_nb       
       
       get_plan_nb = this%m_plan_nb
      end  function get_plan_nb  

!***********************************************************************
!> @brief fixe la masse des planetes 
!!  la fonction set_star_mass doit etre appellee auparavant
!***********************************************************************
      subroutine syspla_set_plan_mass(this, m)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
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

       allocate(this%m_plan_mass0(0:this%m_plan_nb))
       this%m_plan_mass0(0) = m0
       this%m_plan_mass0(1:) = this%m_plan_mass

      end  subroutine syspla_set_plan_mass  

!***********************************************************************
!> @brief fixe les coordonnees des planetes. 
!! Il n'y a aucune conversion des coordonnees.
!! Ces coordoonnees doivent etre compatibles avec l'integrateur
!***********************************************************************
      subroutine syspla_set_plan_coord(this, X, XP)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
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

       call this%energie_with_err(this%m_Hinit)
       call this%mom_cin_with_err(this%m_Cinit)
       
      end  subroutine syspla_set_plan_coord  

!***********************************************************************
!> @brief fixe la masse de l'etoile 
!***********************************************************************
      subroutine set_star_mass(this, cG, m0)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: m0          !< masse de l'etoile
       real(TREAL), intent(in) :: cG          !< constante de Gauss AU**3/an**2
       
       this%m_star_mass = m0
       this%m_cG = cG
      end  subroutine set_star_mass  
      
!***********************************************************************
!> @brief met a jour le temps d'une erreur.
!! indique le comportement du programme qui l'appelle par la suite
!! @return retourne 0 si le programme doit se poursuivre \n
!! retourne 1 si on doit s'arreter d'integrer (mais poursuite de l'execution)
!***********************************************************************
      function syspla_set_error_time(this, t)  result(ret)
       implicit none
       class(t_syspla), intent(inout):: this  !< systeme planetaire
       real(TREAL), intent(in):: t    !< temps a laquelle s'est produite l'erreur
       integer ret ! code de retour
       
       ret = 1 ! en l'absence de particule on s'arrete
       call this%m_plan_arret%set_time(t)
                    
      end function syspla_set_error_time  

!***********************************************************************
!> @brief controle des erreurs du systeme planetaire
!! @return retourne false  si aucune erreur \n
!! retourne true si une erreur s'est produite depuis le dernier appel a cette fonction
!***********************************************************************
      function syspla_check_error(this) result(ret)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       
       logical ret ! code de retour
       ret = .false.
       
         if (this%m_plan_arret%m_stop.ne.0) then
                ret = .true.
         endif
         
      end function  syspla_check_error  

!***********************************************************************
!> @brief recupere le temps de l'erreur dans les donnees
!***********************************************************************
            function syspla_get_error_time(this) result(t)
             class(t_syspla), intent(inout) :: this  !< dummy argument
             real(TREAL) :: t ! temps de l'erreur
             t = this%m_plan_arret%m_time 
            end function syspla_get_error_time

!***********************************************************************
!> @brief recupere le corps ayant genere l'erreur dans les donnees
!***********************************************************************
            function syspla_get_error_body(this) result(body)
             class(t_syspla), intent(inout) :: this  !< dummy argument
             integer :: body ! corps
             body = this%m_plan_arret%m_body 
            end function syspla_get_error_body

!***********************************************************************
!> @brief recupere le code de l'erreur
!***********************************************************************
            function syspla_get_error_code(this) result(code)
             class(t_syspla), intent(inout) :: this  !< dummy argument
             integer :: code ! code
             code = this%m_plan_arret%m_cause 
            end function syspla_get_error_code

!***********************************************************************
!> @brief recupere l'erreur du systeme
!***********************************************************************
      function syspla_get_error(this) result(arret)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       
       type(t_arret) :: arret ! code de retour
       arret = this%m_plan_arret
         
      end function  syspla_get_error  

!***********************************************************************
!> @brief fixe l'erreur du systeme
!***********************************************************************
      subroutine syspla_set_error(this, arret)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       type(t_arret), intent(in) :: arret !< erreur
       
       this%m_plan_arret = arret
       
      end subroutine  syspla_set_error  

!***********************************************************************
!> @brief creation d'une copie pour effectuer les pas de sortie
!! en general, duplique simplement this
!***********************************************************************
      subroutine syspla_create_output(this, sysout)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       class(t_syssymp), allocatable, intent(out) :: sysout  !< systeme cree pour la sortie
       
       allocate(sysout, SOURCE=this)
       select type(sysout)
        class is(t_syspla)
         sysout%m_localstepkind = .true.
         deallocate(sysout%m_plan_mass)
         deallocate(sysout%m_plan_mpsbeta)
         deallocate(sysout%m_plan_coordX)
         deallocate(sysout%m_plan_coordXP)
         deallocate(sysout%m_plan_errX)
         deallocate(sysout%m_plan_errXP)
         allocate(sysout%m_plan_mass(this%m_plan_nb),                    &
     &    SOURCE=this%m_plan_mass)
         allocate(sysout%m_plan_mpsbeta(this%m_plan_nb),                 &
     &     SOURCE=this%m_plan_mpsbeta)
         allocate(sysout%m_plan_coordX(3,this%m_plan_nb),                &
     &    SOURCE=this%m_plan_coordX)
         allocate(sysout%m_plan_coordXP(3,this%m_plan_nb),               &
     &     SOURCE=this%m_plan_coordXP)
         allocate(sysout%m_plan_errX(3,this%m_plan_nb),                  &
     &     SOURCE=this%m_plan_errX)
         allocate(sysout%m_plan_errXP(3,this%m_plan_nb),                 &
     &     SOURCE=this%m_plan_errXP)
        class default
         stop 'Erreur: this ne derive pas de syspla_create_output' 
      end select

      end  subroutine syspla_create_output
        
!***********************************************************************
!> @brief destruction de l'objet cree par create_output
!! en general, detruit  this
!***********************************************************************
      subroutine syspla_delete_output(this)
       implicit none
       class(t_syssymp), allocatable, intent(inout):: this  !< dummy argument
       
       deallocate(this)

      end  subroutine syspla_delete_output
      
!***********************************************************************
!> @brief copie de this dans la copie de sortie (sysout)
!! en general, sysout = this 
!***********************************************************************
      subroutine syspla_copy_output(this, sysout)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       class(t_syssymp), intent(inout) :: sysout  !< systeme cree pour la sortie

       select type(sysout)
        class is(t_syspla)
         sysout%m_plan_coordX(:,:) = this%m_plan_coordX(:,:)
         sysout%m_plan_coordXP(:,:) = this%m_plan_coordXP(:,:)
         sysout%m_plan_errX(:,:) = this%m_plan_errX(:,:)
         sysout%m_plan_errXP(:,:) = this%m_plan_errXP(:,:)
         sysout%m_Hinit = this%m_Hinit
         sysout%m_Cinit = this%m_Cinit
         
        class default
         stop 'Erreur: this ne derive pas de syspla_copy_output' 
      end select
      end  subroutine syspla_copy_output
 

!***********************************************************************
!> @brief copie de this dans dst desinformations transmises 
!! par les buffers de sortie
!! en general, les erreurs
!***********************************************************************
      subroutine syspla_copy_output_feedback(this, dst)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       class(t_syssymp), intent(inout) :: dst  !< systeme distant

      end  subroutine syspla_copy_output_feedback

!***********************************************************************
!> @brief retourne les pas de sortie ou l'integrateur doit appeler la fonction 
!!        write_output
!***********************************************************************
      subroutine syspla_get_outputsteps(this, nstepscall, arstepscall)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       integer, intent(out) :: nstepscall  !< taille de arstepscall
       integer(8), dimension(:), allocatable, intent(out) :: arstepscall !< tableau des pas de sortie
       integer j, k

       nstepscall = 0
       do j=1,size(this%m_stepout)
            if (this%m_localstepkind.eqv.this%m_stepkind(j)) then
            nstepscall = nstepscall+1
            endif
       enddo      
 
       !nstepscall = size(this%m_stepout)
       if (nstepscall.gt.0) then
            allocate(arstepscall(1:nstepscall))
            !arstepscall = this%m_stepout
            k=1
            do j=1,size(this%m_stepout)
                  if (this%m_localstepkind.eqv.this%m_stepkind(j)) then
                  arstepscall(k) = this%m_stepout(j)
                  k=k+1
                  endif
            enddo
       endif
       write(*,*) 'nstepscall',nstepscall, arstepscall

      end  subroutine syspla_get_outputsteps


!***********************************************************************
!> @brief retourne true si c'est un pas de sortie pour ce systeme
!***********************************************************************
      function syspla_isvalidoutput(this, it ,j ) result(r)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: it !< iteration actuelle
       integer, intent(in) :: j !< numero du buffer
       logical :: r

        r = .false.
        if ((mod(it,this%m_stepout(j)).eq.0).and.                            &
     &     (this%m_stepkind(j).eqv.this%m_localstepkind)) then
           r = .true.
        endif  
      end  function syspla_isvalidoutput


!***********************************************************************
!> @brief ecriture de la sortie (appelee par l'integrateur)
!!   a l'iteration donne et au temps donne
!***********************************************************************
      subroutine syspla_write_output(this, it, t)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: it !< iteration actuelle
       real(TREAL), intent(in) :: t  !< temps correspondant a it
       integer :: j, k
       integer nplan
       real(TREAL), dimension(1:5) :: integrales
       type(compensatedsum), dimension(1:3) :: C
       type(compensatedsum) :: E
       
       nplan = this%get_plan_nb()

       if (allocated(this%m_stepout)) then
       do j=1, size(this%m_stepout)
         ! write(*,*) " syspla_write_output", j, it, this%m_stepout(j)
!        if (modX(it,this%m_stepout(j)).eq.0) then
        if (this%isvalidoutput(it,j)) then
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
          call this%energie_with_err(E)
          call this%mom_cin_with_err(C)
          if (it.ne.0) then
          E = E-this%m_Hinit
          do k=1,3
            C(k) = C(k)-this%m_Cinit(k)
          enddo
          end if
          integrales(1) = t
          integrales(2) = E%v*this%m_fact_H ! energie
          do k=1,3
          integrales(2+k) = C(k)%v*this%m_fact_C ! moment cinetique
          enddo
          call this%m_buffer(j)%writedata(integrales)
         endif
         
         ! verifie si une erreur s'est produite dans le successeur
         if (this%m_buffer(j)%check_successor_error()) then
         call this%m_buffer(j)%get_successor_error(this%m_plan_arret)
         endif
        endif
       enddo
       endif
      end  subroutine syspla_write_output

!***********************************************************************
!> @brief force a transmettre le contenu des buffers (appelee par l'integrateur)
!***********************************************************************
      subroutine syspla_flush_output(this)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       integer :: j

       if (allocated(this%m_stepout)) then
        do j=1, size(this%m_stepout)
            if (this%m_stepkind(j).eqv.this%m_localstepkind) then
            call this%m_buffer(j)%flushdata()
            endif
        enddo
       endif
       
      end  subroutine syspla_flush_output

!***********************************************************************
!> @brief routine commune aux planetes et particules pour ajouter un 
!!  nouveau buffer au tableau de ceux existants
!***********************************************************************
      subroutine syspla_add_out_step(this, stepout, buffersize,           &
     &                             buffercs, ncomp, pkind)
       implicit none
       class(t_syspla), intent(inout):: this       !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les planetes
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
       integer, intent (in) :: pkind !< =0 => buffer de planete, = 1 =buffer de particule         
       integer, intent (in) ::  ncomp !< nombre de composantes par pas de sorties
       integer newind
       integer(8), dimension(:), allocatable       :: newstepout 
       type(t_buffer), dimension(:), allocatable :: newbuffer
       integer(8), dimension(:), allocatable       :: newkind
       logical, dimension(:), allocatable       :: newstepkind
       
       newind = 1
       
       if (allocated(this%m_stepout)) then
           newind = size(this%m_stepout)+1
           allocate(newstepout(1:newind))
           allocate(newbuffer(1:newind))
           allocate(newkind(1:newind))
           allocate(newstepkind(1:newind))
           newstepout(1:newind-1) = this%m_stepout
           newbuffer(1:newind-1) = this%m_buffer
           newkind(1:newind-1) = this%m_kind
           newstepkind(1:newind-1) = this%m_stepkind
           call move_alloc(TO=this%m_stepout, FROM=newstepout)
           call move_alloc(TO=this%m_buffer, FROM=newbuffer)
           call move_alloc(TO=this%m_kind, FROM=newkind)
           call move_alloc(TO=this%m_stepkind, FROM=newstepkind)
       else
           allocate(this%m_stepout(1:1))
           allocate(this%m_buffer(1:1))
           allocate(this%m_kind(1:1))
           allocate(this%m_stepkind(1:1))
       endif
      
       this%m_stepout(newind) = stepout
       call this%m_buffer(newind)%init(buffersize, ncomp, buffercs) 
       this%m_kind(newind) = pkind
       this%m_stepkind(newind) = .true.
       
      end  subroutine syspla_add_out_step  


!***********************************************************************
!> @brief routine pour supprimer tous les buffers lies au systeme en sortie
!***********************************************************************
      subroutine syspla_clear_out_step(this)
       implicit none
       class(t_syspla), intent(inout):: this       !< dummy argument
       
       if (allocated(this%m_stepout)) then
        deallocate(this%m_stepout)
        deallocate(this%m_buffer)
        deallocate(this%m_kind)
        deallocate(this%m_stepkind)
       endif
       
      end  subroutine syspla_clear_out_step  



!***********************************************************************
!> @brief fixe le pas de sortie, la taille du buffer des planetes
!! cela cree un nouveau buffer
!! lorsque le buffer est plein, le consommateur buffercs est appelle
!! pour ue sortie normale
!***********************************************************************
      subroutine syspla_add_plan_out_step(this, stepout, buffersize,           &
     &                             buffercs)
       implicit none
       class(t_syspla), intent(inout):: this       !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les planetes
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
       integer nplan
       
       nplan = this%get_plan_nb()
       call this%add_out_step(stepout, buffersize,                      &
     &                               buffercs, 1+nplan*6, 0)
       if (.not.allocated(this%m_plan_dataout)) then
           allocate(this%m_plan_dataout(1:1+nplan*6))
       endif    
       
      end  subroutine syspla_add_plan_out_step  

!***********************************************************************
!> @brief fixe le pas de sortie, la taille du buffer des planetes
!! cela cree un nouveau buffer
!! lorsque le buffer est plein, le consommateur buffercs est appelle
!! si stepout <0 alors on sort sur un pas symplectique
!! si stepout >0 alors on sort sur un pas de sortie normale
!***********************************************************************
      subroutine syspla_add_plan_out_stepmid(this, stepout, buffersize,           &
     &                             buffercs)
       implicit none
       class(t_syspla), intent(inout):: this       !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les planetes
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
       
       call this%add_plan_out_step(abs(stepout),buffersize,buffercs)
       if (stepout.lt.0) then
            this%m_stepkind(size(this%m_stepkind))=.false.
       endif

      end  subroutine syspla_add_plan_out_stepmid  


!***********************************************************************
!> @brief fixe le pas de sortie, la taille du buffer des integrales premieres
!! cela cree un nouveau buffer
!! lorsque le buffer est plein, le consommateur buffercs est appelle
!***********************************************************************
      subroutine syspla_add_int_out_step(this, stepout, buffersize,           &
     &                             buffercs)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les inetgrales premieres
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
               
       call this%add_out_step(stepout, buffersize,                      &
     &             buffercs, 5 , 2) ! 5 = temps+energie +moment cin.
       
      end  subroutine syspla_add_int_out_step  

!***********************************************************************
!> @brief fixe les facteurs multiplicatifs des integrales premieres
!! Ces facteurs multiplicatifs sont par defaut 1
!***********************************************************************
      subroutine syspla_set_factor_int(this, factH, factC)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: factH       !< facteur mutliplicatif pour l'energie
       real(TREAL), intent(in) :: factC       !< facteur mutliplicatif pour le moment cinetique
               
       this%m_fact_H = factH
       this%m_fact_C = factC
       
      end  subroutine syspla_set_factor_int  
      
       
!***********************************************************************
!> @brief  fonction pour calculer l'energie
!***********************************************************************
      subroutine syspla_energie_with_err(this, H)
       implicit none
       class(t_syspla), intent(inout) :: this
       type(compensatedsum), intent(out) :: H
       real(TREAL) :: Hv
       call this%energie(Hv)
       H%v = Hv
       H%err = 0.
      end subroutine syspla_energie_with_err

!***********************************************************************
!> @brief  fonction pour calculer le moment cinetique
!***********************************************************************
      subroutine syspla_momcin_with_err(this, C)
       implicit none
       class(t_syspla), intent(inout) :: this
       type(compensatedsum),dimension(1:3), intent(out) :: C
       real(TREAL),dimension(1:3) :: Cv
       integer j
       call this%mom_cin(Cv)
       do j=1,3
            C(j)%v = Cv(j)
            C(j)%err = 0.
       enddo
      end subroutine syspla_momcin_with_err


!***********************************************************************
!> @brief  fonction pour generer le graph au format dot
!***********************************************************************
       subroutine syspla_graphdot(this, nodesrc,nodecounter,num_file)
        implicit none
        class(t_syspla), intent(inout) :: this  !< dummy argument
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
        
       end subroutine syspla_graphdot

!***********************************************************************
!> @brief specifie le nom pour le graphique DOT utilise par this
!***********************************************************************
       subroutine syspla_set_graphdotname (this, dotname)
        implicit none
        class(t_syspla), intent(inout) :: this  !< dummy argument
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        this%m_dotname = dotname
       end subroutine syspla_set_graphdotname

!***********************************************************************
!> @brief dump binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (in) objet a sauvegarder
!***********************************************************************
      subroutine syspla_dumpbin(this, file) 
       use mod_io_dump
       implicit none
       class(t_syspla), intent(in):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       call syssympbase_dumpbin(this, file)

       write(file%m_nf) this%m_plan_nb 
       write(file%m_nf) this%m_plan_mass 
       write(file%m_nf) this%m_plan_mass0 
       write(file%m_nf) this%m_plan_mpsbeta 
       write(file%m_nf) this%m_plan_coordX 
       write(file%m_nf) this%m_plan_coordXP 
       write(file%m_nf) this%m_plan_errX 
       write(file%m_nf) this%m_plan_errXP 
                  
       write(file%m_nf) this%m_star_mass
       write(file%m_nf) this%m_cG
       write(file%m_nf) this%m_fact_H
       write(file%m_nf) this%m_fact_C
       write(file%m_nf) this%m_Hinit 
       write(file%m_nf) this%m_Cinit 
       
      end subroutine  syspla_dumpbin  

!***********************************************************************
!> @brief restauration binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (inout) objet a restaurer
!***********************************************************************
      subroutine syspla_restorebin(this, file) 
       use mod_io_dump
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       call syssympbase_restorebin(this, file)

       read(file%m_nf) this%m_plan_nb 
       if (.not.allocated(this%m_plan_mass)) then
            allocate(this%m_plan_mass(this%m_plan_nb))
            allocate(this%m_plan_mass0(0:this%m_plan_nb))
            allocate(this%m_plan_mpsbeta(this%m_plan_nb))
            allocate(this%m_plan_coordX(3,this%m_plan_nb))
            allocate(this%m_plan_coordXP(3,this%m_plan_nb))
            allocate(this%m_plan_errX(3,this%m_plan_nb))
            allocate(this%m_plan_errXP(3,this%m_plan_nb))
       endif
       read(file%m_nf) this%m_plan_mass 
       read(file%m_nf) this%m_plan_mass0 
       read(file%m_nf) this%m_plan_mpsbeta 
       read(file%m_nf) this%m_plan_coordX 
       read(file%m_nf) this%m_plan_coordXP 
       read(file%m_nf) this%m_plan_errX 
       read(file%m_nf) this%m_plan_errXP 
                  
       read(file%m_nf) this%m_star_mass
       read(file%m_nf) this%m_cG
       read(file%m_nf) this%m_fact_H
       read(file%m_nf) this%m_fact_C
       read(file%m_nf) this%m_Hinit 
       read(file%m_nf) this%m_Cinit 

      end subroutine  syspla_restorebin  


      end module mod_syspla

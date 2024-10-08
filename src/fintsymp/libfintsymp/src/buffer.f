!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file buffer.f 
!!  \brief gestion des buffers 
!!
!!
!
! history : creation 24/07/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
!> module pour la gestion des buffers
!***********************************************************************
      module mod_buffer
       use mod_arret
       
!***********************************************************************
!> @class t_buffer_consumer
!! classe pour gerer la consommation des donnees d'un buffer.
!!
!! la fonction onbufferfull : est appellee lorsque le buffer est plein.
!! la fonction graphdot : est appellee pour la generation de graphique dot.
!! la fonction onbufferflush : est appellee lorsque le contenu du buffer doit etre propage (fin d'une execution par exemple).
!***********************************************************************
      type t_buffer_consumer 
          !private
          
         procedure(t_procbufferfull),NOPASS,pointer::m_cb_onbufferfull !< fonction appellee lorsque le buffer est plein
         procedure(t_procbufferfull),NOPASS,pointer::m_cb_onbufferflush !< fonction appellee lorsque le contenu du buffer doit etre propage
         procedure(t_procgraphdot), NOPASS, pointer::m_cb_graphdot !< fonction appellee pour generer le graphique au format dot
         class(*), pointer :: m_owner !< proprietaire de ce consommateur         

       contains
          private
                    
          procedure, public:: set_owner => buffer_consumer_setowner ! fixe le proprietaire de ce consommateur
          

      end type t_buffer_consumer  
      
!***********************************************************************
!> @class t_buffer
!! classe de base decrivant 1 buffer
!!
!***********************************************************************
      type t_buffer 
          private
          
          integer(8) :: m_stepmem  !< nombre de pas de sortie a memoriser avant appel du callback 
          integer    :: m_reals    !< nombre de reels par pas de sortie 
          integer(8) :: m_write    !< position de l'ecriture dans le buffer
          real(TREAL), dimension(:), allocatable :: m_array !< donnees du buffer
          procedure(t_procbufferfull), pointer :: m_callbackbufferfull   &
     &               =>NULL() !< fonction appellee lorsque le buffer est plein
          procedure(t_procbufferfull), pointer :: m_callbackbufferflush   &
     &               =>NULL() !< fonction appellee lorsque le contenu du buffer doit etre propage
          procedure(t_procgraphdot), pointer :: m_callbackgraphdot !< fonction appellee pour le graphique dot
          class(*), pointer :: m_userdata !< donnee utilisateur fournie lors de l'appel du callback         

          type(t_arret) :: m_successor_arret !< cause de l'erreur dans les sucesseurs qui provoque l'arret du traitement en cours
          
          type(t_arret), dimension(:), allocatable, public :: m_mce  !< tableau de controles d'erreurs mais ne provoquant pas l'arret (en general 1 par coprs, e.g., 1 pour chaque particule)
          integer :: m_count_multictrlerr !< nombre d'elements de  m_multiplectrlerrors qui ont un stop!=0 ( nombre de quantite en erreurs)

       contains
          procedure :: init => buffer_init ! initialise le buffer 
          procedure :: writedata => buffer_writedata ! ecriture de donnees dans le buffer
          procedure :: readdata => buffer_readdata ! lecture de donnees dans le buffer
          procedure :: flushdata => buffer_flushdata ! force la propagation des donnees dans le buffer
          procedure :: empty => buffer_empty ! vide le buffer (position d'ecriture mis au debut) mais ne libere pas la memoire

          !gestion des erreurs des successeurs : cas de l'erreur globale
          procedure :: notify_successor_error  =>                        &
     &                 buffer_notify_successor_error ! le successeur indique une erreur
          procedure :: check_successor_error =>                          &
     &                 buffer_check_successor_error ! retourne true si une erreur s'est produite dans le successeur
          procedure :: get_successor_error  =>                           &
     &                 buffer_get_successor_error ! recupere l'erreur du successeur
          procedure :: checkget_successor_error =>                          &
     &                 buffer_checkget_successor_error ! recupere l'erreur du successeur uniquement si une erreur s'est produite dans le successeur
       
          !gestion des erreurs de m_mce : cas de l'erreur par corps (e.g. particules)
          procedure :: init_multictrlerr => buffer_init_multictrlerr ! initialise un tableau de gestion d'erreurs (en general, 1 par corps)
          procedure :: set_multictrlerr_buf=>buffer_set_multictrlerr_buf  ! met a jour le tableau des erreurs du buffer  a partir du buffer source 
          procedure :: set_multictrlerr  => buffer_set_multictrlerr  ! met a jour le tableau des erreurs du buffer
          procedure :: get_multictrlerr  => buffer_get_multictrlerr  ! recupere le tableau des erreurs du buffer
          procedure :: tryset_multictrlerr_buf =>                          &
     &                    buffer_tryset_multictrlerr_buf  ! allocation si besoin et  met a jour le tableau des erreurs du buffer a partir du buffer source 
          procedure :: update_count_multictrlerr   =>                        &
     &                 buffer_update_count_multictrlerr  ! met a jour le nombre total d'erreur et retourne true si tous les corps sont en erreur
          
          procedure :: graphdot => buffer_graphdot ! affiche le graphique dot

      end type t_buffer  

      abstract interface
       !> prototype des fonctions appellees lorsque le buffer est plein
       subroutine t_procbufferfull(this, userdata, poswrite)
        import t_buffer
        class(t_buffer), intent(inout) :: this !< buffer ayant genere l'appel
        class(*), intent(inout) :: userdata !< donnee utilisateur
        integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
       end subroutine t_procbufferfull

       !> prototype des fonctions appellees pour generer le graph au format dot
      subroutine t_procgraphdot(this, userdata,nodesrc,nodecounter,       &
     &             num_file)
        import t_buffer
        class(t_buffer), intent(inout) :: this !< buffer ayant genere l'appel
        class(*), intent(inout) :: userdata !< donnee utilisateur
        integer, intent(in) :: nodesrc !< numero du noeud  source 
        integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
        integer, intent(in) :: num_file !<numero du fichier de sortie
      end subroutine t_procgraphdot
      end interface

     
      contains

!***********************************************************************
!> @brief fixe le proprietaire de ce consommateur et ses fonctions callback
!***********************************************************************
      subroutine buffer_consumer_setowner(this, parent,callbackonfull,   &
     &     callbackonflush, callbackgraphdot)
       implicit none
       class(t_buffer_consumer), intent(inout):: this  !< consommateur du buffer : dummy argument
       class(*), target, intent(inout) :: parent !< proprietaire du consommateur
       procedure(t_procbufferfull),pointer,intent(in)::callbackonfull !< fonction a appeler lorsque le buffer est plein
       procedure(t_procbufferfull),pointer,intent(in)::callbackonflush !< fonction a appeler lorsque le contenu du buffer doit etre obilgatoirement propage
       procedure(t_procgraphdot),pointer,intent(in)::callbackgraphdot !< fonction a appeler pour la generation du graph dot
       
       this%m_owner => parent
       this%m_cb_onbufferfull => callbackonfull
       this%m_cb_onbufferflush => callbackonflush
       this%m_cb_graphdot => callbackgraphdot
       
      end subroutine  buffer_consumer_setowner


!***********************************************************************
!> @brief initialise le buffer  et alloue la memoire
!***********************************************************************
      subroutine buffer_init(this, stepmem, nreal, buffercs)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       integer(8), intent(in) :: stepmem   !< nombre de pas de sortie a memoriser avant appel du callback
       integer, intent(in) :: nreal        !< nombre de reels par pas de sortie
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
       
       this%m_stepmem = stepmem
       this%m_reals = nreal
       this%m_callbackbufferfull => buffercs%m_cb_onbufferfull
       this%m_callbackbufferflush => buffercs%m_cb_onbufferflush
       this%m_callbackgraphdot => buffercs%m_cb_graphdot
       this%m_userdata => buffercs%m_owner
       this%m_write = 1
       allocate(this%m_array(1:stepmem*nreal))
      end  subroutine buffer_init  

!***********************************************************************
!> @brief ecriture de donnees dans le buffer. 
!! Si le buffer le est plein, cela declenche le callback
!***********************************************************************
      subroutine buffer_writedata(this, rdata)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       real(TREAL), dimension(:), intent(in) :: rdata   !< "m_reals" donnees a ecrire a la position courante
       
       this%m_array((this%m_write-1)*this%m_reals+1:                     &
     &              this%m_write*this%m_reals)  = rdata
       this%m_write = this%m_write+1
       if (this%m_write.gt.this%m_stepmem) then
        call this%m_callbackbufferfull(this%m_userdata,this%m_write-1)
       endif
      end  subroutine buffer_writedata  

!***********************************************************************
!> @brief force le declenchement du callback si on a au moins une donnee
!***********************************************************************
      subroutine buffer_flushdata(this)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       
       if (associated(this%m_callbackbufferflush)) then
        call this%m_callbackbufferflush(this%m_userdata,this%m_write-1)
       endif
      end  subroutine buffer_flushdata  

!***********************************************************************
!> @brief lecture de donnees dans le buffer a la position j. 
!! le premier indice possible est 1
!***********************************************************************
      subroutine buffer_readdata(this, j, rdata)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       integer(8), intent(in) :: j            !< position des donnees a lire
       real(TREAL), dimension(:), intent(out) :: rdata   !< "m_reals" donnees a ecrire a la position courante
       
       rdata = this%m_array((j-1)*this%m_reals+1:j*this%m_reals)
       end  subroutine buffer_readdata  

!***********************************************************************
!> @brief vide le buffer (position d'ecriture mis au debut) 
!! mais ne libere pas la memoire  
!***********************************************************************
      subroutine buffer_empty(this)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       
       this%m_write = 1
      end  subroutine buffer_empty  


!***********************************************************************
!> @brief indique au buffer courant qu'une erreur s'est produite
!! dans le successeur 
!***********************************************************************
      subroutine buffer_notify_successor_error(this, arret)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       type(t_arret), intent(in):: arret  !< cause de l'arret du successeur
      
       this%m_successor_arret = arret
      end  subroutine buffer_notify_successor_error  

!***********************************************************************
!> @brief verifie si une erreur s'est produite dans le successeur
!! @return retourne false  si aucune erreur \n
!***********************************************************************
      function buffer_check_successor_error(this) result(ret)
       implicit none
       class(t_buffer), intent(inout) :: this  !< buffer : dummy argument
       
       logical ret ! code de retour
       ret = .false.
       
         if (this%m_successor_arret%m_stop.ne.0) then
                ret = .true.
         endif
         
      end function  buffer_check_successor_error  

!***********************************************************************
!> @brief retourne l'erreur du sucesseur 
!***********************************************************************
      subroutine buffer_get_successor_error(this, arret)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       type(t_arret), intent(out):: arret  !< cause de l'arret du successeur
      
       arret = this%m_successor_arret
      end  subroutine buffer_get_successor_error  

!***********************************************************************
!> @brief verifie si une erreur s'est produite dans le successeur
!! et la recupere dans this.
!! s'il n'y a pas d'erreur, on ne fait rien.
!! @return retourne false  si aucune erreur \n
!***********************************************************************
      subroutine buffer_checkget_successor_error(this, sucessor)
       implicit none
       class(t_buffer), intent(inout) :: this  !< buffer : dummy argument
       class(t_buffer), intent(in) :: sucessor  !< buffer sucesseur a verifier
       
         if (sucessor%m_successor_arret%m_stop.ne.0) then
            this%m_successor_arret = sucessor%m_successor_arret
         endif
         
      end subroutine  buffer_checkget_successor_error  

!***********************************************************************
!> @brief alloue "nbstops" controles d'erreur 
!! pour avoir un controle d'erreur par quantite (e.g.,n 1 controle par particule).\n
!! ces erreurs ne sont pas consideres comme critique
!***********************************************************************
      subroutine buffer_init_multictrlerr(this, nbstops)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       integer, intent(in):: nbstops  !< nombre de controle
      
       allocate(this%m_mce(1:nbstops))
       this%m_count_multictrlerr = 0
      end  subroutine buffer_init_multictrlerr  

!***********************************************************************
!> @brief transfere les erreurs (m_mce) de src vers this
!***********************************************************************
      subroutine buffer_set_multictrlerr_buf(this, src)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       class(t_buffer), intent(in):: src  !< buffer source dont on veut transferer les erreurs
      
       this%m_mce(:) = src%m_mce(:)
       this%m_count_multictrlerr = src%m_count_multictrlerr
      end  subroutine buffer_set_multictrlerr_buf  

!***********************************************************************
!> @brief transfere les erreurs (m_mce) de src vers this
!! si le mce a ete alloue sur la source, alloue si besoin this
!***********************************************************************
      subroutine buffer_tryset_multictrlerr_buf(this, src)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       class(t_buffer), intent(in):: src  !< buffer source dont on veut transferer les erreurs
      
       if (allocated(src%m_mce)) then           
        if (.not.allocated(this%m_mce)) then           
            allocate(this%m_mce(1:size(src%m_mce)),source=src%m_mce)
        endif 
        this%m_mce(:) = src%m_mce(:)
        this%m_count_multictrlerr = src%m_count_multictrlerr
       endif
      end  subroutine buffer_tryset_multictrlerr_buf  


!***********************************************************************
!> @brief transfere les erreurs de src vers this
!***********************************************************************
      subroutine buffer_set_multictrlerr(this, src)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       type(t_arret), dimension(:), intent(in):: src  !< tableau des erreurs a copier

       this%m_mce(:) = src(:)
      end  subroutine buffer_set_multictrlerr  

!***********************************************************************
!> @brief transfere les erreurs de this vers dst 
!***********************************************************************
      subroutine buffer_get_multictrlerr(this, dst)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       type(t_arret), dimension(:), intent(out):: dst  !< tableau des erreurs a mettre a jour
      
       dst(:) = this%m_mce(:)
      end  subroutine buffer_get_multictrlerr  

!***********************************************************************
!> @brief met a jour le compteur d'erreurs de m_mce qui ont un stop!=0
!! retourne true si tous les corps sont en erreur
!***********************************************************************
      function buffer_update_count_multictrlerr(this) result(ball)
       implicit none
       class(t_buffer), intent(inout):: this  !< buffer : dummy argument
       integer k, count
       logical ball
       count = 0
       do k=1, size(this%m_mce)
            if (this%m_mce(k)%m_stop.ne.0) then
                count = count +1
            endif
       enddo
       this%m_count_multictrlerr = count
       ball = .false.
       if (count.eq.size(this%m_mce)) then
            !write(*,*) 'count part error',count, '/', size(this%m_mce)
            ball = .true.
       endif
      end  function buffer_update_count_multictrlerr  

!***********************************************************************
!> @brief  fonction pour generer le graph au format dot
!***********************************************************************
       subroutine buffer_graphdot(this,nodesrc,nodecounter,num_file)
         implicit none
         class(t_buffer), intent(inout) :: this  !< buffer : dummy argument
         integer, intent(in) :: nodesrc !< numero du noeud  source 
         integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
         integer, intent(in) :: num_file !<numero du fichier de sortie
             
         call  this%m_callbackgraphdot(this%m_userdata, nodesrc,         &
     &                                nodecounter, num_file)         
       end subroutine buffer_graphdot

!***********************************************************************
!> @brief  ecrit "node"//numsrc -> "node"//numdst dans le fichier num_file.
!! utilise pour la generation du graphique dot
!***********************************************************************
       subroutine buffer_dot_writenode(num_file, numsrc, numdst )
        implicit none
        integer, intent(in) :: num_file !<numero du fichier de sortie
        integer, intent(in) :: numsrc !<numero du noeud source
        integer, intent(in) :: numdst !<numero du noeud destination
        
        character(len=100) str1, str2
        write (str1,fmt="(I0)") numsrc
        write (str2,fmt="(I0)") numdst
        write (num_file,*) 'node_'//trim(adjustl(str1))//' ->'             &
     &    //' node_'//trim(adjustl(str2))//';'
        end subroutine buffer_dot_writenode

!***********************************************************************
!> @brief  ecrit "node"//numsrc =[label ch ]; dans le fichier num_file.
!! utilise pour la generation du graphique dot
!***********************************************************************
       subroutine buffer_dot_writelabel(num_file, numsrc, ch)
        implicit none
        integer, intent(in) :: num_file !<numero du fichier de sortie
        integer, intent(in) :: numsrc !<numero du noeud source
        character(len=*), intent(in) :: ch !< label a ecrire
        
        character(len=100) str1
        write (str1,fmt="(I0)") numsrc
        write (num_file,*) 'node_'//trim(adjustl(str1))//' [ label="'       &
     &    //ch//'"];'
        end subroutine buffer_dot_writelabel

      end module mod_buffer

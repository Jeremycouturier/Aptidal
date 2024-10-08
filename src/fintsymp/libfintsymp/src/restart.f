!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!  \file restart.f 
!  \brief Implementation of :
!    fonction de restart generique : applicable a tous projets
!      -  pour les integrateurs symplectiques
!
!    ce fichier ne doit rien contenir de specifique a un projet 
!
!    File format binary :
!    - 1 integer(4) : which always contains the value 1234
!     (to detect big or little endianness)
!    - 1 integer(4) : version file format
!    - 1 real(8) : torigin : origin of the integration
!    - 1 treal : stepyear : step of integration
!    - 1 integer(8) : currentstep : current step 
!    - 1 integer(8) : stepnumber : number of the final step 
!
!
! history : creation 2010/07/01
!***********************************************************************


!***********************************************************************
! module generique de redemarrage pour tous les integrateurs
!***********************************************************************
       module restart_common
                
!***********************************************************************
!  type pour la gestion commune des fichiers redemarrage 
!***********************************************************************
       type t_restart
                integer :: fileunit   !  unit du fichier utilise
                
                character(len=256) :: filename ! nom du fichier
                integer(4) :: version ! numero de version du fichier de redemarrage
                real(8)    :: torigin ! temps origine des integrations
                real(8)    :: tstep   ! pas des integrations
                
                ! pour auto-increment
                logical    :: autoincrement ! = true => numero de sauvegarde incrementee a chaque sauvegarde
                character(len=256) :: radical ! nom du radical
                character(len=256) :: suffix ! nom du suffix
                integer    :: counter ! numero du fichier (pour autoincrement)
       end type t_restart

!***********************************************************************
! parametres privees au module
!***********************************************************************
       ! value to check endian file format 
       integer(4), parameter, private :: endian=1234 

       contains
       
!***********************************************************************
! initialise un objet "this" gerant les redemarrages
! avec un nom de fichier fixe
!
!  @param this (out) objet initailiser
!  @param fileunit (in) unit du fichier utilise
!  @param filename (in)  nom du fichier
!  @param version  (in) numero de version du fichier de redemarrage
!  @param torigin (in) temps origine des integrations
!  @param tstep   (in) pas des integrations
!***********************************************************************
             subroutine restart_init(this, fileunit, filename,          & 
     &                   version, torigin,  tstep)
                implicit none
                type (t_restart), intent(out) :: this 
                integer, intent(in) :: fileunit  
                character(len=*), intent(in) :: filename 
                integer(4), intent(in) :: version 
                real(8), intent (in)   :: torigin  
                real(8), intent (in)   :: tstep   
               !***********************************
                                     
               ! ouverture
               this%fileunit = fileunit
               this%filename = filename
               this%version = version
               this%torigin = torigin
               this%tstep   = tstep
               this%autoincrement = .false.
               this%counter = 0
            end subroutine restart_init

!***********************************************************************
! initialise un objet "this" gerant les redemarrages
! avec un nom de radical suivi d'un numero qui s'auto-incremente 
! suvi d'un suffixe
!
!
!  @param this (out) objet initailiser
!  @param fileunit (in) unit du fichier utilise
!  @param radical (in)  radical du fichier
!  @param suffix (in)  suffix du fichier
!  @param version  (in) numero de version du fichier de redemarrage
!  @param torigin (in) temps origine des integrations
!  @param tstep   (in) pas des integrations
!***********************************************************************
             subroutine restart_initauto(this, fileunit,                & 
     &                   radical, suffix, version, torigin,  tstep)
                implicit none
                type (t_restart), intent(out) :: this 
                integer, intent(in) :: fileunit  
                character(len=*), intent(in) :: radical 
                character(len=*), intent(in) :: suffix 
                integer(4), intent(in) :: version 
                real(8), intent (in)   :: torigin  
                real(8), intent (in)   :: tstep   
               !***********************************
                                     
               ! ouverture
               this%fileunit = fileunit
               this%filename = ''
               this%version = version
               this%torigin = torigin
               this%tstep   = tstep
               this%autoincrement = .true.
               this%counter = 0
               this%radical = radical
               this%suffix = suffix
            end subroutine restart_initauto


 !***********************************************************************
! met a jour le compteur pour le numero de sauvegarde des fichiers
!
!  @param this   (in) fichier de redemarrage
!  @param counter (in) nouvelle valeur du compteur
!***********************************************************************
             subroutine restart_setautocounter(this,counter)
                implicit none
                type (t_restart), intent(inout) :: this
                integer, intent(in) :: counter 
               !***********************************
                this%counter = counter

            end subroutine restart_setautocounter
           
!***********************************************************************
! ouvre le fichier de redemarrage filename en utilisant l'unit fileunit
! si le fichier existe, il est ecrase
! ecrit l'entete du fichier de redemarrage
!
!  @param this   (in) donnees  a sauver 
!***********************************************************************
             subroutine restart_save_open(this)
                implicit none
                type (t_restart), intent(inout) :: this 
                character(len=256) :: sval                     
               ! construction du nom si auto-incremente
               if (this%autoincrement) then
                    write(sval,'(A,i5.5,A)') '_',this%counter,'_'
                    this%filename = trim(this%radical)                   & 
     &                              //trim(sval)//trim(this%suffix)
                    this%counter = this%counter+1
               endif
               
               ! ouverture
               open(this%fileunit, file=this%filename,                  & 
     &              status="unknown",  form="UNFORMATTED")
     
               ! ecriture entete
               write(this%fileunit) endian, this%version
               write(this%fileunit) this%counter
               write(this%fileunit) this%torigin, this%tstep
               
            end subroutine restart_save_open

!***********************************************************************
!  ferme le fichier de redemarrage ouvert via restart_save_open
!
!  @param this   (in) donnees  a sauver 
!***********************************************************************
             subroutine restart_save_close(this)
                implicit none
                type (t_restart), intent(in) :: this 
                                     
               ! fermeture
               close(this%fileunit)
            end subroutine restart_save_close

!***********************************************************************
! ouvre le fichier de redemarrage filename en utilisant l'unit fileunit
! lit l'entete du fichier de redemarrage et remplit this
!
!  @param this   (inout) donnees lues
!***********************************************************************
             subroutine restart_load_open(this)
                implicit none
                type (t_restart), intent(inout) :: this 
               !***********************************
                integer(4) endianread
                
               ! ouverture
               open(this%fileunit, file=this%filename,                  & 
     &              status="old",  form="UNFORMATTED")
     
               ! ecriture entete
               read(this%fileunit) endianread, this%version
               if (endianread.ne.endian) then
                write(*,*) 'endian n est pas vide : ',endianread,endian
                stop
               endif
               read(this%fileunit) this%counter
               read(this%fileunit) this%torigin, this%tstep
               
            end subroutine restart_load_open

!***********************************************************************
!  ferme le fichier de redemarrage ouvert via restart_load_open
!
!  @param this   (in) donnees  a sauver 
!***********************************************************************
             subroutine restart_load_close(this)
                implicit none
                type (t_restart), intent(in) :: this 
                                     
               ! fermeture
               close(this%fileunit)
            end subroutine restart_load_close

!***********************************************************************
! tronque le fichier de sortie ouvert avec l'unit fileunit jusqu'au pas currentstep 
! (inclus dans la verification).
! Le temps sur chaque ligne est verifiee a "eps" pres.
! Ce fichier commence au temps tinit et les donnees ont ecrites tous les "step_out*tstep". 
!
! 
! Le fichier doit contenir nbcol colonnes dont la 1ere colonne represente le temps
! Les colonnes sont supposees contenir des reels et separees  par des espaces ou tabulations
!
! Il deplace le curseur  apres la ligne corrspoendant a l'iteration currentstep.
! @param fileunit (in) unit du fichier
! @param filename (in) nom du fichier (pour le message d'erreur)
! @param tinit (in) temps initial des fichiers
! @param tstep (in) pas d'integration
! @param currentstep (in) iteration a laquelle on s'arrete 
!               de relire le fichier (iteration incluse)
! @param step_out (in) pas de sortie a laquelle ont ete ecrites les  donnees
!***********************************************************************
      subroutine restart_truncatefile(fileunit,filename, tinit,tstep,      &
     &                  step_out, currentstep, nbcol, eps)
      implicit none
      integer, intent(in) :: fileunit
      character(len=*), intent(in) :: filename
      real(8), intent(in) :: tinit
      real(8), intent(in) :: tstep
      integer(8), intent(in) :: currentstep
      integer(8), intent(in) :: step_out
      integer, intent(in) :: nbcol
      real(8), intent(in) :: eps
      !************************************
      integer(8) j
      real(8) :: tcheck,t 
      real(8),dimension(nbcol-1) :: C
      
      do j=0, currentstep, step_out
       tcheck = tinit+j*tstep
       read(fileunit, *) t, C
         if (abs(t-tcheck)>eps) then
        	write(*,*) 'dans le fichier ', filename
        	write(*,*) ' : temps lu invalide', t, tcheck
        	stop
         endif
      enddo
      end subroutine restart_truncatefile
      
      end module restart_common


!***********************************************************************
!***********************************************************************
!***********************************************************************
! module generique de redemarrage pour tous les integrateurs SYMPLECTIQUES
!***********************************************************************
!***********************************************************************
!***********************************************************************
       module restart_sympl
        use restart_common
                
!***********************************************************************
!  type pour la gestion des fichiers redemarrage
!***********************************************************************
        type t_restart_sympl
                type (t_restart) :: parent ! donnees communes
                integer(8) ::  totalstep  ! nombre d'iteration totale
                integer(8) ::  step_restart  ! pas de creation du fichier de sauvegarde
                                                     ! = 0 => sauvegarde desactive
         end type t_restart_sympl

       contains


!***********************************************************************
! initialise un objet "this" gerant les redemarrages 
!  pour les integrateurs symplectiques
!  @param this   (out) objet  a initialiser 
!  @param fileunit (in) unit du fichier utilise
!  @param filename (in)  nom du fichier
!  @param version  (in) numero de version du fichier de redemarrage
!  @param torigin (in) temps origine des integrations
!  @param tstep   (in) pas des integrations
!  @param totalstep (in) nombre d'iteration totale
!  @param step_restart (in) pas pour la generation des fichiers de sauvegardes
!               =0 => sauvegarde desactive
!***********************************************************************
             subroutine restart_sympl_init(this, fileunit, filename,    & 
     &                   version, torigin,tstep,totalstep,step_restart)
                implicit none
                type (t_restart_sympl), intent(out) :: this 
                integer, intent(in) :: fileunit  
                character(len=*), intent(in) :: filename 
                integer(4), intent(in) :: version 
                real(8), intent (in)   :: torigin 
                real(8), intent (in)   :: tstep  
                integer(8), intent(in) :: totalstep 
                integer(8), intent(in) :: step_restart
               !***********************************
                                     
               call restart_init(this%parent,  fileunit, filename,      & 
     &                   version, torigin,  tstep)
     
               this%totalstep   = totalstep
               this%step_restart = step_restart
            end subroutine restart_sympl_init

!***********************************************************************
! initialise un objet "this" gerant les redemarrages
! avec un nom de radical suivi d'un numero qui s'auto-incremente 
! suvi d'un suffixe
!
!
!  @param this (out) objet initailiser
!  @param fileunit (in) unit du fichier utilise
!  @param radical (in)  radical du fichier
!  @param suffix (in)  suffix du fichier
!  @param version  (in) numero de version du fichier de redemarrage
!  @param torigin (in) temps origine des integrations
!  @param tstep   (in) pas des integrations
!  @param totalstep (in) nombre d'iteration totale
!  @param step_restart (in) pas pour la generation des fichiers de sauvegardes
!***********************************************************************
             subroutine restart_sympl_initauto(this, fileunit,            & 
     &                   radical, suffix, version, torigin,  tstep,       & 
     &                   totalstep, step_restart)
                implicit none
                type (t_restart_sympl), intent(out) :: this 
                integer, intent(in) :: fileunit  
                character(len=*), intent(in) :: radical 
                character(len=*), intent(in) :: suffix 
                integer(4), intent(in) :: version 
                real(8), intent (in)   :: torigin  
                real(8), intent (in)   :: tstep   
                integer(8), intent(in) :: totalstep 
                integer(8), intent(in) :: step_restart
               !***********************************
                                     
               call restart_initauto(this%parent,  fileunit, radical,   & 
     &                   suffix, version, torigin,  tstep)
     
               this%totalstep   = totalstep
               this%step_restart = step_restart
            end subroutine restart_sympl_initauto

!***********************************************************************
! met a jour le compteur pour le numero de sauvegarde des fichiers
!
!  @param this   (in) donnees  a sauver 
!  @param currentstep (in) iteration courante  
!***********************************************************************
             subroutine restart_sympl_setautocounter(this,currentstep)
                implicit none
                type (t_restart_sympl), intent(inout) :: this
                integer(8), intent(in) :: currentstep 
               !***********************************
               integer newcounter
               
                if (this%step_restart/=0) then 
                 newcounter = int((currentstep-1)/this%step_restart,      &
     &                       kind(newcounter))
                 call restart_setautocounter(this%parent, newcounter)
                endif   

            end subroutine restart_sympl_setautocounter

!***********************************************************************
! verifie si currentstep-1 est modulo de this%step_restart
! ouvre le fichier de redemarrage filename en utilisant l'unit fileunit
! si le fichier existe, il est ecrase
! ecrit l'entete du fichier de redemarrage
!
! retourne true si currentstep est modulo de this%step_restart 
!
!  @param this   (in) donnees  a sauver 
!  @param currentstep (in) iteration courante  
!***********************************************************************
             function restart_sympl_save_open(this, currentstep)
                implicit none
                type (t_restart_sympl), intent(inout) :: this
                integer(8), intent(in) :: currentstep 
                integer :: restart_sympl_save_open
               !***********************************
                                     
               if ((this%step_restart/=0).and.                           & 
     &             (mod(currentstep-1,this%step_restart)==0)) then
                    ! ouverture
                    call restart_save_open(this%parent)
                    !ecriture des donnees specifiques aux integrateurs symplectiques
                    write(this%parent%fileunit) this%step_restart
                    write(this%parent%fileunit) this%totalstep
                    write(this%parent%fileunit) currentstep
                    restart_sympl_save_open = this%parent%fileunit
               else
                    restart_sympl_save_open = 0
               endif
            end function restart_sympl_save_open
        
      

!***********************************************************************
!  ferme le fichier de redemarrage ouvert via restart_sympl_save_open
!
!  @param this   (in) donnees  a sauver 
!***********************************************************************
             subroutine restart_sympl_save_close(this)
                implicit none
                type (t_restart_sympl), intent(in) :: this
               !***********************************
                                     
               ! fermeture
               call restart_save_close(this%parent)
            end subroutine restart_sympl_save_close
       
!***********************************************************************
! ouvre le fichier de redemarrage this%filename en utilisant l'unit fileunit
! initialise this sauf fileunit et filename en lisant contenu
!
!  @param this   (inout) donnees lues
!  @param currentstep (out) iteration courante  
!***********************************************************************
             subroutine restart_sympl_load_open(this, currentstep)
                implicit none
                type (t_restart_sympl), intent(inout) :: this
                integer(8), intent(out) :: currentstep 
               !***********************************
                                    
               ! ouverture
               call restart_load_open(this%parent)
                !lecture des donnees specifiques aux integrateurs symplectiques
               read(this%parent%fileunit) this%step_restart
               read(this%parent%fileunit) this%totalstep
               read(this%parent%fileunit) currentstep
               
            end subroutine restart_sympl_load_open

!***********************************************************************
!  ferme le fichier de redemarrage ouvert via restart_sympl_load_open
!
!  @param this   (in) donnees  a sauver 
!***********************************************************************
             subroutine restart_sympl_load_close(this)
                implicit none
                type (t_restart_sympl), intent(in) :: this
               !***********************************
                                     
               ! fermeture
               call restart_load_close(this%parent)
            end subroutine restart_sympl_load_close


       end module restart_sympl



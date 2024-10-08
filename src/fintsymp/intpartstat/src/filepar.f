!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file filepar.f 
!!  \brief gestion du fichier de parametres 
!! les donnees dans le fichier de conditions initiales nf_initext doivent etre en UA et an
!!
!!
! history : creation 15/09/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour gerer le fichier de parametres
!***********************************************************************
      module mod_filepar
        use mod_io_txt_mligncomp
        use mod_io_txt_carpart
        use mod_io_txt_ellpart
       
!***********************************************************************
!> @class t_minmax_par
!! classe de gestion des parametres pour calculer les min/moy/max : 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_minmax_par
      
         !> * = 0 => pas de calcul des min/moy/max
         !! * = 1 => calcul des min/moy/max 
         integer   m_compute 
         integer(8) ::   m_stepcalc !< pas de calcul des min/moy/max 
         integer(8) ::   m_stepout  !< longueur des tranches des min/moy/max 
         integer   m_elltype  !< type des elements elliptiques (memes valeurs que les conditions initiales) 

      end type t_minmax_par

!***********************************************************************
!> @class t_naf_par
!! classe de gestion des parametres pour calculer les analyses en frequence : 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_naf_par

         !> * = 0 => pas de calcul d'analyse en frequence
         !! * = 1 => calcul d'analyse en frequence 
         integer   m_compute 
         integer(8) ::   m_stepcalc !< pas utilise pour l'analyse
         integer(8) ::   m_stepout  !< longueur des tranches pour l'analyse 
         integer   m_elltype  !< type des elements elliptiques  (memes valeurs que les conditions initiales) 
         ! parametre de controle de l'analyse en frequence
         integer   m_naf_nterm !< nombre de termes demandes 
         !> presence de fenetre
         !! *= -1 : fenetre exponentielle PHI(T) = 1/CE*EXP(-1/(1-T^2)) avec CE= 0.22199690808403971891E0
         !! *= 0 : pas de fenetre
         !! *= N > 0 : PHI(T) = CN*(1+COS(PI*T))**N avec CN = 2^N(N!)^2/(2N)!
         integer   m_naf_iw !< nombre de termes demandes 
         !> *= 0 : la methode des secantes n'est pas utilisee.
         !! *= 1 : la methode des secantes est utilisee.
         integer   m_naf_isec 
         real(TREAL) ::   m_naf_dtour !< "longueur d'un tour de cadran"
         real(TREAL) ::   m_naf_tol !< tolerance pour determiner si deux frequences sont identiques
      end type t_naf_par
 

!***********************************************************************
!> @class t_ctrl_diststar_par
!! classe de gestion des parametres pour le controle de distance d'une particule a l'etoile : 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_ctrl_diststar_par
      
         !> * = 0 => pas de controle des distances
         !! * = 1 => controle des distances
         integer          m_compute 
         integer(8)  ::   m_stepcalc !< pas de calcul du controle 
         real(TREAL) ::   m_distmin !< distance minimale a l'etoile. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
         real(TREAL) ::   m_distmax !< distance maximale a l'etoile. si m_distmax =-1, pas de controle, si m_distmax>=0, genere une erreur si distance>m_distmax
      end type t_ctrl_diststar_par

!***********************************************************************
!> @class t_ctrl_distpla_par
!! classe de gestion des parametres pour le controle de distance d'une particule a une planete
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_ctrl_distpla_par
      
         !> * = 0 => pas de controle des distances
         !! * = 1 => controle des distances
         integer          m_compute 
         integer(8)  ::   m_stepcalc !< pas de calcul du controle 
         real(TREAL) ::   m_distmin !< distance minimale d'une particule a une planete. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
      end type t_ctrl_distpla_par

!***********************************************************************
!> @class t_orb_pla_tabulee_par
!! classe de gestion des parametres pour une solution tabulee planetaire fournie par une fonction tabulee
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_orb_pla_tabulee_par
      
         integer          m_coord !< type de coordonnees dans le fichier (memes valeurs que les conditions initiales si m_type=1)  
         character(100)  :: m_nf  !< nom du fichier ou se trouve la solution planetaire
       contains
       procedure :: fread => orb_pla_tabulee_par_fread ! lecture du fichier de parametres de la solution planetaire forcee
       procedure :: fwrite => orb_pla_tabulee_par_fwrite ! ecriture du fichier de parametres au format namelist avec son numero 
       procedure :: fwritetrip => orb_pla_tabulee_par_fwritetrip ! ecriture du fichier de parametres au format trip avec son numero 

      end type t_orb_pla_tabulee_par

!***********************************************************************
!> @class t_filepar
!! classe de gestion du fichier de parametres : 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_filepar

        character(100)  :: chemin,nf_rad,int_type
        character(100) :: nf_initext !< nom du fichier des conditions initiales des planetes
        character(100) :: nf_initpart !< nom du fichier des conditions initiales des particules
        integer :: if_int
        
        integer(8) :: n_iter,n_out
        integer    ::  if_invar
        real(TREAL) :: tinit,dt,dmax
        integer :: part_blocksize !< nombre de particules traitees en simultane
        integer :: ref_gmsun
        
        integer :: if_ell_part !< sortie en elements elliptiques pour les particules
        integer :: if_car_part !< sortie en elements cartesiens pour les particules
        integer :: out_ell_part !< type de sortie des elements elliptiques pour les particules

        ! pour les statistiques des min/moy/max
        type(t_minmax_par) :: minmax_aei  !< min/moy/max en a, e, i
        type(t_minmax_par) :: minmax_diffalp !< min/moy/max en l1-l2,pibar1-pibar2, a1-a2
        integer,dimension(1:1) :: minmax_diffalp_pla !< numero des 1 planete utilisee pour minmax_diffalp


        ! pour les analyses en frequences pour les particules
        type(t_naf_par) :: naf_alkhqp  !< analyse en frequence en a*exp(*il),k+i*h, q+i*p
        type(t_naf_par) :: naf_diffalp  !< analyse en frequence en exp(i*(l1-l2)),exp(i*(pi1-pi2))
        integer,dimension(1:1) :: naf_diffalp_pla !< numero des 1 planete utilisee pour naf_diffalp
        
       ! controle de la distance a l'etoile
        type(t_ctrl_diststar_par) :: ctrl_diststar !< controle de la distance a l'etoile
        
        ! propre aux planetes
        integer :: if_ell_pla !< sortie en elements elliptiques pour les planetes
        integer :: if_car_pla !< sortie en elements cartesiens pour les planetes
        integer :: out_ell_pla !< type de sortie des elements elliptiques pour les planetes

        ! controle de la distance aux planetes
        type(t_ctrl_distpla_par) :: ctrl_distpla !< controle de la distance aux planetes

        ! pour les analyses en frequences sur les planetes
         type(t_naf_par) :: naf_alkhqp_pla  !< analyse en frequence en a*exp(*il),k+i*h, q+i*p
         
        ! pour la solution planetaire forcee
        integer :: if_orb_pla !< solution planetaire integree ou deja precalculee
        type(t_orb_pla_tabulee_par) ::  orb_pla

       contains
       
       procedure :: fread => filepar_fread ! lecture du fichier de parametres
       procedure :: fwrite => filepar_fwrite ! ecriture du fichier de parametres au format namelist avec son numero 
       procedure :: fwritepar => filepar_fwritepar ! ecriture du fichier de parametres au format namelist dont le nom est genere
       procedure :: fwritetrip => filepar_fwritetrip ! ecriture du fichier de parametres au format trip avec son numero 
       procedure :: fwritepartrip => filepar_fwritepartrip ! ecriture du fichier de parametres au format trip dont le nom est genere
                              
      end type t_filepar  

!***********************************************************************
!> @class t_extfilepar
!! classe de gestion du fichier parametres : 
!! contient les donnees du fichier de parametres et les variables derivees
!!  
!***********************************************************************
      type, extends(t_filepar)  :: t_extfilepar

        ! propre aux particules
        character(len=512) :: nf_ci_pla !< nom du fichier des conditions initiales des planetes
        character(len=512) :: nf_ci_part !< nom du fichier des conditions initiales des particules
        character(len=512) :: nf_minmax_aei !< nom du fichier de sortie des min/max a,e,i
        character(len=512) :: nf_minmax_diffalp !< nom du fichier de sortie des min/max a_p1-a_p2,...
        character(len=512) :: nf_naf_alkhqp !< nom du fichier de sortie des analyses en frequence en alkhqp
        character(len=512) :: nf_naf_diffalp !< nom du fichier de sortie des analyses en frequence en exp(i*(l1-l2)), exp(i*(pi1-pi2)) 
        character(len=512) :: nf_car !< nom du fichier de sortie des elements cartesiens des particules
        character(len=512) :: nf_ell !< nom du fichier de sortie des elements elliptiques des particules
        character(len=512) :: nf_control !< nom du fichier de controle de l'integration des particules
        
        integer :: int_var   !< variables de l'integrateur (: 0 : Jacobi  1 :helio canonique)

        type(t_io_txt_car_part) :: io_cart !< sortie des elements cartesiens des particules
        type(t_io_txt_ell_part) :: io_ell !< sortie des elements elliptiques des particules 

        type (t_io_txt_mligncomp) :: io_minmax_aei  !< sortie des min/max en a, e, i
        type (t_io_txt_mligncomp) :: io_minmax_diffalp  !< sortie des min/max des differences a_p1-a_p2,....
        integer :: io_control     !< numero du fichier de controle des integrations des particules

        type (t_io_txt_mligncomp) :: io_naf_alkhqp  !< sortie des analyses en frenquence en a*exp(*il),k+i*h, q+i*p
        type (t_io_txt_mligncomp) :: io_naf_diffalp     !< sortie des analyse en frenquence en exp(i*(l1-l2)),exp(i*(pi1-pi2))

         ! propre aux planetes
        character(len=512) :: nf_int !< nom du fichier de sortie des integrales premieres
        character(len=512) :: nf_car_pla!< nom du fichier de sortie des elements cartesiens des planetes
        character(len=512) :: nf_ell_pla !< nom du fichier de sortie des elements elliptiques des planetes
        character(len=512) :: nf_control_pla !< nom du fichier de controle de l'integration du systeme planetaire
        character(len=512) :: nf_naf_alkhqp_pla !< nom du fichier de sortie des analyses en frequence en alkhqp

        type (t_io_txt_ncomp) :: io_naf_alkhqp_pla  !< sortie des analyses en frenquence en a*exp(*il),k+i*h, q+i*p des planetes
        integer :: io_control_pla     !< numero du fichier de controle des integrations des planetes

       contains
                      
        procedure :: create_filename => extfilepar_create_filename  ! creation des noms de fichiers de sortie  communs a toutes les conditions initiales                     
        procedure :: open_file => extfilepar_open_file ! ouverture en ecriture des fichiers
      end type t_extfilepar  

      contains
      
!***********************************************************************
!> @brief lit le fichier de parametres name et remplit this
!***********************************************************************
      subroutine filepar_fread(this, fname)
       implicit none
       class(t_filepar), intent(inout):: this    !< fichier de parametres
       character(len=*), intent(in) :: fname     !< nom du fichier a lire
       
        character(100)  :: chemin,nf_rad,int_type
        character(100) :: nf_initext, nf_initpart
        integer :: if_int,if_ell_pla,if_car_pla
        integer(8) :: n_iter,n_out
        integer :: out_ell_pla,if_invar, ref_gmsun
        integer :: if_ell_part,if_car_part, out_ell_part
        real(TREAL) :: tinit,dt
        integer :: minmax_aei_compute
        integer(8) :: minmax_aei_stepcalc, minmax_aei_stepout
        integer :: minmax_aei_elltype
        integer :: minmax_diffalp_compute
        integer(8) :: minmax_diffalp_stepcalc, minmax_diffalp_stepout
        integer :: minmax_diffalp_elltype
        integer, dimension(1:1) :: minmax_diffalp_pla
        integer(8) :: rlen
        integer :: naf_alkhqp_compute
        integer(8) :: naf_alkhqp_stepcalc, naf_alkhqp_stepout
        integer :: naf_alkhqp_elltype
        integer :: naf_diffalp_compute
        integer(8) :: naf_diffalp_stepcalc, naf_diffalp_stepout
        integer :: naf_diffalp_elltype
        integer, dimension(1:1) :: naf_diffalp_pla
        integer naf_alkhqp_nterm
        integer naf_alkhqp_isec, naf_alkhqp_iw 
        real(TREAL) naf_alkhqp_dtour, naf_alkhqp_tol 
        integer naf_diffalp_nterm
        integer naf_diffalp_isec, naf_diffalp_iw 
        real(TREAL) naf_diffalp_dtour, naf_diffalp_tol 
        
        integer :: naf_alkhqp_pla_compute
        integer(8) :: naf_alkhqp_pla_stepcalc, naf_alkhqp_pla_stepout
        integer :: naf_alkhqp_pla_elltype
        integer naf_alkhqp_pla_nterm
        integer naf_alkhqp_pla_isec, naf_alkhqp_pla_iw 
        real(TREAL) naf_alkhqp_pla_dtour, naf_alkhqp_pla_tol 

        integer :: ctrl_diststar_compute
        integer(8) :: ctrl_diststar_stepcalc
        real(TREAL) :: ctrl_diststar_distmin
        real(TREAL) :: ctrl_diststar_distmax

        integer :: ctrl_distpla_compute
        integer(8) :: ctrl_distpla_stepcalc
        real(TREAL) :: ctrl_distpla_distmin
        
        integer part_blocksize
        integer if_orb_pla


      namelist/lect/chemin,nf_rad,nf_initext, nf_initpart,
     &   int_type,tinit,dt,n_iter,n_out,
     &   out_ell_pla,int_type,ref_gmsun, if_orb_pla,
     &   if_invar,if_int,if_ell_pla,if_car_pla,
     &   if_ell_part,if_car_part, out_ell_part,
     &   minmax_aei_compute,
     &   minmax_aei_stepcalc, minmax_aei_stepout,
     &   minmax_aei_elltype,
     &   minmax_diffalp_compute,
     &   minmax_diffalp_stepcalc, minmax_diffalp_stepout,
     &   minmax_diffalp_elltype,
     &   minmax_diffalp_pla,
     &   naf_alkhqp_compute,
     &   naf_alkhqp_stepcalc, naf_alkhqp_stepout,
     &   naf_alkhqp_elltype,
     &   naf_alkhqp_nterm,
     &   naf_alkhqp_isec, naf_alkhqp_iw, 
     &   naf_alkhqp_dtour, naf_alkhqp_tol, 
     &   naf_diffalp_compute,
     &   naf_diffalp_stepcalc, naf_diffalp_stepout,
     &   naf_diffalp_elltype,  naf_diffalp_pla,
     &   naf_diffalp_nterm,
     &   naf_diffalp_isec, naf_diffalp_iw, 
     &   naf_diffalp_dtour, naf_diffalp_tol, 
     &   ctrl_diststar_compute, ctrl_diststar_stepcalc, 
     &   ctrl_diststar_distmin, ctrl_diststar_distmax,
     &   ctrl_distpla_compute, ctrl_distpla_stepcalc, 
     &   ctrl_distpla_distmin, part_blocksize,
     &   naf_alkhqp_pla_compute,
     &   naf_alkhqp_pla_stepcalc, naf_alkhqp_pla_stepout,
     &   naf_alkhqp_pla_elltype,
     &   naf_alkhqp_pla_nterm,
     &   naf_alkhqp_pla_isec, naf_alkhqp_pla_iw, 
     &   naf_alkhqp_pla_dtour, naf_alkhqp_pla_tol 
      
      ref_gmsun = 1
      open(10,file=fname,status='old')
      read(10,lect)
      close(10)
      
      this%chemin = chemin
      this%nf_rad = nf_rad
      this%nf_initext = nf_initext
      this%nf_initpart = nf_initpart
      this%int_type = int_type
      this%tinit = tinit
      this%dt = dt
      this%n_iter = n_iter
      this%n_out = n_out
      this%int_type = int_type
      this%if_invar = if_invar
      this%if_int = if_int
      this%if_ell_pla = if_ell_pla
      this%out_ell_pla = out_ell_pla
      this%if_car_pla = if_car_pla
      this%if_ell_part = if_ell_part
      this%out_ell_part = out_ell_part
      this%if_car_part = if_car_part
      this%part_blocksize = part_blocksize
      this%ref_gmsun = ref_gmsun
      this%if_orb_pla = if_orb_pla
      
      ! controle de partblocksize
      if (part_blocksize.le.0) then
        stop "part_blocksize doit etre >0"
      endif

      ! initialisation des min/max a,e,I
      this%minmax_aei%m_compute  = minmax_aei_compute
      this%minmax_aei%m_stepcalc = minmax_aei_stepcalc
      this%minmax_aei%m_stepout  = minmax_aei_stepout
      this%minmax_aei%m_elltype  = minmax_aei_elltype
      ! verifications 
      if (minmax_aei_compute.eq.1) then 
       !verification  de minmax_aei_elltype
       if (minmax_aei_elltype>2) then
        stop "minmax_aei_elltype doit etre 1 ou 2"
       endif
       !verification  de minmax_aei_stepout
       if (mod(minmax_aei_stepout,minmax_aei_stepcalc).ne.0) then
        stop 'minmax_aei_stepout doit etre multiple de'//                  &
     &       ' minmax_aei_stepcalc'
        
       endif
      endif

      ! initialisation des min/max des diff alp
      this%minmax_diffalp%m_compute  = minmax_diffalp_compute
      this%minmax_diffalp%m_stepcalc = minmax_diffalp_stepcalc
      this%minmax_diffalp%m_stepout  = minmax_diffalp_stepout
      this%minmax_diffalp%m_elltype  = minmax_diffalp_elltype
      this%minmax_diffalp_pla = minmax_diffalp_pla
      ! verification de elltype
      if (minmax_diffalp_compute.eq.1) then 
      ! verification de minmax_diffalp_elltype
       if (minmax_diffalp_elltype.ne.6) then
        write(*,*) "minmax_diffalp_elltype doit etre 6"
        stop
       endif
       !verification  de minmax_diffalp_stepout
       rlen = mod(minmax_diffalp_stepout,minmax_diffalp_stepcalc)
       if (rlen.ne.0) then
        stop 'minmax_diffalp_stepout doit etre multiple de'//              &
     &       ' minmax_diffalp_stepcalc'
        
       endif
      endif

      ! initialisation des naf alkhqp
      this%naf_alkhqp%m_compute  = naf_alkhqp_compute
      this%naf_alkhqp%m_stepcalc = naf_alkhqp_stepcalc
      this%naf_alkhqp%m_stepout  = naf_alkhqp_stepout
      this%naf_alkhqp%m_elltype  = naf_alkhqp_elltype
      this%naf_alkhqp%m_naf_nterm  = naf_alkhqp_nterm 
      this%naf_alkhqp%m_naf_iw = naf_alkhqp_iw 
      this%naf_alkhqp%m_naf_isec = naf_alkhqp_isec 
      this%naf_alkhqp%m_naf_dtour = naf_alkhqp_dtour 
      this%naf_alkhqp%m_naf_tol  = naf_alkhqp_tol 
      ! verifications 
      if (naf_alkhqp_compute.ne.0) then 
       !verification  de naf_alkhqp_elltype
       if ((naf_alkhqp_elltype.ne.3).and.(naf_alkhqp_elltype.ne.4))then
        stop "naf_alkhqp_elltype doit etre 3 ou 4"
       endif
       !verification  de naf_alkhqp_stepout
       if (mod(naf_alkhqp_stepout,naf_alkhqp_stepcalc).ne.0) then
        stop 'naf_alkhqp_stepout doit etre multiple de'//                  &
     &       ' naf_alkhqp_stepcalc'
        
       endif
      endif
      
       ! initialisation des naf des diff alp
      this%naf_diffalp%m_compute  = naf_diffalp_compute
      this%naf_diffalp%m_stepcalc = naf_diffalp_stepcalc
      this%naf_diffalp%m_stepout  = naf_diffalp_stepout
      this%naf_diffalp%m_elltype  = naf_diffalp_elltype
      this%naf_diffalp_pla = naf_diffalp_pla
      this%naf_diffalp%m_naf_nterm  = naf_diffalp_nterm 
      this%naf_diffalp%m_naf_iw = naf_diffalp_iw 
      this%naf_diffalp%m_naf_isec = naf_diffalp_isec 
      this%naf_diffalp%m_naf_dtour = naf_diffalp_dtour 
      this%naf_diffalp%m_naf_tol  = naf_diffalp_tol 
      ! verifications 
      if (naf_diffalp_compute.eq.1) then 
       !verification  de naf_diffalp_elltype
       if ((naf_diffalp_elltype.ne.6)) then
        stop "naf_diffalp_elltype doit etre 6"
       endif
       !verification  de naf_diffalp_stepout
       if (mod(naf_diffalp_stepout,naf_diffalp_stepcalc).ne.0) then
        stop 'naf_diffalp_stepout doit etre multiple de'//                  &
     &       ' naf_diffalp_stepcalc'
        
       endif
      endif
      
      ! controle des distances a l'etoile
      this%ctrl_diststar%m_compute = ctrl_diststar_compute
      this%ctrl_diststar%m_stepcalc = ctrl_diststar_stepcalc
      this%ctrl_diststar%m_distmin = ctrl_diststar_distmin
      this%ctrl_diststar%m_distmax = ctrl_diststar_distmax
      
      
      !
      ! pour les planetes
      !
      ! initialisation des naf alkhqp
      this%naf_alkhqp_pla%m_compute  = naf_alkhqp_pla_compute
      this%naf_alkhqp_pla%m_stepcalc = naf_alkhqp_pla_stepcalc
      this%naf_alkhqp_pla%m_stepout  = naf_alkhqp_pla_stepout
      this%naf_alkhqp_pla%m_elltype  = naf_alkhqp_pla_elltype
      this%naf_alkhqp_pla%m_naf_nterm  = naf_alkhqp_pla_nterm 
      this%naf_alkhqp_pla%m_naf_iw = naf_alkhqp_pla_iw 
      this%naf_alkhqp_pla%m_naf_isec = naf_alkhqp_pla_isec 
      this%naf_alkhqp_pla%m_naf_dtour = naf_alkhqp_pla_dtour 
      this%naf_alkhqp_pla%m_naf_tol  = naf_alkhqp_pla_tol 
      ! verifications 
      if (naf_alkhqp_pla_compute.ne.0) then 
       !verification  de naf_alkhqp_elltype
       if ((naf_alkhqp_pla_elltype.ne.3).and.                             &
     &  (naf_alkhqp_pla_elltype.ne.4))then
        stop "naf_alkhqp_pla_elltype doit etre 3 ou 4"
       endif
       !verification  de naf_alkhqp_stepout
       if (mod(naf_alkhqp_pla_stepout,naf_alkhqp_pla_stepcalc).ne.0)       &
     &   then
        stop 'naf_alkhqp_pla_stepout doit etre multiple de'//                  &
     &       ' naf_alkhqp_pla_stepcalc'
        
       endif
      endif

      ! controle des distances aux planetes
      this%ctrl_distpla%m_compute = ctrl_distpla_compute
      this%ctrl_distpla%m_stepcalc = ctrl_distpla_stepcalc
      this%ctrl_distpla%m_distmin = ctrl_distpla_distmin

      ! lecture de la solution orbitale forcee des planetes
      if (this%if_orb_pla.ne.0) then
        call this%orb_pla%fread(fname)
      endif
      end  subroutine filepar_fread  
          

!***********************************************************************
!> @brief ecrit  des parametres de this dans le fichier  nupar
!***********************************************************************
      subroutine filepar_fwrite(this, nupar)
       implicit none
       class(t_filepar), intent(inout):: this    !< fichier de parametres
       integer, intent(in) :: nupar     !< numero de fichier

      write(nupar,*)' chemin     = ',                                    &
     & this%chemin(1:len_trim(this%chemin))
      write(nupar,*)' nf_rad     = ',                                    &
     & this%nf_rad(1:len_trim(this%nf_rad))
      write(nupar,*)' nf_initext = ',                                    &
     & this%nf_initext(1:len_trim(this%nf_initext))
      write(nupar,*)' nf_initpart = ',                                   &
     & this%nf_initpart(1:len_trim(this%nf_initpart))
      write(nupar,*)' int_type   = ',trim(this%int_type)
      write(nupar,*)' if_orb_pla = ',this%if_orb_pla
      write(nupar,*)' ref_gmsun  = ',this%ref_gmsun
      write(nupar,*)' tinit      = ',this%tinit
      write(nupar,*)' dt         = ',this%dt
      write(nupar,*)' n_iter     = ',this%n_iter
      write(nupar,*)' n_out      = ',this%n_out
      write(nupar,*)' if_invar   = ',this%if_invar
      write(nupar,*)' if_int     = ',this%if_int
      write(nupar,*)' out_ell_pla= ',this%out_ell_pla
      write(nupar,*)' if_ell_pla = ',this%if_ell_pla
      write(nupar,*)' if_car_pla = ',this%if_car_pla
      write(nupar,*)' part_blocksize = ',this%part_blocksize
      write(nupar,*)' out_ell_part= ',this%out_ell_part
      write(nupar,*)' if_ell_part = ',this%if_ell_part
      write(nupar,*)' if_car_part = ',this%if_car_part
      

      write(nupar,*)' minmax_aei_compute     = ',                          &
     & this%minmax_aei%m_compute
      write(nupar,*)' minmax_aei_stepcalc     = ',                         &
     & this%minmax_aei%m_stepcalc
      write(nupar,*)' minmax_aei_stepout     = ',                          &
     & this%minmax_aei%m_stepout
      write(nupar,*)' minmax_aei_elltype     = ',                          &
     & this%minmax_aei%m_elltype

      write(nupar,*)' minmax_diffalp_compute     = ',                      &
     & this%minmax_diffalp%m_compute
      write(nupar,*)' minmax_diffalp_stepcalc     = ',                     &
     & this%minmax_diffalp%m_stepcalc
      write(nupar,*)' minmax_diffalp_stepout     = ',                      &
     & this%minmax_diffalp%m_stepout
      write(nupar,*)' minmax_diffalp_elltype     = ',                      &
     & this%minmax_diffalp%m_elltype
      write(nupar,*)' minmax_diffalp_pla(1)     = ',                        &
     & this%minmax_diffalp_pla(1)

      write(nupar,*)' naf_alkhqp_compute     = ',                          &
     & this%naf_alkhqp%m_compute
      write(nupar,*)' naf_alkhqp_stepcalc     = ',                         &
     & this%naf_alkhqp%m_stepcalc
      write(nupar,*)' naf_alkhqp_stepout     = ',                          &
     & this%naf_alkhqp%m_stepout
      write(nupar,*)' naf_alkhqp_elltype     = ',                          &
     & this%naf_alkhqp%m_elltype
      write(nupar,*)' naf_alkhqp_nterm     = ',                            &
     & this%naf_alkhqp%m_naf_nterm
      write(nupar,*)' naf_alkhqp_isec     = ',                             &
     & this%naf_alkhqp%m_naf_isec
      write(nupar,*)' naf_alkhqp_iw     = ',                               &
     & this%naf_alkhqp%m_naf_iw
      write(nupar,fmt='(A,D23.16)')' naf_alkhqp_dtour     = ',             &
     & this%naf_alkhqp%m_naf_dtour
      write(nupar,fmt='(A,D23.16)')' naf_alkhqp_tol     = ',               &
     & this%naf_alkhqp%m_naf_tol

      write(nupar,*)' naf_diffalp_compute     = ',                          &
     & this%naf_diffalp%m_compute
      write(nupar,*)' naf_diffalp_stepcalc     = ',                         &
     & this%naf_diffalp%m_stepcalc
      write(nupar,*)' naf_diffalp_stepout     = ',                          &
     & this%naf_diffalp%m_stepout
      write(nupar,*)' naf_diffalp_elltype     = ',                          &
     & this%naf_diffalp%m_elltype
      write(nupar,*)' naf_diffalp_nterm     = ',                            &
     & this%naf_diffalp%m_naf_nterm
      write(nupar,*)' naf_diffalp_isec     = ',                             &
     & this%naf_diffalp%m_naf_isec
      write(nupar,*)' naf_diffalp_iw     = ',                               &
     & this%naf_diffalp%m_naf_iw
      write(nupar,fmt='(A,D23.16)')' naf_diffalp_dtour     = ',             &
     & this%naf_diffalp%m_naf_dtour
      write(nupar,fmt='(A,D23.16)')' naf_diffalp_tol     = ',               &
     & this%naf_diffalp%m_naf_tol
      write(nupar,*)' naf_diffalp_pla(1)     = ',                           &
     & this%naf_diffalp_pla(1)


      write(nupar,*)' ctrl_diststar_compute     = ',                        &
     & this%ctrl_diststar%m_compute
      write(nupar,*)' ctrl_diststar_stepcalc    = ',                        &
     & this%ctrl_diststar%m_stepcalc
      write(nupar,*)' ctrl_diststar_distmin     = ',                        &
     & this%ctrl_diststar%m_distmin
      write(nupar,*)' ctrl_diststar_distmax     = ',                        &
     & this%ctrl_diststar%m_distmax
      
      write(nupar,*)' ctrl_distpla_compute      = ',                        &
     & this%ctrl_distpla%m_compute
      write(nupar,*)' ctrl_distpla_stepcalc     = ',                        &
     & this%ctrl_distpla%m_stepcalc
      write(nupar,*)' ctrl_distpla_distmin      = ',                        &
     & this%ctrl_distpla%m_distmin

      write(nupar,*)' naf_alkhqp_pla_compute     = ',                          &
     & this%naf_alkhqp_pla%m_compute
      write(nupar,*)' naf_alkhqp_pla_stepcalc     = ',                         &
     & this%naf_alkhqp_pla%m_stepcalc
      write(nupar,*)' naf_alkhqp_pla_stepout     = ',                          &
     & this%naf_alkhqp_pla%m_stepout
      write(nupar,*)' naf_alkhqp_pla_elltype     = ',                          &
     & this%naf_alkhqp_pla%m_elltype
      write(nupar,*)' naf_alkhqp_pla_nterm     = ',                            &
     & this%naf_alkhqp_pla%m_naf_nterm
      write(nupar,*)' naf_alkhqp_pla_isec     = ',                             &
     & this%naf_alkhqp_pla%m_naf_isec
      write(nupar,*)' naf_alkhqp_pla_iw     = ',                               &
     & this%naf_alkhqp_pla%m_naf_iw
      write(nupar,fmt='(A,D23.16)')' naf_alkhqp_pla_dtour     = ',             &
     & this%naf_alkhqp_pla%m_naf_dtour
      write(nupar,fmt='(A,D23.16)')' naf_alkhqp_pla_tol     = ',               &
     & this%naf_alkhqp_pla%m_naf_tol

      ! ecriture de la solution orbitale forcee des planetes
      if( this%if_orb_pla.ne.0) then
        call this%orb_pla%fwrite(nupar)
      endif

      end  subroutine filepar_fwrite  

!***********************************************************************
!> @brief ecrit  des parametres de this dans le fichier 
!! dont le nom est genere a partir de chemin  et nf_rad
!***********************************************************************
      subroutine filepar_fwritepar(this)
       implicit none
       class(t_filepar), intent(inout):: this    !< fichier de parametres

      character(100) nf_par
      integer nu_par

      nf_par = trim(this%chemin)//trim(this%nf_rad)//'.par'
      nu_par = 11
      open(nu_par,file=nf_par,status='unknown')
      call filepar_fwrite(this, nu_par)
      close(nu_par)
      
      end  subroutine filepar_fwritepar        

!***********************************************************************
!> @brief ecrit  des parametres de this dans le fichier  nupar.
!! Le fichier est lible par trip
!***********************************************************************
      subroutine filepar_fwritetrip(this, nupar)
       implicit none
       class(t_filepar), intent(inout):: this    !< fichier de parametres
       integer, intent(in) :: nupar     !< numero de fichier

      write(nupar,*)' chemin     = "',                                    &
     & this%chemin(1:len_trim(this%chemin)), '"$'
      write(nupar,*)' nf_rad     = "',                                    &
     & this%nf_rad(1:len_trim(this%nf_rad)), '"$'
      write(nupar,*)' nf_initext = "',                                    &
     & this%nf_initext(1:len_trim(this%nf_initext)), '"$'
      write(nupar,*)' nf_initpart = "',                                   &
     & this%nf_initpart(1:len_trim(this%nf_initpart)), '"$'
      write(nupar,*)' int_type   = ',trim(this%int_type), '$'
      write(nupar,*)' if_orb_pla = ',this%if_orb_pla, '$'
      write(nupar,*)' tinit      = ',this%tinit, '$'
      write(nupar,*)' dt         = ',this%dt, '$'
      write(nupar,*)' n_iter     = ',this%n_iter, '$'
      write(nupar,*)' n_out      = ',this%n_out, '$'
      write(nupar,*)' if_invar   = ',this%if_invar, '$'
      write(nupar,*)' if_int     = ',this%if_int, '$'
      write(nupar,*)' if_ell_pla = ',this%if_ell_pla, '$'
      write(nupar,*)' if_car_pla = ',this%if_car_pla, '$'
      write(nupar,*)' out_ell_pla= ',this%out_ell_pla, '$'
      write(nupar,*)' part_blocksize = ',this%part_blocksize, '$'
      write(nupar,*)' if_ell_part= ',this%if_ell_part, '$'
      write(nupar,*)' if_car_part= ',this%if_car_part, '$'
      write(nupar,*)' out_ell_part= ',this%out_ell_part, '$'

      write(nupar,*)' minmax_aei_compute     = ',                          &
     & this%minmax_aei%m_compute, '$'
      write(nupar,*)' minmax_aei_stepcalc     = ',                         &
     & this%minmax_aei%m_stepcalc, '$'
      write(nupar,*)' minmax_aei_stepout     = ',                          &
     & this%minmax_aei%m_stepout, '$'
      write(nupar,*)' minmax_aei_elltype     = ',                          &
     & this%minmax_aei%m_elltype, '$'

      write(nupar,*)' minmax_diffalp_compute     = ',                      &
     & this%minmax_diffalp%m_compute, '$'
      write(nupar,*)' minmax_diffalp_stepcalc     = ',                     &
     & this%minmax_diffalp%m_stepcalc, '$'
      write(nupar,*)' minmax_diffalp_stepout     = ',                      &
     & this%minmax_diffalp%m_stepout, '$'
      write(nupar,*)' minmax_diffalp_elltype     = ',                      &
     & this%minmax_diffalp%m_elltype, '$'
      write(nupar,*)' dim minmax_diffalp_pla[1:1]$'
      write(nupar,*)' minmax_diffalp_pla[1]     = ',                       &
     & this%minmax_diffalp_pla(1), '$'

      write(nupar,*)' naf_alkhqp_compute     = ',                          &
     & this%naf_alkhqp%m_compute, '$'
      write(nupar,*)' naf_alkhqp_stepcalc     = ',                         &
     & this%naf_alkhqp%m_stepcalc, '$'
      write(nupar,*)' naf_alkhqp_stepout     = ',                          &
     & this%naf_alkhqp%m_stepout, '$'
      write(nupar,*)' naf_alkhqp_elltype     = ',                          &
     & this%naf_alkhqp%m_elltype, '$'
      write(nupar,*)' naf_alkhqp_nterm     = ',                            &
     & this%naf_alkhqp%m_naf_nterm, '$'
      write(nupar,*)' naf_alkhqp_isec     = ',                             &
     & this%naf_alkhqp%m_naf_isec, '$'
      write(nupar,*)' naf_alkhqp_iw     = ',                               &
     & this%naf_alkhqp%m_naf_iw, '$'
      write(nupar,fmt='(A,D23.16)')' naf_alkhqp_dtour     = ',             &
     & this%naf_alkhqp%m_naf_dtour, '$'
      write(nupar,fmt='(A,D23.16)')' naf_alkhqp_tol     = ',               &
     & this%naf_alkhqp%m_naf_tol, '$'

      write(nupar,*)' naf_diffalp_compute     = ',                          &
     & this%naf_diffalp%m_compute, '$'
      write(nupar,*)' naf_diffalp_stepcalc     = ',                         &
     & this%naf_diffalp%m_stepcalc, '$'
      write(nupar,*)' naf_diffalp_stepout     = ',                          &
     & this%naf_diffalp%m_stepout, '$'
      write(nupar,*)' naf_diffalp_elltype     = ',                          &
     & this%naf_diffalp%m_elltype, '$'
      write(nupar,*)' naf_diffalp_nterm     = ',                            &
     & this%naf_diffalp%m_naf_nterm, '$'
      write(nupar,*)' naf_diffalp_isec     = ',                             &
     & this%naf_diffalp%m_naf_isec, '$'
      write(nupar,*)' naf_diffalp_iw     = ',                               &
     & this%naf_diffalp%m_naf_iw, '$'
      write(nupar,fmt='(A,D23.16,A)')' naf_diffalp_dtour     = ',             &
     & this%naf_diffalp%m_naf_dtour, '$'
      write(nupar,fmt='(A,D23.16,A)')' naf_diffalp_tol     = ',               &
     & this%naf_diffalp%m_naf_tol, '$'
      write(nupar,*)' dim naf_diffalp_pla[1:1]$'
      write(nupar,*)' naf_diffalp_pla[1]     = ',                          &
     & this%naf_diffalp_pla(1), '$'

      write(nupar,*)' ctrl_diststar_compute     = ',                       &
     & this%ctrl_diststar%m_compute, '$'
      write(nupar,*)' ctrl_diststar_stepcalc    = ',                       &
     & this%ctrl_diststar%m_stepcalc, '$'
      write(nupar,*)' ctrl_diststar_distmin     = ',                       &
     & this%ctrl_diststar%m_distmin, '$'
      write(nupar,*)' ctrl_diststar_distmax     = ',                       &
     & this%ctrl_diststar%m_distmax, '$'
      
      write(nupar,*)' ctrl_distpla_compute      = ',                        &
     & this%ctrl_distpla%m_compute, '$'
      write(nupar,*)' ctrl_distpla_stepcalc     = ',                        &
     & this%ctrl_distpla%m_stepcalc, '$'
      write(nupar,*)' ctrl_distpla_distmin      = ',                        &
     & this%ctrl_distpla%m_distmin, '$'

      write(nupar,*)' naf_alkhqp_pla_compute     = ',                          &
     & this%naf_alkhqp_pla%m_compute, '$'
      write(nupar,*)' naf_alkhqp_pla_stepcalc     = ',                         &
     & this%naf_alkhqp_pla%m_stepcalc, '$'
      write(nupar,*)' naf_alkhqp_pla_stepout     = ',                          &
     & this%naf_alkhqp_pla%m_stepout, '$'
      write(nupar,*)' naf_alkhqp_pla_elltype     = ',                          &
     & this%naf_alkhqp_pla%m_elltype, '$'
      write(nupar,*)' naf_alkhqp_pla_nterm     = ',                            &
     & this%naf_alkhqp_pla%m_naf_nterm, '$'
      write(nupar,*)' naf_alkhqp_pla_isec     = ',                             &
     & this%naf_alkhqp_pla%m_naf_isec, '$'
      write(nupar,*)' naf_alkhqp_pla_iw     = ',                               &
     & this%naf_alkhqp_pla%m_naf_iw, '$'
      write(nupar,fmt='(A,D23.16,A)')' naf_alkhqp_pla_dtour     = ',           &
     & this%naf_alkhqp_pla%m_naf_dtour, '$'
      write(nupar,fmt='(A,D23.16,A)')' naf_alkhqp_pla_tol     = ',             &
     & this%naf_alkhqp_pla%m_naf_tol, '$'

      ! ecriture de la solution orbitale forcee des planetes
      if (this%if_orb_pla.ne.0) then
        call this%orb_pla%fwritetrip(nupar)
      endif
      
      end  subroutine filepar_fwritetrip  

!***********************************************************************
!> @brief ecrit  des parametres de this dans le fichier 
!! dont le nom est genere a partir de chemin  et nf_rad\n
!! Le fichier genere est lisible par trip
!***********************************************************************
      subroutine filepar_fwritepartrip(this)
       implicit none
       class(t_filepar), intent(inout):: this    !< fichier de parametres

      character(100) nf_par
      integer nu_par

      nf_par = trim(this%chemin)//trim(this%nf_rad)//'.tpar'
      nu_par = 11
      open(nu_par,file=nf_par,status='unknown')
      call filepar_fwritetrip(this, nu_par)
      close(nu_par)
      
      end  subroutine filepar_fwritepartrip        

!***********************************************************************
!> @brief lit le fichier de parametres name contenant la solution orbitale "forcee"
! et remplit this
!***********************************************************************
      subroutine orb_pla_tabulee_par_fread(this, fname)
       implicit none
       class(t_orb_pla_tabulee_par), intent(inout):: this    !< fichier de parametres
       character(len=*), intent(in) :: fname     !< nom du fichier a lire
       
        character(100)  :: orb_pla_tabulee_nf
        integer orb_pla_tabulee_coord

       namelist/orb_pla_tabulee/ orb_pla_tabulee_nf,                    &
     &    orb_pla_tabulee_coord
      
       open(10,file=fname,status='old')
       read(10,orb_pla_tabulee)
       close(10)
       
       this%m_coord = orb_pla_tabulee_coord
       this%m_nf = orb_pla_tabulee_nf
      
      end  subroutine orb_pla_tabulee_par_fread  

!***********************************************************************
!> @brief ecrit  des parametres de this dans le fichier  nupar
!***********************************************************************
      subroutine orb_pla_tabulee_par_fwrite(this, nupar)
       implicit none
       class(t_orb_pla_tabulee_par), intent(inout):: this    !< fichier de parametres
       integer, intent(in) :: nupar     !< numero de fichier

      write(nupar,*)' orb_pla_tabulee_coord = ', this%m_coord
      write(nupar,*)' orb_pla_tabulee_nf    = ',                             &
     & this%m_nf (1:len_trim(this%m_nf))

      end  subroutine orb_pla_tabulee_par_fwrite  

!***********************************************************************
!> @brief ecrit  des parametres de this dans le fichier  nupar.
!! Le fichier est lible par trip
!***********************************************************************
      subroutine orb_pla_tabulee_par_fwritetrip(this, nupar)
       implicit none
       class(t_orb_pla_tabulee_par), intent(inout):: this    !< fichier de parametres
       integer, intent(in) :: nupar     !< numero de fichier

      write(nupar,*)' orb_pla_tabulee_coord = ', this%m_coord,'$'
      write(nupar,*)' orb_pla_tabulee_nf    = "',                            &
     & this%m_nf (1:len_trim(this%m_nf)), '"$'

      end  subroutine orb_pla_tabulee_par_fwritetrip  

!***********************************************************************
!> creation des noms de fichiers de sortie pour les fichiers communs 
!! a toutes les conditions initiales
!***********************************************************************
      subroutine extfilepar_create_filename(this, rank)
      implicit none
      class(t_extfilepar), intent(inout):: this    !< fichier de parametres
      integer, intent(in) :: rank !< numero du rang MPI
      integer :: l_ch,l_rad,l_nf,ll
      character(len=3) :: supprad ! extension a ajouter au radical
      character(512) nf
    
      write(supprad,fmt='(I3.3)') rank
      l_ch = len_trim(this%chemin)
      l_rad = len_trim(this%nf_rad)
      ll = len_trim(this%int_type)
      l_nf = l_ch+l_rad
      nf = trim(this%chemin)//trim(this%nf_rad)
      nf = trim(nf)//'_proc'//trim(supprad)
      this%nf_ci_pla = trim(nf)//'.ci_pla'
      this%nf_ci_part = trim(nf)//'.ci_part'
      this%nf_control_pla = trim(nf)//'.control'
      this%nf_control = trim(nf)//'.control_part'
      this%nf_minmax_aei = trim(nf)//'.minmax_aei_part'
      this%nf_minmax_diffalp = trim(nf)//'.minmax_diffalp_part'
      select case (this%naf_alkhqp%m_compute)  
       case(0)
         this%nf_naf_alkhqp = "" 
       case(1)
         this%nf_naf_alkhqp = trim(nf)//'.naf_alkhqp_part'
       case (2)
         this%nf_naf_alkhqp = trim(nf)//'.naf_alkh_part'
       case default
         stop "extfilepar_create_filename : naf_alkhqp_compute inconnu"
      end select
      this%nf_naf_diffalp = trim(nf)//'.naf_diffalp_part'

      select case (this%naf_alkhqp_pla%m_compute)  
       case(0)
         this%nf_naf_alkhqp_pla = ""
       case(1)
         this%nf_naf_alkhqp_pla = trim(nf)//'.naf_alkhqp'
       case (2)
         this%nf_naf_alkhqp_pla  = trim(nf)//'.naf_alkh'
       case default
         write(*,*) "extfilepar_create_filename : "
         write(*,*) "naf_alkhqp_pla_compute inconnu"
         stop 1
      end select

      this%nf_ell = trim(nf)//'.ell_part'
      this%nf_car = trim(nf)//'.car_part'

      write(*,*) trim(this%nf_ci_pla)
      write(*,*) trim(this%nf_ci_part)
      write(*,*) trim(this%nf_control_pla)
      write(*,*) trim(this%nf_control)
      write(*,*) trim(this%nf_minmax_aei)
      write(*,*) trim(this%nf_minmax_diffalp)
      write(*,*) trim(this%nf_naf_alkhqp)
      write(*,*) trim(this%nf_naf_diffalp)
      write(*,*) trim(this%nf_naf_alkhqp_pla)
      write(*,*) trim(this%nf_ell)
      write(*,*) trim(this%nf_car)
      
      call extfilepar_create_filename_ci(this,'proc'//trim(supprad))
      
      end  subroutine extfilepar_create_filename    

!***********************************************************************
!> creation des noms de fichiers de sortie pour une seule condition initiale
!***********************************************************************
      subroutine extfilepar_create_filename_ci(this, supprad)
      implicit none
      class(t_extfilepar), intent(inout):: this    !< fichier de parametres
      character(len=*), intent(in) :: supprad !< extension a ajouter au radical
      integer :: l_ch,l_rad,l_nf,ll
      character(512) nf
    
      l_ch = len_trim(this%chemin)
      l_rad = len_trim(this%nf_rad)
      ll = len_trim(this%int_type)
      l_nf = l_ch+l_rad
      nf = trim(this%chemin)//trim(this%nf_rad)
      nf = trim(nf)//'_'//trim(supprad)
      this%nf_int = trim(nf)//'.int'
      this%nf_ell_pla = trim(nf)//'.ell'
      this%nf_car_pla = trim(nf)//'.car'
      write(*,*) trim(this%nf_int)
      write(*,*) trim(this%nf_ell_pla)
      write(*,*) trim(this%nf_car_pla)
     
      end  subroutine extfilepar_create_filename_ci    

!***********************************************************************
!> ouverture en ecriture des fichiers communs 
!! a toutes les conditions initiales
!***********************************************************************
      subroutine extfilepar_open_file(this)
      implicit none
      class(t_extfilepar), intent(inout):: this    !< fichier de parametres
    
      !
      ! pour les planetes
      !
      ! sortie de controle des integrations du systeme planetaire
      this%io_control_pla = 63
      open(this%io_control_pla,file=this%nf_control_pla,                 &
     &     status='unknown')
      ! sortie des analyses en frequence en alkhqp pour les planetes
      if (this%naf_alkhqp_pla%m_compute.ne.0) then
       this%io_naf_alkhqp_pla = t_io_txt_ncomp(                         &
     &    "I/O TXT naf alkhqp_pla")
       call this%io_naf_alkhqp_pla%set_name(trim(                       &
     &    this%nf_naf_alkhqp_pla),70)
       call this%io_naf_alkhqp_pla%fcreate()
      endif
      
      !
      ! pour les particules
      !
      ! sortie de controle des integrations des particules
      this%io_control = 67
      open(this%io_control,file=this%nf_control, status='unknown')

      ! sortie elliptique 
      if (this%if_ell_part.eq.1) then
         this%io_ell = t_io_txt_ell_part("I/O PART ell")
         call this%io_ell%set_name(trim(this%nf_ell),73)
         call this%io_ell%fcreate()
      end if
      
      ! sortie cartesienne
      if (this%if_car_part.eq.1) then
         this%io_cart = t_io_txt_car_part("I/O PART cart")
         call this%io_cart%set_name(trim(this%nf_car),72)
         call this%io_cart%fcreate()
      end if

      ! sortie des min/max en a,e,i
      if (this%minmax_aei%m_compute.eq.1) then
       this%io_minmax_aei = t_io_txt_mligncomp("I/O TXT min/max a,e,i")
       call this%io_minmax_aei%set_name(trim(this%nf_minmax_aei)        &
     &      ,64)
       call this%io_minmax_aei%fcreate()
      endif
      
      ! sortie des min/max en a_particule-a_planete, ....
      if (this%minmax_diffalp%m_compute.eq.1) then
       this%io_minmax_diffalp = t_io_txt_mligncomp                       &
     &      ("I/O min/max a_partk-a_plan,l_partk-l_plan,...")
       call this%io_minmax_diffalp%set_name(                             &
     &       trim(this%nf_minmax_diffalp) ,65)
       call this%io_minmax_diffalp%fcreate()
      endif
     

      ! sortie des analyses en frequence en alkhqp pour les particules
      if (this%naf_alkhqp%m_compute.ne.0) then
       this%io_naf_alkhqp = t_io_txt_mligncomp("I/O TXT naf alkhqp")
       call this%io_naf_alkhqp%set_name(trim(this%nf_naf_alkhqp)        &
     &      ,68)
       call this%io_naf_alkhqp%fcreate()
      endif
      
      ! sortie des  analyses en frequence en exp(i*(l_particule-l_planete)), ....
      if (this%naf_diffalp%m_compute.eq.1) then
       this%io_naf_diffalp = t_io_txt_mligncomp                         &
     &      ("I/O TXT naf exp(i(l_partk-l_plan)),exp(i(pi...-pi...))")
       call this%io_naf_diffalp%set_name(                               &
     &       trim(this%nf_naf_diffalp) ,69)
       call this%io_naf_diffalp%fcreate()
      endif

      end  subroutine extfilepar_open_file    

      end module mod_filepar

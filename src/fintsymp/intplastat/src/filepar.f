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
! history : creation 24/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour gerer le fichier de parametres
!***********************************************************************
      module mod_filepar
        use mod_io_txt_ncomp
        
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
!! classe de gestion des parametres pour le controle de distance a l'etoile : 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_ctrl_diststar_par
      
         !> * = 0 => pas de controle des distances
         !! * = 1 => controle des distances avec notification simple
         !! * = 2 => controle des distances avec dump en coordonnees cartesiennes 
         !! * = 3 => controle des distances avec dump en coordonnees elliptiques 
         integer          m_compute 
         integer(8)  ::   m_stepcalc !< pas de calcul du controle 
         real(TREAL) ::   m_distmin !< distance minimale a l'etoile. si m_distmin =-1, pas de controle, si m_distmin>=0, genere une erreur si distance<m_distmin
         real(TREAL) ::   m_distmax !< distance maximale a l'etoile. si m_distmax =-1, pas de controle, si m_distmax>=0, genere une erreur si distance>m_distmax
      end type t_ctrl_diststar_par

!***********************************************************************
!> @class t_ctrl_distpla_par
!! classe de gestion des parametres pour le controle de distance entre planetes : 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_ctrl_distpla_par
      
         !> * = 0 => pas de controle des distances
         !! * = 1 => controle des distances avec notification simple
         !! * = 2 => controle des distances avec dump en coordonnees cartesiennes 
         !! * = 3 => controle des distances avec dump en coordonnees elliptiques 
         integer          m_compute 
         integer(8)  ::   m_stepcalc !< pas de calcul du controle 
         character(100) ::  m_nfdistmin !< fichier de distances minimales entre planetes.
         ! Ce fichier doit avoir le meme nombre de lignes que celui des conditions initiales 
      end type t_ctrl_distpla_par

!***********************************************************************
!> @class t_ctrl_energie_par
!! classe de gestion des parametres pour controle de l'energie 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_ctrl_energie_par
      
         !> * = 0 => pas de controle de l'energie
         !! * = 1 => controle de l'energie
         integer          m_compute 
         integer(8)  ::   m_stepcalc !< pas de calcul du controle 
         real(TREAL) ::   m_relenermax !< difference maximale de l'erreur relative sur l'energie. si l'erreur relative de l'energie >m_relenermax =-1, alors on genere une erreur 
      end type t_ctrl_energie_par

!***********************************************************************
!> @class t_filepar
!! classe de gestion du fichier de parametres : 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_filepar

        character(100)  :: chemin,nf_rad,int_type
        character(100) :: nf_initext
        integer :: if_int,if_ell,if_car, if_dump, if_restart
        integer(8) :: n_iter,n_out, n_dump
        integer ::    out_ell,if_invar
        real(TREAL) :: tinit,dt,dmax
        integer :: ref_gmsun
        integer :: type_pas
        
        ! pour les statistiques des min/moy/max
        type(t_minmax_par) :: minmax_aei  !< min/moy/max en a, e, i
        type(t_minmax_par) :: minmax_diffalp !< min/moy/max en l1-l2,pibar1-pibar2, a1-a2 entre [-pi,pi] et [0, 2*pi]
        integer,dimension(1:2) :: minmax_diffalp_pla !< numero des 2 planetes utilisees pour minmax_diffalp
        type(t_minmax_par) :: minmax_diffalc !< min/moy/max en l1-l2,pibar1-pibar2, a1-a2 avec angle continu
        integer,dimension(1:2) :: minmax_diffalc_pla !< numero des 2 planetes utilisees pour minmax_diffalpc
        type(t_minmax_par) :: minmax_diffae2 !< min/moy/max en (a1-a2)**2, (e1-e2)**2, (a1-a2)**2+(e1-e2)**2
        integer,dimension(1:2) :: minmax_diffae2_pla !< numero des 2 planetes utilisees pour minmax_diffae2


        ! pour les analyses en frequences
        type(t_naf_par) :: naf_alkhqp  !< analyse en frequence en a*exp(*il),k+i*h, q+i*p
        type(t_naf_par) :: naf_diffalp  !< analyse en frequence en exp(i*(l1-l2)),exp(i*(pi1-pi2))
        integer,dimension(1:2) :: naf_diffalp_pla !< numero des 2 planetes utilisees pour naf_diffalp
        type(t_naf_par) :: naf_difflpm  !< analyse en frequence en ((l1-l2),0) et ((pi1-pi2),0) modulo [0,2pi] et [-pi,pi]
        integer,dimension(1:2) :: naf_difflpm_pla !< numero des 2 planetes utilisees pour naf_difflpm
        
        ! controle de la distance a l'etoile
        type(t_ctrl_diststar_par) :: ctrl_diststar !< controle de la distance a l'etoile
        ! controle de la distance entre planetes
        type(t_ctrl_distpla_par) :: ctrl_distpla   !< controle de la distance entre planetes
        ! controle de l'energie
        type(t_ctrl_energie_par) :: ctrl_energie !< controle de l'energie

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

        character(len=512) :: nf_int !< nom du fichier de sortie des integrales premieres
        character(len=512) :: nf_car !< nom du fichier de sortie des elements cartesiens
        character(len=512) :: nf_ell !< nom du fichier de sortie des elements elliptiques
        character(len=512) :: nf_sch !< nom du fichier de sortie des donnees du schema adaptatif
        character(len=512) :: nf_dump !< nom du fichier de sortie des donnees des dump bnaires
        character(len=512) :: nf_ci !< nom du fichier des conditions initiales
        character(len=512) :: nf_control !< nom du fichier de controle de l'integration
        character(len=512) :: nf_ctrlstar_dump !< nom du fichier du dump sur la distance a l'etoile lors d'un arret
        character(len=512) :: nf_ctrlpla_dump  !< nom du fichier du dump sur les distances entre planetes lors d'un arret
        character(len=512) :: nf_cictrlpladistmin !< nom du fichier des conditions initiales pour les distances minimales entre planetes
        character(len=512) :: nf_minmax_aei !< nom du fichier de sortie des min/max a,e,i
        character(len=512) :: nf_minmax_diffalp !< nom du fichier de sortie des min/max a_p1-a_p2,... entre [-pi,pi] et [0, 2*pi]
        character(len=512) :: nf_minmax_diffalc !< nom du fichier de sortie des min/max a_p1-a_p2,... avec angle continu
        character(len=512) :: nf_minmax_diffae2 !< nom du fichier de sortie des min/max (a_p1-a_p2)**2,...
        character(len=512) :: nf_naf_alkhqp !< nom du fichier de sortie des analyses en frequence en alkhqp
        character(len=512) :: nf_naf_diffalp !< nom du fichier de sortie des analyses en frequence en exp(i*(l1-l2)), exp(i*(pi1-pi2)) 
        character(len=512) :: nf_naf_difflpm !< nom du fichier de sortie des analyses en frequence en (l1-l2) 
        
        integer :: int_var   !< variables de l'integrateur (: 0 : Jacobi  1 :helio canonique)

         type (t_io_txt_ncomp) :: io_minmax_aei  !< sortie des min/max en a, e, i
         type (t_io_txt_ncomp) :: io_minmax_diffalp  !< sortie des min/max des differences a_p1-a_p2,.... entre [-pi,pi] et [0, 2*pi]
         type (t_io_txt_ncomp) :: io_minmax_diffalc  !< sortie des min/max des differences a_p1-a_p2,.... avec angle continu
         type (t_io_txt_ncomp) :: io_minmax_diffae2  !< sortie des min/max des differences (a_p1-a_p2)**2,....

         integer :: io_control     !< numero du fichier du controle des integrations
         type (t_io_txt_ncomp) :: io_ctrlstar_dump !< sortie du dump sur la distance a l'etoile lors d'un arret
         type (t_io_txt_ncomp) :: io_ctrlpla_dump  !< sortie du dump sur les distances entre planetes lors d'un arret
         

         type (t_io_txt_ncomp) :: io_naf_alkhqp  !< sortie des analyses en frenquence en a*exp(*il),k+i*h, q+i*p
         type (t_io_txt_ncomp) :: io_naf_diffalp !< sortie des analyses en frenquence en exp(i*(l1-l2)),exp(i*(pi1-pi2))
         type (t_io_txt_ncomp) :: io_naf_difflpm !< sortie des analyses en frenquence en (l1-l2)

       contains
                      
        procedure :: create_filename => extfilepar_create_filename  ! creation des noms de fichiers de sortie  communs a toutes les conditions initiales                     
        procedure ::create_filename_ci=>extfilepar_create_filename_ci  ! creation des noms de fichiers de sortie propre a une condition initiale                     
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
        character(100) :: nf_initext
        integer :: if_int,if_ell,if_car, if_dump, if_restart
        integer(8) :: n_iter,n_out, n_dump
        integer :: out_ell,if_invar, ref_gmsun
        real(TREAL) :: tinit,dt
        integer :: type_pas
        integer :: minmax_aei_compute
        integer(8) :: minmax_aei_stepcalc, minmax_aei_stepout
        integer :: minmax_aei_elltype
        integer :: minmax_diffalp_compute
        integer(8) :: minmax_diffalp_stepcalc, minmax_diffalp_stepout
        integer :: minmax_diffalp_elltype
        integer, dimension(1:2) :: minmax_diffalp_pla
        integer :: minmax_diffalc_compute
        integer(8) :: minmax_diffalc_stepcalc, minmax_diffalc_stepout
        integer :: minmax_diffalc_elltype
        integer, dimension(1:2) :: minmax_diffalc_pla
        integer :: minmax_diffae2_compute
        integer(8) :: minmax_diffae2_stepcalc, minmax_diffae2_stepout
        integer :: minmax_diffae2_elltype
        integer, dimension(1:2) :: minmax_diffae2_pla
        integer(8) :: rlen
        integer :: naf_alkhqp_compute
        integer(8) :: naf_alkhqp_stepcalc, naf_alkhqp_stepout
        integer :: naf_alkhqp_elltype
        integer naf_alkhqp_nterm
        integer naf_alkhqp_isec, naf_alkhqp_iw 
        real(TREAL) naf_alkhqp_dtour, naf_alkhqp_tol 
        integer :: naf_diffalp_compute
        integer(8) :: naf_diffalp_stepcalc, naf_diffalp_stepout
        integer :: naf_diffalp_elltype
        integer, dimension(1:2) :: naf_diffalp_pla
        integer naf_diffalp_nterm
        integer naf_diffalp_isec, naf_diffalp_iw 
        real(TREAL) naf_diffalp_dtour, naf_diffalp_tol 
        
        integer :: naf_difflpm_compute
        integer(8) :: naf_difflpm_stepcalc, naf_difflpm_stepout
        integer :: naf_difflpm_elltype
        integer, dimension(1:2) :: naf_difflpm_pla
        integer naf_difflpm_nterm
        integer naf_difflpm_isec, naf_difflpm_iw 
        real(TREAL) naf_difflpm_dtour, naf_difflpm_tol 

        integer :: ctrl_diststar_compute
        integer(8) :: ctrl_diststar_stepcalc
        real(TREAL) :: ctrl_diststar_distmin
        real(TREAL) :: ctrl_diststar_distmax

        integer :: ctrl_distpla_compute
        integer(8) :: ctrl_distpla_stepcalc
        character(100) :: ctrl_distpla_nfdistmin
        
        integer :: ctrl_energie_compute
        integer(8) :: ctrl_energie_stepcalc
        real(TREAL) :: ctrl_energie_relenermax


      namelist/lect/chemin,nf_rad,nf_initext,
     &   int_type,tinit,dt,n_iter,n_out, n_dump,
     &   out_ell,int_type, ref_gmsun,type_pas,
     &   if_invar,if_int,if_ell,if_car, if_dump, if_restart,
     &   minmax_aei_compute,
     &   minmax_aei_stepcalc, minmax_aei_stepout,
     &   minmax_aei_elltype,
     &   minmax_diffalp_compute,
     &   minmax_diffalp_stepcalc, minmax_diffalp_stepout,
     &   minmax_diffalp_elltype,
     &   minmax_diffalp_pla,
     &   minmax_diffalc_compute,
     &   minmax_diffalc_stepcalc, minmax_diffalc_stepout,
     &   minmax_diffalc_elltype,
     &   minmax_diffalc_pla,
     &   minmax_diffae2_compute,
     &   minmax_diffae2_stepcalc, minmax_diffae2_stepout,
     &   minmax_diffae2_elltype,
     &   minmax_diffae2_pla,
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
     &   naf_difflpm_compute,
     &   naf_difflpm_stepcalc, naf_difflpm_stepout,
     &   naf_difflpm_elltype,  naf_difflpm_pla,
     &   naf_difflpm_nterm,
     &   naf_difflpm_isec, naf_difflpm_iw, 
     &   naf_difflpm_dtour, naf_difflpm_tol, 
     &   ctrl_diststar_compute, ctrl_diststar_stepcalc, 
     &   ctrl_diststar_distmin, ctrl_diststar_distmax,
     &   ctrl_distpla_compute, ctrl_distpla_stepcalc, 
     &   ctrl_distpla_nfdistmin, 
     &   ctrl_energie_compute, ctrl_energie_stepcalc,
     &   ctrl_energie_relenermax
      
      ref_gmsun = 1
      type_pas = 0
      n_dump = 0
      if_restart = 0
      if_dump = 0
      open(10,file=fname,status='old')
      read(10,lect)
      close(10)
      
      this%chemin = chemin
      this%nf_rad = nf_rad
      this%nf_initext = nf_initext
      this%int_type = int_type
      this%tinit = tinit
      this%dt = dt
      this%n_iter = n_iter
      this%n_out = n_out
      this%n_dump = n_dump
      this%out_ell = out_ell
      this%int_type = int_type
      this%if_invar = if_invar
      this%if_int = if_int
      this%if_ell = if_ell
      this%if_car = if_car
      this%if_dump = if_dump
      this%if_restart = if_restart
      this%ref_gmsun = ref_gmsun
      this%type_pas = type_pas

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

      ! initialisation des min/max des diff alc
      this%minmax_diffalc%m_compute  = minmax_diffalc_compute
      this%minmax_diffalc%m_stepcalc = minmax_diffalc_stepcalc
      this%minmax_diffalc%m_stepout  = minmax_diffalc_stepout
      this%minmax_diffalc%m_elltype  = minmax_diffalc_elltype
      this%minmax_diffalc_pla = minmax_diffalc_pla
      ! verification de elltype
      if (minmax_diffalc_compute.eq.1) then 
      ! verification de minmax_diffalp_elltype
       if (minmax_diffalc_elltype.ne.6) then
        write(*,*) "minmax_diffalc_elltype doit etre 6"
        stop
       endif
       !verification  de minmax_diffalp_stepout
       rlen = mod(minmax_diffalc_stepout,minmax_diffalc_stepcalc)
       if (rlen.ne.0) then
        stop 'minmax_diffalc_stepout doit etre multiple de'//              &
     &       ' minmax_diffalc_stepcalc'
        
       endif
      endif


      ! initialisation des min/max des diff ae2
      this%minmax_diffae2%m_compute  = minmax_diffae2_compute
      this%minmax_diffae2%m_stepcalc = minmax_diffae2_stepcalc
      this%minmax_diffae2%m_stepout  = minmax_diffae2_stepout
      this%minmax_diffae2%m_elltype  = minmax_diffae2_elltype
      this%minmax_diffae2_pla = minmax_diffae2_pla
      ! verification de elltype
      if (minmax_diffae2_compute.eq.1) then 
      ! verification de minmax_diffae2_elltype
       if (minmax_diffae2_elltype.ne.6) then
        write(*,*) "minmax_diffae2_elltype doit etre 6"
        stop
       endif
       !verification  de minmax_diffae2_stepout
       rlen = mod(minmax_diffae2_stepout,minmax_diffae2_stepcalc)
       if (rlen.ne.0) then
        stop 'minmax_diffae2_stepout doit etre multiple de'//              &
     &       ' minmax_diffae2_stepcalc'
        
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
      
       ! initialisation des naf des diff lam
      this%naf_difflpm%m_compute  = naf_difflpm_compute
      this%naf_difflpm%m_stepcalc = naf_difflpm_stepcalc
      this%naf_difflpm%m_stepout  = naf_difflpm_stepout
      this%naf_difflpm%m_elltype  = naf_difflpm_elltype
      this%naf_difflpm_pla = naf_difflpm_pla
      this%naf_difflpm%m_naf_nterm  = naf_difflpm_nterm 
      this%naf_difflpm%m_naf_iw = naf_difflpm_iw 
      this%naf_difflpm%m_naf_isec = naf_difflpm_isec 
      this%naf_difflpm%m_naf_dtour = naf_difflpm_dtour 
      this%naf_difflpm%m_naf_tol  = naf_difflpm_tol 
      ! verifications 
      if (naf_difflpm_compute.eq.1) then 
       !verification  de naf_difflpm_elltype
       if ((naf_difflpm_elltype.ne.6)) then
        stop "naf_difflpm_elltype doit etre 6"
       endif
       !verification  de naf_difflpm_stepout
       if (mod(naf_difflpm_stepout,naf_difflpm_stepcalc).ne.0) then
        stop 'naf_difflpm_stepout doit etre multiple de'//                  &
     &       ' naf_difflpm_stepcalc'
        
       endif
      endif

      ! controle des distances a l'etoile
      this%ctrl_diststar%m_compute = ctrl_diststar_compute
      this%ctrl_diststar%m_stepcalc = ctrl_diststar_stepcalc
      this%ctrl_diststar%m_distmin = ctrl_diststar_distmin
      this%ctrl_diststar%m_distmax = ctrl_diststar_distmax
      
      ! controle des distances a la planete
      this%ctrl_distpla%m_compute = ctrl_distpla_compute
      this%ctrl_distpla%m_stepcalc = ctrl_distpla_stepcalc
      this%ctrl_distpla%m_nfdistmin = ctrl_distpla_nfdistmin

      ! controle de l'energie
      this%ctrl_energie%m_compute = ctrl_energie_compute
      this%ctrl_energie%m_stepcalc = ctrl_energie_stepcalc
      this%ctrl_energie%m_relenermax = ctrl_energie_relenermax

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
      write(nupar,*)' int_type   = ',trim(this%int_type)
      write(nupar,*)' ref_gmsun  = ',this%ref_gmsun
      write(nupar,*)' tinit      = ',this%tinit
      write(nupar,*)' dt         = ',this%dt
      write(nupar,*)' n_iter     = ',this%n_iter
      write(nupar,*)' n_out      = ',this%n_out
      write(nupar,*)' n_dump     = ',this%n_dump
      write(nupar,*)' out_ell    = ',this%out_ell
      write(nupar,*)' if_invar   = ',this%if_invar
      write(nupar,*)' if_int     = ',this%if_int
      write(nupar,*)' if_ell     = ',this%if_ell
      write(nupar,*)' if_car     = ',this%if_car
      write(nupar,*)' if_dump    = ',this%if_dump
      write(nupar,*)' if_restart = ',this%if_restart
      write(nupar,*)' type_pas   = ',this%type_pas

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
      write(nupar,*)' minmax_diffalp_pla(1)     = ',                       &
     & this%minmax_diffalp_pla(1)
      write(nupar,*)' minmax_diffalp_pla(2)     = ',                       &
     & this%minmax_diffalp_pla(2)

      write(nupar,*)' minmax_diffalc_compute     = ',                      &
     & this%minmax_diffalc%m_compute
      write(nupar,*)' minmax_diffalc_stepcalc     = ',                     &
     & this%minmax_diffalc%m_stepcalc
      write(nupar,*)' minmax_diffalc_stepout     = ',                      &
     & this%minmax_diffalc%m_stepout
      write(nupar,*)' minmax_diffalc_elltype     = ',                      &
     & this%minmax_diffalc%m_elltype
      write(nupar,*)' minmax_diffalc_pla(1)     = ',                       &
     & this%minmax_diffalc_pla(1)
      write(nupar,*)' minmax_diffalc_pla(2)     = ',                       &
     & this%minmax_diffalc_pla(2)

      write(nupar,*)' minmax_diffae2_compute     = ',                      &
     & this%minmax_diffae2%m_compute
      write(nupar,*)' minmax_diffae2_stepcalc     = ',                     &
     & this%minmax_diffae2%m_stepcalc
      write(nupar,*)' minmax_diffae2_stepout     = ',                      &
     & this%minmax_diffae2%m_stepout
      write(nupar,*)' minmax_diffae2_elltype     = ',                      &
     & this%minmax_diffae2%m_elltype
      write(nupar,*)' minmax_diffae2_pla(1)     = ',                       &
     & this%minmax_diffae2_pla(1)
      write(nupar,*)' minmax_diffae2_pla(2)     = ',                       &
     & this%minmax_diffae2_pla(2)

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
      write(nupar,*)' naf_diffalp_pla(2)     = ',                           &
     & this%naf_diffalp_pla(2)

      write(nupar,*)' naf_difflpm_compute     = ',                          &
     & this%naf_difflpm%m_compute
      write(nupar,*)' naf_difflpm_stepcalc     = ',                         &
     & this%naf_difflpm%m_stepcalc
      write(nupar,*)' naf_difflpm_stepout     = ',                          &
     & this%naf_difflpm%m_stepout
      write(nupar,*)' naf_difflpm_elltype     = ',                          &
     & this%naf_difflpm%m_elltype
      write(nupar,*)' naf_difflpm_nterm     = ',                            &
     & this%naf_difflpm%m_naf_nterm
      write(nupar,*)' naf_difflpm_isec     = ',                             &
     & this%naf_difflpm%m_naf_isec
      write(nupar,*)' naf_difflpm_iw     = ',                               &
     & this%naf_difflpm%m_naf_iw
      write(nupar,fmt='(A,D23.16)')' naf_difflpm_dtour     = ',             &
     & this%naf_difflpm%m_naf_dtour
      write(nupar,fmt='(A,D23.16)')' naf_difflpm_tol     = ',               &
     & this%naf_difflpm%m_naf_tol
      write(nupar,*)' naf_difflpm_pla(1)     = ',                           &
     & this%naf_difflpm_pla(1)
      write(nupar,*)' naf_difflpm_pla(2)     = ',                           &
     & this%naf_difflpm_pla(2)

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
      write(nupar,*)' ctrl_distpla_nfdistmin    = ',                        &
     & this%ctrl_distpla%m_nfdistmin

      write(nupar,*)' ctrl_energie_compute      = ',                        &
     & this%ctrl_energie%m_compute
      write(nupar,*)' ctrl_energie_stepcalc     = ',                        &
     & this%ctrl_energie%m_stepcalc
      write(nupar,*)' ctrl_energie_relenermax   = ',                        &
     & this%ctrl_energie%m_relenermax

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
      write(nupar,*)' int_type   = ',trim(this%int_type), '$'
      write(nupar,*)' type_pas   = ',this%type_pas, '$'
      write(nupar,*)' tinit      = ',this%tinit, '$'
      write(nupar,*)' dt         = ',this%dt, '$'
      write(nupar,*)' n_iter     = ',this%n_iter, '$'
      write(nupar,*)' n_out      = ',this%n_out, '$'
      write(nupar,*)' n_dump     = ',this%n_dump, '$'
      write(nupar,*)' out_ell    = ',this%out_ell, '$'
      write(nupar,*)' if_invar   = ',this%if_invar, '$'
      write(nupar,*)' if_int     = ',this%if_int, '$'
      write(nupar,*)' if_ell     = ',this%if_ell, '$'
      write(nupar,*)' if_car     = ',this%if_car, '$'
      write(nupar,*)' if_dump    = ',this%if_dump, '$'
      write(nupar,*)' if_restart = ',this%if_restart, '$'

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
      write(nupar,*)' dim minmax_diffalp_pla[1:2]$'
      write(nupar,*)' minmax_diffalp_pla[1]     = ',                        &
     & this%minmax_diffalp_pla(1), '$'
      write(nupar,*)' minmax_diffalp_pla[2]     = ',                        &
     & this%minmax_diffalp_pla(2), '$'

      write(nupar,*)' minmax_diffalc_compute     = ',                      &
     & this%minmax_diffalc%m_compute, '$'
      write(nupar,*)' minmax_diffalc_stepcalc     = ',                     &
     & this%minmax_diffalc%m_stepcalc, '$'
      write(nupar,*)' minmax_diffalc_stepout     = ',                      &
     & this%minmax_diffalc%m_stepout, '$'
      write(nupar,*)' minmax_diffalc_elltype     = ',                      &
     & this%minmax_diffalc%m_elltype, '$'
      write(nupar,*)' dim minmax_diffalc_pla[1:2]$'
      write(nupar,*)' minmax_diffalc_pla[1]     = ',                        &
     & this%minmax_diffalc_pla(1), '$'
      write(nupar,*)' minmax_diffalc_pla[2]     = ',                        &
     & this%minmax_diffalc_pla(2), '$'

      write(nupar,*)' minmax_diffae2_compute     = ',                      &
     & this%minmax_diffae2%m_compute, '$'
      write(nupar,*)' minmax_diffae2_stepcalc     = ',                     &
     & this%minmax_diffae2%m_stepcalc, '$'
      write(nupar,*)' minmax_diffae2_stepout     = ',                      &
     & this%minmax_diffae2%m_stepout, '$'
      write(nupar,*)' minmax_diffae2_elltype     = ',                      &
     & this%minmax_diffae2%m_elltype, '$'
      write(nupar,*)' dim minmax_diffae2_pla[1:2]$'
      write(nupar,*)' minmax_diffae2_pla[1]     = ',                       &
     & this%minmax_diffae2_pla(1), '$'
      write(nupar,*)' minmax_diffae2_pla[2]     = ',                       &
     & this%minmax_diffae2_pla(2), '$'

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
      write(nupar,fmt='(A,D23.16)')' naf_diffalp_dtour     = ',             &
     & this%naf_diffalp%m_naf_dtour, '$'
      write(nupar,fmt='(A,D23.16)')' naf_diffalp_tol     = ',               &
     & this%naf_diffalp%m_naf_tol, '$'
      write(nupar,*)' dim naf_diffalp_pla[1:2]$'
      write(nupar,*)' naf_diffalp_pla[1]     = ',                          &
     & this%naf_diffalp_pla(1), '$'
      write(nupar,*)' naf_diffalp_pla[2]     = ',                          &
     & this%naf_diffalp_pla(2), '$'

      write(nupar,*)' naf_difflpm_compute     = ',                          &
     & this%naf_difflpm%m_compute, '$'
      write(nupar,*)' naf_difflpm_stepcalc     = ',                         &
     & this%naf_difflpm%m_stepcalc, '$'
      write(nupar,*)' naf_difflpm_stepout     = ',                          &
     & this%naf_difflpm%m_stepout, '$'
      write(nupar,*)' naf_difflpm_elltype     = ',                          &
     & this%naf_difflpm%m_elltype, '$'
      write(nupar,*)' naf_difflpm_nterm     = ',                            &
     & this%naf_difflpm%m_naf_nterm, '$'
      write(nupar,*)' naf_difflpm_isec     = ',                             &
     & this%naf_difflpm%m_naf_isec, '$'
      write(nupar,*)' naf_difflpm_iw     = ',                               &
     & this%naf_difflpm%m_naf_iw, '$'
      write(nupar,fmt='(A,D23.16)')' naf_difflpm_dtour     = ',             &
     & this%naf_difflpm%m_naf_dtour, '$'
      write(nupar,fmt='(A,D23.16)')' naf_difflpm_tol     = ',               &
     & this%naf_difflpm%m_naf_tol, '$'
      write(nupar,*)' dim naf_difflpm_pla[1:2]$'
      write(nupar,*)' naf_difflpm_pla[1]     = ',                          &
     & this%naf_difflpm_pla(1), '$'
      write(nupar,*)' naf_difflpm_pla[2]     = ',                          &
     & this%naf_difflpm_pla(2), '$'

      write(nupar,*)' ctrl_diststar_compute     = ',                       &
     & this%ctrl_diststar%m_compute, '$'
      write(nupar,*)' ctrl_diststar_stepcalc    = ',                       &
     & this%ctrl_diststar%m_stepcalc, '$'
      write(nupar,*)' ctrl_diststar_distmin     = ',                       &
     & this%ctrl_diststar%m_distmin, '$'
      write(nupar,*)' ctrl_diststar_distmax     = ',                       &
     & this%ctrl_diststar%m_distmax, '$'
      
      write(nupar,*)' ctrl_distpla_compute      = ',                       &
     & this%ctrl_distpla%m_compute, '$'
      write(nupar,*)' ctrl_distpla_stepcalc     = ',                       &
     & this%ctrl_distpla%m_stepcalc, '$'
      write(nupar,*)' ctrl_distpla_nfdistmin    = "',                      &
     & trim(this%ctrl_distpla%m_nfdistmin), '"$'


      write(nupar,*)' ctrl_energie_compute      = ',                       &
     & this%ctrl_energie%m_compute, '$'
      write(nupar,*)' ctrl_energie_stepcalc     = ',                       &
     & this%ctrl_energie%m_stepcalc, '$'
      write(nupar,*)' ctrl_energie_relenermax   = ',                       &
     & this%ctrl_energie%m_relenermax, '$'

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
      this%nf_ci = trim(nf)//'.ci'
      
      this%nf_control = trim(nf)//'.control'
      ! controle a l'etoile
      select case (this%ctrl_diststar%m_compute)  
       case (0,1)
         this%nf_ctrlstar_dump = ''
       case (2)
         this%nf_ctrlstar_dump = trim(nf)//'.ctrlstar_car'
       case (3)
         this%nf_ctrlstar_dump = trim(nf)//'.ctrlstar_ell'
       case default
         stop "extfilepar_create_filename : "                            &
     &        //"ctrl_diststar_compute inconnu"
      end select
      
      ! controle entre planetes
      this%nf_cictrlpladistmin=trim(nf)//'.ci_ctrlpladistmin'
      select case (this%ctrl_distpla%m_compute)  
       case (0,1)
         this%nf_ctrlpla_dump = ''
       case (2)
         this%nf_ctrlpla_dump = trim(nf)//'.ctrlpla_car'
       case (3)
         this%nf_ctrlpla_dump = trim(nf)//'.ctrlpla_ell'
       case default
         stop "extfilepar_create_filename : "                            &
     &        //"ctrl_distpla_compute inconnu"
      end select

      
      this%nf_minmax_aei = trim(nf)//'.minmax_aei'
      this%nf_minmax_diffalp = trim(nf)//'.minmax_diffalp'
      this%nf_minmax_diffalc = trim(nf)//'.minmax_diffalc'
      this%nf_minmax_diffae2 = trim(nf)//'.minmax_diffae2'
      select case (this%naf_alkhqp%m_compute)  
       case(0)
       case(1)
         this%nf_naf_alkhqp = trim(nf)//'.naf_alkhqp'
       case (2)
         this%nf_naf_alkhqp = trim(nf)//'.naf_alkh'
       case default
         stop "extfilepar_create_filename : naf_alkhqp_compute inconnu"
      end select
      this%nf_naf_diffalp = trim(nf)//'.naf_diffalp'
      this%nf_naf_difflpm = trim(nf)//'.naf_difflpm'

      write(*,*) trim(this%nf_ci)
      write(*,*) trim(this%nf_control)
      write(*,*) trim(this%nf_ctrlstar_dump)
      write(*,*) trim(this%nf_ctrlpla_dump)
      write(*,*) trim(this%nf_cictrlpladistmin)
      write(*,*) trim(this%nf_minmax_aei)
      write(*,*) trim(this%nf_minmax_diffalp)
      write(*,*) trim(this%nf_minmax_diffalc)
      write(*,*) trim(this%nf_minmax_diffae2)
      write(*,*) trim(this%nf_naf_alkhqp)
      write(*,*) trim(this%nf_naf_diffalp)
      write(*,*) trim(this%nf_naf_difflpm)
      
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
      this%nf_ell = trim(nf)//'.ell'
      this%nf_car = trim(nf)//'.car'
      this%nf_sch = trim(nf)//'.sch'  
      this%nf_dump = trim(nf)//'.dump'  
      write(*,*) trim(this%nf_int)
      write(*,*) trim(this%nf_ell)
      write(*,*) trim(this%nf_car)
      write(*,*) trim(this%nf_sch)
      write(*,*) trim(this%nf_dump)
    
      end  subroutine extfilepar_create_filename_ci    

!***********************************************************************
!> ouverture en ecriture des fichiers communs 
!! a toutes les conditions initiales
!! valeur du dernier fichier ouvert = 72
!***********************************************************************
      subroutine extfilepar_open_file(this)
      implicit none
      class(t_extfilepar), intent(inout):: this    !< fichier de parametres
    
      ! sortie de controle des integrations
      this%io_control = 63
      open(this%io_control,file=this%nf_control,status='unknown')

      ! sortie des min/max en a,e,i
      if (this%minmax_aei%m_compute.eq.1) then
       this%io_minmax_aei = t_io_txt_ncomp("I/O TXT min/max a,e,i")
       call this%io_minmax_aei%set_name(trim(this%nf_minmax_aei)        &
     &      ,64)
       call this%io_minmax_aei%fcreate()
      endif
      
      ! sortie des min/max en a_p1-a_p2, .... entre [-pi,pi] et [0, 2*pi]
      if (this%minmax_diffalp%m_compute.eq.1) then
       this%io_minmax_diffalp = t_io_txt_ncomp                           &
     &      ("I/O min/max a_1-a_2,l_1-l_2,... [-pi,pi] et [0, 2*pi]")
       call this%io_minmax_diffalp%set_name(                             &
     &       trim(this%nf_minmax_diffalp) ,65)
       call this%io_minmax_diffalp%fcreate()
      endif
     
      ! sortie des min/max en a_p1-a_p2, .... angle continu
      if (this%minmax_diffalc%m_compute.eq.1) then
       this%io_minmax_diffalc = t_io_txt_ncomp                           &
     &      ("I/O min/max a_1-a_2,l_1-l_2,... angle continu")
       call this%io_minmax_diffalc%set_name(                             &
     &       trim(this%nf_minmax_diffalc) ,72)
       call this%io_minmax_diffalc%fcreate()
      endif
     
      ! sortie des min/max en (a_p1-a_p2)**2, ....
      if (this%minmax_diffae2%m_compute.eq.1) then
       this%io_minmax_diffae2 = t_io_txt_ncomp                           &
     &      ("I/O min/max (a_1-a_2)**2,(e_1-e_2)**2,...")
       call this%io_minmax_diffae2%set_name(                             &
     &       trim(this%nf_minmax_diffae2) ,71)
       call this%io_minmax_diffae2%fcreate()
      endif

      ! sortie des analyses en frequence en alkhqp
      if (this%naf_alkhqp%m_compute.ne.0) then
       this%io_naf_alkhqp = t_io_txt_ncomp("I/O TXT naf alkhqp")
       call this%io_naf_alkhqp%set_name(trim(this%nf_naf_alkhqp)        &
     &      ,68)
       call this%io_naf_alkhqp%fcreate()
      endif
      
      ! sortie des  analyses en frequence en exp(i*(l1-l2)), ....
      if (this%naf_diffalp%m_compute.eq.1) then
       this%io_naf_diffalp = t_io_txt_ncomp                             &
     &      ("I/O TXT naf exp(i(l1-l2)),exp(i(pi1-pi2))")
       call this%io_naf_diffalp%set_name(                               &
     &       trim(this%nf_naf_diffalp) ,69)
       call this%io_naf_diffalp%fcreate()
      endif

      ! sortie des  analyses en frequence en ((l1-l2),0), ... modulo 2*pi et +-pi
      if (this%naf_difflpm%m_compute.eq.1) then
       this%io_naf_difflpm = t_io_txt_ncomp                             &
     &      ("I/O TXT naf (l1-l2) mod 2*PI ou +-PI")
       call this%io_naf_difflpm%set_name(                               &
     &       trim(this%nf_naf_difflpm) ,70)
       call this%io_naf_difflpm%fcreate()
      endif

      ! sortie du controle de distance a l'etoile
      if (this%ctrl_diststar%m_compute.gt.1) then
       this%io_ctrlstar_dump = t_io_txt_ncomp("I/O TXT ctrl diststar")
       call this%io_ctrlstar_dump%set_name(                             &
     &       trim(this%nf_ctrlstar_dump) ,73)
       call this%io_ctrlstar_dump%fcreate()
      endif

      ! sortie du controle de distance entre planetes
      if (abs(this%ctrl_distpla%m_compute).gt.1) then
       this%io_ctrlpla_dump = t_io_txt_ncomp("I/O TXT ctrl distpla")
       call this%io_ctrlpla_dump%set_name(                              &
     &       trim(this%nf_ctrlpla_dump) ,74)
       call this%io_ctrlpla_dump%fcreate()
      endif

      end  subroutine extfilepar_open_file    

      end module mod_filepar

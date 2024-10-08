!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file naf_base.f 
!!  \brief Calcule l'analyse en frequence de P composantes complexes parmi M composantes complexes
!! sur une tranche de temps donne
!!
!!    --------------        naf         --------------\n
!!    | buffer src |  ----------------> | buffer dst |\n
!!    --------------                    --------------\n
!!
! history : creation 28/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le calcul de l'analyse en frequence
!***********************************************************************
      module mod_naf_base
       use mod_converter
       
!***********************************************************************
!> @class t_naf_base
!! Calcule l'analyse en frequence de P composantes complexes 
!! parmi M composantes complexes sur une tranche de temps donne. \n
!! Le buffer d'entree, contient le temps + M composantes composantes.\n
!! Le buffer de sortie contient le temps + les P*NTERM*(ampl complexe, freq)\n
!! Le stockage dans le buffer de sortie se fait de maniere suivante :\n
!!    temps , (ampl1, freq1, ampl2, freq2, ..., ampln, freqn) de composante_1, (ampl1, freq1,...) de composante_2,  ....
!!
!! La configuration des composantes a analyser se fait par la fonction
!! set_component.
!!
!! P sera nomme : m_nbcalccomp
!! M sera nomme : m_nbinputcomp
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_naf_base
          private

          integer :: m_nbcalccomp = 0 !< nombre de composantes complexes dont on calcule l'analyse en frequence
          integer, public :: m_nbinputcomp = 0 !< nombre de composantes dont on extrait les m_nbcalccomp
          !> tableau de "m_nbcalccomp" indices indiquant quelles composantes sont extraites parmi les "m_nbinputcomp" composantes d'entree
          integer, dimension(:), allocatable :: m_indexcomp 
          complex(TREAL),dimension(:,:),allocatable,public :: m_zcomp !< tableau des valeurs d'entree des donnees (1:m_nbcalccomp,1:m_stepout)
          real(TREAL), dimension(:), allocatable, public :: m_tcomp !< tableau des temps d'entree (1:m_stepout)
          integer(8), public :: m_stepout !< nombre de pas pour calculer l'analyse en frequence et ecrire le resultat
          integer(8), public :: m_counter = 0 !< compteur de pas accumule
          integer, dimension(:), allocatable :: m_indexmce !< tableau de  "m_nbcalccomp" des indices des arrets assoccies aux composantes analysees
          

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

          ! resultats
          real(TREAL), dimension(:,:), allocatable :: m_freq !< tableau des frequences trouvees (1:m_nbcalccomp,1:m_naf_nterm)
          complex(TREAL), dimension(:,:), allocatable :: m_zamp !< tableau des amplitudes trouvees (1:m_nbcalccomp,1:m_naf_nterm)
          integer, dimension(:), allocatable :: m_nbfreq !< tableau du nombre de frequences trouvees (1:m_nbcalccomp)


      contains
          procedure :: set_param_naf => naf_base_set_param_naf ! fixe les parametres de naf
          procedure :: set_output => naf_base_set_output ! fixe le buffer de sortie et le consommateur
          procedure :: set_stepout => naf_base_set_stepout ! fixe la longueur d'une tranche pour les calculs (en nombre de pas d'entree)
          procedure :: set_component => naf_base_set_component ! definit les composantes a analyser
          procedure :: set_multictrlerr => naf_base_set_multictrlerr ! associe aux composantes a analyser
          
          procedure :: oninputbufferfull => naf_base_oninputbufferfull ! procedure appellee lorsque le buffer d'entree est plein

      end type t_naf_base  


      private naf_base_resetval
      
      contains

!***********************************************************************
!> @brief fixe la taille du buffer de sortie et son consommateur
!***********************************************************************
         subroutine naf_base_set_param_naf(this, nterm, iw, isec,dtour,    &
     &                     tol)
          implicit none
          class(t_naf_base), intent(inout):: this  !< dummy argument
          integer, intent(in) :: nterm !< nombre de termes demandes
         !> presence de fenetre
         !! *= -1 : fenetre exponentielle PHI(T) = 1/CE*EXP(-1/(1-T^2)) avec CE= 0.22199690808403971891E0
         !! *= 0 : pas de fenetre
         !! *= N > 0 : PHI(T) = CN*(1+COS(PI*T))**N avec CN = 2^N(N!)^2/(2N)!
          integer, intent(in) :: iw
         !> *= 0 : la methode des secantes n'est pas utilisee.
         !! *= 1 : la methode des secantes est utilisee.
          integer, intent(in) :: isec
          real(TREAL), intent(in) ::   dtour !< "longueur d'un tour de cadran"
          real(TREAL), intent(in) ::   tol !< tolerance pour determiner si deux frequences sont identiques
         
          this%m_naf_dtour = dtour
          this%m_naf_nterm = nterm
          this%m_naf_iw = iw
          this%m_naf_isec = isec
          this%m_naf_tol = tol
         end  subroutine naf_base_set_param_naf

!***********************************************************************
!> @brief fixe la taille du buffer de sortie et son consommateur
!***********************************************************************
         subroutine naf_base_set_output(this,buffersize,buffercs)
          implicit none
          class(t_naf_base), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%setconsumer("NAF")
          call this%m_buffer%init(buffersize,                              &
     &              1+this%m_naf_nterm*this%m_nbcalccomp*3, buffercs) 
         end  subroutine naf_base_set_output

!***********************************************************************
!> @brief definit les composantes a analyser
!! a partir du tableau d'entiers fournis(indice par rapport au buffer d'entree).\n
!! Remarque : les entiers sont la position en composantes reelles du debut 
!! de chaque composante complexe.
!***********************************************************************
         subroutine naf_base_set_component(this,arPcomp, Mcomp)
          implicit none
          class(t_naf_base), intent(inout):: this  !< dummy argument
          integer, dimension(:), intent(in) :: arPcomp  !< tableau d'indices des composantes a analyser (commence a 1)
          integer, intent(in) :: Mcomp  !< nombre de composantes du buffer d'entree
          
          this%m_nbcalccomp = size(arPcomp)
          this%m_nbinputcomp = Mcomp
          allocate(this%m_indexcomp(1:this%m_nbcalccomp))
         
          this%m_indexcomp(:) = arPcomp(1:this%m_nbcalccomp)
          call naf_base_resetval(this)
          this%m_counter = 0 
         
         end  subroutine naf_base_set_component

!***********************************************************************
!> @brief associe aux composantes a analyser les indices dans le tableau des arrets
!! arindcomp doit avoir autant d'elements que arPcomp de la fonction precedente
!***********************************************************************
         subroutine naf_base_set_multictrlerr(this,arindcomp)
          implicit none
          class(t_naf_base), intent(inout):: this  !< dummy argument
          integer, dimension(:), intent(in) :: arindcomp  !< tableau des indices des arrets associes aux composantes a analyser 
          
          allocate(this%m_indexmce(1:this%m_nbcalccomp))
          this%m_indexmce(:)=arindcomp(:)
         
         end  subroutine naf_base_set_multictrlerr


!***********************************************************************
!> @brief fixe la longueur d'une tranche d'analyses (en nombre de pas d'entree)
!***********************************************************************
         subroutine naf_base_set_stepout(this,stepout)
          implicit none
          class(t_naf_base), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: stepout   !< nombre de pas de la tranche (en pas d'entree)
          
          this%m_stepout =  stepout
         end  subroutine naf_base_set_stepout

!***********************************************************************
!> @brief reinitialise l'analyse pour commencer la tanche
!***********************************************************************
         subroutine naf_base_resetval(this)
          implicit none
          class(t_naf_base), intent(inout):: this  !< dummy argument
          
          this%m_counter = 0
          
         end  subroutine naf_base_resetval
     
!***********************************************************************
!> @brief initialise les tableaux d'entree et sortie pour les analyses en frequences
!***********************************************************************
         subroutine naf_base_init_zcomp(this) 
          implicit none
          class(t_naf_base), intent(inout) :: this   !< donnee utilisateur 
          
          integer nbcalccomp
          
          nbcalccomp = this%m_nbcalccomp
          
          if (.not.allocated(this%m_zcomp)) then
           allocate(this%m_zcomp(1:nbcalccomp,1:this%m_stepout))
           allocate(this%m_tcomp(1:this%m_stepout))
           allocate(this%m_freq(1:nbcalccomp,1:this%m_naf_nterm))
           allocate(this%m_zamp(1:nbcalccomp,1:this%m_naf_nterm))
           allocate(this%m_nbfreq(1:nbcalccomp))
          endif
        end subroutine naf_base_init_zcomp
        
!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree est plein
!! stocke les composantes dans le tableau intermediaire
!! des composantes selectionnees
!***********************************************************************
         subroutine naf_base_oninputbufferfull(this, buffer, poswrite) 
          implicit none
          class(t_naf_base), intent(inout) :: this   !< donnee utilisateur 
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL), dimension(1+this%m_nbinputcomp) :: R
          real(TREAL) t, v1,v2
          integer k
          complex(TREAL) v
          
          integer nbcalccomp
          
          nbcalccomp = this%m_nbcalccomp
          
          call naf_base_init_zcomp(this)
          
          call this%m_buffer%tryset_multictrlerr_buf(buffer)
          
            do j=1, poswrite
             ! lecture de l'entree
             call buffer%readdata(j, R)            
             t = R(1)
             this%m_tcomp(this%m_counter+1) = t
             ! remplissage de m_zcomp
             do k=1, this%m_nbcalccomp
              v1 = R(1+this%m_indexcomp(k))
              v2 = R(1+this%m_indexcomp(k)+1)
              v = CMPLX(v1,v2, kind(REALC(1.0)))
              this%m_zcomp(k, this%m_counter+1) = v
             enddo
             
             this%m_counter = this%m_counter+1
             if (this%m_counter==this%m_stepout) then
              ! fin de la tranche 
              call naf_base_analyzeall(this, buffer)
             endif
            enddo
          
          call buffer%tryset_multictrlerr_buf(this%m_buffer)
          call buffer%empty()
         
         end subroutine naf_base_oninputbufferfull

!***********************************************************************
!> @brief fonction qui realise les analyses en frequence pour toutes les composantes et stocke le resultat dans le buffer de sortie
!! Elle suppose que les tableaux tcomp/zcomp sont completement remplies.
!! 
!***********************************************************************
         subroutine naf_base_analyzeall(this, buffer) 
          implicit none
          class(t_naf_base), intent(inout) :: this   !< donnee utilisateur 
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie

          real(TREAL), dimension(:), allocatable :: W
          real(TREAL) zr, zi
          integer nbcalccomp
          integer k, l
          
          nbcalccomp = this%m_nbcalccomp
          
           if (allocated(this%m_indexmce)) then
              ! il y a un controle des erreurs par corps
              do k=1, this%m_nbcalccomp
               if (buffer%m_mce(this%m_indexmce(k))%m_stop.eq.0) then
                call naf_base_analyze(this, k)
               endif
              enddo
           else
              ! => effectuer les calculs
              ! pas de controle des erreurs par corps
              do k=1, this%m_nbcalccomp
               call naf_base_analyze(this, k)
              enddo
          endif
          
          ! envoi des donnees au buffer de sortie
          allocate(W(1:1+this%m_naf_nterm*nbcalccomp*3))
          W(1) = this%m_tcomp(1)
          do k=1, this%m_nbcalccomp
           do l=1, this%m_naf_nterm
           W(1+((k-1)*this%m_naf_nterm+(l-1))*3+1)=this%m_freq(k,l)
            zr = real(this%m_zamp(k,l))
            zi = imag(this%m_zamp(k,l))
            W(1+((k-1)*this%m_naf_nterm+(l-1))*3+2)=zr
            W(1+((k-1)*this%m_naf_nterm+(l-1))*3+3)=zi
           enddo
          enddo
          call this%m_buffer%writedata(W)   
          deallocate(W)
          
          ! initialisation de la prochaine tranche
          call naf_base_resetval(this)

         end subroutine naf_base_analyzeall
         
!***********************************************************************
!> @brief fonction qui realise l'analyse en frequence pour la composante m_zcomp(kcomp)
!! 
!***********************************************************************
         subroutine naf_base_analyze(this,kcomp) 
          use NAFF
          implicit none
          class(t_naf_base), intent(inout) :: this   !< ensemble des donnees pour l'analyse en frequence 
          integer, intent(in) :: kcomp !< composante de m_zcomp a traiter
          
          integer j
          real(TREAL) nafeps
          real(TREAL) nofr
          
! INITIALISER LES VARIABLES :
           NAF_DTOUR = this%m_naf_dtour
           NAF_KTABS = ((this%m_stepout-1)/6)*6
           NAF_XH = this%m_tcomp(2)-this%m_tcomp(1)
           NAF_NTERM = this%m_naf_nterm
           NAF_IW = this%m_naf_iw
           NAF_T0 = this%m_tcomp(1)
           NAF_ICPLX = 1
           NAF_ISEC = this%m_naf_isec
           NAF_NFPRT = 6
           NAF_IPRT = -1
           NAF_TOL = this%m_naf_tol
           
!       CALL INITNAF
           call initnaf
!      1)  REMPLIR LE TABLEAU NAF_ZTABS
             do j=0, NAF_KTABS
             NAF_ZTABS(j) = this%m_zcomp(kcomp, j+1)
             enddo
!      2)  CALL MFTNAF(NBTERM,EPS)
           nafeps = ABSR(NAF_FREFON)/REALC(1E100)
           call mftnaf(NAF_NTERM, real(nafeps, kind=8)) !TREAL))
            
           !write(*,*) NAF_KTABS, NAF_NTERM, NAF_NFS 
           !write(*,*) (NAF_TFS(j),j=1,NAF_NFS)
           !write(*,*) (NAF_ZAMP(j),j=1,NAF_NFS)
           this%m_nbfreq(kcomp) = NAF_NFS
           do j=1, NAF_NFS
            this%m_freq(kcomp, j) = NAF_TFS(j)
            this%m_zamp(kcomp, j) = NAF_ZAMP(j)
           enddo
           nofr = REALC(1E98)
           do j=NAF_NFS+1, NAF_NTERM 
            this%m_freq(kcomp, j) = nofr
            this%m_zamp(kcomp, j) = CMPLX(nofr, nofr, kind(nofr))
           enddo
           
!      3) CALL CLEANNAF
           call cleannaf

         end subroutine naf_base_analyze

      end module
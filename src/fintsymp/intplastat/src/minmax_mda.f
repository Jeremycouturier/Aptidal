!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file minmax_mda.f 
!!  \brief Calcule le min/moy/max des a_p(1)-a_p(2),  la_p(1)-la_p(2), 
!! et pi_p(1)-pi_p(2) sur une tranche de temps donne avec p(1) et p(2) les indices de 2 planetes
!!
!!
!! Une double determination est faite pour les angles : entre [0, 2PI] et [-PI,PI].
!!
!!    --------------   min/moy/max      --------------\n
!!    | buffer src |  ----------------> | buffer dst |\n
!!    --------------                    --------------\n
!!
! history : creation 25/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le calcul des min/moy/max avec des angles entre [0, 2PI] et [-PI,PI]
!***********************************************************************
      module mod_minmax_mda
       use mod_minmax
       
              
!***********************************************************************
!> @class t_minmax_mda
!! Calcule le min/moy/max des a_p(1)-a(p2),  la_p(1)-la(p2), 
!! et pi_p(1)-pi_(p2) sur une tranche de temps donne avec p(1) et p(2) les indices de 2 planetes
!! Le buffer d'entree, contient le temps + M composantes.\n
!! les M composants doivent etre des elements elliptiques dont les 3 dernieres composantes sont des angles.
!! Le buffer de sortie contient le temps + les ???*P min/moy/max.\n
!! Le stockage dans le buffer de sortie se fait de maniere suivante :\n
!!    temps , min/moy/max de composante_1, 
!!    min/moy/max de composante_2 entre [0,2*PI],  
!!    min/moy/max de composante_2 entre [-PI,PI],
!!    min/moy/max de composante_3 entre [0,2*PI],  
!!    min/moy/max de composante_3 entre [-PI,PI]
!!
!! La configuration des composantes a analyser se fait par la fonction
!! set_component.
!!
!! P sera nomme : m_nbcalccomp
!! M sera nomme : m_nbinputcomp
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_minmax_mda
          private

          integer :: m_nbcalccomp = 0 !< nombre de composantes dont on calcule le min/moy/max
          integer :: m_nbstorecomp = 0 !< nombre de composantes dont on stocke le min/moy/max (car 2 valeurs pour les angles)
          integer :: m_nbinputcomp = 0 !< nombre de composantes dont on extrait les m_nbcalccomp
          !> tableau de "m_nbcalccomp" indices indiquant quelles composantes sont extraites parmi les "m_nbinputcomp" composantes d'entree 
          !! pour la premiere planete
          integer, dimension(:), allocatable :: m_indexcomp1 
          !> tableau de "m_nbcalccomp" indices indiquant quelles composantes sont extraites parmi les "m_nbinputcomp" composantes d'entree 
          !! pour la deuxieme planete
          integer, dimension(:), allocatable :: m_indexcomp2 
          !> tableau de "m_nbcalccomp" booleens indiquant si la composante extraite est un angle ou pas
          logical, dimension(:), allocatable :: m_anglecomp 
          !> tableau de "m_nbcalccomp" indices indiquant ou on stocke les donnees dans le buffer de sortie
          integer, dimension(:), allocatable :: m_storecomp 
          real(TREAL), dimension(:), allocatable :: m_minval !< valeur minimale courante des "m_nbstorecomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_maxval !< valeur maximale courante des "m_nbstorecomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_sumval !< somme courante des "m_nbstorecomp" composantes 
          integer(8) :: m_stepout !< nombre de pas pour calculer le min/max et ecrire le resultat
          integer(8) :: m_counter = 0 !< compteur de pas accumule
          real(TREAL) :: m_pi !< valeur de pi

      contains
          procedure :: set_output => minmaxmda_set_output ! fixe le buffer de sortie et le consommateur
          procedure :: set_component => minmaxmda_set_component ! definit les composantes a analyser
          procedure :: set_stepout => minmaxmda_set_stepout ! fixe la longueur d'une tranche pour les calculs (en nombre de pas d'entree)
          
          procedure :: oninputbufferfull => minmaxmda_oninputbufferfull ! procedure appellee lorsque le buffer d'entree est plein

      end type t_minmax_mda  


      private minmaxmda_resetval
      
      contains

!***********************************************************************
!> @brief fixe la taille du buffer de sortie et son consommateur
!***********************************************************************
         subroutine minmaxmda_set_output(this,buffersize,buffercs)
          implicit none
          class(t_minmax_mda), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%setconsumer("MIN/MAX")
          call this%m_buffer%init(buffersize, 1+3*this%m_nbstorecomp,       &
     &              buffercs) 
         end  subroutine minmaxmda_set_output

!***********************************************************************
!> @brief definit les composantes a analyser
!! a partir du tableau d'entiers fournis(indice par rapport au buffer d'entree)
!***********************************************************************
         subroutine minmaxmda_set_component(this,arPcomp1, arPcomp2,        &
     &                       anglecomp, Mcomp)
          implicit none
          class(t_minmax_mda), intent(inout):: this  !< dummy argument
          integer, dimension(:), intent(in) :: arPcomp1  !< tableau d'indices des composantes a analyser pour la premiere planete
          integer, dimension(:), intent(in) :: arPcomp2  !< tableau d'indices des composantes a analyser pour la deuxieme planete
          logical, dimension(:), intent(in) :: anglecomp  !< tableau de booleen =true => c'est un angle
          integer, intent(in) :: Mcomp  !< nombre de composantes du buffer d'entree
          
          integer k,s
          
          this%m_nbcalccomp =  size(arPcomp1)
          this%m_nbinputcomp =  Mcomp
          allocate(this%m_indexcomp1(1:this%m_nbcalccomp))
          allocate(this%m_indexcomp2(1:this%m_nbcalccomp))
          allocate(this%m_anglecomp(1:this%m_nbcalccomp))
          allocate(this%m_storecomp(1:this%m_nbcalccomp))
         
          this%m_indexcomp1(:) = arPcomp1(1:this%m_nbcalccomp)
          this%m_indexcomp2(:) = arPcomp2(1:this%m_nbcalccomp)
          this%m_anglecomp(:) = anglecomp(1:this%m_nbcalccomp)
          !initialise storecomp
          s=1
          do k=1, this%m_nbcalccomp
           this%m_storecomp(k) = s
           s=s+1
           if (anglecomp(k)) s=s+1
          enddo
          
          this%m_nbstorecomp = s-1
          allocate(this%m_minval(1:this%m_nbstorecomp))
          allocate(this%m_maxval(1:this%m_nbstorecomp))
          allocate(this%m_sumval(1:this%m_nbstorecomp))

          call minmaxmda_resetval(this)
          this%m_counter = -1 ! pour forcer la 1ere valeur a afficher
          
          this%m_pi = atan2(REALC(0.), REALC(-1.))
          !write(*,*) this%m_pi
          !write(*,*) "storecomp", this%m_storecomp
          !write(*,*) "m_indexcomp1", this%m_indexcomp1
          !write(*,*) "m_indexcomp2", this%m_indexcomp2
          !write(*,*) "m_anglecomp", this%m_anglecomp
           
         end  subroutine minmaxmda_set_component

!***********************************************************************
!> @brief fixe la longueur d'une tranche d'analyses (en nombre de pas d'entree)
!***********************************************************************
         subroutine minmaxmda_set_stepout(this,stepout)
          implicit none
          class(t_minmax_mda), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: stepout   !< nombre de pas de la tranche (en pas d'entree)
          
          this%m_stepout =  stepout
         end  subroutine minmaxmda_set_stepout

!***********************************************************************
!> @brief reinitialise les min/sum/max pour commencer la tranche
!***********************************************************************
         subroutine minmaxmda_resetval(this)
          implicit none
          class(t_minmax_mda), intent(inout):: this  !< dummy argument
          
          integer j

          do j=1, this%m_nbstorecomp
            this%m_minval(j) = 1D300
            this%m_maxval(j) = -1D300
            this%m_sumval(j) = 0
          enddo
          this%m_counter = 0
          
         end  subroutine minmaxmda_resetval
     
!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree est plein
!! calcule la somme (pour la moyenne), le minium et le maximum 
!! des composantes selectionnees
!***********************************************************************
       subroutine minmaxmda_oninputbufferfull(this, buffer, poswrite) 
        implicit none
        class(t_minmax_mda), intent(inout) :: this   !< donnee utilisateur 
        class(t_buffer), intent(inout) :: buffer !< buffer de sortie
        integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
        integer(8) :: j
        real(TREAL), dimension(1+this%m_nbinputcomp) :: R
        real(TREAL), dimension(1+3*this%m_nbstorecomp) :: W
        real(TREAL) t, v1, v2, v, v_0_2pi, v_npi_pi
        integer k,s
        
        
         do j=1, poswrite
           ! lecteure de l'entree
           call buffer%readdata(j, R)            
           t = R(1)
           ! calcule le min/max/sum
           do k=1, this%m_nbcalccomp
            v1 = R(1+this%m_indexcomp1(k))
            v2 = R(1+this%m_indexcomp2(k))
            v = v1-v2
            
            s = this%m_storecomp(k)
            if (this%m_anglecomp(k)) then
            ! cas d'un angle => double determination
             ! entre 0 et 2*pi
             v_0_2pi = modulo(v, 2*this%m_pi)
             if (v_0_2pi<REALC(0.)) v_0_2pi=v_0_2pi+2*this%m_pi
             if (v_0_2pi<this%m_minval(s)) this%m_minval(s) = v_0_2pi
             if (v_0_2pi>this%m_maxval(s)) this%m_maxval(s) = v_0_2pi
             this%m_sumval(s) = this%m_sumval(s)+v_0_2pi
             s = s+1
             ! entre -pi et pi
             v_npi_pi = mod(v, -2*this%m_pi)
             if (v_npi_pi>this%m_pi) v_npi_pi = v_npi_pi-2*this%m_pi
             if (v_npi_pi<-this%m_pi) v_npi_pi = v_npi_pi+2*this%m_pi
             if (v_npi_pi<this%m_minval(s)) this%m_minval(s)=v_npi_pi
             if (v_npi_pi>this%m_maxval(s)) this%m_maxval(s)=v_npi_pi
             this%m_sumval(s) = this%m_sumval(s)+v_npi_pi
            else
             ! cas d'une variable ordinaire => une seule determination
             if (v<this%m_minval(s)) this%m_minval(s) = v
             if (v>this%m_maxval(s)) this%m_maxval(s) = v
             this%m_sumval(s) = this%m_sumval(s)+v
            endif 
           enddo
           
           this%m_counter = this%m_counter+1
           if (this%m_counter==this%m_stepout) then
            ! fin de la tranche 
            ! la moyenne se fait sur m_counter+1 valeurs (lie a la presence pour le temps t=0)
            W(1) = t
            do k=1, this%m_nbstorecomp
             !write(*,*) "w", k, this%m_minval(k), this%m_maxval(k) 
             W(1+3*(k-1)+1) = this%m_minval(k) 
             W(1+3*(k-1)+3) = this%m_maxval(k) 
             W(1+3*(k-1)+2) = this%m_sumval(k)/(this%m_counter+1) 
            enddo
            
            call this%m_buffer%writedata(W)   
            
            ! reinitialise la tranche
            call minmaxmda_resetval(this)        
            ! intialise le min/max/sum avec le dernier point
            do k=1, this%m_nbcalccomp
                v1 = R(1+this%m_indexcomp1(k))
                v2 = R(1+this%m_indexcomp2(k))
                v = v1-v2
                s = this%m_storecomp(k)
                if (this%m_anglecomp(k)) then
                 ! cas d'un angle => double determination
                 ! entre 0 et 2*pi
                 v_0_2pi = modulo(v, 2*this%m_pi)
                 if (v_0_2pi<REALC(0.)) v_0_2pi=v_0_2pi+2*this%m_pi
                 this%m_minval(s) = v_0_2pi
                 this%m_maxval(s) = v_0_2pi
                 this%m_sumval(s) = v_0_2pi
                 s = s+1
                 ! entre -pi et pi
                 v_npi_pi = modulo(v, -2*this%m_pi)
                 if (v_npi_pi>this%m_pi) v_npi_pi=v_npi_pi-2*this%m_pi
                 if (v_npi_pi<-this%m_pi) v_npi_pi=v_npi_pi+2*this%m_pi
                 this%m_minval(s) = v_npi_pi
                 this%m_maxval(s) = v_npi_pi
                 this%m_sumval(s) = v_npi_pi
                else
                 ! cas d'une variable ordinaire => une seule determination
                 this%m_minval(s) = v
                 this%m_maxval(s) = v
                 this%m_sumval(s) = v
                endif 
            enddo

           endif
          enddo
        
        call buffer%empty()
       
       end subroutine minmaxmda_oninputbufferfull


      end module
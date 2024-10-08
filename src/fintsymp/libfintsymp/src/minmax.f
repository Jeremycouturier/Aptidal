!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file minmax.f 
!!  \brief Calcule le min/moy/max de P composantes parmi M composantes
!! sur une tranche de temps donne
!!
!!    --------------   min/moy/max      --------------\n
!!    | buffer src |  ----------------> | buffer dst |\n
!!    --------------                    --------------\n
!!
! history : creation 21/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le calcul des min/moy/max
!***********************************************************************
      module mod_minmax
       use mod_converter
       
              
!***********************************************************************
!> @class t_minmax
!! Calcule le min, moy et max de P composantes parmi M composantes
!! sur une tranche de temps donne. 
!! Le buffer d'entree, contient le temps + M composantes.\n
!! Le buffer de sortie contient le temps + les 3*P min/moy/max.\n
!! Le stockage dans le buffer de sortie se fait de maniere suivante :\n
!!    temps , min/moy/max de composante_1, min/moy/max de composante_2,  ....
!!
!! La configuration des composantes a analyser se fait par la fonction
!! set_component.
!!
!! P sera nomme : m_nbcalccomp
!! M sera nomme : m_nbinputcomp
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_minmax
          private

          integer :: m_nbcalccomp = 0 !< nombre de composantes dont on calcule le min/moy/max
          integer :: m_nbinputcomp = 0 !< nombre de composantes dont on extrait les m_nbcalccomp
          !> tableau de "m_nbcalccomp" indices indiquant quelles composantes sont extraites parmi les "m_nbinputcomp" composantes d'entree
          integer, dimension(:), allocatable :: m_indexcomp 
          real(TREAL), dimension(:), allocatable :: m_minval !< valeur minimale courante des "m_nbcalccomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_maxval !< valeur maximale courante des "m_nbcalccomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_sumval !< somme courante des "m_nbcalccomp" composantes 
          integer(8) :: m_stepout !< nombre de pas pour calculer le min/max et ecrire le resultat
          integer(8) :: m_counter = 0 !< compteur de pas accumule

      contains
          procedure :: set_output => minmax_set_output ! fixe le buffer de sortie et le consommateur
          procedure :: set_component => minmax_set_component ! definit les composantes a analyser
          procedure :: set_stepout => minmax_set_stepout ! fixe la longueur d'une tranche pour les calculs (en nombre de pas d'entree)
          
          procedure :: oninputbufferfull => minmax_oninputbufferfull ! procedure appellee lorsque le buffer d'entree est plein

      end type t_minmax  


      private minmax_resetval
      
      contains

!***********************************************************************
!> @brief fixe la taille du buffer de sortie et son consommateur
!***********************************************************************
         subroutine minmax_set_output(this,buffersize,buffercs)
          implicit none
          class(t_minmax), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%setconsumer("MIN/MAX")
          call this%m_buffer%init(buffersize, 1+3*this%m_nbcalccomp,        &
     &              buffercs) 
         end  subroutine minmax_set_output

!***********************************************************************
!> @brief definit les composantes a analyser
!! a partir du tableau d'entiers fournis(indice par rapport au buffer d'entree)
!***********************************************************************
         subroutine minmax_set_component(this,arPcomp, Mcomp)
          implicit none
          class(t_minmax), intent(inout):: this  !< dummy argument
          integer, dimension(:), intent(in) :: arPcomp  !< tableau d'indices des composantes a analyser (commence a 1)
          integer, intent(in) :: Mcomp  !< nombre de composantes du buffer d'entree
          
          
          this%m_nbcalccomp =  size(arPcomp)
          this%m_nbinputcomp =  Mcomp
          allocate(this%m_indexcomp(1:this%m_nbcalccomp))
          allocate(this%m_minval(1:this%m_nbcalccomp))
          allocate(this%m_maxval(1:this%m_nbcalccomp))
          allocate(this%m_sumval(1:this%m_nbcalccomp))
         
          this%m_indexcomp(:) = arPcomp(1:this%m_nbcalccomp)
          call minmax_resetval(this)
          this%m_counter = -1 ! pour forcer la 1ere valeur a afficher
         
         end  subroutine minmax_set_component

!***********************************************************************
!> @brief fixe la longueur d'une tranche d'analyses (en nombre de pas d'entree)
!***********************************************************************
         subroutine minmax_set_stepout(this,stepout)
          implicit none
          class(t_minmax), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: stepout   !< nombre de pas de la tranche (en pas d'entree)
          
          this%m_stepout =  stepout
         end  subroutine minmax_set_stepout

!***********************************************************************
!> @brief reinitialise les min/sum/max pour commencer la tanche
!***********************************************************************
         subroutine minmax_resetval(this)
          implicit none
          class(t_minmax), intent(inout):: this  !< dummy argument
          
          integer j

          do j=1, this%m_nbcalccomp
            this%m_minval(j) = 1D300
            this%m_maxval(j) = -1D300
            this%m_sumval(j) = 0
          enddo
          this%m_counter = 0
          
         end  subroutine minmax_resetval
     
!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree est plein
!! calcule la somme (pour la moyenne), le minium et le maximum 
!! des composantes selectionnees
!***********************************************************************
         subroutine minmax_oninputbufferfull(this, buffer, poswrite) 
          implicit none
          class(t_minmax), intent(inout) :: this   !< donnee utilisateur 
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL), dimension(1+this%m_nbinputcomp) :: R
          real(TREAL), dimension(1+3*this%m_nbcalccomp) :: W
          real(TREAL) t, v
          integer k
          
            call this%m_buffer%tryset_multictrlerr_buf(buffer)

            do j=1, poswrite
             ! lecteure de l'entree
             call buffer%readdata(j, R)            
             t = R(1)
             ! calcule le min/max/sum
             do k=1, this%m_nbcalccomp
              v = R(1+this%m_indexcomp(k))
              if (v<this%m_minval(k)) this%m_minval(k) = v
              if (v>this%m_maxval(k)) this%m_maxval(k) = v
              this%m_sumval(k) = this%m_sumval(k)+v
             enddo
             
             this%m_counter = this%m_counter+1
             if (this%m_counter==this%m_stepout) then
              ! fin de la tranche 
              ! la moyenne se fait sur m_counter+1 valeurs (lie a la presence pour le temps t=0)
              W(1) = t
              do k=1, this%m_nbcalccomp
               W(1+3*(k-1)+1) = this%m_minval(k) 
               W(1+3*(k-1)+3) = this%m_maxval(k) 
               W(1+3*(k-1)+2) = this%m_sumval(k)/(this%m_counter+1) 
              enddo
              call this%m_buffer%writedata(W)   
              ! reinitialise la tranche
              call minmax_resetval(this)        
              ! intialise le min/max/sum avec le dernier point
              do k=1, this%m_nbcalccomp
               v = R(1+this%m_indexcomp(k))
               this%m_minval(k) = v
               this%m_maxval(k) = v
               this%m_sumval(k) = v
              enddo

             endif
            enddo
            
            call buffer%tryset_multictrlerr_buf(this%m_buffer)
          
          call buffer%empty()
         
         end subroutine minmax_oninputbufferfull


      end module
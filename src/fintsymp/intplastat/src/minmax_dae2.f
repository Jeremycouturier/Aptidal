!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file minmax_dae2.f 
!!  \brief Calcule le min/moy/max des (a_p(1)-a_p(2))**2, (e_p(1)-e_p(2))**2, 
!! et (a_p(1)-a_p(2))**2+(e_p(1)-e_p(2))**2 sur une tranche de temps donne 
!! avec p(1) et p(2) les indices de 2 planetes
!!
!!
!!
!!    --------------   min/moy/max      --------------\n
!!    | buffer src |  ----------------> | buffer dst |\n
!!    --------------                    --------------\n
!!
! history : creation 24/03/2016
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le calcul des min/moy/max
!***********************************************************************
      module mod_minmax_dae2
       use mod_minmax
       
              
!***********************************************************************
!> @class t_minmax_dae2
!! Calcule le min/moy/max des (a_p(1)-a_p(2))**2, (e_p(1)-e_p(2))**2, 
!! et (a_p(1)-a_p(2))**2+(e_p(1)-e_p(2))**2 sur une tranche de temps donne avec p(1) et p(2) les indices de 2 planetes
!! Le buffer d'entree, contient le temps + M composantes.\n
!! les M composants doivent etre des elements elliptiques.
!! Le buffer de sortie contient le temps + les ???*P min/moy/max.\n
!! Le stockage dans le buffer de sortie se fait de maniere suivante :\n
!!    temps , min/moy/max de composante_1, 
!!    min/moy/max de composante_2,  
!!    min/moy/max de (composante_1+composante_2)
!!
!! La configuration des composantes a analyser se fait par la fonction
!! set_component.
!!
!! P sera nomme : m_nbcalccomp
!! M sera nomme : m_nbinputcomp
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_minmax_dae2
          private

          integer :: m_nbcalccomp = 0 !< nombre de composantes dont on calcule le min/moy/max
          integer :: m_nbstorecomp = 0 !< nombre de composantes dont on stocke le min/moy/max
          integer :: m_nbinputcomp = 0 !< nombre de composantes dont on extrait les m_nbcalccomp
          !> tableau de "m_nbcalccomp" indices indiquant quelles composantes sont extraites parmi les "m_nbinputcomp" composantes d'entree 
          !! pour la premiere planete
          integer, dimension(:), allocatable :: m_indexcomp1 
          !> tableau de "m_nbcalccomp" indices indiquant quelles composantes sont extraites parmi les "m_nbinputcomp" composantes d'entree 
          !! pour la deuxieme planete
          integer, dimension(:), allocatable :: m_indexcomp2 
          !> tableau de "m_nbstorecomp" indices indiquant ou on stocke les donnees dans le buffer de sortie
          integer, dimension(:), allocatable :: m_storecomp 
          real(TREAL), dimension(:), allocatable :: m_minval !< valeur minimale courante des "m_nbstorecomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_maxval !< valeur maximale courante des "m_nbstorecomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_sumval !< somme courante des "m_nbstorecomp" composantes 
          integer(8) :: m_stepout !< nombre de pas pour calculer le min/max et ecrire le resultat
          integer(8) :: m_counter = 0 !< compteur de pas accumule

      contains
          procedure :: set_output => minmaxdae2_set_output ! fixe le buffer de sortie et le consommateur
          procedure :: set_component => minmaxdae2_set_component ! definit les composantes a analyser
          procedure :: set_stepout => minmaxdae2_set_stepout ! fixe la longueur d'une tranche pour les calculs (en nombre de pas d'entree)
          
          procedure:: oninputbufferfull => minmaxdae2_oninputbufferfull ! procedure appellee lorsque le buffer d'entree est plein

      end type t_minmax_dae2  


      private minmaxdae2_resetval
      
      contains

!***********************************************************************
!> @brief fixe la taille du buffer de sortie et son consommateur
!***********************************************************************
         subroutine minmaxdae2_set_output(this,buffersize,buffercs)
          implicit none
          class(t_minmax_dae2), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%setconsumer("MIN/MAX DAE2")
          call this%m_buffer%init(buffersize, 1+3*this%m_nbstorecomp,       &
     &              buffercs) 
         end  subroutine minmaxdae2_set_output

!***********************************************************************
!> @brief definit les composantes a analyser
!! a partir du tableau d'entiers fournis(indice par rapport au buffer d'entree)
!***********************************************************************
         subroutine minmaxdae2_set_component(this,arPcomp1, arPcomp2,        &
     &                       Mcomp)
          implicit none
          class(t_minmax_dae2), intent(inout):: this  !< dummy argument
          integer, dimension(:), intent(in) :: arPcomp1  !< tableau d'indices des composantes a analyser pour la premiere planete
          integer, dimension(:), intent(in) :: arPcomp2  !< tableau d'indices des composantes a analyser pour la deuxieme planete
          integer, intent(in) :: Mcomp  !< nombre de composantes du buffer d'entree
          
          integer k
          
          this%m_nbcalccomp =  size(arPcomp1)
          this%m_nbstorecomp = this%m_nbcalccomp+1 ! +1 pour le (a_p(1)-a_p(2))**2+(e_p(1)-e_p(2))**2
          if (this%m_nbcalccomp.ne.2) then 
           write(*,*) 'minmaxdae2_set_component requiert 2 composantes'
           stop
          endif
          this%m_nbinputcomp =  Mcomp
          allocate(this%m_indexcomp1(1:this%m_nbcalccomp))
          allocate(this%m_indexcomp2(1:this%m_nbcalccomp))
          allocate(this%m_storecomp(1:this%m_nbstorecomp))
         
          this%m_indexcomp1(:) = arPcomp1(1:this%m_nbcalccomp)
          this%m_indexcomp2(:) = arPcomp2(1:this%m_nbcalccomp)
          !initialise storecomp
          do k=1, this%m_nbcalccomp
           this%m_storecomp(k) = k
          enddo
          ! pour le (a_p(1)-a_p(2))**2+(e_p(1)-e_p(2))**2
          this%m_storecomp(this%m_nbcalccomp+1) = this%m_nbcalccomp+1
          
          allocate(this%m_minval(1:this%m_nbstorecomp))
          allocate(this%m_maxval(1:this%m_nbstorecomp))
          allocate(this%m_sumval(1:this%m_nbstorecomp))

          call minmaxdae2_resetval(this)
          this%m_counter = -1 ! pour forcer la 1ere valeur a afficher
          
          !write(*,*) "storecomp", this%m_storecomp
          !write(*,*) "m_indexcomp1", this%m_indexcomp1
          !write(*,*) "m_indexcomp2", this%m_indexcomp2
           
         end  subroutine minmaxdae2_set_component

!***********************************************************************
!> @brief fixe la longueur d'une tranche d'analyses (en nombre de pas d'entree)
!***********************************************************************
         subroutine minmaxdae2_set_stepout(this,stepout)
          implicit none
          class(t_minmax_dae2), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: stepout   !< nombre de pas de la tranche (en pas d'entree)
          
          this%m_stepout =  stepout
         end  subroutine minmaxdae2_set_stepout

!***********************************************************************
!> @brief reinitialise les min/sum/max pour commencer la tanche
!***********************************************************************
         subroutine minmaxdae2_resetval(this)
          implicit none
          class(t_minmax_dae2), intent(inout):: this  !< dummy argument
          
          integer j

          do j=1, this%m_nbstorecomp
            this%m_minval(j) = 1D300
            this%m_maxval(j) = -1D300
            this%m_sumval(j) = 0
          enddo
          this%m_counter = 0
          
         end  subroutine minmaxdae2_resetval
     
!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree est plein
!! calcule la somme (pour la moyenne), le minium et le maximum 
!! des composantes selectionnees
!***********************************************************************
       subroutine minmaxdae2_oninputbufferfull(this, buffer, poswrite) 
        implicit none
        class(t_minmax_dae2), intent(inout) :: this   !< donnee utilisateur 
        class(t_buffer), intent(inout) :: buffer !< buffer de sortie
        integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
        integer(8) :: j
        real(TREAL), dimension(1+this%m_nbinputcomp) :: R
        real(TREAL), dimension(1+3*this%m_nbstorecomp) :: W
        real(TREAL) t, v1, v2, v, sumv
        integer k,s
        
        
         do j=1, poswrite
           ! lecteure de l'entree
           call buffer%readdata(j, R)            
           t = R(1)
           ! calcule le min/max/sum
           sumv = 0.D0
           do k=1, this%m_nbcalccomp
            v1 = R(1+this%m_indexcomp1(k))
            v2 = R(1+this%m_indexcomp2(k))
            v = (v1-v2)**2
            sumv = sumv+v
            
            s = this%m_storecomp(k)
            ! cas d'une variable ordinaire => une seule determination
            if (v<this%m_minval(s)) this%m_minval(s) = v
            if (v>this%m_maxval(s)) this%m_maxval(s) = v
            this%m_sumval(s) = this%m_sumval(s)+v
           enddo
           ! traitement particulier pour (a_p(1)-a_p(2))**2+(e_p(1)-e_p(2))**2
            s = this%m_storecomp(this%m_nbcalccomp+1)
            ! cas d'une variable ordinaire => une seule determination
            if (sumv<this%m_minval(s)) this%m_minval(s) = sumv
            if (sumv>this%m_maxval(s)) this%m_maxval(s) = sumv
            this%m_sumval(s) = this%m_sumval(s)+sumv
           
           
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
            call minmaxdae2_resetval(this)        
            ! initialise le min/max/sum avec le dernier point
            sumv = 0.D0
            do k=1, this%m_nbcalccomp
                v1 = R(1+this%m_indexcomp1(k))
                v2 = R(1+this%m_indexcomp2(k))
                v = (v1-v2)**2
                sumv = sumv+v
                s = this%m_storecomp(k)
                this%m_minval(s) = v
                this%m_maxval(s) = v
                this%m_sumval(s) = v
            enddo
           ! traitement particulier pour (a_p(1)-a_p(2))**2+(e_p(1)-e_p(2))**2
            s = this%m_storecomp(this%m_nbcalccomp+1)
            ! cas d'une variable ordinaire => une seule determination
            this%m_minval(s) = sumv
            this%m_maxval(s) = sumv
            this%m_sumval(s) = sumv

           endif
          enddo
        
        call buffer%empty()
       
       end subroutine minmaxdae2_oninputbufferfull


      end module
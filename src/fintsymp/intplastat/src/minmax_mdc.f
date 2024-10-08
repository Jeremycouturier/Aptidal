!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file minmax_mdc.f 
!!  \brief Calcule le min/moy/max des a_p(1)-a_p(2),  la_p(1)-la_p(2), 
!! et pi_p(1)-pi_p(2) sur une tranche de temps donne avec p(1) et p(2) les indices de 2 planetes
!!
!!
!! Les angles sont redresses pour etre continu.
!! Les angles sont redresses lorsque les sauts sont superieurs a 3.
!!
!!    --------------   min/moy/max      --------------\n
!!    | buffer src |  ----------------> | buffer dst |\n
!!    --------------                    --------------\n
!!
! history : creation 26/03/2016
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le calcul des min/moy/max avec des angles continus et redresses
!***********************************************************************
      module mod_minmax_mdc
       use mod_minmax
       
              
!***********************************************************************
!> @class t_minmax_mdc
!! Calcule le min/moy/max des a_p(1)-a(p2),  la_p(1)-la(p2), 
!! et pi_p(1)-pi_(p2) sur une tranche de temps donne avec p(1) et p(2) les indices de 2 planetes
!! Le buffer d'entree, contient le temps + M composantes.\n
!! les M composants doivent etre des elements elliptiques dont les 3 dernieres composantes sont des angles.
!! Le buffer de sortie contient le temps + les ???*P min/moy/max.\n
!! Le stockage dans le buffer de sortie se fait de maniere suivante :\n
!!    temps , min/moy/max de composante_1, 
!!    min/moy/max de composante_2 redresses,
!!    min/moy/max de composante_3 redresses
!!
!! La configuration des composantes a analyser se fait par la fonction
!! set_component.
!!
!! P sera nomme : m_nbcalccomp
!! M sera nomme : m_nbinputcomp
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_minmax_mdc
          private

          integer :: m_nbcalccomp = 0 !< nombre de composantes dont on calcule le min/moy/max
          integer :: m_nbinputcomp = 0 !< nombre de composantes dont on extrait les m_nbcalccomp
          !> tableau de "m_nbcalccomp" indices indiquant quelles composantes sont extraites parmi les "m_nbinputcomp" composantes d'entree 
          !! pour la premiere planete
          integer, dimension(:), allocatable :: m_indexcomp1 
          !> tableau de "m_nbcalccomp" indices indiquant quelles composantes sont extraites parmi les "m_nbinputcomp" composantes d'entree 
          !! pour la deuxieme planete
          integer, dimension(:), allocatable :: m_indexcomp2 
          !> tableau de "m_nbcalccomp" booleens indiquant si la composante extraite est un angle ou pas
          logical, dimension(:), allocatable :: m_anglecomp 
          real(TREAL), dimension(:), allocatable :: m_oldangle !< derniere valeur des "m_nbcalccomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_accumsaut !< entier correspodnant a la somme cumulee des sauts de +-2*pi a effectuer 
          real(TREAL), dimension(:), allocatable :: m_minval !< valeur minimale courante des "m_nbcalccomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_maxval !< valeur maximale courante des "m_nbcalccomp" composantes 
          real(TREAL), dimension(:), allocatable :: m_sumval !< somme courante des "m_nbcalccomp" composantes 
          integer(8) :: m_stepout !< nombre de pas pour calculer le min/max et ecrire le resultat
          integer(8) :: m_counter = 0 !< compteur de pas accumule
          logical :: m_firstcall = .true. !< = true => au 1er appel, on ne redresse pas les angles. =false => on redresse les angles
          real(TREAL) :: m_pi !< valeur de pi

      contains
          procedure :: set_output => minmaxmdc_set_output ! fixe le buffer de sortie et le consommateur
          procedure :: set_component => minmaxmdc_set_component ! definit les composantes a analyser
          procedure :: set_stepout => minmaxmdc_set_stepout ! fixe la longueur d'une tranche pour les calculs (en nombre de pas d'entree)
          
          procedure :: oninputbufferfull => minmaxmdc_oninputbufferfull ! procedure appellee lorsque le buffer d'entree est plein

      end type t_minmax_mdc  


      private minmaxmdc_resetval
      
      contains

!***********************************************************************
!> @brief fixe la taille du buffer de sortie et son consommateur
!***********************************************************************
         subroutine minmaxmdc_set_output(this,buffersize,buffercs)
          implicit none
          class(t_minmax_mdc), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%setconsumer("MIN/MAX MDC")
          call this%m_buffer%init(buffersize, 1+3*this%m_nbcalccomp,       &
     &              buffercs) 
         end  subroutine minmaxmdc_set_output

!***********************************************************************
!> @brief definit les composantes a analyser
!! a partir du tableau d'entiers fournis(indice par rapport au buffer d'entree)
!***********************************************************************
         subroutine minmaxmdc_set_component(this,arPcomp1, arPcomp2,        &
     &                       anglecomp, Mcomp)
          implicit none
          class(t_minmax_mdc), intent(inout):: this  !< dummy argument
          integer, dimension(:), intent(in) :: arPcomp1  !< tableau d'indices des composantes a analyser pour la premiere planete
          integer, dimension(:), intent(in) :: arPcomp2  !< tableau d'indices des composantes a analyser pour la deuxieme planete
          logical, dimension(:), intent(in) :: anglecomp  !< tableau de booleen =true => c'est un angle
          integer, intent(in) :: Mcomp  !< nombre de composantes du buffer d'entree
          
          this%m_nbcalccomp =  size(arPcomp1)
          this%m_nbinputcomp =  Mcomp
          allocate(this%m_indexcomp1(1:this%m_nbcalccomp))
          allocate(this%m_indexcomp2(1:this%m_nbcalccomp))
          allocate(this%m_anglecomp(1:this%m_nbcalccomp))
         
          this%m_indexcomp1(:) = arPcomp1(1:this%m_nbcalccomp)
          this%m_indexcomp2(:) = arPcomp2(1:this%m_nbcalccomp)
          this%m_anglecomp(:) = anglecomp(1:this%m_nbcalccomp)
          
          allocate(this%m_minval(1:this%m_nbcalccomp))
          allocate(this%m_maxval(1:this%m_nbcalccomp))
          allocate(this%m_sumval(1:this%m_nbcalccomp))
          allocate(this%m_oldangle(1:this%m_nbcalccomp))
          allocate(this%m_accumsaut(1:this%m_nbcalccomp))
          this%m_accumsaut(:) = REALC(0.)

          call minmaxmdc_resetval(this)
          this%m_counter = -1 ! pour forcer la 1ere valeur a afficher
          
          this%m_pi = atan2(REALC(0.), REALC(-1.))
          !write(*,*) this%m_pi
          !write(*,*) "storecomp", this%m_storecomp
          !write(*,*) "m_indexcomp1", this%m_indexcomp1
          !write(*,*) "m_indexcomp2", this%m_indexcomp2
          !write(*,*) "m_anglecomp", this%m_anglecomp
           
         end  subroutine minmaxmdc_set_component

!***********************************************************************
!> @brief fixe la longueur d'une tranche d'analyses (en nombre de pas d'entree)
!***********************************************************************
         subroutine minmaxmdc_set_stepout(this,stepout)
          implicit none
          class(t_minmax_mdc), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: stepout   !< nombre de pas de la tranche (en pas d'entree)
          
          this%m_stepout =  stepout
         end  subroutine minmaxmdc_set_stepout

!***********************************************************************
!> @brief reinitialise les min/sum/max pour commencer la tranche
!***********************************************************************
         subroutine minmaxmdc_resetval(this)
          implicit none
          class(t_minmax_mdc), intent(inout):: this  !< dummy argument
          
          integer j

          do j=1, this%m_nbcalccomp
            this%m_minval(j) = 1D300
            this%m_maxval(j) = -1D300
            this%m_sumval(j) = 0
          enddo
          this%m_counter = 0
          
         end  subroutine minmaxmdc_resetval
     
!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree est plein
!! calcule la somme (pour la moyenne), le minium et le maximum 
!! des composantes selectionnees
!***********************************************************************
       subroutine minmaxmdc_oninputbufferfull(this, buffer, poswrite) 
        implicit none
        class(t_minmax_mdc), intent(inout) :: this   !< donnee utilisateur 
        class(t_buffer), intent(inout) :: buffer !< buffer de sortie
        integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
        integer(8) :: j
        real(TREAL), dimension(1+this%m_nbinputcomp) :: R
        real(TREAL), dimension(1+3*this%m_nbcalccomp) :: W
        real(TREAL) t, v1, v2, v, vprev
        integer k,s
        logical :: firstcall
        
        
          
         do j=1, poswrite
           firstcall = this%m_firstcall
           ! lecture de l'entree
           call buffer%readdata(j, R)            
           t = R(1)
           ! calcule le min/max/sum
           do k=1, this%m_nbcalccomp
            v1 = R(1+this%m_indexcomp1(k))
            v2 = R(1+this%m_indexcomp2(k))
            v = v1-v2
            
            s = k
            if (this%m_anglecomp(k)) then
            ! cas d'un angle => redressement si necessaire
             vprev = this%m_oldangle(s)
             this%m_oldangle(s) = v
             if (firstcall) then
              this%m_firstcall = .false.
             else
              if (v-vprev>REALC(3.)) then
               this%m_accumsaut(s)=this%m_accumsaut(s)-1
              else if (v-vprev<-REALC(3.)) then
               this%m_accumsaut(s)=this%m_accumsaut(s)+1
              endif
              v= v+this%m_accumsaut(s)*REALC(2.)*this%m_pi
             endif
            endif 
            if (v<this%m_minval(s)) this%m_minval(s) = v
            if (v>this%m_maxval(s)) this%m_maxval(s) = v
            this%m_sumval(s) = this%m_sumval(s)+v
           enddo
           
           this%m_counter = this%m_counter+1
           if (this%m_counter==this%m_stepout) then
            ! fin de la tranche 
            ! la moyenne se fait sur m_counter+1 valeurs (lie a la presence pour le temps t=0)
            W(1) = t
            do k=1, this%m_nbcalccomp
             !write(*,*) "w", k, this%m_minval(k), this%m_maxval(k) 
             W(1+3*(k-1)+1) = this%m_minval(k) 
             W(1+3*(k-1)+3) = this%m_maxval(k) 
             W(1+3*(k-1)+2) = this%m_sumval(k)/(this%m_counter+1) 
            enddo
            
            call this%m_buffer%writedata(W)   
            
            ! reinitialise la tranche
            call minmaxmdc_resetval(this)        
            ! intialise le min/max/sum avec le dernier point
            do k=1, this%m_nbcalccomp
                v1 = R(1+this%m_indexcomp1(k))
                v2 = R(1+this%m_indexcomp2(k))
                v = v1-v2
                s = k
                if (this%m_anglecomp(k)) then
                 ! cas d'un angle => redressement si necessaire
                 v = v+this%m_accumsaut(s)*REALC(2.)*this%m_pi
                endif
                this%m_minval(s) = v
                this%m_maxval(s) = v
                this%m_sumval(s) = v
            enddo

           endif
          enddo
        
        call buffer%empty()
       
       end subroutine minmaxmdc_oninputbufferfull


      end module
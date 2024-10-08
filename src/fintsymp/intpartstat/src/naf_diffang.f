!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file naf_diffang.f 
!!  \brief Calcule l'analyse en frequence de exp(i*(angle_p1(k)-angle_p2(k)))
!! sur une tranche de temps donne
!!
!!    --------------        naf exp(i*...)        --------------\n
!!    | buffer src |  --------------------------> | buffer dst |\n
!!    --------------                              --------------\n
!!
! history : creation 17/09/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le calcul de l'analyse en frequence
!***********************************************************************
      module mod_naf_diffang
       use mod_naf_base
       
!***********************************************************************
!> @class t_naf_diffang
!! Calcule l'analyse en frequence de exp(i*(angle_p1(k)-angle_p2(k)))
!! parmi M composantes sur une tranche de temps donne. \n
!! Le buffer d'entree, contient le temps + M composantes.\n
!!
!! Le buffer de sortie contient le temps + les P*NTERM*(ampl complexe, freq)\n
!! Le stockage dans le buffer de sortie se fait de maniere suivante :\n
!!    temps , (ampl1, freq1, ampl2, freq2, ..., ampln, freqn) de exp(i*(angle_p1(1)-angle_p2(1)))....
!!   repete ensuite pour les angle_p1(i)-angle_p2(i) suivants 
!! 
!! La configuration des composantes a analyser se fait par la fonction
!! set_component.
!!
!! M sera nomme : m_nbinputcomp
!!  
!***********************************************************************
      type, extends(t_naf_base) :: t_naf_diffang
          private

          integer, dimension(:), allocatable :: m_arPcomp1  !< tableau d'indices des composantes a analyser pour le premier corps
          integer, dimension(:), allocatable :: m_arPcomp2  !< tableau d'indices des composantes a analyser pour le deuxieme corps

      contains
          procedure :: set_component12 => naf_diffang_set_component ! definit les composantes a analyser
          
          procedure :: oninputbufferfull=>naf_diffang_oninputbufferfull ! procedure appellee lorsque le buffer d'entree est plein

      end type t_naf_diffang  

      contains



!***********************************************************************
!> @brief definit les composantes a analyser angle_p1(), angle_p2()
!! a partir du tableau d'entiers fournis(indice par rapport au buffer d'entree).\n
!***********************************************************************
         subroutine naf_diffang_set_component(this,arPcomp1, arPcomp2,        &
     &                       Mcomp)
          implicit none
          class(t_naf_diffang), intent(inout):: this  !< dummy argument
          integer, dimension(:), intent(in) :: arPcomp1  !< tableau d'indices des composantes a analyser pour la premiere planete
          integer, dimension(:), intent(in) :: arPcomp2  !< tableau d'indices des composantes a analyser pour la deuxieme planete
          integer, intent(in) :: Mcomp  !< nombre de composantes du buffer d'entree
          
          integer, dimension(1:size(arPcomp1)) :: nafzindice 
          integer k
          integer ncomp
          
          ncomp = size(arPcomp1)
          
          allocate(this%m_arPcomp1(1:ncomp))
          allocate(this%m_arPcomp2(1:ncomp))
          this%m_arPcomp1(:) = arPcomp1(:)
          this%m_arPcomp2(:) = arPcomp2(:)
          
          do k=1, ncomp
            nafzindice(k) = 2*(k-1)+1  ! 2 car composante complexe
          enddo
             
          call naf_base_set_component(this,nafzindice, Mcomp)
          
         end  subroutine naf_diffang_set_component

!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree est plein
!! stocke les composantes dans le tableau intermediaire
!! des composantes selectionnees
!***********************************************************************
         subroutine naf_diffang_oninputbufferfull(this,buffer,poswrite)
          implicit none
          class(t_naf_diffang), intent(inout) :: this   !< donnee utilisateur 
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL), dimension(1+this%m_nbinputcomp) :: R
          real(TREAL) t, l1,l2
          complex(2*TREAL) zl
          integer k
          
          call naf_base_init_zcomp(this)
          
          call this%m_buffer%tryset_multictrlerr_buf(buffer)
          
          do j=1, poswrite
             ! lecture de l'entree
             call buffer%readdata(j, R)            
             t = R(1)
             this%m_tcomp(this%m_counter+1) = t
             ! calcul des  exp(i*(l_p1(k)-l_p2(k))) 
             do k=1, size(this%m_arPcomp1)
                l1 = R(1+this%m_arPcomp1(k))
                l2 = R(1+this%m_arPcomp2(k))
                zl = CMPLX(cos(l1-l2),sin(l1-l2), kind(REALC(1.0)))
                this%m_zcomp(k, this%m_counter+1) = zl
             enddo

             this%m_counter = this%m_counter+1
             if (this%m_counter==this%m_stepout) then
              ! tranche finie => analyse toutes le composantes
              call naf_base_analyzeall(this, buffer)
             endif
          enddo
            
          call buffer%tryset_multictrlerr_buf(this%m_buffer)
          call buffer%empty()
         
         end subroutine naf_diffang_oninputbufferfull

      end module
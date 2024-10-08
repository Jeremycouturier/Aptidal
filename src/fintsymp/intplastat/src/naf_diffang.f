!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file naf_diffang.f 
!!  \brief Calcule l'analyse en frequence de exp(i*(l_p(1)-l_p(2))) et  exp(i*(pi_p(1)-pi_p(2)))
!! sur une tranche de temps donne
!!
!!    --------------        naf exp(i*...)        --------------\n
!!    | buffer src |  --------------------------> | buffer dst |\n
!!    --------------                              --------------\n
!!
! history : creation 02/04/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le calcul de l'analyse en frequence
!***********************************************************************
      module mod_naf_diffang
       use mod_naf_base
       
!***********************************************************************
!> @class t_naf_diffang
!! Calcule l'analyse en frequence de exp(i*(l_p(1)-l_p(2))) et  exp(i*(pi_p(1)-pi_p(2))) 
!! parmi M composantes sur une tranche de temps donne. \n
!! Le buffer d'entree, contient le temps + M composantes.\n
!! Les M composants doivent etre des elements elliptiques dont les 3 dernieres composantes sont des angles (la, pi, Omega).\n
!!
!! Le buffer de sortie contient le temps + les P*NTERM*(ampl complexe, freq)\n
!! Le stockage dans le buffer de sortie se fait de maniere suivante :\n
!!    temps , (ampl1, freq1, ampl2, freq2, ..., ampln, freqn) de exp(i*(l_p(1)-l_p(2))), (ampl1, freq1,...) de exp(i*(pi_p(1)-pi_p(2)))
!!
!! La configuration des composantes a analyser se fait par la fonction
!! set_component.
!!
!! M sera nomme : m_nbinputcomp
!!  
!***********************************************************************
      type, extends(t_naf_base) :: t_naf_diffang
          private

          integer, dimension(1:2) :: m_arPcomp1  !< tableau d'indices des composantes a analyser pour la premiere planete
          integer, dimension(1:2) :: m_arPcomp2  !< tableau d'indices des composantes a analyser pour la deuxieme planete

      contains
          procedure :: set_component12 => naf_diffang_set_component ! definit les composantes a analyser
          
          procedure :: oninputbufferfull=>naf_diffang_oninputbufferfull ! procedure appellee lorsque le buffer d'entree est plein

      end type t_naf_diffang  

      contains



!***********************************************************************
!> @brief definit les composantes a analyser l_p(1),l_p(2), pi_p(1) et pi_p(2)
!! a partir du tableau d'entiers fournis(indice par rapport au buffer d'entree).\n
!***********************************************************************
         subroutine naf_diffang_set_component(this,arPcomp1, arPcomp2,        &
     &                       Mcomp)
          implicit none
          class(t_naf_diffang), intent(inout):: this  !< dummy argument
          integer, dimension(1:2), intent(in) :: arPcomp1  !< tableau d'indices des composantes a analyser pour la premiere planete
          integer, dimension(1:2), intent(in) :: arPcomp2  !< tableau d'indices des composantes a analyser pour la deuxieme planete
          integer, intent(in) :: Mcomp  !< nombre de composantes du buffer d'entree
          
          integer, dimension(1:2) :: nafzindice 
          
          this%m_arPcomp1 =  arPcomp1
          this%m_arPcomp2 =  arPcomp2
          
          nafzindice(1) = 1
          nafzindice(2) = 3
          
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
          real(TREAL) t, l1,l2, pi1,pi2
          complex(TREAL) zl, zp
          
          
          call naf_base_init_zcomp(this)
          
          do j=1, poswrite
             ! lecture de l'entree
             call buffer%readdata(j, R)            
             t = R(1)
             this%m_tcomp(this%m_counter+1) = t
             ! calcul des  exp(i*(l_p(1)-l_p(2))) et  exp(i*(pi_p(1)-pi_p(2)))
             l1 = R(1+this%m_arPcomp1(1))
             l2 = R(1+this%m_arPcomp2(1))
             pi1 = R(1+this%m_arPcomp1(2))
             pi2 = R(1+this%m_arPcomp2(2))
             zl = CMPLX(cos(l1-l2),sin(l1-l2), kind(REALC(1.0)))
             zp = CMPLX(cos(pi1-pi2),sin(pi1-pi2), kind(REALC(1.0)))
             this%m_zcomp(1, this%m_counter+1) = zl
             this%m_zcomp(2, this%m_counter+1) = zp

             this%m_counter = this%m_counter+1
             if (this%m_counter==this%m_stepout) then
              ! tranche finie => analyse toutes le composantes
              call naf_base_analyzeall(this, buffer)
             endif
          enddo
            
          call buffer%empty()
         
         end subroutine naf_diffang_oninputbufferfull

      end module
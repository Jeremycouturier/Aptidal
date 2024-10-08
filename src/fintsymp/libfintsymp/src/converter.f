!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file converter.f 
!!  \brief conversion de coordonnees entre 2 buffers pour les planetes
!!
!!    --------------   conversion   --------------\n
!!    | buffer src |  ------------> | buffer dst |\n
!!    --------------                --------------\n
!!
! history : creation 03/08/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la conversion de coordonnees entre 2 buffers
!***********************************************************************
      module mod_converter
       use mod_buffer
       use mod_arret
       
!***********************************************************************
!> @class t_plan_converter
!! classe de base pour assurer les conversions de  coordonnees.
!! Cette classe ne doit pas etre appellee directement.
!! Il faut utiliser une des classes derivees
!!  
!***********************************************************************
      type :: t_plan_converter   !// , abstract
          private

          type(t_buffer), public :: m_buffer !< buffer de sortie contenant la sortie
          type(t_buffer_consumer), pointer :: m_input_buffer =>NULL() !< consommateur du buffer en entree

          type(t_arret),public :: m_arret !< cause de l'erreur dans les conversions
          
          character(len=50) :: m_dotname !< nom du convertisseur pour le graphique dot

      contains
          procedure :: getconsumer => plan_converter_getconsumer   ! retourne le consommateur du buffer  en entree
          procedure :: setconsumer => plan_converter_setconsumer  ! fixe le proprietaire du consommateur du buffer en entree

          procedure, NON_OVERRIDABLE :: set_graphdotname =>             &
     &              plan_converter_set_graphdotname ! affiche le graphique dot
        
          procedure :: oninputbufferflush =>                              &
     &           plan_converter_oninputbufferflush ! procedure appellee lorsque le buffer d'entree doit etre propage

          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          procedure :: oninputbufferfull =>                              &
     &           plan_converter_oninputbufferfull! procedure appellee lorsque le buffer d'entree est plein

          final ::  plan_converter_destructor ! destructor       

      end type t_plan_converter  

      abstract interface 
!***********************************************************************
!>  type de fonction pour generer le graph au format dot
!***********************************************************************
       subroutine t_oninputbuffer(this,buffer,poswrite)
             import t_plan_converter
             import t_buffer
             class(t_plan_converter), intent(inout) :: this  !< dummy argument
             class(t_buffer), intent(inout) :: buffer !< buffer de sortie
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
        
       end subroutine t_oninputbuffer
      end interface 
      
!***********************************************************************
!> @class t_plan_PhVb2PhVh
!! conversion de position heliocentrique/vitesse barycentrique 
!!  en position heliocentrique/vitesse heliocentrique
!!  pour des planetes
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_plan_PhVb2PhVh
          ! private : following data should only be used by derived type
          real(TREAL), dimension(:), allocatable :: m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          integer :: m_plan_nb !< nombre de planetes
          
       contains
          
          procedure :: set_mass =>plan_PhVb2PhVh_set_mass ! fixe les masses a utiliser
          procedure :: set_output => plan_PhVb2PhVh_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein
         
          procedure::oninputbufferfull => PhVb2PhVh_oninputbufferfull ! procedure appellee lorsque le buffer est plein
                    
      end type t_plan_PhVb2PhVh  

!***********************************************************************
!> @class t_plan_PjVj2PhVh
!! conversion de position jacobi/vitesse jacobi 
!!  en position heliocentrique/vitesse heliocentrique
!!  pour des planetes
!!
!! m_buffer  de t_plan_converter :  buffer de sortie contenant les positions heliocentriques/vitesses heliocentriques
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_plan_PjVj2PhVh
          private
          real(TREAL), dimension(:), allocatable :: m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          real(TREAL), dimension(:), allocatable :: m_plan_eta   !< eta jacobi
          integer :: m_plan_nb !< nombre de planetes
          
       contains
          
          procedure :: set_mass =>plan_PjVj2PhVh_set_mass ! fixe les masses a utiliser
          procedure :: set_output => plan_PjVj2PhVh_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein
         
          procedure::oninputbufferfull => PjVj2PhVh_oninputbufferfull ! procedure appellee lorsque le buffer est plein
                    
      end type t_plan_PjVj2PhVh  

!***********************************************************************
!> @class t_plan_PjVj2PhVb
!! conversion de position jacobi/vitesse jacobi 
!!  en position heliocentrique/vitesse varycentrique
!!  pour des planetes
!!
!! m_buffer  de t_plan_converter :  buffer de sortie contenant les positions heliocentriques/vitesses barycentriques
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_plan_PjVj2PhVb
          private
          real(TREAL), dimension(:), allocatable :: m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          real(TREAL), dimension(:), allocatable :: m_plan_eta   !< eta jacobi
          integer :: m_plan_nb !< nombre de planetes

       contains
          
          procedure :: set_mass =>plan_PjVj2PhVb_set_mass ! fixe les masses a utiliser
          procedure :: set_output => plan_PjVj2PhVb_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein
         
          procedure::oninputbufferfull => PjVj2PhVb_oninputbufferfull ! procedure appellee lorsque le buffer est plein
                    
      end type t_plan_PjVj2PhVb  


!***********************************************************************
!> @class t_plan_PhVb2ell
!! conversion de position heliocentrique/vitesse barycentrique 
!!  en elements elliptiques canoniques ou non canoniques du type
!!  (a,e,I,M,om,Om) ou (a,la,k,h,q,p)  \n
!!  pour des planetes \n
!! le type de coordonnees en sortie est fixe par set_kind
!!
!! m_buffer  de t_plan_converter :  buffer de sortie contenant les elements elliptiques
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_plan_PhVb2ell
          !private
          real(TREAL), dimension(:), allocatable :: m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          real(TREAL), dimension(:), allocatable :: m_plan_mpsbeta !< (m_p+m_star)/m_star pour chaque  planete
          real(TREAL), dimension(:), allocatable ::   m_plan_mu_helio !< mu heliocentrique
          integer :: m_plan_nb !< nombre de planetes
          !> type de coordonnee du buffer en sortie 
          !! 1 = (a,e,I,M,om,Om) elliptiques canoniques
          !! 2 = (a,e,I,M,om,Om) elliptiques non canoniques
          !! 3 = (a,la,k,h,q,p) elliptiques canoniques
          !! 4 = (a,la,k,h,q,p) elliptiques non canonique
          integer :: m_type_output = 1
          
       contains
          
          procedure :: set_kind =>plan_PhVb2ell_set_kind ! fixe le type de coordonnees en sortie
          procedure :: set_mass =>plan_PhVb2ell_set_mass ! fixe les masses a utiliser
          procedure :: set_output => plan_PhVb2ell_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein
         
          procedure::oninputbufferfull => PhVb2ell_oninputbufferfull ! procedure appellee lorsque le buffer est plein
          
      end type t_plan_PhVb2ell  
    
!***********************************************************************
!> @class t_plan_PhVh2ell
!! conversion de position heliocentrique/vitesse heliocentrique 
!!  en elements elliptiques canoniques ou non canoniques du type
!!  (a,e,I,M,om,Om) ou (a,la,k,h,q,p)  \n
!!  pour des planetes \n
!! le type de coordonnees en sortie est fixe par set_kind
!!
!! m_buffer  de t_plan_converter :  buffer de sortie contenant les elements elliptiques
!!  
!***********************************************************************
      type, extends(t_plan_PhVb2ell) :: t_plan_PhVh2ell
          
       contains
          
          procedure::oninputbufferfull => PhVh2ell_oninputbufferfull ! procedure appellee lorsque le buffer est plein
          
      end type t_plan_PhVh2ell  

      contains

!***********************************************************************
!***********************************************************************
!***********************************************************************
! implementation
!***********************************************************************
!***********************************************************************
!***********************************************************************

!***********************************************************************
!> @brief destructor
!***********************************************************************
       subroutine plan_converter_destructor (this)
        implicit none
        type ( t_plan_converter ) :: this  !< dummy argument

             if (associated(this%m_input_buffer)) then
              deallocate(this%m_input_buffer)
             endif
       end subroutine plan_converter_destructor


!***********************************************************************
!> @brief retourne le consommateur du buffer d'entree
!***********************************************************************
            function plan_converter_getconsumer(this) result(cs)
             implicit none
             class(t_plan_converter), intent(inout):: this  !< dummy argument
#if __INTEL_COMPILER <= 1600
             type(t_buffer_consumer), pointer :: cs 
#else
             class(t_buffer_consumer), pointer :: cs 
#endif
             cs => this%m_input_buffer
            end  function plan_converter_getconsumer

!***********************************************************************
!> @brief fixe le proprietaire de ce consommateur et le nom de this
!***********************************************************************
      subroutine plan_converter_setconsumer(this, dotname) 
       implicit none
       class(t_plan_converter), intent(inout):: this  !< dummy argument
       character(len=*), intent(in):: dotname !< nom du convertisseur
       procedure(t_procbufferfull), pointer :: pfcnfull
       procedure(t_procbufferfull), pointer :: pfcnflush
       procedure(t_procgraphdot), pointer :: pfcndot
       
        allocate(this%m_input_buffer)
        pfcnfull => plan_converter_onbufferfull_cb
        pfcnflush => plan_converter_onbufferflush_cb
        pfcndot => plan_converter_ongraphdot
        call this%m_input_buffer%set_owner(this, pfcnfull, pfcnflush,       &
     &   pfcndot)
     
        this%m_dotname = dotname
      end subroutine  plan_converter_setconsumer

!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
         subroutine plan_converter_oninputbufferfull(this, buffer,           &
     &        poswrite) 
          implicit none
          class(t_plan_converter), intent(inout) :: this   !< donnee utilisateur 
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
            
            ! cette fonction doit etre surchargee par le type derive
            stop 'plan_converter_oninputbufferfull'
         
         end subroutine plan_converter_oninputbufferfull

!***********************************************************************
! @brief fonction appellee lorsque le buffer d'entree est plein
!! fonction callback
!***********************************************************************
         subroutine plan_converter_onbufferfull_cb(this, userdata,               &
     &                   poswrite) 
          implicit none
          class(t_buffer), intent(inout) :: this !< buffer de sortie
          class(*), intent(inout) :: userdata    !< de type t_plan_PhVb2PhVh
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)

          select type(userdata)
           class is(t_plan_converter)
            call userdata%oninputbufferfull(this, poswrite)
           class default
            stop 'plan_converter_onbufferfull_cb : bad class'
          end select 
                   
         end subroutine plan_converter_onbufferfull_cb


!***********************************************************************
! @brief fonction appellee lorsque le buffer d'entree doit etre propage (par exemple a la fin)
!***********************************************************************
         subroutine plan_converter_oninputbufferflush(this, buffer,           &
     &        poswrite) 
          implicit none
          class(t_plan_converter), intent(inout) :: this   !< donnee utilisateur 
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
            
            call this%m_buffer%flushdata()
         
         end subroutine plan_converter_oninputbufferflush

!***********************************************************************
! @brief fonction appellee lorsque le buffer d'entree doit etre propage (par exemple a la fin)
!! fonction callback
!***********************************************************************
         subroutine plan_converter_onbufferflush_cb(this, userdata,               &
     &                   poswrite) 
          implicit none
          class(t_buffer), intent(inout) :: this !< buffer de sortie
          class(*), intent(inout) :: userdata    !< de type t_plan_PhVb2PhVh
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)

          select type(userdata)
           class is(t_plan_converter)
            call userdata%oninputbufferfull(this, poswrite)
            call userdata%oninputbufferflush(this, poswrite)
           class default
          end select 
                   
         end subroutine plan_converter_onbufferflush_cb

         
!***********************************************************************
!> @brief  fonction appellee pour generer le graph au format dot
!***********************************************************************
       subroutine plan_converter_ongraphdot(this,userdata, nodesrc,          &
     &                   nodecounter, num_file)
        implicit none
        class(t_buffer), intent(inout) :: this !< buffer ayant genere l'appel
        class(*), intent(inout) :: userdata !< donnee utilisateur de type t_plan_converter
        integer, intent(in) :: nodesrc !< numero du noeud  source 
        integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
        integer, intent(in) :: num_file !<numero du fichier de sortie
             
        integer curnode

        curnode = nodecounter+1

          select type(userdata)
           class is(t_plan_converter)
             call buffer_dot_writelabel(num_file,curnode,                 &
     &          trim(userdata%m_dotname))
             call buffer_dot_writenode(num_file, nodesrc,curnode)
             nodecounter = curnode
             call userdata%m_buffer%graphdot(curnode,nodecounter,          &
     &          num_file)
           class default
            stop 'plan_converter_ongraphdot : bad class'
          end select 

       end subroutine plan_converter_ongraphdot

!***********************************************************************
!> @brief specifie le nom pour le graphique DOT utilise par this
!***********************************************************************
       subroutine plan_converter_set_graphdotname (this, dotname)
        implicit none
        class(t_plan_converter), intent(inout) :: this  !< dummy argument
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        this%m_dotname = dotname
       end subroutine plan_converter_set_graphdotname
       

!***********************************************************************
!> @brief fixe la masse de l'etoile et des planetes a utiliser pour la conversion 
!***********************************************************************
         subroutine plan_PhVb2PhVh_set_mass(this, m0, mpl)
          implicit none
          class(t_plan_PhVb2PhVh), intent(inout):: this  !< dummy argument
          real(TREAL), intent(in) :: m0  !< masse de l'etoile
          real(TREAL), dimension(:), intent(in) :: mpl  !< masse des planetes
          
          this%m_plan_nb = size(mpl)
          allocate(this%m_plan_mass0(0:this%m_plan_nb))
          this%m_plan_mass0(0) = m0
          this%m_plan_mass0(1:) = mpl
          call this%setconsumer("PhVb2PhVh")
          
         end  subroutine plan_PhVb2PhVh_set_mass

!***********************************************************************
!> @brief fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein 
!***********************************************************************
         subroutine plan_PhVb2PhVh_set_output(this,buffersize,           &
     &                          buffercs)
          implicit none
          class(t_plan_PhVb2PhVh), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%m_buffer%init(buffersize, 1+6*this%m_plan_nb,        &
     &              buffercs) 
          
         end  subroutine plan_PhVb2PhVh_set_output


     
!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
         subroutine PhVb2PhVh_oninputbufferfull(this,buffer,poswrite) 
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_plan_PhVb2PhVh), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL), dimension(1+6*this%m_plan_nb) :: R
          real(TREAL), dimension(3,this%m_plan_nb) :: Vb
          real(TREAL), dimension(3,this%m_plan_nb) :: Vh
          integer nplan
          
            nplan = this%m_plan_nb
            do j=1, poswrite
             call buffer%readdata(j, R)
             Vb = reshape(R(2+3*nplan:1+6*nplan), shape(Vb))
             call coord_vb2vh(nplan,this%m_plan_mass0, Vb,Vh)
             R(2+3*nplan:1+6*nplan)=reshape(Vh,                            &
     &               shape(R(2+3*nplan:1+6*nplan)))
             call this%m_buffer%writedata(R)           
            enddo
          
          ! verifie si une erreur s'est produite dans le successeur et la recupere
          call buffer%checkget_successor_error(this%m_buffer)

          call buffer%empty()
         
         end subroutine PhVb2PhVh_oninputbufferfull

!***********************************************************************
!> @brief fixe la masse de l'etoile et des planetes a utiliser pour la conversion 
!***********************************************************************
         subroutine plan_PjVj2PhVh_set_mass(this, m0, mpl)
          implicit none
          class(t_plan_PjVj2PhVh), intent(inout):: this  !< dummy argument
          real(TREAL), intent(in) :: m0  !< masse de l'etoile
          real(TREAL), dimension(:), intent(in) :: mpl  !< masse des planetes
          integer i
          
          this%m_plan_nb = size(mpl)
          allocate(this%m_plan_mass0(0:this%m_plan_nb))
          this%m_plan_mass0(0) = m0
          this%m_plan_mass0(1:) = mpl
          
          allocate(this%m_plan_eta(0:this%m_plan_nb))
          this%m_plan_eta(0) = m0 
          do i=1,this%m_plan_nb
            this%m_plan_eta(i) = this%m_plan_eta(i-1) + mpl(i)
          end do 

           call this%setconsumer("PjVj2PhVh")

        end  subroutine plan_PjVj2PhVh_set_mass

!***********************************************************************
!> @brief fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein 
!***********************************************************************
         subroutine plan_PjVj2PhVh_set_output(this,buffersize,           &
     &                          buffercs)
          implicit none
          class(t_plan_PjVj2PhVh), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%m_buffer%init(buffersize, 1+6*this%m_plan_nb,        &
     &              buffercs) 
          
         end  subroutine plan_PjVj2PhVh_set_output


     
!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
         subroutine PjVj2PhVh_oninputbufferfull(this, buffer,             &
     &                   poswrite) 
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_plan_PjVj2PhVh), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL), dimension(1+6*this%m_plan_nb) :: R
          real(TREAL), dimension(3,this%m_plan_nb) :: Pj,Vj
          real(TREAL), dimension(3,this%m_plan_nb) :: Ph,Vh
          integer nplan
          
            nplan = this%m_plan_nb
            do j=1, poswrite
             call buffer%readdata(j, R)
             Pj = reshape(R(2:1+3*nplan), shape(Pj))
             Vj = reshape(R(2+3*nplan:1+6*nplan), shape(Vj))
             call coord_j2h(nplan,this%m_plan_mass0,                      &
     &               this%m_plan_eta,Pj,Ph)
             call coord_j2h(nplan,this%m_plan_mass0,                      &
     &               this%m_plan_eta,Vj,Vh)
             R(2:1+3*nplan)=reshape(Ph, shape(R(2:1+3*nplan)))
             R(2+3*nplan:1+6*nplan)=reshape(Vh,                           &
     &               shape(R(2+3*nplan:1+6*nplan)))
             call this%m_buffer%writedata(R)           
            enddo
          
          ! verifie si une erreur s'est produite dans le successeur et la recupere
          call buffer%checkget_successor_error(this%m_buffer)

          call buffer%empty()
         
         end subroutine PjVj2PhVh_oninputbufferfull

!***********************************************************************
!> @brief fixe la masse de l'etoile et des planetes a utiliser pour la conversion 
!***********************************************************************
         subroutine plan_PjVj2PhVb_set_mass(this, m0, mpl)
          implicit none
          class(t_plan_PjVj2PhVb), intent(inout):: this  !< dummy argument
          real(TREAL), intent(in) :: m0  !< masse de l'etoile
          real(TREAL), dimension(:), intent(in) :: mpl  !< masse des planetes
          integer i
          
          this%m_plan_nb = size(mpl)
          allocate(this%m_plan_mass0(0:this%m_plan_nb))
          this%m_plan_mass0(0) = m0
          this%m_plan_mass0(1:) = mpl
          
          allocate(this%m_plan_eta(0:this%m_plan_nb))
          this%m_plan_eta(0) = m0 
          do i=1,this%m_plan_nb
            this%m_plan_eta(i) = this%m_plan_eta(i-1) + mpl(i)
          end do 

          call this%setconsumer("PjVj2PhVb")

         end  subroutine plan_PjVj2PhVb_set_mass

!***********************************************************************
!> @brief fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein 
!***********************************************************************
         subroutine plan_PjVj2PhVb_set_output(this,buffersize,           &
     &                          buffercs)
          implicit none
          class(t_plan_PjVj2PhVb), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%m_buffer%init(buffersize, 1+6*this%m_plan_nb,        &
     &              buffercs) 
          
         end  subroutine plan_PjVj2PhVb_set_output


     
!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
         subroutine PjVj2PhVb_oninputbufferfull(this,  buffer,           &
     &                   poswrite) 
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_plan_PjVj2PhVb), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL), dimension(1+6*this%m_plan_nb) :: R
          real(TREAL), dimension(3,this%m_plan_nb) :: Pj,Vj
          real(TREAL), dimension(3,this%m_plan_nb) :: Ph,Vh, Vb
          integer nplan
          
            nplan = this%m_plan_nb
            do j=1, poswrite
             call buffer%readdata(j, R)
             Pj = reshape(R(2:1+3*nplan), shape(Pj))
             Vj = reshape(R(2+3*nplan:1+6*nplan), shape(Vj))
             call coord_j2h(nplan,this%m_plan_mass0,                     &
     &               this%m_plan_eta,Pj,Ph)
             call coord_j2h(nplan,this%m_plan_mass0,                     &
     &               this%m_plan_eta,Vj,Vh)
             call coord_h2b(nplan,this%m_plan_mass0,Vh,Vb)
             R(2:1+3*nplan)=reshape(Ph, shape(R(2:1+3*nplan)))
             R(2+3*nplan:1+6*nplan)=reshape(Vb,                            &
     &               shape(R(2+3*nplan:1+6*nplan)))
             call this%m_buffer%writedata(R)           
            enddo
            
          ! verifie si une erreur s'est produite dans le successeur et la recupere
          call buffer%checkget_successor_error(this%m_buffer)
          
          call buffer%empty()
         
         end subroutine PjVj2PhVb_oninputbufferfull


!***********************************************************************
!> @brief fixe le nom pour le graphique dot 
!***********************************************************************
         subroutine plan_PhVb2ell_set_graphdotname(this)
          implicit none
          class(t_plan_PhVb2ell), intent(inout):: this  !< dummy argument

          character (len=50), dimension(8) :: dotname
          dotname(1) = "PhVb2(a,e,I,M,om,Om)_canon"
          dotname(2) = "PhVb2(a,e,I,M,om,Om)_non_canon"
          dotname(3) = "PhVb2(a,la,k,h,q,p)_canon"
          dotname(4) = "PhVb2(a,la,k,h,q,p)_non_canon"
          dotname(5) = "??"
          dotname(6) = "PhVb2(a,e,I,la,pi,Omega)_non_canon"
          dotname(7) = "PhVb2(a*cos(la),a*sin(la),k,h,q,p)_canon"
          dotname(8) = "PhVb2(a*cos(la),a*sin(la),k,h,q,p)_non_canon"
          
           call this%set_graphdotname(dotname(this%m_type_output))
           
         end  subroutine plan_PhVb2ell_set_graphdotname
         
         
!***********************************************************************
!> @brief fixe le type de coordonnee voulue en sortie 
!***********************************************************************
         subroutine plan_PhVb2ell_set_kind(this, itype)
          implicit none
          class(t_plan_PhVb2ell), intent(inout):: this  !< dummy argument
          !> type de coordonnee en sortie du buffer
          !! 1 = (a,e,I,M,om,Om) elliptiques canoniques
          !! 2 = (a,e,I,M,om,Om) elliptiques non canoniques
          !! 3 = (a,la,k,h,q,p) elliptiques canoniques
          !! 4 = (a,la,k,h,q,p) elliptiques non canonique
          !! 6 = (a,e,I,la,pi,Omega) elliptiques non canoniques
          !! 7 = (a*cos(la),a*sin(la),k,h,q,p) elliptiques canoniques
          !! 8 = (a*cos(la),a*sin(la),k,h,q,p) elliptiques non canonique
          !! 11...18 = reserve pour les circum-binaires
          integer, intent(in) :: itype  
          
          this%m_type_output = itype
           call plan_PhVb2ell_set_graphdotname(this)
         
         end  subroutine plan_PhVb2ell_set_kind

!***********************************************************************
!> @brief fixe la constante de gauss, masse de l'etoile et des planetes a utiliser pour la conversion 
!***********************************************************************
         subroutine plan_PhVb2ell_set_mass(this, cG, m0, mpl)
          implicit none
          class(t_plan_PhVb2ell), intent(inout):: this  !< dummy argument
          real(TREAL), intent(in) :: cG          !< constante de Gauss AU**3/an**2
          real(TREAL), intent(in) :: m0  !< masse de l'etoile
          real(TREAL), dimension(:), intent(in) :: mpl  !< masse des planetes
          
          this%m_plan_nb = size(mpl)
          allocate(this%m_plan_mass0(0:this%m_plan_nb))
          allocate(this%m_plan_mu_helio(this%m_plan_nb))
          allocate(this%m_plan_mpsbeta(this%m_plan_nb))
          this%m_plan_mass0(0) = m0
          this%m_plan_mass0(1:) = mpl
          this%m_plan_mu_helio = cG*(mpl+m0)
          this%m_plan_mpsbeta = (mpl+m0)/m0
         
           call this%setconsumer("")
           call plan_PhVb2ell_set_graphdotname(this)

        end  subroutine plan_PhVb2ell_set_mass

!***********************************************************************
!> @brief fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein 
!***********************************************************************
         subroutine plan_PhVb2ell_set_output(this,buffersize,           &
     &                          buffercs)
          implicit none
          class(t_plan_PhVb2ell), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%m_buffer%init(buffersize, 1+6*this%m_plan_nb,        &
     &              buffercs) 
          
         end  subroutine plan_PhVb2ell_set_output


!***********************************************************************
! @brief fonction appellee lorsque le buffer d'entree est plein
!! 
!***********************************************************************
         subroutine PhVb2ell_oninputbufferfull(this,buffer,                &
     &                   poswrite) 
          use mod_elliptid
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_plan_PhVb2ell), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          integer i
          real(TREAL), dimension(1+6*this%m_plan_nb) :: R
          real(TREAL), dimension(3,this%m_plan_nb) :: Vb,Vbc
          real(TREAL), dimension(3,this%m_plan_nb) :: Vh
          real(TREAL), dimension(3,this%m_plan_nb) :: Ph
          real(TREAL), dimension(6,this%m_plan_nb) :: ell,ell1
          real(TREAL) :: t, mu_helio, mpsbeta
          real(TREAL) :: ac, as
          integer nplan
          
            nplan = this%m_plan_nb
            do j=1, poswrite
             ! lecture dans le buffer d'entree
             call buffer%readdata(j, R)
             t = R(1)
             Ph = reshape(R(2:1+3*nplan), shape(Ph))
             Vb = reshape(R(2+3*nplan:1+6*nplan), shape(Vb))
             
             i = 0
             ! conversion selon le type
             select case (this%m_type_output)             
          case(1)
!           (a,e,I,M,om,Om) elliptiques canoniques'
             do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               mpsbeta = this%m_plan_mpsbeta(i)
               Vbc(:,i) = mpsbeta*Vb(:,i)
               call xyzkhqp1(Ph(:,i),Vbc(:,i),mu_helio,ell1(:,i),        &
     &              this%m_arret)
             end do
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(nplan,ell1,ell)
             end if 
          case(2)
!           (a,e,I,M,om,Om) elliptiques non canoniques'
            call coord_vb2vh(nplan,this%m_plan_mass0, Vb,Vh)
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell1(:,i),         &
     &              this%m_arret)
            end do
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(nplan,ell1,ell)
            end if 
         case(3)
!           (a,la,k,h,q,p) elliptiques canoniques'
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               mpsbeta = this%m_plan_mpsbeta(i)
               Vbc(:,i) = mpsbeta*Vb(:,i)
               call xyzkhqp1(Ph(:,i),Vbc(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
            end do
         case(4)
!           (a,la,k,h,q,p) elliptiques non canoniques'
            call coord_vb2vh(nplan,this%m_plan_mass0, Vb,Vh)
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
            end do

          case(6)
!           (a,e,I,la,pi,Omega) elliptiques non canoniques'
            call coord_vb2vh(nplan,this%m_plan_mass0, Vb,Vh)
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell1(:,i),         &
     &              this%m_arret)
            end do
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_pi(nplan,ell1,ell)
            end if 
            
         case(7)
!           (a*cos(la), a*sin(la),k,h,q,p) elliptiques canoniques'
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               mpsbeta = this%m_plan_mpsbeta(i)
               Vbc(:,i) = mpsbeta*Vb(:,i)
               call xyzkhqp1(Ph(:,i),Vbc(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
               ac = ell(1,i)*cos(ell(2,i))
               as = ell(1,i)*sin(ell(2,i))
               ell(1,i) = ac
               ell(2,i) = as
            end do
         case(8)
!           (a*cos(la), a*sin(la),k,h,q,p) elliptiques non canoniques'
            call coord_vb2vh(nplan,this%m_plan_mass0, Vb,Vh)
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
               ac = ell(1,i)*cos(ell(2,i))
               as = ell(1,i)*sin(ell(2,i))
               ell(1,i) = ac
               ell(2,i) = as
            end do

         case default
            write(*,*) "type de conversion", this%m_type_output
            stop " conversion non supportee vers ell "
         end select
             
             if (this%m_arret%m_stop.ne.0) then 
               call this%m_arret%set_body(i)
               call this%m_arret%set_time(t)
               call buffer%notify_successor_error(this%m_arret)
             endif 

             R(2:1+6*nplan)=reshape(ell,shape(R(2:1+6*nplan)))
             call this%m_buffer%writedata(R)           
            enddo
          
          call buffer%empty()
         
         end subroutine PhVb2ell_oninputbufferfull

!***********************************************************************
! @brief fonction appellee lorsque le buffer d'entree est plein
!! 
!***********************************************************************
         subroutine PhVh2ell_oninputbufferfull(this,buffer,                &
     &                   poswrite) 
          use mod_elliptid
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_plan_PhVh2ell), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          integer i
          real(TREAL), dimension(1+6*this%m_plan_nb) :: R
          real(TREAL), dimension(3,this%m_plan_nb) :: Vb,Vbc
          real(TREAL), dimension(3,this%m_plan_nb) :: Vh
          real(TREAL), dimension(3,this%m_plan_nb) :: Ph
          real(TREAL), dimension(6,this%m_plan_nb) :: ell,ell1
          real(TREAL) :: t, mu_helio, mpsbeta
          real(TREAL) :: ac, as
          integer nplan
          
            nplan = this%m_plan_nb
            do j=1, poswrite
             ! lecture dans le buffer d'entree
             call buffer%readdata(j, R)
             t = R(1)
             Ph = reshape(R(2:1+3*nplan), shape(Ph))
             Vh = reshape(R(2+3*nplan:1+6*nplan), shape(Vh))
             
             i = 0
             ! conversion selon le type
             select case (this%m_type_output)             
          case(1)
!           (a,e,I,M,om,Om) elliptiques canoniques'
             call coord_vh2vb(nplan,this%m_plan_mass0, Vh,Vb)
             do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               mpsbeta = this%m_plan_mpsbeta(i)
               Vbc(:,i) = mpsbeta*Vb(:,i)
               call xyzkhqp1(Ph(:,i),Vbc(:,i),mu_helio,ell1(:,i),        &
     &              this%m_arret)
             end do
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(nplan,ell1,ell)
             end if 
          case(2)
!           (a,e,I,M,om,Om) elliptiques non canoniques'
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell1(:,i),         &
     &              this%m_arret)
            end do
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(nplan,ell1,ell)
            end if 
         case(3)
!           (a,la,k,h,q,p) elliptiques canoniques'
            call coord_vh2vb(nplan,this%m_plan_mass0, Vh,Vb)
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               mpsbeta = this%m_plan_mpsbeta(i)
               Vbc(:,i) = mpsbeta*Vb(:,i)
               call xyzkhqp1(Ph(:,i),Vbc(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
            end do
         case(4)
!           (a,la,k,h,q,p) elliptiques non canoniques'
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
            end do

          case(6)
!           (a,e,I,la,pi,Omega) elliptiques non canoniques'
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell1(:,i),         &
     &              this%m_arret)
            end do
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_pi(nplan,ell1,ell)
            end if 
            
         case(7)
!           (a*cos(la), a*sin(la),k,h,q,p) elliptiques canoniques'
            call coord_vh2vb(nplan,this%m_plan_mass0, Vh,Vb)
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               mpsbeta = this%m_plan_mpsbeta(i)
               Vbc(:,i) = mpsbeta*Vb(:,i)
               call xyzkhqp1(Ph(:,i),Vbc(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
               ac = ell(1,i)*cos(ell(2,i))
               as = ell(1,i)*sin(ell(2,i))
               ell(1,i) = ac
               ell(2,i) = as
            end do
         case(8)
!           (a*cos(la), a*sin(la),k,h,q,p) elliptiques non canoniques'
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
               ac = ell(1,i)*cos(ell(2,i))
               as = ell(1,i)*sin(ell(2,i))
               ell(1,i) = ac
               ell(2,i) = as
            end do

         case default
            write(*,*) "type de conversion", this%m_type_output
            stop " conversion non supportee vers ell "
         end select
             
             if (this%m_arret%m_stop.ne.0) then 
               call this%m_arret%set_body(i)
               call this%m_arret%set_time(t)
               call buffer%notify_successor_error(this%m_arret)
             endif 

             R(2:1+6*nplan)=reshape(ell,shape(R(2:1+6*nplan)))
             call this%m_buffer%writedata(R)           
            enddo
          
          call buffer%empty()
         
         end subroutine PhVh2ell_oninputbufferfull
     
!***********************************************************************
!> @brief PhVb2inv
!! conversion de position heliocentrique/vitesse barycentrique 
!!  vers son propre repere invariant
!!  
!***********************************************************************
         subroutine PhVb2inv(m_plan_nb,m_plan_mass0,Ph,Vb,Phinv,Vbinv)
          implicit none
          integer, intent(in) :: m_plan_nb !< nombre de planetes
          real(TREAL),dimension(0:m_plan_nb),intent(in)::m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          real(TREAL),dimension(3,m_plan_nb), intent(in) :: Ph !< position heliocentrique exprimee dans un repere quelconque
          real(TREAL),dimension(3,m_plan_nb), intent(in) :: Vb !< vitesse barycentrique exprimee dans un repere quelconque
          real(TREAL),dimension(3,m_plan_nb),intent(out):: Phinv !< position heliocentrique exprimee dans le repere inariant du systeme
          real(TREAL),dimension(3,m_plan_nb),intent(out):: Vbinv !< vitesse barycentrique exprimee dans le repere inariant du systeme
          
          real(TREAL), dimension(3) :: C
          call mon_cin2(m_plan_nb,m_plan_mass0,Ph,Vb,C)
          write(*,*)'PhVb2inv C'
          write(*,*) C
          call pass_invar(m_plan_nb,Ph,Vb,C,Phinv,Vbinv)
          write(*,*)'PhVb2inv Phinv'
          write(*,*) Phinv
          write(*,*)'PhVb2inv Vbinv'
          write(*,*) Vbinv

         end subroutine PhVb2inv

!***********************************************************************
!> @brief CI2PhVh
!! conversion de c.i. vers position heliocentrique/vitesse helio 
!!  
!***********************************************************************
         subroutine CI2PhVh(nplan,m_plan_mass0,cG, ci_type,CI,Ph,Vh)
          use mod_elliptid
          implicit none
          integer, intent(in) :: nplan !< nombre de planetes
          real(TREAL),dimension(0:nplan),intent(in)::m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          !> type de donnees dans CI
          !! 1 = (a,e,I,M,om,Om) elliptiques canoniques (angles en radians)
          !! 2 = (a,e,I,M,om,Om) elliptiques non canoniques (angles en radians)
          !! 3 = (a,la,k,h,q,p) elliptiques canoniques (angles en radians)
          !! 4 = (a,la,k,h,q,p) elliptiques non canoniques (angles en radians)
          !! 5 = (Ph, Vh) positions-vitesses helicentriques cartesiennes         
          integer, intent(in) :: ci_type
          real(TREAL), intent(in) :: cG !< constante de gauss en UA, an
          real(TREAL), dimension(:,:), intent(in) :: CI !< exprime dans un repere quelconque
          real(TREAL),dimension(3,nplan),intent(out):: Ph !< position heliocentrique
          real(TREAL),dimension(3,nplan),intent(out):: Vh !< vitesse heliocentrique
          
          real(TREAL) ::mu_helio_i
          real(TREAL),dimension(6,nplan) :: ell_om,ell_kh
           
          integer i 

      write(*,*)  ' entree CI2PhVh'
      select case (ci_type) 
      case(1)
         write(*,*) '(a,e,I,M,om,Om) elliptiques canoniques'
         ell_om(1:2,:) = CI(1:2,:)
         ell_om(3:6,:) = CI(3:6,:)
         call ell_om2ell_kh(nplan,ell_om,ell_kh)
         call ell_hcan2PhVh(nplan,m_plan_mass0,cG,ell_kh,Ph,Vh)
      case(2)
         write(*,*) '(a,e,I,M,om,Om) elliptiques non canoniques'
         ell_om(1:2,:) = CI(1:2,:)
         ell_om(3:6,:) = CI(3:6,:)
         call ell_om2ell_kh(nplan,ell_om,ell_kh)
         do i=1,nplan
            mu_helio_i = cG*(m_plan_mass0(i)+m_plan_mass0(0))
            call ellipx1(ell_kh(:,i),mu_helio_i,Ph(:,i),Vh(:,i))
         end do
      case(3)
          write(*,*) '(a,la,k,h,q,p) elliptiques canoniques'
          ell_kh(:,:) = CI(:,:)
          call ell_hcan2PhVh(nplan,m_plan_mass0,cG,ell_kh,Ph,Vh)
      case(4)
          write(*,*) '(a,la,k,h,q,p) elliptiques non canoniques'
          ell_kh(:,:) = CI(:,:)
         do i=1,nplan
            mu_helio_i = cG*(m_plan_mass0(i)+m_plan_mass0(0))
            call ellipx1(ell_kh(:,i),mu_helio_i,Ph(:,i),Vh(:,i))
         end do
      case(5)
         write(*,*) 'Conditions initiales cartesiennes'
         Ph(:,:) =  CI(1:3,:)
         Vh(:,:) =  CI(4:6,:)
      case default
         write(*,*) 'ci_type =',ci_type 
         stop
      end select
      
      write(*,*) ' sortie CI2PhVh'

         end subroutine CI2PhVh

!***********************************************************************
!> @brief PhVh2PVj
!! conversion de position heliocentrique/vitesse heliocentrique
!!  en position/vitesse jacobi 
!!  
!***********************************************************************
         subroutine PhVh2PjVj(nplan,m_plan_mass0, Ph,Vh, Pj, Vj)
          implicit none
          integer, intent(in) :: nplan !< nombre de planetes
          real(TREAL),dimension(0:nplan),intent(in)::m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          real(TREAL),dimension(3,nplan),intent(in):: Ph !< position helio
          real(TREAL),dimension(3,nplan),intent(in):: Vh !< vitesse helio
          real(TREAL),dimension(3,nplan),intent(out):: Pj !< position jacobi
           real(TREAL),dimension(3,nplan),intent(out):: Vj !< vitesse jacobi
         real(TREAL),dimension(0:nplan) :: eta 
          integer i 
          
          eta(0) = m_plan_mass0(0) 
          do i=1,nplan
            eta(i) = eta(i-1) + m_plan_mass0(i)
          end do 
          
          call coord_h2j(nplan,m_plan_mass0,eta,Ph,Pj)
          call coord_h2j(nplan,m_plan_mass0,eta,Vh,Vj)

         end subroutine PhVh2PjVj

      end module mod_converter


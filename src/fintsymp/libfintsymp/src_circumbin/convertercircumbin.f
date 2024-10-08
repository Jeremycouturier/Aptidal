!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file convertercircumbin.f 
!!  \brief conversion de coordonnees entre 2 buffers pour les systemes circum-binaires
!!
!!    --------------   conversion   --------------\n
!!    | buffer src |  ------------> | buffer dst |\n
!!    --------------                --------------\n
!!
! history : creation 19/11/2015
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la conversion de coordonnees entre 2 buffers
!***********************************************************************
      module mod_convertercircumbin
       use mod_converter
       
!***********************************************************************
!> @class t_circum_VVhat2PhVh
!! conversion de position/vitesse circum-binaire 
!!  en position heliocentrique/vitesse heliocentrique
!!  pour des planetes
!!
!! m_buffer  de t_plan_converter :  buffer de sortie contenant les positions heliocentriques/vitesses heliocentriques
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_circum_VVhat2PhVh
          private
          real(TREAL), dimension(:), allocatable :: m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          integer :: m_plan_nb !< nombre de planetes
          
       contains
          
          procedure :: set_mass =>circum_VVhat2PhVh_set_mass ! fixe les masses a utiliser
          procedure :: set_output => circum_VVhat2PhVh_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein
         
          procedure::oninputbufferfull => VVhat2PhVh_oninputbufferfull ! procedure appellee lorsque le buffer est plein
                    
      end type t_circum_VVhat2PhVh  

!***********************************************************************
!> @class t_circum_VVhat2ell
!! conversion de coordonnees circum-binaire (position, vitesse) = (V, V^) 
!!  en elements elliptiques circum-binaire canoniques ou non canoniques du type
!!  (a,e,I,M,om,Om) ou (a,la,k,h,q,p)  \n
!!  pour la 2nde etoile et les planetes \n
!! le type de coordonnees en sortie est fixe par set_kind
!!
!! m_buffer  de t_plan_converter :  buffer de sortie contenant les elements elliptiques
!!  
!***********************************************************************
      type, extends(t_plan_converter) :: t_circum_VVhat2ell
          !private
          real(TREAL), dimension(:), allocatable :: m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          real(TREAL), dimension(:), allocatable :: m_plan_mpsbeta !< (m_p+m_star)/m_star pour chaque  planete
          real(TREAL), dimension(:), allocatable :: m_mu !< mu(1) = G(somme masses 2 etoiles) , mu(2:) = G*(masses 2 etoiles + masse de la planete)
          integer :: m_plan_nb !< nombre de planetes +1 (2nde etoile incluse)
          !> type de coordonnee du buffer en sortie 
          !! 11 = (a,e,I,M,om,Om) elliptiques canoniques
          !! 13 = (a,la,k,h,q,p) elliptiques canoniques
          integer :: m_type_output = 1
          
       contains
          
          procedure :: set_kind =>circum_VVhat2ell_set_kind ! fixe le type de coordonnees en sortie
          procedure :: set_mass =>circum_VVhat2ell_set_mass ! fixe les masses a utiliser
          procedure :: set_output => circum_VVhat2ell_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein
         
          procedure::oninputbufferfull => circum_VVhat2ell_onibufferfull ! procedure appellee lorsque le buffer est plein
                    
      end type t_circum_VVhat2ell  

      contains

!***********************************************************************
!***********************************************************************
!***********************************************************************
! implementation
!***********************************************************************
!***********************************************************************
!***********************************************************************

!***********************************************************************
!> @brief fixe la masse de l'etoile et des planetes a utiliser pour la conversion 
!***********************************************************************
         subroutine circum_VVhat2PhVh_set_mass(this, m0, mpl)
          implicit none
          class(t_circum_VVhat2PhVh), intent(inout):: this  !< dummy argument
          real(TREAL), intent(in) :: m0  !< masse de la 1ere etoile
          real(TREAL), dimension(:), intent(in) :: mpl  !< masse de la 2nde etoile (en (1)) et des planetes (en (2:))
          integer nplan1
          
          nplan1 = size(mpl)
          this%m_plan_nb =nplan1
          allocate(this%m_plan_mass0(0:nplan1))
          this%m_plan_mass0(0) = m0
          this%m_plan_mass0(1:) = mpl

           call this%setconsumer("vv^2PhVh")

        end  subroutine circum_VVhat2PhVh_set_mass

!***********************************************************************
!> @brief fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein 
!***********************************************************************
         subroutine circum_VVhat2PhVh_set_output(this,buffersize,           &
     &                          buffercs)
          implicit none
          class(t_circum_VVhat2PhVh), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%m_buffer%init(buffersize, 1+6*this%m_plan_nb,        &
     &              buffercs) 
          
         end  subroutine circum_VVhat2PhVh_set_output


     
!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
         subroutine VVhat2PhVh_oninputbufferfull(this, buffer,             &
     &                   poswrite) 
          use mod_coordcircumbin
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_circum_VVhat2PhVh), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL), dimension(1+6*this%m_plan_nb) :: R
          real(TREAL), dimension(3,this%m_plan_nb) :: V,Vhat
          real(TREAL), dimension(3,this%m_plan_nb) :: Ph,Vh
          integer nplan
          
            nplan = this%m_plan_nb
            do j=1, poswrite
             call buffer%readdata(j, R)
             V = reshape(R(2:1+3*nplan), shape(V))
             Vhat = reshape(R(2+3*nplan:1+6*nplan), shape(Vh))
             call circum_vvhat2PhVh(nplan,this%m_plan_mass0,V,Vhat,       &
     &               Ph,Vh)
             R(2:1+3*nplan)=reshape(Ph, shape(R(2:1+3*nplan)))
             R(2+3*nplan:1+6*nplan)=reshape(Vh,                           &
     &               shape(R(2+3*nplan:1+6*nplan)))
             call this%m_buffer%writedata(R)           
            enddo
          
          call buffer%empty()
         
         end subroutine VVhat2PhVh_oninputbufferfull



!***********************************************************************
!> @brief fixe le nom pour le graphique dot 
!***********************************************************************
         subroutine circum_VVhat2ell_set_graphdotname(this)
          implicit none
          class(t_circum_VVhat2ell), intent(inout):: this  !< dummy argument

          character (len=40), dimension(11:18) :: dotname
          dotname(11) ="(v,v^)_circum2(a,e,I,M,om,Om)_canon"
          dotname(12) ="(v,v^)_circum2(a,e,I,M,om,Om)_non_canon"
          dotname(13) ="(v,v^)_circum2(a,la,k,h,q,p)_canon"
          dotname(14) ="(v,v^)_circum2(a,la,k,h,q,p)_non_canon"
          dotname(15) ="??"
          dotname(16) ="??"
          dotname(17) ="PhVb2(a*cos(la),a*sin(la),k,h,q,p)_canon"
          dotname(18) ="??"
          
          call this%set_graphdotname(dotname(this%m_type_output))
           
         end  subroutine circum_VVhat2ell_set_graphdotname
         
         
!***********************************************************************
!> @brief fixe le type de coordonnee voulue en sortie 
!***********************************************************************
         subroutine circum_VVhat2ell_set_kind(this, itype)
          implicit none
          class(t_circum_VVhat2ell), intent(inout):: this  !< dummy argument
          !> type de coordonnee en sortie du buffer
          !! 11 = (a,e,I,M,om,Om) elliptiques canoniques : associee a (V,Vhat)
          !! 12 = (a,e,I,M,om,Om) elliptiques non canoniques : associee a (V,Vdot)
          !! 13 = (a,la,k,h,q,p) elliptiques canoniques : associee a (V,Vhat)
          !! 14 = (a,la,k,h,q,p) elliptiques non canoniques : associee a (V,Vdot)
          !! 17 = (a*cos(la),a*sin(la),k,h,q,p) elliptiques canoniques
          integer, intent(in) :: itype  
          
          this%m_type_output = itype
           call circum_VVhat2ell_set_graphdotname(this)
         
         end  subroutine circum_VVhat2ell_set_kind

!***********************************************************************
!> @brief fixe la constante de gauss, masse des etoiles et des planetes a utiliser pour la conversion 
!***********************************************************************
         subroutine circum_VVhat2ell_set_mass(this, cG, m0, mpl)
          use mod_coordcircumbin
          implicit none
          class(t_circum_VVhat2ell), intent(inout):: this  !< dummy argument
          real(TREAL), intent(in) :: cG          !< constante de Gauss AU**3/an**2
          real(TREAL), intent(in) :: m0  !< masse de la 1ere etoile
          real(TREAL), dimension(:), intent(in) :: mpl  !< masse de la 2nde etoile (en (1)) et des planetes (en (2:))
          integer nplan1
          nplan1 = size(mpl)
          this%m_plan_nb =nplan1
          allocate(this%m_plan_mass0(0:nplan1))
          allocate(this%m_mu(nplan1))
          allocate(this%m_plan_mpsbeta(nplan1))
          this%m_plan_mass0(0) = m0
          this%m_plan_mass0(1:) = mpl
          call circum_calc_mu_mpsbeta(nplan1,cG, this%m_plan_mass0,           &
     &           this%m_mu,this%m_plan_mpsbeta)
                   
           call this%setconsumer("")
           call circum_VVhat2ell_set_graphdotname(this)

        end  subroutine circum_VVhat2ell_set_mass

!***********************************************************************
!> @brief fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein 
!***********************************************************************
         subroutine circum_VVhat2ell_set_output(this,buffersize,           &
     &                          buffercs)
          implicit none
          class(t_circum_VVhat2ell), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          call this%m_buffer%init(buffersize, 1+6*this%m_plan_nb,        &
     &              buffercs) 
          
         end  subroutine circum_VVhat2ell_set_output


!***********************************************************************
! @brief fonction appellee lorsque le buffer d'entree est plein
!! 
!***********************************************************************
         subroutine circum_VVhat2ell_onibufferfull(this,buffer,                &
     &                   poswrite) 
          use mod_elliptid
          use mod_coordcircumbin
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_circum_VVhat2ell), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL), dimension(1+6*this%m_plan_nb) :: R
          real(TREAL), dimension(3,this%m_plan_nb) :: V,Vhat
          real(TREAL), dimension(6,this%m_plan_nb) :: ell,ell1
          real(TREAL) :: t
          real(TREAL) :: ac, as
          integer nplan, i

            nplan = this%m_plan_nb
            do j=1, poswrite
             ! lecture dans le buffer d'entree
             call buffer%readdata(j, R)
             t = R(1)
             V = reshape(R(2:1+3*nplan), shape(V))
             Vhat = reshape(R(2+3*nplan:1+6*nplan), shape(Vhat))
             
             ! conversion selon le type
             select case (this%m_type_output)             
          case(11)
!           (a,e,I,M,om,Om) elliptiques canoniques'
             call circum_vvhat2ell_hcan(nplan, this%m_mu,                &
     &         this%m_plan_mpsbeta, V, Vhat,ell1, this%m_arret)
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(nplan,ell1,ell)
             end if 
             
          case(12)
!           (a,e,I,M,om,Om) elliptiques non canoniques'
             call circum_vvhat2ell_h(nplan,this%m_plan_mass0,this%m_mu, &
     &         this%m_plan_mpsbeta, V, Vhat,ell1, this%m_arret)
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(nplan,ell1,ell)
             end if 
             
         case(13)
!           (a,la,k,h,q,p) elliptiques canoniques'
             call circum_vvhat2ell_hcan(nplan, this%m_mu,                &
     &          this%m_plan_mpsbeta, V, Vhat,ell, this%m_arret)

         case(14)
!           (a,la,k,h,q,p) elliptiques non canoniques'
             call circum_vvhat2ell_h(nplan,this%m_plan_mass0,this%m_mu,  &
     &          this%m_plan_mpsbeta, V, Vhat,ell, this%m_arret)

         case(17)
!           (a*cos(la), a*sin(la),k,h,q,p) elliptiques canoniques'
             call circum_vvhat2ell_hcan(nplan, this%m_mu,                &
     &          this%m_plan_mpsbeta, V, Vhat,ell, this%m_arret)
             if (this%m_arret%m_stop.eq.0) then 
            do i=1,nplan
               ac = ell(1,i)*cos(ell(2,i))
               as = ell(1,i)*sin(ell(2,i))
               ell(1,i) = ac
               ell(2,i) = as
            end do
             end if 

         case default
            write(*,*) "type de conversion", this%m_type_output
            stop " conversion non supportee vers ell "
         end select
             
             if (this%m_arret%m_stop.ne.0) then 
               call this%m_arret%set_time(t)
               call buffer%notify_successor_error(this%m_arret)
             endif 

             R(2:1+6*nplan)=reshape(ell,shape(R(2:1+6*nplan)))
             call this%m_buffer%writedata(R)           
            enddo

          call buffer%empty()
        
         end subroutine circum_VVhat2ell_onibufferfull

!***********************************************************************
!> @brief circum_VVhat2inv
!! conversion de  coordonnees circum-binaires 
!!  vers son propre repere invariant pour les planetes et etoiles
!!  
!***********************************************************************
         subroutine circum_VVhat2inv(m_plan_nb,m_plan_mass0,V,Vhat,Vinv,   &
     &                Vhatinv)
          implicit none
          integer, intent(in) :: m_plan_nb !< nombre de planetes+1
          real(TREAL),dimension(0:m_plan_nb),intent(in)::m_plan_mass0 !< masses des planetes en (2:) + masse des etoiles en (0:1)
          real(TREAL),dimension(3,m_plan_nb), intent(in) :: V !< position de la 2nde etoile en (1) et des  planetes en (2:) en coordonnees circum-binaires exprimee dans un repere quelconque
          real(TREAL),dimension(3,m_plan_nb), intent(in) :: Vhat !< vitesse de la 2nde etoile en (1) et des  planetes en (2:) en coordonnees circum-binaires exprimee dans un repere quelconque
          real(TREAL),dimension(3,m_plan_nb), intent(out) :: Vinv !< position de la 2nde etoile en (1) et des  planetes en (2:) en coordonnees circum-binaires exprimee dans le repere inariant du systeme
          real(TREAL),dimension(3,m_plan_nb), intent(out) :: Vhatinv !< vitesse de la 2nde etoile en (1) et des  planetes en (2:) en coordonnees circum-binaires exprimee dans le repere inariant du systeme
          
          real(TREAL), dimension(3) :: C
          call mon_cin2(m_plan_nb,m_plan_mass0,V,Vhat,C)
          write(*,*)'circum_VVhat2inv C'
          write(*,*) C
          call pass_invar(m_plan_nb,V,Vhat,C,Vinv,Vhatinv)
          write(*,*)'circum_VVhat2inv Vinv'
          write(*,*) Vinv
          write(*,*)'circum_VVhat2inv Vhatinv'
          write(*,*) Vhatinv

         end subroutine circum_VVhat2inv

!***********************************************************************
!> @brief circum_CI2VVhat
!! conversion de c.i. des 2 etoiles et des planetes 
!! vers (V,V^)
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!!  
!***********************************************************************
         subroutine circum_CI2VVhat(nplan,m_plan_mass0,cG,ci_type,          &
     &    CI,V,Vhat)
          use mod_coordcircumbin
          use mod_elliptid
          implicit none
          integer, intent(in) :: nplan !< nombre de planetes + 1 (car 2nde etoile en position 1)
          real(TREAL),dimension(0:nplan),intent(in)::m_plan_mass0 !< masses des planetes en (2:) + masse des etoiles en (0:1)
          !> type de donnees dans CI et CIpart
          !! 11 = (a,e,I,M,om,Om) elliptiques canoniques (angles en radians)
          !! 13 = (a,la,k,h,q,p) elliptiques canoniques (angles en radians)
          integer, intent(in) :: ci_type
          real(TREAL), intent(in) :: cG !< constante de gauss en UA, an
          real(TREAL), dimension(:,:), intent(in) :: CI !< planetes exprime dans un repere quelconque
          real(TREAL),dimension(3,nplan),intent(out):: V !< position de la 2nde etoile et des planetes en coordonnees circum-binaires
          real(TREAL),dimension(3,nplan),intent(out):: Vhat !< vitesse de la 2nde etoile et des planetes en coordonnees circum-binaires
          
          real(TREAL),dimension(6,nplan) :: ell_om,ell_kh
          real(TREAL),dimension(3,nplan) :: Ph, Vh
          
      write(*,*)  ' entree circum_CI2VVhat'
      select case (ci_type) 
      case(11)
         write(*,*) '(a,e,I,M,om,Om) circum elliptiques canoniques'
         ell_om(1:6,:) = CI(1:6,:)
         call ell_om2ell_kh(nplan,ell_om,ell_kh)
         call circum_ell_hcan2VVhat(nplan,m_plan_mass0,cG,ell_kh,V,Vhat)
      case(12)
         write(*,*) '(a,e,I,M,om,Om) circum elliptiques non canoniques'
         ell_om(1:6,:) = CI(1:6,:)
         call ell_om2ell_kh(nplan,ell_om,ell_kh)
         call circum_ell_h2VVhat(nplan,m_plan_mass0,cG,ell_kh,V,Vhat)
         
      case(13)
         write(*,*) '(a,la,k,h,q,p) circum elliptiques canoniques'
         ell_kh(:,:) = CI(:,:)
         call circum_ell_hcan2VVhat(nplan,m_plan_mass0,cG,ell_kh,V,Vhat)
         
      case(14)
         write(*,*) '(a,la,k,h,q,p) circum elliptiques non canoniques'
         ell_kh(:,:) = CI(:,:)
         call circum_ell_h2VVhat(nplan,m_plan_mass0,cG,ell_kh,V,Vhat)
         
      case(5)
         write(*,*) 'Conditions initiales cartesiennes Ph/Vh'
         Ph(:,:) =  CI(1:3,:)
         Vh(:,:) =  CI(4:6,:)
         call circum_PhVh2vvhat(nplan,m_plan_mass0,Ph,Vh,V,Vhat)
       case default
         write(*,*) 'ci_type =',ci_type 
         stop
      end select
      
      write(*,*) ' sortie circum_CI2VVhat'

         end subroutine circum_CI2VVhat

      end module mod_convertercircumbin


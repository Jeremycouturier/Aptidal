!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file converterpart.f 
!!  \brief conversion de coordonnees entre 2 buffers pour les particules
!!
!!    --------------   conversion   --------------\n
!!    | buffer src |  ------------> | buffer dst |\n
!!    --------------                --------------\n
!!
! history : creation 26/08/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la conversion de coordonnees entre 2 buffers
!***********************************************************************
      module mod_converterpart
       use mod_converter
       use mod_coordpart
       
!***********************************************************************
!> @class t_part_PhVb2PhVh
!! conversion de position heliocentrique/vitesse barycentrique 
!!  en position heliocentrique/vitesse heliocentrique
!!  pour des particules
!!  
!***********************************************************************
      type, extends(t_plan_PhVb2PhVh) :: t_part_PhVb2PhVh
          private
          integer :: m_part_nb !< nombre de particules
          
       contains
          
          procedure:: set_nb => part_PhVb2PhVh_set_nb ! fixe le nombre de particules 

          procedure :: set_output => part_PhVb2PhVh_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein
         
          procedure::oninputbufferfull =>                                &
     &              part_PhVb2PhVh_oninputbufferfull ! procedure appellee lorsque le buffer est plein
                    
      end type t_part_PhVb2PhVh  
   
!***********************************************************************
!> @class t_part_PhVb2ell
!! conversion de position heliocentrique/vitesse barycentrique 
!!  en elements elliptiques canoniques ou non canoniques du type
!!  (a,e,I,M,om,Om) ou (a,la,k,h,q,p)  \n
!!  pour des particules \n
!! le type de coordonnees en sortie est fixe par set_kind
!!
!!  
!***********************************************************************
      type, extends(t_plan_PhVb2ell) :: t_part_PhVb2ell
          private
          integer :: m_part_nb !< nombre de particules
          real(TREAL) ::   m_part_mu_helio !< mu heliocentrique de particules
          
       contains
          
          procedure:: set_nb => part_PhVb2ell_set_nb    ! fixe le nombre de particules 
          procedure:: set_mass =>part_PhVb2ell_set_mass ! fixe les masses a utiliser
          procedure:: set_output => part_PhVb2ell_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein
          procedure:: oninputbufferfull =>                               &
     &               part_PhVb2ell_oninputbufferfull    ! procedure appellee lorsque le buffer est plein
          
          
      end type t_part_PhVb2ell  

      contains

!***********************************************************************
!***********************************************************************
!***********************************************************************
! implementation
!***********************************************************************
!***********************************************************************
!***********************************************************************

!***********************************************************************
!> @brief fixe le nombre de particules 
!***********************************************************************
      subroutine part_PhVb2PhVh_set_nb(this, npart)
       implicit none
       class(t_part_PhVb2PhVh), intent(inout):: this  !< dummy argument
       integer, intent(in) :: npart                !< nombre de particules
       
       this%m_part_nb = npart
       call this%m_buffer%init_multictrlerr(npart)

      end  subroutine part_PhVb2PhVh_set_nb  

!***********************************************************************
!> @brief fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein 
!***********************************************************************
         subroutine part_PhVb2PhVh_set_output(this,buffersize,           &
     &                          buffercs)
          implicit none
          class(t_part_PhVb2PhVh), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          integer ntot6
          
          ntot6=6*this%m_plan_nb+6*this%m_part_nb
          call this%m_buffer%init(buffersize, 1+ntot6, buffercs) 

         end  subroutine part_PhVb2PhVh_set_output


     
!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
         subroutine part_PhVb2PhVh_oninputbufferfull(this,buffer,        &
     &         poswrite) 
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_part_PhVb2PhVh), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          real(TREAL),dimension(1+6*this%m_plan_nb+6*this%m_part_nb)::R
          real(TREAL), dimension(3,this%m_plan_nb) :: Vb
          real(TREAL), dimension(3,this%m_plan_nb) :: Vh
          real(TREAL), dimension(3,this%m_part_nb) :: Vbpart
          real(TREAL), dimension(3,this%m_part_nb) :: Vhpart
          integer nplan, npart, ntot
          integer ideb, ifin
          
            nplan = this%m_plan_nb
            npart = this%m_part_nb
            ntot = nplan+npart
            
            call this%m_buffer%set_multictrlerr_buf(buffer)
            do j=1, poswrite
             call buffer%readdata(j, R)
             
             ! pour les planetes
             ideb = 2+3*nplan
             ifin = 1+6*nplan
             Vb = reshape(R(ideb:ifin), shape(Vb))
             call coord_vb2vh(nplan,this%m_plan_mass0, Vb,Vh)
             R(ideb:ifin) = reshape(Vh, shape(R(ideb:ifin)))
             
             ! pour les particules
             ideb = 2+3*npart+6*nplan
             ifin = 1+6*npart+6*nplan
             Vbpart = reshape(R(ideb:ifin), shape(Vbpart))
             call coord_vb2vh_part(nplan,this%m_plan_mass0,Vb,          &
     &               npart, Vbpart, Vhpart)
             R(ideb:ifin) = reshape(Vhpart, shape(R(ideb:ifin)))

             call this%m_buffer%writedata(R)           
            enddo
            call buffer%set_multictrlerr_buf(this%m_buffer)
          
          call buffer%empty()
         
         end subroutine part_PhVb2PhVh_oninputbufferfull

!***********************************************************************
!> @brief fixe le nombre de particules 
!***********************************************************************
      subroutine part_PhVb2ell_set_nb(this, npart)
       implicit none
       class(t_part_PhVb2ell), intent(inout):: this  !< dummy argument
       integer, intent(in) :: npart                !< nombre de particules
       
       this%m_part_nb = npart
       call this%m_buffer%init_multictrlerr(npart)

      end  subroutine part_PhVb2ell_set_nb  

!***********************************************************************
!> @brief fixe la constante de gauss, masse de l'etoile et des planetes a utiliser pour la conversion 
!***********************************************************************
        subroutine part_PhVb2ell_set_mass(this, cG, m0, mpl)
          implicit none
          class(t_part_PhVb2ell), intent(inout):: this  !< dummy argument
          real(TREAL), intent(in) :: cG          !< constante de Gauss AU**3/an**2
          real(TREAL), intent(in) :: m0  !< masse de l'etoile
          real(TREAL), dimension(:), intent(in) :: mpl  !< masse des planetes
          
          call plan_PhVb2ell_set_mass(this, Cg, m0, mpl)
          this%m_part_mu_helio = cG*m0

        end  subroutine part_PhVb2ell_set_mass

!***********************************************************************
!> @brief fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein 
!***********************************************************************
         subroutine part_PhVb2ell_set_output(this,buffersize,           &
     &                          buffercs)
          implicit none
          class(t_part_PhVb2ell), intent(inout):: this  !< dummy argument
          integer(8), intent(in) :: buffersize        !< nombre de pas (pouvant stocker le buffer de sortie) avant d'appeller la fonction usercallback
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          integer ntot6
          
          ntot6=6*this%m_plan_nb+6*this%m_part_nb
          call this%m_buffer%init(buffersize, 1+ntot6, buffercs) 
          
         end  subroutine part_PhVb2ell_set_output


     
!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
         subroutine part_PhVb2ell_oninputbufferfull(this,buffer,        &
     &         poswrite) 
          use mod_elliptid
          implicit none
          class(t_buffer), intent(inout) :: buffer !< buffer de sortie
          class(t_part_PhVb2ell), intent(inout) :: this    !< donnee utilisateur 
          integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
          integer(8) :: j
          integer i
          real(TREAL),dimension(1+6*this%m_plan_nb+6*this%m_part_nb)::R
          real(TREAL), dimension(6,this%m_plan_nb) :: ell,ell1
          real(TREAL), dimension(3,this%m_plan_nb) :: Vb, Vbc
          real(TREAL), dimension(3,this%m_plan_nb) :: Ph
          real(TREAL), dimension(3,this%m_plan_nb) :: Vh
          real(TREAL), dimension(6,this%m_part_nb) :: ellpart,ell1part
          real(TREAL), dimension(3,this%m_part_nb) :: Vbpart
          real(TREAL), dimension(3,this%m_part_nb) :: Phpart
          real(TREAL), dimension(3,this%m_part_nb) :: Vhpart
          real(TREAL) :: t, mu_helio, mpsbeta, mu_helio_part
          real(TREAL) :: ac, as
          integer nplan, npart, ntot
          integer ideb, ifin
          
            nplan = this%m_plan_nb
            npart = this%m_part_nb
            ntot = nplan+npart
            mu_helio_part = this%m_part_mu_helio

            call this%m_buffer%set_multictrlerr_buf(buffer)
            do j=1, poswrite
             call buffer%readdata(j, R)
             
             t=R(1)
             ! pour les planetes
             ideb = 2
             ifin = 1+3*nplan
             Ph = reshape(R(ideb:ifin), shape(Ph))
             ideb = 2+3*nplan
             ifin = 1+6*nplan
             Vb = reshape(R(ideb:ifin), shape(Vb))
             
             ! pour les particules
             ideb = 2+6*nplan
             ifin = 1+3*npart+6*nplan
             Phpart = reshape(R(ideb:ifin), shape(Phpart))
             ideb = 2+3*npart+6*nplan
             ifin = 1+6*npart+6*nplan
             Vbpart = reshape(R(ideb:ifin), shape(Vbpart))

             i = 0
             ! conversion selon le type
             select case (this%m_type_output)             
          case(1)
!           (a,e,I,M,om,Om) elliptiques canoniques'
             ! planetes
             do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               mpsbeta = this%m_plan_mpsbeta(i)
               Vbc(:,i) = mpsbeta*Vb(:,i)
               call xyzkhqp1(Ph(:,i),Vbc(:,i),mu_helio,ell1(:,i),         &
     &              this%m_arret)
             end do
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(nplan,ell1,ell)
             end if 
             ! particules
             do i=1, npart
               call xyzkhqp1(Phpart(:,i),Vbpart(:,i),mu_helio_part,     &
     &              ell1part(:,i), this%m_buffer%m_mce(i))
             enddo  
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(npart,ell1part,ellpart)
             end if 

          case(2)
!           (a,e,I,M,om,Om) elliptiques non canoniques'
             ! planetes
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
            
            ! particules
            call coord_vb2vh_part(nplan,this%m_plan_mass0,Vb,           &
     &               npart, Vbpart, Vhpart)
            do i=1, npart
               call xyzkhqp1(Phpart(:,i),Vhpart(:,i),mu_helio_part,     &
     &              ell1part(:,i), this%m_buffer%m_mce(i))
            enddo  
             if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_om(npart,ell1part,ellpart)
            end if
            
         case(3)
!           (a,la,k,h,q,p) elliptiques canoniques'
             ! planetes
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               mpsbeta = this%m_plan_mpsbeta(i)
               Vbc(:,i) = mpsbeta*Vb(:,i)
               call xyzkhqp1(Ph(:,i),Vbc(:,i),mu_helio,ell(:,i),        &
     &              this%m_arret)
            end do
            
            ! particules
            do i=1, npart
               call xyzkhqp1(Phpart(:,i),Vbpart(:,i),mu_helio_part,     &
     &              ellpart(:,i), this%m_buffer%m_mce(i))
            enddo  

         case(4)
!           (a,la,k,h,q,p) elliptiques non canoniques'
            ! planetes
            call coord_vb2vh(nplan,this%m_plan_mass0, Vb,Vh)
            do while ((i.lt.nplan).and.(this%m_arret%m_stop.eq.0))
               i = i + 1
               mu_helio = this%m_plan_mu_helio(i)
               call xyzkhqp1(Ph(:,i),Vh(:,i),mu_helio,ell(:,i),         &
     &              this%m_arret)
            end do
            ! particules
            call coord_vb2vh_part(nplan,this%m_plan_mass0,Vb,           &
     &               npart, Vbpart, Vhpart)
            do i=1, npart
               call xyzkhqp1(Phpart(:,i),Vhpart(:,i),mu_helio_part,     &
     &              ellpart(:,i), this%m_buffer%m_mce(i))
            enddo  
            
          case(6)
!           (a,e,I,la,pi,Omega) elliptiques non canoniques'
            ! planetes
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
            ! particules
            call coord_vb2vh_part(nplan,this%m_plan_mass0,Vb,           &
     &               npart, Vbpart, Vhpart)
            do i=1, npart
               call xyzkhqp1(Phpart(:,i),Vhpart(:,i),mu_helio_part,     &
     &              ell1part(:,i), this%m_buffer%m_mce(i))
            enddo  
            if (this%m_arret%m_stop.eq.0) then 
               call ell_kh2ell_pi(npart,ell1part,ellpart)
            end if 
            
         case(7)
!           (a*cos(la), a*sin(la),k,h,q,p) elliptiques canoniques'
            ! planetes
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
            ! particules
            do i=1, npart
               call xyzkhqp1(Phpart(:,i),Vbpart(:,i),mu_helio_part,     &
     &              ellpart(:,i), this%m_buffer%m_mce(i))
               ac = ellpart(1,i)*cos(ellpart(2,i))
               as = ellpart(1,i)*sin(ellpart(2,i))
               ellpart(1,i) = ac
               ellpart(2,i) = as
            enddo  
         case(8)
!           (a*cos(la), a*sin(la),k,h,q,p) elliptiques non canoniques'
            ! planetes
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
            ! particules
            call coord_vb2vh_part(nplan,this%m_plan_mass0,Vb,           &
     &               npart, Vbpart, Vhpart)
            ! particules
            do i=1, npart
               call xyzkhqp1(Phpart(:,i),Vhpart(:,i),mu_helio_part,     &
     &              ellpart(:,i), this%m_buffer%m_mce(i))
               ac = ellpart(1,i)*cos(ellpart(2,i))
               as = ellpart(1,i)*sin(ellpart(2,i))
               ellpart(1,i) = ac
               ellpart(2,i) = as
            enddo  
         case default
            stop  " conversion non supportee vers ell-part "
            stop 1
         end select
                  
             if (this%m_arret%m_stop.ne.0) then 
               call this%m_arret%set_body(i)
               call this%m_arret%set_time(t)
               call buffer%notify_successor_error(this%m_arret)
             endif 
            
             ! pour les planetes
             ideb = 2
             ifin = 1+6*nplan
             R(ideb:ifin)=reshape(ell,shape(R(ideb:ifin)))
             ! pour les particules
             ideb = 2+6*nplan
             ifin = 1+6*npart+6*nplan
             R(ideb:ifin)=reshape(ellpart,shape(R(ideb:ifin)))

             call this%m_buffer%writedata(R)           
            enddo
            
            call buffer%set_multictrlerr_buf(this%m_buffer)
          
          call buffer%empty()
         
         end subroutine part_PhVb2ell_oninputbufferfull

!***********************************************************************
!> @brief PhVb2inv_part
!! conversion de position heliocentrique/vitesse barycentrique 
!!  vers son propre repere invariant  pour les planetes et particules
!!  
!***********************************************************************
         subroutine PhVb2inv_part(m_plan_nb,m_plan_mass0,Ph,Vb,Phinv,   &
     &                Vbinv,npart, Phpart,Vbpart,Phinvpart,Vbinvpart)
          implicit none
          integer, intent(in) :: m_plan_nb !< nombre de planetes
          real(TREAL),dimension(0:m_plan_nb),intent(in)::m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          real(TREAL),dimension(3,m_plan_nb), intent(in) :: Ph !< position heliocentrique des planetes exprimee dans un repere quelconque
          real(TREAL),dimension(3,m_plan_nb), intent(in) :: Vb !< vitesse barycentrique des planetes exprimee dans un repere quelconque
          real(TREAL),dimension(3,m_plan_nb),intent(out):: Phinv !< position heliocentrique des planetes exprimee dans le repere inariant du systeme
          real(TREAL),dimension(3,m_plan_nb),intent(out):: Vbinv !< vitesse barycentrique des planetes exprimee dans le repere inariant du systeme
          integer, intent(in) :: npart !< nombre de particules
          real(TREAL),dimension(3,npart), intent(in) :: Phpart !< position heliocentrique des particules exprimee dans un repere quelconque
          real(TREAL),dimension(3,npart), intent(in) :: Vbpart !< vitesse barycentrique des particules exprimee dans un repere quelconque
          real(TREAL),dimension(3,npart),intent(out):: Phinvpart !< position heliocentrique des particules exprimee dans le repere inariant du systeme
          real(TREAL),dimension(3,npart),intent(out):: Vbinvpart !< vitesse barycentrique des particules exprimee dans le repere inariant du systeme
          
          real(TREAL), dimension(3) :: C
          call mon_cin2(m_plan_nb,m_plan_mass0,Ph,Vb,C)
          write(*,*)'PhVb2inv_part C'
          write(*,*) C
          call pass_invar(m_plan_nb,Ph,Vb,C,Phinv,Vbinv)
          call pass_invar(npart,Phpart,Vbpart,C,Phinvpart,Vbinvpart)
          write(*,*)'PhVb2inv_part Phinv'
          write(*,*) Phinv
          write(*,*)'PhVb2inv_part Vbinv'
          write(*,*) Vbinv
          write(*,*)'PhVb2inv_part Phinvpart'
          write(*,*) Phinvpart
          write(*,*)'PhVb2inv_part Vbinvpart'
          write(*,*) Vbinvpart

         end subroutine PhVb2inv_part

!***********************************************************************
!> @brief CI2PhVh_part
!! conversion de c.i. des planetes et des particules
!! vers position heliocentrique/vitesse helio pour les planetes et particules
!!  
!***********************************************************************
         subroutine CI2PhVh_part(nplan,m_plan_mass0,cG,ci_type,          &
     &    CI,Ph,Vh, npart, CIpart, Phpart, Vhpart)
          use mod_elliptidpart
          implicit none
          integer, intent(in) :: nplan !< nombre de planetes
          real(TREAL),dimension(0:nplan),intent(in)::m_plan_mass0 !< masses des planetes + masse de l'etoile en 0
          !> type de donnees dans CI et CIpart
          !! 1 = (a,e,I,M,om,Om) elliptiques canoniques (angles en radians)
          !! 2 = (a,e,I,M,om,Om) elliptiques non canoniques (angles en radians)
          !! 3 = (a,la,k,h,q,p) elliptiques canoniques (angles en radians)
          !! 4 = (a,la,k,h,q,p) elliptiques non canoniques (angles en radians)
          !! 5 = (Ph, Vh) positions-vitesses helicentriques cartesiennes         
          integer, intent(in) :: ci_type
          real(TREAL), intent(in) :: cG !< constante de gauss en UA, an
          real(TREAL), dimension(:,:), intent(in) :: CI !< planetes exprime dans un repere quelconque
          real(TREAL),dimension(3,nplan),intent(out):: Ph !< position heliocentrique des planetes
          real(TREAL),dimension(3,nplan),intent(out):: Vh !< vitesse heliocentrique des planetes
          integer, intent(in) :: npart !< nombre de particules
          real(TREAL), dimension(:,:), intent(in) :: CIpart !< particules exprime dans un repere quelconque
          real(TREAL),dimension(3,npart),intent(out):: Phpart !< position heliocentrique des particules
          real(TREAL),dimension(3,npart),intent(out):: Vhpart !< vitesse heliocentrique des particules
          
          real(TREAL) ::mu_helio_i
          real(TREAL),dimension(6,npart) :: ell_om,ell_kh
          real(TREAL),dimension(6,nplan) :: ell_ompla,ell_khpla
           
          integer i 

          call CI2PhVh(nplan,m_plan_mass0,cG,ci_type,CI,Ph,Vh)
          
      write(*,*)  ' entree CI2PhVh_part'
      select case (ci_type) 
      case(1)
         write(*,*) '(a,e,I,M,om,Om) elliptiques canoniques'
         ell_ompla(1:6,:) = CI(1:6,:)
         call ell_om2ell_kh(nplan,ell_ompla,ell_khpla)
         ell_om(1:6,:) = CIpart(1:6,:)
         call ell_om2ell_kh(npart,ell_om,ell_kh)
         call ell_hcan2PhVh_part(nplan,m_plan_mass0,cG,ell_khpla, npart,    &
     &     ell_kh, Phpart,Vhpart)
      case(2)
         write(*,*) '(a,e,I,M,om,Om) elliptiques non canoniques'
         ell_om(1:6,:) = CIpart(1:6,:)
         call ell_om2ell_kh(npart,ell_om,ell_kh)
         mu_helio_i = cG*m_plan_mass0(0)
         do i=1,npart
            call ellipx1(ell_kh(:,i),mu_helio_i,Phpart(:,i),Vhpart(:,i))
         end do
      case(3)
         write(*,*) '(a,la,k,h,q,p) elliptiques canoniques'
         ell_khpla(:,:) = CI(:,:)
         ell_kh(:,:) = CIpart(:,:)
         call ell_hcan2PhVh_part(nplan,m_plan_mass0,cG,ell_khpla, npart,    &
     &     ell_kh, Phpart,Vhpart)
      case(4)
         write(*,*) '(a,la,k,h,q,p) elliptiques non canoniques'
         ell_kh(:,:) = CIpart(:,:)
         mu_helio_i = cG*m_plan_mass0(0)
         do i=1,npart
            call ellipx1(ell_kh(:,i),mu_helio_i,Phpart(:,i),Vhpart(:,i))
         end do
      case(5)
         write(*,*) 'Conditions initiales cartesiennes'
         Phpart(:,:) =  CIpart(1:3,:)
         Vhpart(:,:) =  CIpart(4:6,:)
      case default
         write(*,*) 'ci_type =',ci_type 
         stop
      end select
      
      write(*,*) ' sortie CI2PhVh_part'

         end subroutine CI2PhVh_part

      end module mod_converterpart


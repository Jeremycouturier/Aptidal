!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file sysplaode.f 
!!  \brief systeme planetaires avec 1 etoile et N planetes
!!         pouvant etre integre par un integrateur numerique du type odex ou dopri.
!!
!
! history : creation 25/01/2016
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme a 1 etoile avec N planetes 
! exprime en variables heliocentriques
! les interactions newtonniennes sont prises en compte
!***********************************************************************
      module mod_sysplaodehelio
       use mod_sysplaode
       private energie_phvh, moment

!***********************************************************************
!> @class t_sysplaodehelio
!! classe decrivant 1 systeme planetaire avec 1 etoile et N planetes
!!  les interactions newtonniennes sont prises en compte.
!!  
!***********************************************************************
      type, extends(t_sysplaode) :: t_sysplaodehelio 
     
          real(TREAL), dimension(:), allocatable ::   m_plan_gm !<  = G*(masse de la planete)
          
       contains
          
          procedure :: set_plan_mass => sysplaodehelio_set_plan_mass ! fixe la masse de l' etoiles et des planetes

          procedure :: energie => sysplaodehelio_energie 
          procedure :: mom_cin => sysplaodehelio_mom_cin 
          
          procedure :: allocate_y  => sysplaodehelio_allocate_y ! alloue le vecteur d'etat
          procedure :: fill_from_y => sysplaodehelio_fill_from_y ! recupere les donnees depuis le vecteur d'etat
          procedure :: fill_to_y   => sysplaodehelio_fill_to_y ! envoie les donnees vers le vecteur d'etat
          procedure :: fcn         => sysplaodehelio_fcn ! evalue le second membre pour l'ODE
          
      end type t_sysplaodehelio  
      
      contains
         
!***********************************************************************
!> @brief fixe la masse des planetes 
!!  la fonction set_star_mass doit etre appellee auparavant
!***********************************************************************
      subroutine sysplaodehelio_set_plan_mass(this, m)
       implicit none
       class(t_sysplaodehelio), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:), intent(in) :: m  !< masse des planetes
       real(TREAL) :: m0
       
       call sysplaode_set_plan_mass(this, m)
       call this%set_graphdotname("SYSPLA Helio")

       allocate(this%m_plan_gm(this%m_plan_nb))
       if( any( shape(this%m_plan_gm) /= shape(m) ) ) then
          stop "Erreur : set_plan_mass : dimensions differentes"
       endif

       m0 = this%m_star_mass
       this%m_plan_gm =this%m_cG*m
       
      end  subroutine sysplaodehelio_set_plan_mass  

!***********************************************************************
!> @brief allouer le vecteur d'etat du second membre
!***********************************************************************
      subroutine sysplaodehelio_allocate_y(this, y)
       implicit none
       class(t_sysplaodehelio), intent(inout) :: this
       real(TREAL), dimension(:), allocatable, intent(inout) :: y   !< vecteur d'etat au temps t
       integer nplan, ndim
       
       nplan = this%get_plan_nb()
       ndim = 6*nplan
       allocate(y(1:ndim))
      end subroutine sysplaodehelio_allocate_y

!***********************************************************************
! initialiser le vecteur d'etat du second membre
!***********************************************************************
      subroutine sysplaodehelio_fill_to_y(this, y, t)
       implicit none
       class(t_sysplaodehelio), intent(inout) :: this
       real(TREAL), dimension(:), intent(inout) :: y   !< vecteur d'etat au temps t
       real(TREAL), intent(in) :: t                 !< temps ou est appelle t_odefill
       integer nplan, j, indv
       
       nplan = this%get_plan_nb()
       indv = 3*nplan
       do j=1, nplan
        y(3*(j-1)+1:3*(j-1)+3) = this%m_plan_coordX(:,j)
        y(indv+3*(j-1)+1:indv+3*(j-1)+3) = this%m_plan_coordXP(:,j)
       enddo
       
      end subroutine sysplaodehelio_fill_to_y

!***********************************************************************
! recupere le vecteur d'etat du second membre pour remplir this
!***********************************************************************
      subroutine sysplaodehelio_fill_from_y(this, y)
       implicit none
       class(t_sysplaodehelio), intent(inout) :: this
       real(TREAL), dimension(:), intent(in) :: y   !< vecteur d'etat au temps t
       integer nplan, j, indv
       
       nplan = this%get_plan_nb()
       indv = 3*nplan
       do j=1, nplan
        this%m_plan_coordX(:,j) = y(3*(j-1)+1:3*(j-1)+3)
        this%m_plan_coordXP(:,j) = y(indv+3*(j-1)+1:indv+3*(j-1)+3)
       enddo

      end subroutine sysplaodehelio_fill_from_y


!***********************************************************************
! fonction pour calculer le second membre de l'ODE
! avec y=(r,rdot) et retourne  dy=(dr/dt=rdot, Drdot/dt)
!***********************************************************************
       subroutine sysplaodehelio_fcn(this, n, t, y, dy)
        use mod_coordcircumbin
        implicit none
        class(t_sysplaodehelio), intent(inout) :: this
        integer, intent(in) :: n                     !< nombre de composante de y
        real(TREAL), intent(in) :: t                 !< temps ou est appelle t_odefcn
        real(TREAL), dimension(:), intent(in) :: y   !< vecteur d'etat au temps t
        real(TREAL), dimension(:), intent(out) :: dy !< derivee du vecteur d'etat calcule au temps t
        real(TREAL) :: mu, temp,r, r3, cmusol, mui, muj
        real(TREAL), dimension(3)  :: sol
        integer indv, nplan, ip, k, indi, jp, indj
        
        nplan = this%get_plan_nb()

        !call sysplaodehelio_fill_from_y(this, y)
        !dy(indv+1:n)=REALC(0.)
        if (n.eq.6*nplan) then
            indv = 3*nplan
            dy(1:indv)=y(indv+1:n)
        else
            indv = 0
        endif    
        
!--------!  cmusol : masse du soleil = 1 x G
        cmusol = this%m_cG
!--------! s(1..3) acceleration du soleil       
        sol = 0
        do ip = 1, nplan
           mu = this%m_plan_gm(ip)
           indi=3*(ip-1)
           temp = y(indi+1)**2+y(indi+2)**2+y(indi+3)**2
           r = sqrt(temp)
           r3 = REALC(1.0)/(temp*r)
           do k=1,3
              dy(indv+indi+k) = -cmusol*y(indi+k)*r3
              sol(k) = sol(k)-mu*y(indi+k)*r3
           end do
        enddo
        do ip = 1, nplan
           mui = this%m_plan_gm(ip)
           indi = 3*(ip-1)
           do jp = ip+1, nplan
              muj = this%m_plan_gm(jp)
              indj = 3*(jp-1)
              temp = REALC(0.)
              do k = 1,3
                 temp = temp + (y(indi+k)-y(indj+k))**2
              end do
              r = sqrt(temp)
              r3 = REALC(1.)/(temp*r)
              do k = 1,3
                 dy(indv+indi+k)=dy(indv+indi+k)                            &
     &            +muj*(y(indj+k)-y(indi+k))*r3
                 dy(indv+indj+k)=dy(indv+indj+k)                            &
     &            +mui*(y(indi+k)-y(indj+k))*r3
              end do 
           end do
           do k = 1,3
               dy(indv+indi+k) = dy(indv+indi+k)+sol(k)
           end do
        end do
        
         
       end subroutine sysplaodehelio_fcn

       subroutine energie_phvh(np, aG, cmu, y,yp, etot)
*---------------------------------------------------
*  Calcul de l'energie du systeme
*
*   cmu : tableau des (masses/masse du soleil)*aG
*
* ASD 4/8/95 (*m/4)
*---------------------------------------------------
        implicit none
        integer, intent(in)  :: np
        real(TREAL), dimension(1:3*np), intent(in)  :: y, yp
        real(TREAL), intent(out)  :: etot
        real(TREAL), intent(in)  :: aG
        real(TREAL), dimension(1:3*np), intent(in)  :: cmu
        real(TREAL), dimension(1:3) :: v
        real(TREAL) :: tot, ep,r, ec
        integer ip, i, indi, jp, indj, k

*--------! v :  vitesse du barycentre
*--------! tot : masse totale ( masse du soleil =1)
       do i=1,3
          v(i)=REALC(0.)
       enddo
       ep = REALC(0.)
       tot = aG
       do ip = 1, np
          do i=1,3
             v(i) = v(i) + cmu(ip)*yp(3*(ip-1)+i)
         end do
         tot = tot + cmu(ip)
         r=0.D0
         do i=1,3
             r=r+y(3*(ip-1)+i)**2
         end do
         r = sqrt(r) 
         ep = ep - aG*cmu(ip)/r
      end do
      do i=1,3
        v(i) = v(i)/tot
      end do
*--------! energie cinetique du soleil
      ec = REALC(0.5)*aG*(v(1)**2 + v(2)**2 + v(3)**2)
      do  ip = 1, np
        do i=1,3
        ec = ec+REALC(0.5)*cmu(ip)*(yp(3*(ip-1)+i)-v(i))**2
        end do
        do jp = ip+1, np
           r=REALC(0.)
           indi=3*(ip-1)
           indj=3*(jp-1)
           do k=1,3
             r=r+(y(indi+k)-y(indj+k))**2 
           end do
           r =sqrt(r)
           ep = ep - cmu(ip)*cmu(jp)/r
        end do
       end do
       etot = (ec + ep)/aG
       return
       end subroutine energie_phvh

!***********************************************************************
!> @brief calcule l'energie de this
!***********************************************************************
      subroutine sysplaodehelio_energie(this, H)
       implicit none
       class(t_sysplaodehelio), intent(inout):: this  !< dummy argument
       real(TREAL), intent(out) :: H  !< energie calculee
       
       call energie_phvh(this%m_plan_nb, this%m_cG , this%m_plan_gm,             &
     &    this%m_plan_coordX, this%m_plan_coordXP,H)
       
      end  subroutine sysplaodehelio_energie  

        subroutine moment(np, aG, cmu, y,yp, c)
*-------------------------------------------
*   Calcule le moment cinetique global
*
*   cmu : tableau des (masses/masse du soleil)*aG
*
* ASD 4/8/95 (*m/4)
*-------------------------------------------   
        implicit none
        integer, intent(in)  :: np
        real(TREAL), dimension(1:3*np), intent(in)  :: y, yp
        real(TREAL), dimension(1:3), intent(out)  :: c
        real(TREAL), intent(in)  :: aG
        real(TREAL), dimension(1:3*np), intent(in)  :: cmu
        real(TREAL), dimension(1:3) :: p, v
        real(TREAL) :: tot
        integer ip, i, ind
*---------! p,v : position et vitesse du soleil par rapport
*---------!         au barycentre
*---------! tot : masse totale ( masse du soleil =1)
        p=0
        v=0
        tot=aG
        do ip=1,np
          ind = 3*(ip-1)
          tot =tot+cmu(ip) 
          do i=1,3
            p(i)=p(i)+cmu(ip)*y(ind+i)
            v(i)=v(i)+cmu(ip)*yp(ind+i)
          end do
        end do
        c(1)=-(p(2)*v(3)-p(3)*v(2))/tot
        c(2)=-(p(3)*v(1)-p(1)*v(3))/tot
        c(3)=-(p(1)*v(2)-p(2)*v(1))/tot
        do ip=1,np
          ind=3*(ip-1)
          c(1)=c(1)+cmu(ip)*(y(ind+2)*yp(ind+3)-y(ind+3)*yp(ind+2))  
          c(2)=c(2)+cmu(ip)*(y(ind+3)*yp(ind+1)-y(ind+1)*yp(ind+3))
          c(3)=c(3)+cmu(ip)*(y(ind+1)*yp(ind+2)-y(ind+2)*yp(ind+1))
        end do
        do i=1,3
           c(i)=c(i)/aG
         end do
        return
        end subroutine moment


!***********************************************************************
!> @brief calcule le moment cinetique de this
!***********************************************************************
      subroutine sysplaodehelio_mom_cin(this, C)
       implicit none
       class(t_sysplaodehelio), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(1:3), intent(out) :: C  !< moment cinetique

       call moment(this%m_plan_nb, this%m_cG , this%m_plan_gm,             &
     &    this%m_plan_coordX, this%m_plan_coordXP,C)

      end  subroutine sysplaodehelio_mom_cin  

      end module mod_sysplaodehelio

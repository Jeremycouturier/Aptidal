!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file syscircummaree.f 
!!  \brief systeme planetaires avec 2 etoiles et N planetes
!!         pouvant etre integre par un schema symplectique.
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!!         les interactions newtonniennes sont prises en compte.
!!         un effet de maree est applique sur la 2nde etoile mais pas sur les planetes
!!
!! le systeme integre est (dV/dt, DVhat/dt)=f(V,Vhat)
!!
!
! history : creation 04/12/2015
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme a 2 etoiles avec N planetes 
! exprime en variables heliocentriques
! les interactions newtonniennes sont prises en compte
!  un effet de maree est applique sur la 2nde etoile mais pas sur les planetes
!***********************************************************************
      module mod_syscircummaree
       use mod_sysplaode
       use mod_syscircumnewt

!***********************************************************************
! type de coordonnee complexe (v0, v+) selon equation 29
!***********************************************************************
       type :: coordvarshalovich
         real(TREAL) :: v0 ! composante v0
         complex(TREAL) :: vp ! composante v+
        contains
         procedure:: get_xyz => coordvarshalovich_get_xyz !  composantes xyz
         procedure:: get_norm => coordvarshalovich_get_norm !  norme du vecteur complexe
       end type coordvarshalovich


!***********************************************************************
! type de coordonnee complexe pour les operateurs (v0, v+) avec v0 complex aussi
!***********************************************************************
       type :: operavarshalovich
        complex(TREAL) :: v0 ! composante v0
        complex(TREAL) :: vp ! composante v+
       end type operavarshalovich

       
!***********************************************************************
!> @class t_syscircummaree
!! classe decrivant 1 systeme planetaire avec 2 etoiles et N planetes
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!!  les interactions newtonniennes sont prises en compte.
!!  un effet de maree est applique sur la 2nde etoile mais pas sur les planetes
!!  
!! attention : t_syspla::m_star_mass = masse de la 1ere etoile
!!             t_syspla::m_plan_nb = nombre de planetes +1 (2nde etoile a l'indice 1)
!!             t_syspla::m_plan_mass(1) = masse de la 2nde etoile
!!             t_syspla::m_plan_mass(2:) = masse des planetes
!!             t_syspla::m_plan_mpsbeta(1) = m_etoile1/beta_1 avec beta_1=(m_etoile0*m_etoile1)/(m_etoile0+m_etoile1) pour la 2nde etoile
!!             t_syspla::m_plan_mpsbeta(2:) = m_p/beta_p avec beta_p=((m_etoile+m_etoile1)*m_p)/((m_etoile0+m_etoile1)+m_p) pour les planetes 
!!  
!***********************************************************************
      type, extends(t_sysplaode) :: t_syscircummaree 
     
          real(TREAL), dimension(:), allocatable ::   m_mu !< mu(1) = G(somme masses 2 etoiles) , mu(2:) = G*(masses 2 etoiles + masse de la planete)
          real(TREAL) ::m_summass_star !< somme de la masse des 2 etoiles
          real(TREAL), dimension(:), allocatable :: m_plan_mass0 !< masse des 2 etoiles (en 0:1) et masse des planetes (en 2:)
          
          ! pour les effets de marees
          integer :: m_if_maree = 0 !< =1 => application des effet de maree
          integer :: m_ind_maree    !< premier indice des elements de maree dans le vecteur d'etat
          real(TREAL) :: m_maree_Rpla     !< rayon de la planete (=2nde etoile : ici)
          real(TREAL) :: m_maree_xi       !< moment d'inertie reduit de la planete (C = xi*m*R**2)
          real(TREAL) :: m_maree_k20      !< second nombre de Love fluide de la planete
          real(TREAL) :: m_maree_taue     !< temps de relaxation elastique de Maxwell de la planete
          real(TREAL) :: m_maree_tau2     !< temps de relaxation global de la planete
          type(coordvarshalovich) :: m_maree_lcin !< moment cinetique de la planete
          complex(TREAL), dimension(0:2) :: m_maree_znu2 !< coefficients de deformation de la planete
                    
       contains
          
          procedure :: set_plan_nb => syscircummaree_set_plan_nb ! fixe le nombre d'etoiles et de planetes 
          procedure :: set_mass    => syscircummaree_set_mass ! fixe la masse des etoiles et des planetes
          procedure :: set_tides   => syscircummaree_set_tides ! fixe les effets de marees

          procedure :: set_lcin    => syscircummaree_set_lcin !  fixe le moment cinetique initial de la premiere planete
          procedure :: set_vrot    => syscircummaree_set_vrot !  fixe le vecteur instantane de rotation de la 1ere planete
          procedure :: set_znu2    => syscircummaree_set_znu2 !  fixe la deformation initiale de la planete
          
          procedure :: write_output=>syscircummaree_write_output ! ecriture de la sortie (appelee par l'integrateur)
          procedure :: add_maree_out_step =>                            &
     &                             syscircummaree_add_maree_out_step !  fixe le pas de sortie, la taille du buffer des marees

          procedure :: energie => syscircummaree_energie 
          procedure :: mom_cin => syscircummaree_mom_cin 
          
          procedure :: allocate_y  => syscircummaree_allocate_y ! alloue le vecteur d'etat
          procedure :: fill_from_y => syscircummaree_fill_from_y ! recupere les donnees depuis le vecteur d'etat
          procedure :: fill_to_y   => syscircummaree_fill_to_y ! envoie les donnees vers le vecteur d'etat
          procedure :: fcn         => syscircummaree_fcn ! evalue le second membre pour l'ODE
          
          procedure, private :: tides => syscircummaree_tides ! calcule les effets de maree
      end type t_syscircummaree  
      
      contains
         
!***********************************************************************
!> @brief retourne les coordonnees xyz d'un vecteur complex
!***********************************************************************
      function coordvarshalovich_get_xyz(this) result(xyz)
       implicit none
       class(coordvarshalovich), intent(inout):: this !< dummy argument
       real(TREAL), dimension(1:3) :: xyz

       xyz(1) = -real (this%vp)*sqrt(REALC(2.))
       xyz(2) = -aimag(this%vp)*sqrt(REALC(2.))
       xyz(3) =        this%v0
      end function coordvarshalovich_get_xyz
         
!***********************************************************************
!> @brief retourne la norme d'un vecteur complex
!***********************************************************************
      function coordvarshalovich_get_norm(this)
       implicit none
       class(coordvarshalovich), intent(in):: this !< dummy argument
       real(TREAL) :: coordvarshalovich_get_norm
       real(TREAL) :: v0, vp2
       complex(TREAL) :: vp

       v0  = this%v0
       vp  = this%vp
       vp2 = real(vp*conjg(vp))   !< = |vp|**2
       coordvarshalovich_get_norm = sqrt(v0**2+REALC(2.)*vp2)
      end function coordvarshalovich_get_norm

!***********************************************************************
!> @brief fixe le nombre de planetes
!***********************************************************************
      subroutine syscircummaree_set_plan_nb(this, nplan)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       integer,  intent(in) :: nplan  !<nombre de planete
       
       call sysplaode_set_plan_nb(this, nplan+1) ! +1 utile pour stocker la 2nde etoile a l'indice 1

      end  subroutine syscircummaree_set_plan_nb  

!***********************************************************************
!> @brief fixe la masse des 2 etoiles et des planetes
!***********************************************************************
      subroutine syscircummaree_set_mass(this, cG, m)
       use mod_coordcircumbin
       implicit none
       class(t_syscircummaree), intent(inout) :: this  !< dummy argument
       real(TREAL), intent(in) :: cG          !< constante de Gauss AU**3/an**2
       real(TREAL), dimension(0:), intent(in) :: m  !< masse des etoiles et planetes (0:1) = etoiles, (2:) : planetes
       real(TREAL) :: mstar ! somme des masses des etoiles
       integer nplan
       
       nplan = this%get_plan_nb()
       call sysplaode_set_star_mass(this, cG, m(0))
       mstar = m(0)+m(1) 
       this%m_summass_star = mstar
       call this%set_plan_mass(m(1:nplan))
       call this%set_graphdotname("SYS CIRCUMBINAIRE MAREE")

       allocate(this%m_mu(1:nplan))
       call circum_calc_mu_mpsbeta(nplan,cG, m(0:nplan),                 &
     &           this%m_mu,this%m_plan_mpsbeta)

       allocate(this%m_plan_mass0(0:nplan)) 
       this%m_plan_mass0(0:nplan) = m(0:nplan)
        
      end  subroutine syscircummaree_set_mass  

!***********************************************************************
!> @brief fixe les effets de marees
!***********************************************************************
      subroutine syscircummaree_set_tides(this, if_maree, Rpla, xi,      &
     &           k20, taue, tau2)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       integer, intent(in) :: if_maree          !< =1 => application des effets de marees
       real(TREAL), intent(in) :: Rpla          !< rayon de la planete en UA
       real(TREAL), intent(in) :: xi            !< moment d'inertie de la planete tel que C = xi*m*R**2
       real(TREAL), intent(in) :: k20           !< second nombre de Love fluide de la planete
       real(TREAL), intent(in) :: taue          !< temps de relaxation elastique (Maxwell) de la planete
       real(TREAL), intent(in) :: tau2          !< temps de relaxation global de la planete
       integer nplan
       
       nplan = this%get_plan_nb()
       this%m_if_maree   = if_maree
       this%m_maree_Rpla = Rpla
       this%m_maree_xi   = xi
       this%m_maree_k20  = k20
       this%m_maree_taue = taue
       this%m_maree_tau2 = tau2
       this%m_ind_maree  = 6*nplan
      end  subroutine syscircummaree_set_tides

!***********************************************************************
!> @brief fixe le moment cinetique initial de rotation de la premiere planete
!***********************************************************************
      subroutine syscircummaree_set_lcin(this, lcin_ci)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(3), intent(in) :: lcin_ci  !< moment cinetique de rotation
       type(coordvarshalovich) :: lcin

       call compx(lcin_ci, lcin)
       this%m_maree_lcin = lcin
      end subroutine syscircummaree_set_lcin

!***********************************************************************
!> @brief fixe le vecteur instantane de rotation de la planete interne
!***********************************************************************
      subroutine syscircummaree_set_vrot(this, vrot_ci)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(3), intent(in) :: vrot_ci  !< vecteur instantane de rotation
       real(TREAL), dimension(3)             :: lcin_ci  !< moment cinetique de rotation
       type(coordvarshalovich) :: lcin
       real(TREAL) :: xi, m, R

       xi = this%m_maree_xi
       m  = this%m_plan_mass0(1)
       R  = this%m_maree_Rpla
       lcin_ci = xi*m*R**2*vrot_ci
       call compx(lcin_ci, lcin)
       this%m_maree_lcin = lcin
      end subroutine syscircummaree_set_vrot

!***********************************************************************
!> @brief fixe la deformation initiale de la premiere planete
!***********************************************************************
      subroutine syscircummaree_set_znu2(this)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       type (coordvarshalovich) :: x ! coordonnee complexe de la 2nde etoile
       real(TREAL), dimension(1:3) :: v1 ! coordonnee de la 2nde etoile
       complex(TREAL), dimension(2:3,-3:3) :: Yx  !< Y_{l,m}(v/||v||)
       complex(TREAL), dimension(2:3,-3:3) :: Rom  !< om**l*Y_{l,m}(om/||om||)
       type (coordvarshalovich) :: om    !< omega
       type (coordvarshalovich) :: lcin  !< moment cinetique l
       complex(TREAL), dimension(0:2) :: ze2  !< z^e(2,m=0..2)
       real(TREAL) :: cG, m, m0, R, xi, k20
       complex(TREAL) Z0
       
       v1   = this%m_plan_coordX(:,1)
       cG   = this%m_cG
       m0   = this%m_plan_mass0(0)
       m    = this%m_plan_mass0(1)
       R    = this%m_maree_Rpla
       xi   = this%m_maree_xi
       k20  = this%m_maree_k20
       lcin = this%m_maree_lcin

       call compx(v1, x)
       call compY(x, Yx)
       call compom(lcin, xi, m, R, om)
       call compR(om, Rom)
       call compze2(x, R, cG, m0, m, Rom, Yx, k20, ze2)

       !DEBUG
       !this%m_maree_znu2(0:2) = ze2(0:2)
       Z0 = CMPLX(REALC(0.), REALC(0.0), kind(REALC(1.)))
       this%m_maree_znu2(0:2) = Z0
      end subroutine syscircummaree_set_znu2

!***********************************************************************
!> @brief allouer le vecteur d'etat du second membre
!***********************************************************************
      subroutine syscircummaree_allocate_y(this, y)
       implicit none
       class(t_syscircummaree), intent(inout) :: this
       real(TREAL), dimension(:), allocatable, intent(inout) :: y   !< vecteur d'etat au temps t
       integer nplan, ndim
       
       nplan = this%get_plan_nb()
       ndim = 6*nplan
       if (this%m_if_maree.eq.1) then
           ndim = ndim + 8
       endif
       allocate(y(1:ndim))
      end subroutine syscircummaree_allocate_y

!***********************************************************************
! initialiser le vecteur d'etat du second membre
!***********************************************************************
      subroutine syscircummaree_fill_to_y(this, y, t)
       implicit none
       class(t_syscircummaree), intent(inout) :: this
       real(TREAL), dimension(:), intent(inout) :: y   !< vecteur d'etat au temps t
       real(TREAL), intent(in) :: t                 !< temps ou est appelle t_odefill
       type (coordvarshalovich) :: lcin !< moment cinetique
       complex(TREAL), dimension(0:2) :: znu2 !< moment cinetique
       integer nplan, j, indv, m, im
       
       nplan = this%get_plan_nb()
       indv = 3*nplan
       do j=1, nplan
        y(3*(j-1)+1:3*(j-1)+3) = this%m_plan_coordX(:,j)
        y(indv+3*(j-1)+1:indv+3*(j-1)+3) = this%m_plan_coordXP(:,j)
       enddo
       
       if (this%m_if_maree.eq.1) then
        im = this%m_ind_maree
        lcin = this%m_maree_lcin
        znu2 = this%m_maree_znu2
        y(im+1:im+3) = lcin%get_xyz()
        y(im+4) = real(znu2(0))
        do m=1, 2
         y(im+5+2*(m-1)) = real(znu2(m))
         y(im+6+2*(m-1)) = aimag(znu2(m))
        enddo
       endif
      end subroutine syscircummaree_fill_to_y

!***********************************************************************
! recupere le vecteur d'etat du second membre pour remplir this
!***********************************************************************
      subroutine syscircummaree_fill_from_y(this, y)
       implicit none
       class(t_syscircummaree), intent(inout) :: this
       real(TREAL), dimension(:), intent(in) :: y   !< vecteur d'etat au temps t
       real(TREAL), dimension(1:3) :: lcin !< coordonnees xyz du moment cinetique
       complex(TREAL) :: znu2
       integer nplan, j, indv, im, m
       
       nplan = this%get_plan_nb()
       indv = 3*nplan
       do j=1, nplan
        this%m_plan_coordX(:,j) = y(3*(j-1)+1:3*(j-1)+3)
        this%m_plan_coordXP(:,j) = y(indv+3*(j-1)+1:indv+3*(j-1)+3)
       enddo

       if (this%m_if_maree.eq.1) then 
        im = this%m_ind_maree
        lcin = y(im+1: im+3)
        call compx(lcin, this%m_maree_lcin)
        znu2 = CMPLX(y(im+4), REALC(0.), kind(REALC(1.0)))
        this%m_maree_znu2(0) = znu2
        do m=1, 2
         znu2 = CMPLX(y(im+5+2*(m-1)),y(im+6+2*(m-1)),kind(REALC(1.0)))
         this%m_maree_znu2(m) = znu2
        enddo
       endif
      end subroutine syscircummaree_fill_from_y


!***********************************************************************
! fonction pour calculer le second membre de l'ODE
! avec y=(V,Vhat) et retourne  dy=(dV/dt, DVhat/dt)
!***********************************************************************
       subroutine syscircummaree_fcn(this, n, t, y, dy)
        use mod_coordcircumbin
        implicit none
        class(t_syscircummaree), intent(inout) :: this
        integer, intent(in) :: n                     !< nombre de composante de y
        real(TREAL), intent(in) :: t                 !< temps ou est appelle t_odefcn
        real(TREAL), dimension(:), intent(in) :: y   !< vecteur d'etat au temps t
        real(TREAL), dimension(:), intent(out) :: dy !< derivee du vecteur d'etat calcule au temps t
        real(TREAL), dimension(3) :: vi
        real(TREAL) :: vi3, mpsbeta, mu
        real(TREAL),dimension(3,this%m_plan_nb)  :: Vdot
        integer indv, j, nplan, im, m
        type(coordvarshalovich) :: pdot, ldot
        real(TREAL), dimension(3) :: dVhat
        complex(TREAL), dimension(0:2)  :: znu2dot 
        
        nplan = this%get_plan_nb()
        indv = 3*nplan

        call syscircummaree_fill_from_y(this, y)
        dy(1:n)=0.d0
        
        call circum_vhat2vdot(nplan,this%m_plan_mass0,                   &
     &    this%m_plan_mpsbeta, this%m_plan_coordXP, Vdot)
       
        ! derive de V
        dy(1:3*nplan)=reshape(Vdot,shape(dy))
        
        ! derivee de Vhat
        do j=1, nplan
           mpsbeta = this%m_plan_mpsbeta(j)
           mu = this%m_mu(j)
           vi = this%m_plan_coordX(:,j)
           vi3 = sqrt(dot_product(vi,vi))**3 
           dy(indv+3*(j-1)+1:indv+3*(j-1)+3)=-mu*vi/vi3/mpsbeta
        enddo   
        
        ! interaction mutuelle
        call syscircummaree_pasB2(this, dy)
        
        ! effet de maree sur la 2nde etoile
        if (this%m_if_maree.eq.1) then
            call this%tides(pdot, ldot, znu2dot)
            mpsbeta = this%m_plan_mpsbeta(1)
            dVhat = pdot%get_xyz()/ mpsbeta
            dy(indv+1:indv+3) = dy(indv+1:indv+3) + dVhat  !< vitesse de la planete
            im = this%m_ind_maree
            dy(im+1:im+3) = ldot%get_xyz()  !< moment cinetique de rotation de la planete
            dy(im+4) = real(znu2dot(0))             !< coefficients de
            do m=1, 2                               !< deformation de la
               dy(im+5+2*(m-1)) = real(znu2dot(m))  !< planete
               dy(im+6+2*(m-1)) = aimag(znu2dot(m)) !< znu2(0:2)
            enddo
        endif

       end subroutine syscircummaree_fcn


!***********************************************************************
!> @brief interaction mutuelle entre les corps
!***********************************************************************
      subroutine syscircummaree_pasB2(this, dy)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:), intent(inout) :: dy !< derivee du vecteur d'etat calcule au temps t
       
      real(TREAL),dimension(3,this%m_plan_nb,this%m_plan_nb) :: TD
      real(TREAL),dimension(3,this%m_plan_nb)  :: Acc
      real(TREAL),dimension(3)                 :: Vect
      real(TREAL) :: mstar ! somme des masses des 2 etoiles
      real(TREAL) :: mplj, factmass, n0, n1, m0,m1, vi3
      real(TREAL),dimension(2:this%m_plan_nb) :: vi_n0_v1,vi_n1_v1
      real(TREAL),dimension(3) :: v1, vi
      integer i, j, nplan, indv, indj1

      mstar = this%m_summass_star
      m0 = this%m_star_mass
      m1 = this%m_plan_mass(1)
      n0 = -m1/mstar !=-m1/(m0+m1)
      n1 = m0/mstar !=m0/(m0+m1)
      nplan = this%m_plan_nb
      v1(:) = this%m_plan_coordX(:,1) ! 2nde etoile

      Acc = REALC(0.E0)
      
      do i=2,nplan
        vi = this%m_plan_coordX(:,i)
        vi_n0_v1(i)=REALC(1.E0)/sqrt(dot_product(vi-n0*v1,vi-n0*v1))**3 
        vi_n1_v1(i)=REALC(1.E0)/sqrt(dot_product(vi-n1*v1,vi-n1*v1))**3 
      enddo

!
!----- pour la 2nde etoile----------------
!      
      do i=2,nplan
        factmass = m0*this%m_plan_mass(i)/mstar
        vi = this%m_plan_coordX(:,i)
        Acc(:,1) = Acc(:,1) + factmass*(vi_n0_v1(i)-vi_n1_v1(i))*vi      &
     &         +factmass*(-n0*vi_n0_v1(i)+n1*vi_n1_v1(i))*v1
      enddo
!
!----- pour les planetes----------------
!
      factmass = m0*m1/mstar
      do i=2,nplan
        vi = this%m_plan_coordX(:,i)
        vi3 = REALC(1.E0)/sqrt(dot_product(vi,vi))**3 
        Acc(:,i) = Acc(:,i)                                              &
     &             + (m0*vi_n0_v1(i)+m1*vi_n1_v1(i)-mstar*vi3)*vi        &
     &             +factmass*(vi_n0_v1(i)-vi_n1_v1(i))*v1
      enddo


!      
!----- Calcul du tableaux TD 
      do i=2,nplan
         do j=i+1,nplan
            Vect= this%m_plan_coordX(:,i)-this%m_plan_coordX(:,j)
            TD(:,i,j) = Vect/sqrt(dot_product(Vect,Vect))**3
            TD(:,j,i) = -TD(:,i,j)
         end do
      end do
!----- Calcul des acceleration Acc
      do i=2,nplan
         do j=2,i-1
            mplj = this%m_plan_mass(j)
            Acc(:,i) = Acc(:,i) + mplj*TD(:,i,j)
         enddo
         do j=i+1,nplan
            mplj = this%m_plan_mass(j)
            Acc(:,i) = Acc(:,i) + mplj*TD(:,i,j)
         enddo
      enddo 
      Acc = this%m_cG*Acc

      indv = 3*nplan
      do j=1,nplan
         indj1 = indv+3*(j-1)
         dy(indj1+1:indj1+3)= dy(indj1+1:indj1+3) - Acc(:,j)
      end do
       
      end  subroutine syscircummaree_pasB2  

!***********************************************************************
!> @brief type de fonction pour calculer le second membre de l'ODE
! y=(V,Vhat) et  dy=(dV/dt, DVhat/dt)
!***********************************************************************
       subroutine syscircummaree_tides(this, pdot, ldot,znu2dot)
        implicit none
        class(t_syscircummaree), intent(inout):: this  !< dummy argument
        type (coordvarshalovich), intent(out) :: pdot !< partie due a la maree de l'equation 24b 
        type (coordvarshalovich), intent(out) :: ldot !< quantite de l'equation 24c
        complex(TREAL), dimension(0:2), intent(out)  :: znu2dot !< z^nu_{2,m=0..2}  : quantite de l'equation 24d
        type (coordvarshalovich) :: x ! coordonnee complexe de la 2nde etoile
        real(TREAL), dimension(1:3) :: v1 ! coordonnee de la 2nde etoile
        complex(TREAL), dimension(2:3,-3:3) :: Yx  !< Y_{l,m}(v/||v||)
        complex(TREAL), dimension(2:3,-3:3) :: Rom  !< om**l * Yom_{l,m}(om/||om||)
        type (operavarshalovich), dimension(-2:2):: JY    !< J(Y_{2,m=-2..2}(v/||v||)) 
        type (operavarshalovich), dimension(-2:2) :: nablaY    !< nablaY_{2,m=-2..2}(v/||v||) selon eq 34
        complex(TREAL), dimension(0:2,-2:2) :: J2om    !< J^2(omega) selon eq 44
        type (coordvarshalovich) :: om    !< omega
        type (coordvarshalovich) :: lcin  !< moment cinetique l
        complex(TREAL), dimension(-2:2) :: z2   !< z_{2,m=-2..2}  (equation 17)
        complex(TREAL), dimension(0:2) :: znu2 !< z^nu(2,m=0..2)
        complex(TREAL), dimension(0:2) :: ze2  !< z^e(2,m=0..2)
        complex(TREAL), dimension(-2:2) :: znu2all !< z^nu(2,m=-2..2)
        real(TREAL) :: nv ! =||v||
        real(TREAL) :: fact, cG, m, m0, R, xi, mu, taue, tau2, k20
        complex(TREAL) :: l0,p0, ZI, Z0
        integer :: mi, mp
        
        v1   = this%m_plan_coordX(:,1)
        nv   = sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
        cG   = this%m_cG
        m0   = this%m_plan_mass0(0)
        m    = this%m_plan_mass0(1)
        mu   = this%m_mu(1)
        R    = this%m_maree_Rpla
        xi   = this%m_maree_xi
        k20  = this%m_maree_k20
        taue = this%m_maree_taue
        tau2 = this%m_maree_tau2
        znu2 = this%m_maree_znu2
        lcin = this%m_maree_lcin
        Z0 = CMPLX(REALC(0.), REALC(0.0), kind(REALC(1.)))
        ZI = CMPLX(REALC(0.), REALC(1.0), kind(REALC(1.)))

        znu2all(0:2) = znu2
        znu2all(-1) = -conjg(znu2(1))
        znu2all(-2) = +conjg(znu2(2))
        
        call compx(v1, x)
        call compY(x, Yx)
        call compom(lcin, xi, m, R, om)
        call compR(om, Rom)
        call compJY(Yx, JY)
        call compnablaY(nv,Yx,nablaY)
        call compze2(x, R, cG, m0, m, Rom, Yx, k20, ze2)
        call compz2(tau2, taue, znu2, ze2, z2)
        call compJ2om(om,J2om)
                
        ! partie due la maree de l'equation 24b
        p0 = Z0
        pdot%vp = Z0
        do mi=-2,2
            p0 = p0+z2(mi)*nablaY(mi)%v0
            pdot%vp = pdot%vp+z2(mi)*nablaY(mi)%vp
        enddo
        pdot%v0=real(mu*R**2*p0)
        pdot%vp=     mu*R**2*pdot%vp

        ! equation 24c
        fact = cG*m*m0*R**2/nv**3
        l0 = Z0
        ldot%vp = Z0
        do mi=-2,2
            l0 = l0+z2(mi)*JY(mi)%v0
            ldot%vp = ldot%vp+z2(mi)*JY(mi)%vp
        enddo
        ldot%v0=real(-ZI*fact*l0)
        ldot%vp=     -ZI*fact*ldot%vp

        ! equation 24d
        do mi=0,2
            znu2dot(mi) = (ze2(mi)-znu2(mi))/tau2
            do mp=-2,2 
                znu2dot(mi) = znu2dot(mi)-ZI*J2om(mi,mp)*znu2all(mp)
            enddo
        enddo
        
      end subroutine syscircummaree_tides  

!***********************************************************************
!> @brief calcule om a partir de lcin = C*om avec C = xi*m*R**2
!***********************************************************************
      subroutine compom(lcin, xi, m, R, om)
       implicit none
       type (coordvarshalovich), intent(in) :: lcin  !< moment cinetique de rotation de la planete
       real(TREAL), intent(in) :: xi  !< moment d'inertie de la planete C/(mR**2)
       real(TREAL), intent(in) :: m   !< masse de la planete
       real(TREAL), intent(in) :: R   !< rayon de la planete
       type (coordvarshalovich), intent(out) :: om   !< vecteur instantane de rotation de la planete
       real(TREAL) :: C

       C = xi*m*R**2
       om%v0 = lcin%v0/C
       om%vp = lcin%vp/C
      
      end subroutine compom

!***********************************************************************
!> @brief calcule (x0, x+) selon equation 29 
!***********************************************************************
      subroutine compx(v, xc)
       implicit none
       real(TREAL), dimension(1:3), intent(in) :: v   !< vecteur position de la planete
       type (coordvarshalovich), intent(out) :: xc    !< coordonnee (x0,x+)  de x/||x|| selon eq 29
     
       xc%v0 = v(3)
       xc%vp = -CMPLX(v(1),v(2), kind(REALC(1.)))/sqrt(REALC(2.))
      
      end subroutine compx

!***********************************************************************
!> @brief calcule Y_{l,m}(v/||v||) selon equation 28 et 31
!***********************************************************************
      subroutine compY(xc, Y)
       implicit none
       type (coordvarshalovich), intent(in) :: xc    !< coordonnee (x0,x+)  de x selon eq 29
       real(TREAL) :: nx !< norme de xc
       complex(TREAL), dimension(2:3,-3:3), intent(out) :: Y  !< quantite Y_{l,m}(v/||v||)
       real(TREAL) :: x0        !< coordonnee x0  de x/||x|| selon eq 29
       complex(TREAL)  :: xp   !< coordonnee x+  de x/||x|| selon eq 29
       type (coordvarshalovich) :: xn    !< coordonnee (x0,x+) normalise
      
       x0 = xc%v0
       xp = xc%vp

       ! normalisation de x
       nx = xc%get_norm()

       xn%v0 = x0/nx
       xn%vp = xp/nx

       call compR(xn, Y)
       
      end subroutine compY

!***********************************************************************
!> @brief calcule R_{l,m}(v) (harmonique solide reguliere)
!***********************************************************************
      subroutine compR(xc, Y)
       implicit none
       type (coordvarshalovich), intent(in) :: xc    !< coordonnee (x0,x+)  de x selon eq 29
       real(TREAL) :: nx  !< norme de xc
       real(TREAL) :: nx2 !< norme**2 de xc
       complex(TREAL), dimension(2:3,-3:3), intent(out) :: Y  !< quantite ||v||^l*Y_{l,m}(v/||v||)
       real(TREAL) :: x0        !< coordonnee x0  de x selon eq 29
       complex(TREAL)  :: xp   !< coordonnee x+  de x selon eq 29
      
       x0  = xc%v0
       xp  = xc%vp
       nx  = xc%get_norm()
       nx2 = nx**2

       Y(:,:)=1D100 ! invalid value
       ! equation 31
       Y(2,0) = REALC(0.5)*(REALC(3.)*x0**2-nx2)
       Y(2,1) = sqrt(REALC(3.))*x0*xp
       Y(2,2) = sqrt(REALC(6.))*xp**2/REALC(2.)
       
       Y(3,0) = REALC(5.)/REALC(2.)*x0**3-REALC(3.)/REALC(2.)*x0*nx2
       Y(3,1) = sqrt(REALC(6.))/REALC(4.)*(REALC(5.)*x0**2-nx2)*xp
       Y(3,2) = sqrt(REALC(30.))/REALC(2.)*x0*xp**2
       Y(3,3) = sqrt(REALC(10.))/REALC(2.)*xp**3
       
       ! equation 28
       Y(2,-1) = -conjg(Y(2,1))
       Y(2,-2) = +conjg(Y(2,2))
       
       Y(3,-1) = -conjg(Y(3,1))
       Y(3,-2) = +conjg(Y(3,2))
       Y(3,-3) = -conjg(Y(3,3))
      
      end subroutine compR


!***********************************************************************
!> @brief calcule J(Y) selon equation 35
!***********************************************************************
      subroutine compJY(Y,JY)
       implicit none
       complex(TREAL), dimension(2:3,-3:3), intent(in) :: Y  !< quantite Y_{l,m}(v/||v||)
       type (operavarshalovich), dimension(-2:2), intent(out) :: JY    !< quantite J(Y_{2,m=-2..2}(v/||v||)) selon eq 35
       integer m, l
       real(TREAL) :: num
       
       l=2
       do m=-2,1
         JY(m)%v0=m*Y(2,m)
         num = l*(l+1)-m*(m+1)
         JY(m)%vp = -sqrt(num/REALC(2.))*Y(l,m+1)
       enddo
       m=2
       JY(m)%v0=m*Y(2,m)
       JY(m)%vp=CMPLX(REALC(0.), REALC(0.), kind(REALC(1.0)))
      
      end subroutine compJY

!***********************************************************************
!> @brief calcule nabla(Y) selon equation 34
!***********************************************************************
      subroutine compnablaY(nx,Y,nablaY)
       implicit none
       real(TREAL), intent(in) :: nx  !< quantite ||v||)
       complex(TREAL), dimension(2:3,-3:3), intent(in) :: Y  !< quantite Y_{l,m}(v/||v||)
       type (operavarshalovich), dimension(-2:2), intent(out) :: nablaY !< quantite nablaY_{2,m=-2..2}(v/||v||) selon eq 34
       integer m, l
       real(TREAL) :: num, denom

       l=2
       denom = nx**(l+2)
       do m=-2,2
         num = (l+m+1)*(l-m+1)
         nablaY(m)%v0=-sqrt(num)*Y(l+1,m)/denom
         num = (l+m+1)*(l+m+2)
         nablaY(m)%vp = -sqrt(num/REALC(2.))*Y(l+1,m+1)/denom
       enddo
      
      end subroutine compnablaY

!***********************************************************************
!> @brief calcule J^2(w) selon equation 44
!***********************************************************************
      subroutine compJ2om(om,J2om)
       implicit none
       type (coordvarshalovich), intent(in) :: om  !< omega
       complex(TREAL), dimension(0:2,-2:2), intent(out) :: J2om    !< quantite J2(omega) selon eq 44
       integer m, l
       real(TREAL) :: num
       complex(TREAL) :: omp, omm
       
       omp = om%vp
       omm = -conjg(omp)
       J2om(:,:)=CMPLX(REALC(0.), REALC(0.0), kind(REALC(1.)))
       
       l=2
       do m=-2,2
         if (m-1.ge.0) then
            num = l*(l+1)-m*(m-1)
            J2om(m-1,m) = -sqrt(num/REALC(2.))*omp
         endif
         if (m.ge.0) then
            J2om(m,m) = m*om%v0
         endif
         if ((m+1.ge.0).and.(m+1.le.2)) then
            num = l*(l+1)-m*(m+1)
            J2om(m+1,m) = +sqrt(num/REALC(2.))*omm
         endif  
       enddo
      
      end subroutine compJ2om


!***********************************************************************
!> @brief calcule ze2 selon l'equation 10
!***********************************************************************
      subroutine compze2(v1, R, cG, m0, m, Rom, Yx, k20, ze2)
       implicit none
       real(TREAL), intent(in) :: R !< rayon de la planete
       real(TREAL), intent(in) :: cG 
       real(TREAL), intent(in) :: m0 !< masse de l'etoile 1 
       real(TREAL), intent(in) :: m !< masse de la seconde etoile
       real(TREAL), intent(in) :: k20 !< contante de love
       type(coordvarshalovich), intent(in) :: v1 !< veteur position de la 2nde etoile
       complex(TREAL), dimension(2:3,-3:3), intent(in) :: Rom  !< quantite omega**l*Y_{l,m}(omega/||omega||)
       complex(TREAL), dimension(2:3,-3:3), intent(in) :: Yx  !< quantite Y_{l,m}(V/||V||)
       complex(TREAL), dimension(0:2), intent(out) :: ze2 !< z^e(2,m=-2..2)
       real(TREAL) nv ! norme de v1
       real(TREAL) fact1, fact2
       integer mp
       
       nv=v1%get_norm()
       
       fact1 = R**3/(REALC(3.)*cG*m)
       fact2 = m0/m*(R/nv)**3
       do mp=0,2
            ze2(mp) = conjg( - fact1*Rom(2,mp) + fact2*Yx(2,mp) )
       enddo
       
       ze2=k20*ze2
  
      end subroutine compze2

!***********************************************************************
!> @brief calcule z2 selon l'equation 17
!***********************************************************************
      subroutine compz2(tau2, taue, znu2, ze2, z2)
       implicit none
        real(TREAL), intent(in) :: tau2, taue !< temps de dephasage
        complex(TREAL), dimension(0:2), intent(in) :: znu2 !< z^nu_(2,m=0..2)
        complex(TREAL), dimension(0:2), intent(in) :: ze2 !< z^e_(2,m=0..2)
        complex(TREAL), dimension(-2:2), intent(out) :: z2 !< z_{2,m=-2..2} 
        integer m
        
        do m=0,2
         z2(m) = (REALC(1.)-taue/tau2)*znu2(m)+(taue/tau2)*ze2(m)
        enddo 
        z2(-1) = -conjg(z2(1))
        z2(-2) = +conjg(z2(2))
      end subroutine compz2

!***********************************************************************
!> @brief ecriture de la sortie (appelee par l'integrateur)
!!   a l'iteration donne et au temps donne
!***********************************************************************
      subroutine syscircummaree_write_output(this, it, t)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: it !< iteration actuelle
       real(TREAL), intent(in) :: t  !< temps correspondant a it
       real(TREAL), dimension(1:9) :: marees
       type(coordvarshalovich) :: lcin
       complex(TREAL), dimension(0:2) :: znu2
       integer :: j, m

       ! appel de la subroutine par defaut
       call sysplaode_write_output(this, it, t)

       do j=1, size(this%m_stepout)
        if (mod(it,this%m_stepout(j)).eq.0) then

         ! cas d'un bufer de maree
         if (this%m_kind(j).eq.3) then
          ! write(*,*) " sortie maree ",it,t
          lcin = this%m_maree_lcin
          znu2 = this%m_maree_znu2
          marees(1)   = t
          marees(2:4) = lcin%get_xyz()
          marees(5)   = real(znu2(0))
          do m=1, 2
           marees(6+2*(m-1)) = real(znu2(m))
           marees(7+2*(m-1)) = aimag(znu2(m))
          enddo
          call this%m_buffer(j)%writedata(marees)
         endif

         ! verifie si une erreur s'est produite dans le successeur
         if (this%m_buffer(j)%check_successor_error()) then
         call this%m_buffer(j)%get_successor_error(this%m_plan_arret)
         endif
        endif
       enddo


      end subroutine syscircummaree_write_output

!***********************************************************************
!> @brief fixe le pas de sortie, la taille du buffer des marees
!! cela cree un nouveau buffer
!! lorsque le buffer est plein, le consommateur buffercs est appelle
!***********************************************************************
      subroutine syscircummaree_add_maree_out_step(this, stepout,          &
     &                             buffersize,  buffercs)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les inetgrales premieres
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 

       call sysplaode_add_out_step(this, stepout, buffersize,                &
     &             buffercs, 9 , 3) ! 9 = temps+moment cin+deformation.

      end  subroutine syscircummaree_add_maree_out_step

!***********************************************************************
!> @brief calcule l'energie de this
!***********************************************************************
      subroutine syscircummaree_energie(this, H)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       real(TREAL), intent(out) :: H  !< energie calculee
       
       H = 0
       write(*,*) "syscircummaree_energie : non supporte car ode"
       
       stop      
      end  subroutine syscircummaree_energie  

!***********************************************************************
!> @brief calcule le moment cinetique de this
!***********************************************************************
      subroutine syscircummaree_mom_cin(this, C)
       implicit none
       class(t_syscircummaree), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(1:3), intent(out) :: C  !< moment cinetique
       type(coordvarshalovich) :: lcin !< moment cinetique de rotation de la planete

       call mom_cin_circum(this%m_plan_nb,this%m_plan_mass,              &
     &    this%m_plan_coordX, this%m_plan_coordXP,C)

       if (this%m_if_maree.eq.1) then
         lcin = this%m_maree_lcin
         C = C + lcin%get_xyz()   !< on ajoute le moment cinetique de rotation de la planete
       endif
       !write(*,*) "syscircummaree_mom_cin",C

      end  subroutine syscircummaree_mom_cin  

      end module mod_syscircummaree

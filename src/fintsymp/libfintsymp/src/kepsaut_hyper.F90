!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Petit, Boué, Gastineau
! 
!  \file kepsaut.f 
!  \brief Implementation of :
!     kepsaut with hyperbolic motion taken into account
!      inspired by Rein and Tamayo 2015 (package REBOUND : https://github.com/hannorein/rebound/ )
!
!    ce fichier ne doit rien contenir de specifique a un projet 
!     Contains the function kepsaut that take position and velocity as well as a dt
!     and return the variation in positions and velocity 
!     Using Stumpf series and Gauss functions.
!     

!***********************************************************************

#include "realc.f"

module mod_kepsaut_hyper
  use mod_arret
  
  private stumpff_cs, stiefel_Gs
  
  real(TREAL),parameter :: TWO_PI = REALC(2.E0)*ATAN2(REALC(0.E0),REALC(-1.E0))
  real(TREAL), dimension(0:34) :: invfact =[ REALC(1.E0), REALC(1.E0), REALC(1.E0)/REALC(2.E0), REALC(1.E0)/REALC(6.E0), REALC(1.E0)/REALC(24.E0), REALC(1.E0)/REALC(120.E0), REALC(1.E0)/REALC(720.E0), REALC(1.E0)/REALC(5040.E0), REALC(1.E0)/REALC(40320.E0), REALC(1.E0)/REALC(362880.E0), REALC(1.E0)/REALC(3628800.E0), REALC(1.E0)/REALC(39916800.E0), REALC(1.E0)/REALC(479001600.E0), REALC(1.E0)/REALC(6227020800.E0), REALC(1.E0)/REALC(87178291200.E0), REALC(1.E0)/REALC(1307674368000.E0), REALC(1.E0)/REALC(20922789888000.E0), REALC(1.E0)/REALC(355687428096000.E0), REALC(1.E0)/REALC(6402373705728000.E0), REALC(1.E0)/REALC(121645100408832000.E0), REALC(1.E0)/REALC(2432902008176640000.E0), REALC(1.E0)/REALC(51090942171709440000.E0), REALC(1.E0)/REALC(1124000727777607680000.E0), REALC(1.E0)/REALC(25852016738884976640000.E0), REALC(1.E0)/REALC(620448401733239439360000.E0), REALC(1.E0)/REALC(15511210043330985984000000.E0),REALC(1.E0)/REALC(403291461126605635584000000.E0), REALC(1.E0)/REALC(10888869450418352160768000000.E0), REALC(1.E0)/REALC(304888344611713860501504000000.E0), REALC(1.E0)/REALC(8841761993739701954543616000000.E0), REALC(1.E0)/REALC(265252859812191058636308480000000.E0), REALC(1.E0)/REALC(8222838654177922817725562880000000.E0), REALC(1.E0)/REALC(263130836933693530167218012160000000.E0), REALC(1.E0)/REALC(8683317618811886495518194401280000000.E0), REALC(1.E0)/REALC(295232799039604140847618609643520000000.E0) ]
  
  type :: t_kepler_solver
     private
     real(TREAL) :: m_beta, m_eta0, m_zeta0
     real(TREAL) :: m_invper, m_r0, m_dt, m_Ms,norm2_v0
     real(TREAL) :: X,X_per_period,r0i,ri
     
     real(TREAL) , dimension(0:3):: Gs
     
   contains
     procedure :: init_solver
     procedure :: solve_kepler_eq
     
  end type t_kepler_solver
contains
  
  !Main routine that should be called from outside
  subroutine kepsaut(r0,v0,M,dt,dr, dv, arret)
    implicit none     
    type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
    real(TREAL),intent(in) :: M,dt
    real(TREAL),intent(in) :: r0(3),v0(3)
    real(TREAL),intent(out) :: dr(3),dv(3)
    type(t_kepler_solver) :: solver
    real(TREAL) :: f,g,df,dg,cs(0:3)
    
    call init_solver(solver, r0, v0, M, dt)
    call solve_kepler_eq(solver,arret)
    
    call stumpff_cs(solver%m_beta*solver%X*solver%X,cs)
    call stiefel_Gs(solver%X, solver%m_beta, solver%Gs)
    
    f= -solver%m_Ms*solver%Gs(2)*solver%r0i
    g= solver%m_dt - solver%m_Ms*solver%Gs(3)
    df= -solver%m_Ms*solver%Gs(1)*solver%r0i*solver%ri
    dg= -solver%m_Ms*solver%Gs(2)*solver%ri
    
    dr = f*r0 + g*v0
    dv = df*r0 + dg*v0

  end subroutine kepsaut
  
  subroutine stumpff_cs(zin,cs)
    implicit none
    real(TREAL), intent(in) :: zin
    real(TREAL), intent(out) :: cs(0:3)
    real(TREAL):: z, c_odd, c_even
    integer n,nmax,k
    cs = REALC(0.E0)
    n=0
    z=zin
    do while (abs(z) .gt. REALC(0.1E0))
       z = z/REALC(4.E0)
       n = n+1
    end do
    ! write(*,*) "z= ", z,"zin= ", zin, "n= ", n
    !     To compute the cs series at machine precision, it is sufficient to stop at 13 for double and 31 for quadruple
    if (TREAL .eq.8) then 
       nmax = 13
    else if (TREAL.eq.16) then
       nmax = 33
    else
       write(*,*)  'The precision is not implemented in kepsaut_hyper'
       write(*,*) "TREAL", TREAL
       stop
    end if
    c_odd = invfact(nmax)
    c_even = invfact(nmax-1)
    ! write (*,*) "niter", nmax
    ! write (*,*) "c_odd", c_odd
    ! write (*,*) "c_even", c_even
    do k=nmax-2,3,-2
       c_odd  = invfact(k) - z*c_odd
       c_even = invfact(k-1) - z*c_even
       !   write (*,*) "niter", k
       !  write (*,*) "c_odd", c_odd
       ! write (*,*) "c_even", c_even
    end do
    
    cs(3) = c_odd
    cs(2) = c_even
    cs(1) = invfact(1)-z*c_odd
    cs(0) = invfact(0)-z*c_even
    ! write (*,*) "niter", 1
    ! write (*,*) "c_odd", cs(1)
    ! write (*,*) "c_even", cs(0)
    do while (n .gt. 0) 
       cs(3) = (cs(2)+cs(0)*cs(3))*REALC(0.25E0)
       cs(2) = cs(1)*cs(1)*REALC(0.5E0)
       cs(1) = cs(0)*cs(1)
       cs(0) = REALC(2.E0)*cs(0)*cs(0)-REALC(1.E0)
       n = n-1
       !   write(*,*) "Going back n", n
       !   write(*,*), cs
    end do
  end subroutine stumpff_cs
  
  subroutine stiefel_Gs(X,beta,Gs)
    implicit none
    real(TREAL), intent(in) :: X, beta
    real(TREAL), intent(out) :: Gs(0:3)
    real(TREAL) :: X2
    X2 = X*X      
    call stumpff_cs(beta*X2,Gs)
    Gs(1) = X*Gs(1) 
    Gs(2) = X2*Gs(2)
    Gs(3) = X2*X*Gs(3)
  end subroutine stiefel_Gs
  
  subroutine init_solver(this, r0, v0, Ms, dt)
    use IEEE_ARITHMETIC
    implicit none
    
    class(t_kepler_solver), intent(inout) :: this
    real(TREAL), dimension(:), intent(in) :: r0,v0
    real(TREAL), intent(in) :: Ms, dt
    real(TREAL) :: norm2_r0, norm2_v0, norm_r0,sbeta,beta,dt_r0i,r0i
    
    norm2_r0 = sum(r0*r0)
    norm_r0 = sqrt(norm2_r0)
    norm2_v0 =sum(v0*v0)
    r0i = REALC(1.0E0)/norm_r0
    beta =  REALC(2.0E0)*Ms*r0i-norm2_v0
    dt_r0i = dt*r0i
    
    this%m_Ms = Ms
    this%m_dt = dt
    this%m_r0 = norm_r0
    this%r0i = r0i
    this%norm2_v0 = norm2_v0
    this%m_beta = beta
    this%m_eta0 = sum(r0*v0)
    this%m_zeta0 = Ms - beta*norm_r0
    
    if (beta .gt. REALC(0.E0)) then
       sbeta = sqrt(beta)
       this%m_invper = beta*sbeta/(TWO_PI*Ms)
       this%X_per_period = TWO_PI/sbeta
       !Second order initial guess
       this%X = dt_r0i*(REALC(1.E0)-dt_r0i*this%m_eta0*REALC(0.5E0)*r0i)
    else
       !For hyperbolic orbits
       this%m_invper = REALC(0.E0)
       this%X = REALC(0.E0)
       this%X_per_period = IEEE_VALUE(this%X_per_period, IEEE_QUIET_NAN)
    end if
  end subroutine init_solver
  
  
  subroutine solve_kepler_eq(this,arret)
    !Should have been initialized
    implicit none
    class(t_kepler_solver), intent(inout) :: this
    type(t_arret), intent(inout) :: arret
    integer, parameter :: N_MAX_NEWT=32
    integer, parameter  :: N_MAX_QUART=64
    integer, parameter  :: N_MAX_SIZE=MAX(N_MAX_NEWT, N_MAX_QUART)
    integer converged,niter
    integer ntest
    real(TREAL), dimension(1:N_MAX_SIZE) :: prevX
    real(TREAL) :: f, fp, fpp,denom
    real(TREAL) oldX,X,eta0Gs1zeta0Gs2, ri
    real(TREAL) X_min,X_max
    real(TREAL) h2,q,vq,Ms2
    real(TREAL) temp,s
    real(TREAL) Xprevn1

    converged = 0
    X =this%X
    oldX = X
    
    !     One Newton step
    call stiefel_Gs(X, this%m_beta, this%Gs)
    eta0Gs1zeta0Gs2 = this%m_eta0*this%Gs(1)+this%m_zeta0*this%Gs(2)
    ri = REALC(1.E0)/(this%m_r0+eta0Gs1zeta0Gs2)
    X = ri*(X*eta0Gs1zeta0Gs2-this%m_eta0*this%Gs(2)- this%m_zeta0*this%Gs(3)+this%m_dt)
    
    !     Test if the first step is small enough
    if (abs(X-oldX).gt.REALC(0.01E0)*this%X_per_period) then
       !Second order method
       !New initial guess
       X = this%m_beta*this%m_dt/this%m_Ms
       Xprevn1 = X

       do niter=1,N_MAX_QUART
          call stiefel_Gs(X, this%m_beta, this%Gs)
          f = this%m_r0*X + this%m_eta0*this%Gs(2) + this%m_zeta0*this%Gs(3) - this%m_dt
          fp= this%m_r0+this%m_eta0*this%Gs(1)+this%m_zeta0*this%Gs(2)
          fpp = this%m_eta0*this%Gs(0) + this%m_zeta0*this%Gs(1)
          denom = fp + sqrt(abs(REALC(16.E0)*fp*fp - REALC(20.E0)*f*fpp))
          X = (X*denom - REALC(5.E0)*f)/denom
          do ntest =1,niter-1
             if (X .eq. prevX(ntest) ) then
                !     Converged. Exit.
                converged = 1
                exit
             end if
             ! write (*,*) "Test",ntest,  X-prevX(ntest)
          end do
          prevX(niter) = X
          ! write (*,*) X-prevX(niter-1)
          if (converged .eq.1) then
             exit
          end if
       end do
       eta0Gs1zeta0Gs2 =this%m_eta0*this%Gs(1)+this%m_zeta0*this%Gs(2)
       ri = REALC(1.E0)/(this%m_r0+eta0Gs1zeta0Gs2)
       !write (*,*) "Quartic method ", niter
    else
       !     Newton method (also for hyperbolic)
       Xprevn1 = X
       do niter = 1, N_MAX_NEWT
          call stiefel_Gs(X, this%m_beta, this%Gs)
          eta0Gs1zeta0Gs2 = this%m_eta0*this%Gs(1)+this%m_zeta0*this%Gs(2)
          ri =  REALC(1.E0)/(this%m_r0+eta0Gs1zeta0Gs2)
          X = ri*(X*eta0Gs1zeta0Gs2-this%m_eta0*this%Gs(2)- this%m_zeta0*this%Gs(3)+this%m_dt)
          Xprevn1 = X
          do ntest =niter-1,1,-1
             if (X .eq. prevX(ntest) ) then
                ! Converged. Exit.
                converged = 1
                exit
             end if
          end do
          prevX(niter) = X
          ! write (*,*) X-prevX(niter-1)
          if (converged .eq.1) then
             exit
          end if
       end do
       !write (*,*) "Newton ", niter
    end if
    if (converged.eq.0) then
       !     Bisection
       write (*,*) "Bissection", X-prevX(niter-1)
       if (this%m_beta .gt. REALC(0.E0)) then
          !     Elliptic
          X_min = this%X_per_period *floor(this%m_dt*this%m_invper)
          X_max = X_min +  this%X_per_period
       else
          !     Hyperbolic
          h2=this%m_r0*this%m_r0*this%norm2_v0-this%m_eta0*this%m_eta0
          Ms2=this%m_Ms*this%m_Ms
          q=h2/(this%m_Ms*( REALC(1.E0)+sqrt( REALC(1.E0)-h2*this%m_beta/(Ms2))))
          vq = sign( sqrt(h2)/q, this%m_dt)
          ! X_max and X_min correspond to dt/r_min and dt/r_max
          !     which are reachable in this timestep
          !     r_max = vq*_dt+r0
          !     r_min = pericenter
          X_min = this%m_dt/(abs(vq*this%m_dt)+this%m_r0) 
          X_max = this%m_dt/q
          if (this%m_dt< REALC(0.E0)) then
             temp= X_min
             X_min = X_max
             X_max = temp
          end if
       end if
       X = (X_min+X_max)/REALC(2.E0)
       do while (abs(X_max-X_min) .gt. abs(X_min+X_max)*REALC(1.E-15) )
          call stiefel_Gs(X, this%m_beta, this%Gs) 
          s= this%m_r0*X + this%m_eta0*this%Gs(2) + this%m_zeta0*this%Gs(3)-this%m_dt
          if (s.ge.REALC(0.E0))then
             X_max = X
          else
             X_min =X
          end if
          X = (X_min+X_max)/REALC(2.E0)
       end do
       eta0Gs1zeta0Gs2=this%m_eta0*this%Gs(1)+this%m_zeta0*this%Gs(2)
       ri = REALC(1.E0)/(this%m_r0+eta0Gs1zeta0Gs2)
       converged = 1
    end if
    if (isnan(ri)) then
       ri =0.D0
       this%Gs(1:3) = REALC(0.E0)
       converged = 1
    end if
    
    if (converged .eq. 1) then
       call stiefel_Gs(X, this%m_beta, this%Gs)
       eta0Gs1zeta0Gs2=this%m_eta0*this%Gs(1)+this%m_zeta0*this%Gs(2)
       ri = REALC(1.E0)/(this%m_r0+eta0Gs1zeta0Gs2)
       this%X = X
       this%ri = ri
    else
       WRITE(*,*) 'ERREUR FATALE DANS KEPSAUT:solve_kepler_eq'
       WRITE(*,*) 'X,DIFF',X,X_min-X_max
       !-------------------gestion de l'erreur 
       call arret%set_error(-3)
    end if
  end subroutine solve_kepler_eq
end  module mod_kepsaut_hyper


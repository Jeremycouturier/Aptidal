!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!  \file kahansum.f 
!  \brief Implementation of :
!    somme compensee (kahan summation)
!
!    ce fichier ne doit rien contenir de specifique a un projet 
!
!
!
! history : creation 2012/03/07
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la sommation compensee
!***********************************************************************
      module kahansum
      
!***********************************************************************
! compensated sum type
!***********************************************************************
        type compensatedsum
             real(TREAL) :: v = REALC(0.)
             real(TREAL) :: err = REALC(0.)
        end type 
        
!***********************************************************************
! operator compensatedsum = real
!***********************************************************************
        interface assignment(=)
        module procedure compensatedsum_operator_assign
        end interface         

!***********************************************************************
! operator compensatedsum + real
!***********************************************************************
        interface operator(+)
        module procedure compensatedsum_operator_pluss
        end interface         
                 
!***********************************************************************
! operator compensatedsum + compensatedsum
!***********************************************************************
        interface operator(+)
        module procedure compensatedsum_operator_plus
        end interface         

!***********************************************************************
! operator compensatedsum - real
!***********************************************************************
        interface operator(-)
        module procedure compensatedsum_operator_minuss
        end interface         
                 
!***********************************************************************
! operator compensatedsum - compensatedsum
!***********************************************************************
        interface operator(-)
        module procedure compensatedsum_operator_minus
        end interface         

!***********************************************************************
! operator  - compensatedsum
!***********************************************************************
        interface operator(-)
        module procedure compensatedsum_operator_minusunary
        end interface         

!***********************************************************************
! operator compensatedsum / real
!***********************************************************************
        interface operator(/)
        module procedure compensatedsum_operator_div
        end interface         

!***********************************************************************
! operator compensatedsum * real
!***********************************************************************
        interface operator(*)
        module procedure compensatedsum_operator_muls
        end interface         

!***********************************************************************
! operator real * compensatedsum
!***********************************************************************
        interface operator(*)
        module procedure compensatedsum_operator_smul
        end interface         

      contains
!***********************************************************************
! calcule (X+err)= (X+err)+v en utilisant la methode de Kahan
! avec X et err a 3 composantes
!
!  @param X (inout) vecteurs 3 composantes
!  @param err (inout) erreur sur X (3 composantes)
!  @param v (in) vecteur a ajouter a X
!***********************************************************************
         subroutine kahansum3(X,err, v)   
          implicit none
            real(TREAL) , intent(inout), dimension(1:3) :: X
            real(TREAL) , intent(inout), dimension(1:3) :: err
            real(TREAL) , intent(in), dimension(1:3) :: v
!            real(TREAL), dimension(1:3) :: AUX
            real(TREAL) :: sum_, news, temp
            integer :: j
#if 0
            AUX = X
            err =  err +v
            X = AUX + err
            err = err + (AUX - X)
#else
            do j=1, 3
            if (v(j).ne.REALC(0.E0)) then
                if (ABS(X(j))>=ABS(v(j))) then 
                sum_ = X(j)
                news = v(j) - err(j)
                temp = sum_ + news;
                err(j) = (temp - sum_) - news
                X(j) = temp;
                else
                sum_ = V(j)
                news = X(j)
                temp = sum_ + news
                err(j)  = err(j)+ ((temp - sum_) - news)
                X(j) = temp
                endif
            endif
            enddo
            
#endif
            end subroutine kahansum3
      
!***********************************************************************
! calcule (X+err)= (X+err)+v en utilisant la methode de Kahan
! avec X et err a 3 composantes
!
!  @param X (inout) vecteurs 3 composantes
!  @param err (inout) erreur sur X (3 composantes)
!  @param v (in) vecteur a ajouter a X
!***********************************************************************
         subroutine kahansum1(X,err, v)   
          implicit none
            real(TREAL) , intent(inout) :: X
            real(TREAL) , intent(inout) :: err
            real(TREAL) , intent(in)    :: v
!            real(TREAL)                 :: AUX
            real(TREAL)                 :: sum_, news, temp
#if 0
            AUX = X
            err =  err +v
            X = AUX + err
            err = err + (AUX - X)
#else
            if (v.ne.REALC(0.E0)) then
                if (ABS(X)>=ABS(v)) then 
                sum_ = X
                news = v - err
                temp = sum_ + news;
                err = (temp - sum_) - news
                X = temp;
                else
                sum_ = V
                news = X
                temp = sum_ + news
                err  = err+ ((temp - sum_) - news)
                X = temp
                endif
            endif
#endif
         end subroutine kahansum1    
         
!***********************************************************************
! operator  compensatedsum = real -> compensatedsum
!***********************************************************************
        subroutine compensatedsum_operator_assign(v1,s2)
          implicit none
            type(compensatedsum), intent(out) :: v1
            real(TREAL), intent(in) :: s2
            v1%v = s2
            v1%err = s2
        end subroutine compensatedsum_operator_assign

!***********************************************************************
! operator  compensatedsum + real -> compensatedsum
!***********************************************************************
        function compensatedsum_operator_pluss(v1,s2) result(v3)
          implicit none
            type(compensatedsum), intent(in) :: v1
            real(TREAL), intent(in) :: s2
            type(compensatedsum) :: v3
            v3 = v1
            call kahansum1(v3%v, v3%err, s2) 
        end function compensatedsum_operator_pluss

!***********************************************************************
! operator  compensatedsum + compensatedsum -> compensatedsum
!***********************************************************************
        function compensatedsum_operator_plus(v1,v2) result(v3)
          implicit none
            type(compensatedsum), intent(in) :: v1, v2
            type(compensatedsum) :: v3
            v3 = v1
            call kahansum1(v3%v, v3%err, v2%v) 
            call kahansum1(v3%v, v3%err, v2%err) 
        end function compensatedsum_operator_plus
             
!***********************************************************************
! operator  compensatedsum - real -> compensatedsum
!***********************************************************************
        function compensatedsum_operator_minuss(v1,s2) result(v3)
          implicit none
            type(compensatedsum), intent(in) :: v1
            real(TREAL), intent(in) :: s2
            type(compensatedsum) :: v3
            v3 = v1
            call kahansum1(v3%v, v3%err, -s2) 
        end function compensatedsum_operator_minuss

!***********************************************************************
! operator  compensatedsum - compensatedsum -> compensatedsum
!***********************************************************************
        function compensatedsum_operator_minus(v1,v2) result(v3)
          implicit none
            type(compensatedsum), intent(in) :: v1, v2
            type(compensatedsum) :: v3
            v3 = v1
            call kahansum1(v3%v, v3%err, -v2%v) 
            call kahansum1(v3%v, v3%err, -v2%err) 
        end function compensatedsum_operator_minus

!***********************************************************************
! operator  - compensatedsum -> compensatedsum
!***********************************************************************
        function compensatedsum_operator_minusunary(v1) result(v3)
          implicit none
            type(compensatedsum), intent(in) :: v1
            type(compensatedsum) :: v3
            v3%v =  -v1%v
            v3%err =  -v1%err
        end function compensatedsum_operator_minusunary

!***********************************************************************
! operator  compensatedsum / real -> compensatedsum
!***********************************************************************
        function compensatedsum_operator_div(v1,s2) result(v3)
          implicit none
            type(compensatedsum), intent(in) :: v1
            real(TREAL), intent(in) :: s2
            type(compensatedsum) :: v3
            v3%v = v1%v/s2
            v3%err = v1%err/s2
        end function compensatedsum_operator_div
             
!***********************************************************************
! operator  compensatedsum * real -> compensatedsum
!***********************************************************************
        function compensatedsum_operator_muls(v1,s2) result(v3)
          implicit none
            type(compensatedsum), intent(in) :: v1
            real(TREAL), intent(in) :: s2
            type(compensatedsum) :: v3
            v3%v = v1%v*s2
            v3%err = v1%err*s2
        end function compensatedsum_operator_muls
             
!***********************************************************************
! operator  real * compensatedsum -> compensatedsum
!***********************************************************************
        function compensatedsum_operator_smul(s2, v1) result(v3)
          implicit none
            type(compensatedsum), intent(in) :: v1
            real(TREAL), intent(in) :: s2
            type(compensatedsum) :: v3
            v3%v = v1%v*s2
            v3%err = v1%err*s2
        end function compensatedsum_operator_smul
       end module kahansum
      

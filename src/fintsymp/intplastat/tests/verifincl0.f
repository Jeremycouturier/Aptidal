!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
! programme de verification de  XYZKHQP1 pour les inclinaisons i=0 ou i=PI
!***********************************************************************
#include "realc.f"
        program verifinc0
          use mod_elliptid
          use mod_arret
          implicit none

          real(TREAL),dimension(6) :: ell1, ell2
          real(TREAL) e, varpi, lenpi, mu
          integer cpta, cptl, cpte, cptpi, k
          real(TREAL),dimension(3) :: P, V

          type(t_arret) :: arret 
          
          
          lenpi = 4.D0*ATAN(1.D0)
          do cpta=0,10
           do cptl=0,10
            do cpte=0,10
            do cptpi=0,10
         
            ell1(1) = 0.2D0+cpta*0.2D0
            ell1(2) = cptl/10.D0*2.D0*lenpi
            e = cpte/11.D0
            varpi = cptpi/10.D0*2.D0*lenpi
            ell1(3) = e*cos(varpi)
            ell1(4) = e*sin(varpi)
             mu =  39.5164033477944D0  
           
            !cas i = 0
            ell1(5) = 0.D0
            ell1(6) = 0.D0
            call ellipx1(ell1, mu,P,V)
            call XYZKHQP1(P,V,mu,ell2, arret)
             do k=1,6
              if (k.eq.2) then
              if ((abs(sin(ell1(k))-sin(ell2(k))).gt.1D-12).or.         &
     &            (abs(cos(ell1(k))-cos(ell2(k))).gt.1D-12)) then
               write(*,*) "ell1", ell1
               write(*,*) "ell2", ell2
               write(*,*) "diff",  ell1-ell2
               write(*,*) "P",  P
               write(*,*) "V",  V
               stop 1
              endif
              else
              if (abs(ell1(k)-ell2(k)).gt.1D-12) then
               write(*,*) "ell1", ell1
               write(*,*) "ell2", ell2
               write(*,*) "diff",  ell1-ell2
               write(*,*) "P",  P
               write(*,*) "V",  V
               stop 1
              endif
              endif
             enddo   

            ! cas i=pi
            ell1(5) = 1.D0
            ell1(6) = 0.D0
            if (e.eq.0) then
             ell1(3) = 0.1D0*cos(varpi)
             ell1(4) = 0.1D0*sin(varpi)
            endif
            call ellipx1(ell1, mu,P,V)
            call XYZKHQP1(P,V,mu,ell2, arret)
             do k=1,6
              if (k.eq.2) then
              if ((abs(sin(ell1(k))-sin(ell2(k))).gt.1D-12).or.         &
     &            (abs(cos(ell1(k))-cos(ell2(k))).gt.1D-12)) then
               write(*,*) "ell1", ell1
               write(*,*) "ell2", ell2
               write(*,*) "diff",  ell1-ell2
               write(*,*) "P",  P
               write(*,*) "V",  V
               stop 1
              endif
              else
              if (abs(ell1(k)-ell2(k)).gt.1D-12) then
               write(*,*) "ell1", ell1
               write(*,*) "ell2", ell2
               write(*,*) "diff",  ell1-ell2
               write(*,*) "P",  P
               write(*,*) "V",  V
               stop 1
              endif
              endif
             enddo   


            enddo
            enddo
           enddo
          enddo
          
        end  
          
!            write(*,*) "ell2=", ELL2
!            write(*,*) "arret=", arret%m_cause
!            write(*,*) "i=", i
!            write(*,*) "ell=", ell_kh(:,i)
!            write(*,*) "PH=", Ph
!            write(*,*) "Vh=", Vh

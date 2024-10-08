      subroutine interherm(m,x,y,a)
c*****************************************************************************
c     Methode d'interpolation d'Hermite
c     calcule les coefficients a(i) du polynome 
c        P(x)= a(0)+a(1)[[x-x_0]]+...+a(m)[[x-x_0]]^m
c     qui verifie P(x_i)=y_i,...,P^(j)(x_i)=y_i^j).
c     (convention de notation du livre de Stoer et Burlish:
c                 Introduction to Numerical Analysis
c                 Chap 2, pages 52-57.)
c ---> A utiliser couple avec une routine d'evaluation de ce polynome
c     Frederic Joutel (*m/4)
c     (c) ASD/BdL - 11/1991 -
c     M. GASTINEAU 14/09/2016 : suppression des warnings et suppression du common
c*****************************************************************************
c
c***** PARAMETRE
c        nmax : degre maximal du polynome cherche
c***** ENTREE 
c        m : le degre du polynome (m<= nmax)
c        x(0:m) les abscisses rangees en une suite monotone:
c                       x_0,x_0,.....,x_1,x_1,...
c        y(0:m) les ordonnees et les derivees:
c                       y_0,y_0^1,...,y_1,y_1^1,...
c 
c***** SORTIE 
c        a(0:m) les coefficients du polynome d'interpolation d'Hermite
c
c*****************************************************************************
      !implicit double precision (a-h,q-y)
      implicit none
      integer, intent(in) ::  m
      real(8), intent(in) :: x,y
      real(8), intent(out) :: a
      dimension x(0:m),y(0:m),a(0:m)
      integer nmax, k, i
      real(8) :: travail, nr, t, fact
      parameter (nmax=10) 
      dimension travail(0:nmax) 
      dimension nr(0:nmax)
!      COMMON/ASD/IPTR,NFPTR,NERROR
c Recherche de la place des abscisses distinctes
c et initialisation du tableau de travail 
      k=0
      t=y(0)
      travail(0)=t
      nr(k)=0
      do i=1,m
         if (x(i).eq.x(i-1)) then
            nr(i)=k
            travail(i)=t
         else
            k=i
            nr(i)=i
            t=y(i)
            travail(i)=t
         endif
      end do
c Calcul des coefficients
      a(0)=travail(0)
      fact=1.D0
      do k=1,m
         fact=fact*k 
         do i=0,m-k
            if (x(i).eq.x(i+k)) then
               travail(i)=y(nr(i)+k)/fact
            else
               travail(i)=(travail(i)-travail(i+1))/(x(i)-x(i+k))
            endif
         end do
         a(k)=travail(0)
      end do
!      IF (IPTR.GE.1) THEN
!         WRITE(NFPTR,*) 'Routine interherm'
!         IF (IPTR.GE.2) THEN
!            WRITE(NFPTR,*) ' Resultats :'
!            DO I=0,M
!               WRITE(NFPTR,*) I,' : ',a(I)
!            END DO
!         ENDIF
!      ENDIF 
      return
      end



       


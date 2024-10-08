!-------------------------------------------------------------------
! driver de l'interpolation
!-------------------------------------------------------------------
#include "realc.f"
         module mod_interdriver
          
          type t_interdriver
          ! donnees a interpoler
          real(TREAL), allocatable, dimension(:,:) :: mie_y
             ! 1ere dimension = mie_col : nombre de colonnes de y 
             ! 2eme dimension = mie_lig : nombre de lignes de y 
          ! nombre de colonnes
          integer                :: mie_col
          ! nombre de lignes
          integer               :: mie_lig
          
          ! pas de temps entre 2 donnees
          real(TREAL)                :: mie_tpas 
          
          ! intervalle courant   
          integer                :: mie_interv
          
          ! nombre d'intervalle courant   
          integer                :: mie_nbinterv

          
          ! coefficients des polynomes d'interpolation
          real(TREAL), allocatable, dimension(:,:,:) :: mie_cp
             ! 1ere dimension = mie_col : nombre de colonnes de y 
             ! 2eme dimension = 0:mie_nbinterv : nombre d'intervalles de y
             ! 3eme dimension = 0:mie_ni : coefficients du polynome
             
          ! coefficients des polynomes d'interpolation
          real(TREAL), allocatable, dimension(:) :: mie_tps
             ! 1ere dimension = mie_lig : nombre de lignes de y 

         contains
          procedure :: load => interdriver_load ! charge le fichier a interpoler
          procedure :: calc => interdriver_calc ! interpole les donnees en un point
        
         end type t_interdriver
         
          ! nombre de points utilises pour l'interpolation
          integer,private, parameter            :: mie_ni=8 

         contains 
         
!-------------------------------------------------------------------
! initialisation de l'interpolation
! charge le fichier nomfich et utilise ncol+1 colonnes dans ce fichier
! on interpolera ncol colonnes (la 1ere colonne est le temps (=>+1)).
!-------------------------------------------------------------------
         subroutine interdriver_load(this, nomfich, ncol)
          use modreadfile
          implicit none
          class(t_interdriver), intent(inout):: this  !< dummy argument
          character(len=*) , intent(in) :: nomfich
          integer, intent(in) :: ncol
          integer j, k, l, imin, cptligread
          real(TREAL), dimension(0:mie_ni-1):: xa, ya, ca
          integer IPTR,NFPTR,NERROR
          COMMON/ASD/IPTR,NFPTR,NERROR   
          
          IPTR = -1
          NFPTR = 6
          NERROR  = 0     
          ! lecture du fichier
          this%mie_col = ncol+1
          !call this%readfile(nomfich, mie_col)
          call readfile(nomfich, this%mie_col, this%mie_y, cptligread)
         
          ! determination du pas des donnees
          this%mie_tpas = this%mie_y(1,2)-this%mie_y(1,1) 
          this%mie_lig = size(this%mie_y,2)      
          this%mie_nbinterv = abs((this%mie_y(1,this%mie_lig)             &
     &        -this%mie_y(1,1))/this%mie_tpas)
                    
          ! intervalle courant = -1
          ! => recharger un nouvel intervalle a la fois suivante
          this%mie_interv = -1
          
          ! donnees lues
!          write(*,*)'pas de temps (tpas)=', mie_tpas
!          write(*,*)'nombre de lignes lues (mie_lig)=', mie_lig
!          write(*,*)'nombre d intervalle (mie_nbinterv)=',mie_nbinterv
          
          allocate(this%mie_tps(this%mie_lig))
          do j=1, this%mie_lig
           this%mie_tps(j) = this%mie_y(1,j)
          enddo
          
          allocate(this%mie_cp(ncol,0:this%mie_nbinterv,0:mie_ni))

          
          ! Calcul des coefficients des polynomes
          ! pour chaque intervalle
          do j=0, this%mie_nbinterv
           ! determiner les points a utiliser
           imin = MIN(MAX(1,j+2-mie_ni/2),this%mie_lig-mie_ni+1)
           do l=0, mie_ni-1
            xa(l) = this%mie_tps(imin+l)
           enddo
           ! pour chaque colonne
           do k=1, ncol
            do l=0, mie_ni-1
             ya(l) = this%mie_y(k+1,imin+l)
            enddo
            
            call interherm(mie_ni-1,xa, ya, ca)
            
            do l=0, mie_ni-1
             this%mie_cp(k,j, l) = ca(l)
            enddo
           enddo
          enddo
          
          deallocate(this%mie_y)
          write(*,*) 'init de l interpolation terminee'
         return
         end subroutine
         
!-------------------------------------------------------------------
! calcul des elements y(1:n) au temps t 
! en utilisant comme temps mie_y(1,:) 
! et elements mie_cp(:)
!-------------------------------------------------------------------
         subroutine interdriver_calc(this, t,y)
          implicit none
          class(t_interdriver), intent(inout):: this  !< dummy argument
          real(TREAL), intent(in) :: t
          real(TREAL), intent(out), dimension(this%mie_col-1) :: y
          integer curinterv, imin
          integer j, l
          real(TREAL), dimension(0:mie_ni-1) :: xa
          
          ! calcul de l'intervalle a utiliser
          curinterv = (t-this%mie_tps(1))/this%mie_tpas
         
          ! verifie la validite de curinterv
          if ((curinterv.lt.0).or.(curinterv.gt.this%mie_nbinterv)) then
           write(*,*) 'intervalle de temps non disponible'
           write(*,*) 't=',t
           write(*,*) 'temps initial=',this%mie_tps(1)
           write(*,*) 'temps final=', this%mie_tps(this%mie_lig)
           stop
          endif
          
          !recupere le temps
          imin = MIN(MAX(1,curinterv+2-mie_ni/2),                          &
     &     this%mie_lig-mie_ni+1)
          do j=0, mie_ni-1
           xa(j) = this%mie_tps(imin+j)
          enddo

           
           do l=1,this%mie_col-1
             y(l) = this%mie_cp(l,curinterv,mie_ni-1)
             do j=mie_ni-2, 0, -1
              y(l) = y(l)*(t-xa(j) )+ this%mie_cp(l,curinterv,j)
             enddo           
          enddo
                            
         return
         end subroutine
         
#if 0
!-------------------------------------------------------------------
! lecture d'un fichier ascii contenant n colonnes
! stockage d'un tableau contenant n colonnes
!-------------------------------------------------------------------
       subroutine intell_readfile(nomfich, ncol)
       use modintell
       implicit none
        character(len=*) , intent(in) :: nomfich
        integer, intent(in) :: ncol
        real(8), allocatable, dimension(:,:) :: ytemp, ytemp2
       ! variable locale
        integer, parameter :: indfich=20
        integer  maxline
        integer cptlig, j
        integer cptligread
       
       ! ouverture du fichier
       open(indfich, file=nomfich,form='FORMATTED',status='old')
       
       maxline=150000
       allocate(ytemp(ncol,maxline))
       
       ! lecture des donnees
       cptligread=0
       do cptlig=1, maxline
        read(indfich,end=1000, fmt=*) (ytemp(j,cptlig),j=1,ncol)
        cptligread = cptligread+1
       enddo
       
       write(*,*) 'fichier trop long augmenter la constante maxline'
       stop
       
       
       ! lecture terminee
1000   allocate(mie_y(ncol,cptligread))
       do cptlig=1, cptligread
        mie_y(:,cptlig) = ytemp(:,cptlig)
       enddo
       
       deallocate(ytemp)
       
       !fermeture du fichier
       close(indfich)
       
       return 
       end
#endif
       
       end module mod_interdriver
       

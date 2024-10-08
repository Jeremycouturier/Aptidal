!********************************************************
!>  \file readfile.f90
!!
!!  \brief lecture de fichiers  ascii        
!!           
!! 
!!  
!! Author : M. GASTINEAU
!!
!! Laboratoire : Astronomie et Systemes Dynamiques, IMCCE, Paris
!!
!! Version : 2.0 
!!
!! Date de derniere modification :  29/01/2018

!! Modification MG 19/05/2011 : suppression de maxline 
!! Modification HM 19/08/2011 : suppression des interfaces, integration des routines dans le module
!!
!<
!********************************************************
#include "realc.f"
       module modreadfile


       contains
      

!-------------------------------------------------------------------
!> lecture d'un fichier ascii contenant n colonnes
!! stockage d'un tableau contenant n colonnes
!! retourne dans cptligread le nombre de lignes lues
!! @param nomfich (in) nom du fichier
!! @param ncol (in) nombre de colonnes a lire
!! @param mie_y (out) tableau lu
!! @param cptligread (out) nombre de lignes lues
!<
!-------------------------------------------------------------------
       subroutine readfile(nomfich, ncol, mie_y, cptligread)
       implicit none
        character(len=*) , intent(in) :: nomfich
        integer, intent(in) :: ncol
        real(TREAL), allocatable, dimension(:,:), intent(inout) :: mie_y
        integer, intent(out) :: cptligread
       ! variable locale
        integer, parameter :: indfich=20
        integer cptlig, j
       

       ! type definissant la liste des lignes
       type LineData
        type(LineData), pointer :: next 
        real(kind=TREAL), dimension(:),allocatable :: value
       end type LineData
       
       type(LineData), pointer :: list, pnew
       
       nullify(list)

       ! ouverture du fichier
       open(indfich, file=nomfich,form='FORMATTED',status='old')
       
       
       ! lecture des donnees
       cptligread=0
       do while(1.eq.1)
        allocate(pnew)
        allocate(pnew%value(1:ncol))
        read(indfich,end=1000, fmt=*) (pnew%value(j),j=1,ncol)
        cptligread = cptligread+1
        pnew%next => list
        list => pnew
        if (mod(cptligread,1000000).eq.0) then
         write(*,*) 'ligne lue : ', cptligread
        endif
       enddo
       
       ! lecture terminee
1000   deallocate(pnew%value)
       deallocate(pnew)
       allocate(mie_y(ncol,cptligread))
       do cptlig=1, cptligread
        pnew=>list
        mie_y(:,cptligread-cptlig+1) = pnew%value(:)
        list=>list%next
        deallocate(pnew%value)
        deallocate(pnew)
       enddo
       
       
       !fermeture du fichier
       close(indfich)
       
       return 
       end subroutine readfile
      
      
!-------------------------------------------------------------------
!> lecture d'un fichier ascii contenant n colonnes en quadruple precision 
!! stockage d'un tableau contenant n colonnes
!! retourne dans cptligread le nombre de lignes lues
!! @param nomfich (in) nom du fichier
!! @param ncol (in) nombre de colonnes a lire
!! @param mie_y (out) tableau lu
!! @param cptligread (out) nombre de lignes lues
!<
!-------------------------------------------------------------------
       subroutine readfile16(nomfich, ncol, mie_y, cptligread)
       implicit none
        character(len=*) , intent(in) :: nomfich
        integer, intent(in) :: ncol
        real(16), allocatable, dimension(:,:), intent(inout) :: mie_y
        integer, intent(out) :: cptligread
       ! variable locale
        integer, parameter :: indfich=20
        integer cptlig, j
       
       ! type definissant la liste des lignes
       type LineData
        type(LineData), pointer :: next 
        real(kind=16), dimension(:),allocatable :: value
       end type LineData
       
       type(LineData), pointer :: list, pnew
       
       nullify(list)

       ! ouverture du fichier
       open(indfich, file=nomfich,form='FORMATTED',status='old')
       
       
       ! lecture des donnees
       cptligread=0
       do while(1.eq.1)
        allocate(pnew)
        allocate(pnew%value(1:ncol))
        read(indfich,end=1000, fmt=*) (pnew%value(j),j=1,ncol)
        cptligread = cptligread+1
        pnew%next => list
        list => pnew
        if (mod(cptligread,1000000).eq.0) then
         write(*,*) 'ligne lue : ', cptligread
        endif
       enddo
       
       
       ! lecture terminee
1000   deallocate(pnew%value)
       deallocate(pnew)
       allocate(mie_y(ncol,cptligread))
       do cptlig=1, cptligread
        pnew=>list
        mie_y(:,cptligread-cptlig+1) = pnew%value(:)
        list=>list%next
        deallocate(pnew%value)
        deallocate(pnew)
       enddo
       
       !fermeture du fichier
       close(indfich)
       
       write(*,*) 'total de lignes lues : ', cptligread
       return 
       end subroutine readfile16
      
!-------------------------------------------------------------------
!> lecture d'un fichier ascii sous forme d'un tableau de chaine
!! stockage d'un tableau contenant n lignes. Chaque ligne est une ligne de texte
!! retourne dans cptligread le nombre de lignes lues
!! @param nomfich (in) nom du fichier
!! @param lenline (in) longueur maximale admissible d'une ligne
!! @param mie_y (out) tableau lu
!! @param cptligread (out) nombre de lignes lues
!<
!-------------------------------------------------------------------
       subroutine readfiletxtlines(nomfich,lenline, mie_y, cptligread)
       implicit none
        character(len=*) , intent(in) :: nomfich
        integer, intent(in) :: lenline
        character(len=lenline), allocatable, dimension(:),               &
     &       intent(inout) :: mie_y
        integer, intent(out) :: cptligread
       ! variable locale
        integer, parameter :: indfich=20
        integer cptlig
       

       ! type definissant la liste des lignes
       type LineData
        type(LineData), pointer :: next 
        character(len=:), allocatable :: value
       end type LineData
       
       type(LineData), pointer :: list, pnew
       
       nullify(list)

       ! ouverture du fichier
       open(indfich, file=nomfich,form='FORMATTED',status='old')
       
       
       ! lecture des donnees
       cptligread=0
       do while(1.eq.1)
        allocate(pnew)
        allocate(character(len=lenline) :: pnew%value)
        read(indfich,end=1000, fmt='(A)') pnew%value
        cptligread = cptligread+1
        pnew%next => list
        list => pnew
        if (mod(cptligread,1000000).eq.0) then
         write(*,*) 'ligne lue : ', cptligread
        endif
       enddo
       
       ! lecture terminee
1000   deallocate(pnew%value)
       deallocate(pnew)
       allocate(mie_y(cptligread))
       do cptlig=1, cptligread
        pnew=>list
        mie_y(cptligread-cptlig+1) = pnew%value
        list=>list%next
        deallocate(pnew)
       enddo
       
       
       !fermeture du fichier
       close(indfich)
       
       return 
       end subroutine readfiletxtlines
      
!-------------------------------------------------------------------
!> lecture d'un fichier ascii sous forme de 2 tableau de chaines
!! le premier tableau est la premiere collone du fichier.
!! Le second tableau represente la ligne complete (y compris la 1ere colonne).
!! retourne dans cptligread le nombre de lignes lues
!! @param nomfich (in) nom du fichier
!! @param lencol1 (in) longueur maximale admissible de la 1ere colonne
!! @param lenline (in) longueur maximale admissible de la ligne
!! @param mie_key (out) tableau lu pour la 1ere colonne
!! @param mie_y (out) tableau lu pour le reste des lignes
!! @param cptligread (out) nombre de lignes lues
!<
!-------------------------------------------------------------------
       subroutine readfiletxtkeylines(nomfich,lencol1,lenline,mie_key,     &
     &         mie_y,   cptligread)
       implicit none
        character(len=*) , intent(in) :: nomfich
        integer, intent(in) :: lenline, lencol1
        character(len=lenline), allocatable, dimension(:),               &
     &       intent(inout) :: mie_y
        character(len=lencol1), allocatable, dimension(:),               &
     &       intent(inout) :: mie_key
        integer, intent(out) :: cptligread
       ! variable locale
        integer, parameter :: indfich=20
        integer cptlig

       ! type definissant la liste des lignes
       type LineData
        type(LineData), pointer :: next 
        character(len=:), allocatable :: value
       end type LineData
       
       type(LineData), pointer :: list, pnew
       
       nullify(list)

       ! ouverture du fichier
       open(indfich, file=nomfich,form='FORMATTED',status='old')
       
       
       ! lecture des donnees
       cptligread=0
       do while(1.eq.1)
        allocate(pnew)
        allocate(character(len=lenline) :: pnew%value)
        read(indfich,end=1000, fmt='(A)') pnew%value
        cptligread = cptligread+1
        pnew%next => list
        list => pnew
        if (mod(cptligread,1000000).eq.0) then
         write(*,*) 'ligne lue : ', cptligread
        endif
       enddo
       
       ! lecture terminee
1000   deallocate(pnew%value)
       deallocate(pnew)
       allocate(mie_y(cptligread))
       allocate(mie_key(cptligread))
       do cptlig=1, cptligread
        pnew=>list
        mie_y(cptligread-cptlig+1) = pnew%value
        read(pnew%value, *) mie_key(cptligread-cptlig+1)
        write(*,*) 'ligne ', cptlig
        write(*,*) '  key ', mie_key(cptligread-cptlig+1)
        write(*,*) '  complete ', mie_y(cptligread-cptlig+1)
        list=>list%next
        deallocate(pnew)
       enddo
       
       
       !fermeture du fichier
       close(indfich)
       
       return 
       end subroutine readfiletxtkeylines

      end module
        



!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file io_txt_mligncomp.f 
!!  \brief ecriture en ascii sur fichier de n composantes sur m lignes
!!
!!     cette classe t_io_txt_mligncomp gere la sortie de n composantes sous forme de n colonnes sur m lignes.\n
!!  Chaque ligne peut etre prefixee par un radical fourni.
!!     
!!
! history : creation 29/08/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture en ascii sur fichier de n composantes sur m lignes
!***********************************************************************
      module mod_io_txt_mligncomp
       use mod_io_txt_ncomp
       
!***********************************************************************
!> @class t_prefixline
!! prefixe de chaque ligne (pour avoir des tableaux de chaine en fortran)
!!  
!***********************************************************************
      type :: t_prefixline

          character(len=:), allocatable :: m_prefixe !< prefixe

       end type t_prefixline

!***********************************************************************
!> @class t_io_txt_mligncomp
!! classe d'ecriture en ascii sur fichier de n composantes sur m lignes
!! il faut donc m*n composantes en entree.
!! Exception: si le buffer a une taille de 1+(n-1)*m en entree, 
!! alors la 1ere composante est repete sur chaque ligne
!!  
!***********************************************************************
      type, extends(t_io_txt_ncomp) :: t_io_txt_mligncomp

          integer, private :: m_comp_input !< nombre de composantes en entree dans le buffer
          integer, private :: m_col_nb !< nombre de composantes par ligne
          integer, private :: m_lig_nb !< nombre de lignes a ecrire a chaque fois 
          logical, private :: m_repeatfirstcomp !< repete ou la 1ere composante
          
          type(t_prefixline),dimension(:),allocatable :: m_pref2dline !< prefixe au debut de chaque ligne
          integer, dimension(:,:), allocatable :: m_ind2dreorder !< indice pour le reordonnancement des colonnes

       contains
          
          
          procedure :: set_col_nb =>io_txt_mligncomp_set_comp_nb ! fixe le nombre de composantes
          
          procedure :: set_pref2dline =>io_txt_mligncomp_set_prefli ! fixe le prefixe sur chaque ligne
          procedure :: getconsumer => io_txt_mligncomp_getconsumer   ! retourne le consommateur du buffer  en entree
        
          
      end type t_io_txt_mligncomp  

      interface t_io_txt_mligncomp
            module procedure t_io_txt_mligncomp_dotname 
       end interface t_io_txt_mligncomp 
          
      contains
         
!***********************************************************************
!> @brief constructeur qui specifie le nom pour le graphique DOT
!***********************************************************************
       function t_io_txt_mligncomp_dotname (dotname) result (iotxt)
        implicit none
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        type ( t_io_txt_mligncomp ) :: iotxt
        iotxt%m_dotname = dotname
       end function t_io_txt_mligncomp_dotname
         
!***********************************************************************
!> @brief retourne le consommateur du buffer d'entree
!***********************************************************************
            function io_txt_mligncomp_getconsumer(this) result(cs)
             implicit none
             class(t_io_txt_mligncomp), intent(inout):: this  !< dummy argument
#if __INTEL_COMPILER <= 1600
             type(t_buffer_consumer), pointer :: cs 
#else
             class(t_buffer_consumer), pointer :: cs 
#endif
             procedure(t_procbufferfull), pointer :: pfcnfull
             procedure(t_procgraphdot), pointer :: pfcndot
             
             if (.not.associated(this%m_input_buffer)) then
              allocate(this%m_input_buffer)
              pfcnfull=>io_txt_mligncomp_onbufferfull
              pfcndot=>io_txt_ncomp_ongraphdot
              call this%m_input_buffer%set_owner(this,                    &
     &         pfcnfull, pfcnfull, pfcndot)
             endif
             
             cs => this%m_input_buffer
            end  function io_txt_mligncomp_getconsumer

!***********************************************************************
!> @brief fixe le nombre de composantes par ligne et nombre de lignes
!! indique egalement le nombre de composantes en entree 
!! attention, il faut : nbinput >=nbcomp*nblig
!***********************************************************************
            subroutine io_txt_mligncomp_set_comp_nb(this,nbcomp,nblig,    &
     &              nbinput)
             implicit none
             class(t_io_txt_mligncomp), intent(inout):: this  !< dummy argument
             integer, intent(in) :: nbcomp  !< nombre de composantes par ligne
             integer, intent(in) :: nblig  !< nombre de lignes
             integer, intent(in) :: nbinput !< nombre de composantes dans le buffer d'entree
             integer ndiff
             procedure(t_procbufferfull), pointer :: pfcnfull
             procedure(t_procgraphdot), pointer :: pfcndot
             
             this%m_col_nb = nbcomp
             this%m_lig_nb = nblig
             this%m_comp_input = nbinput
             allocate(this%m_pref2dline(1:nblig))
             
             ! determine si la 1ere composante (e.g., le temps) est repete sur les lignes
             ndiff= nbcomp*nblig-nbinput
             this%m_repeatfirstcomp = .false.
             if (ndiff.ne.0)  this%m_repeatfirstcomp = .true.
             
             if (.not.associated(this%m_input_buffer)) then
              allocate(this%m_input_buffer)
              pfcnfull=>io_txt_mligncomp_onbufferfull
              pfcndot=>io_txt_ncomp_ongraphdot
              call this%m_input_buffer%set_owner(this,                    &
     &         pfcnfull, pfcnfull, pfcndot)
             endif
            end  subroutine io_txt_mligncomp_set_comp_nb


!***********************************************************************
!> @brief indique le prefixe qui sera ajoute au debut de chaque ligne
!***********************************************************************
            subroutine io_txt_mligncomp_set_prefli(this, ilig, prefix)
             implicit none
             class(t_io_txt_mligncomp), intent(inout):: this  !< dummy argument
             integer, intent (in) :: ilig !< numero de la ligne
             character(len=*), intent(in) :: prefix  !< prefix sur la ligne ilig
             integer l
             
             l = len_trim(prefix)
             if (allocated(this%m_pref2dline(ilig)%m_prefixe)) then
              deallocate(this%m_pref2dline(ilig)%m_prefixe)
             endif
             allocate(character(len=l) ::                                    &
     &                 this%m_pref2dline(ilig)%m_prefixe)
             this%m_pref2dline(ilig)%m_prefixe = trim(prefix)
             
            end  subroutine io_txt_mligncomp_set_prefli

!***********************************************************************
!> @brief fixe l'ordre des colonnes pour l'ecriture
!! newindex(i,j) indique que la colonne i de la ligne j contient les valeurs de l'entree newindex(i,j)
!***********************************************************************
            subroutine io_txt_mligncomp_reorder(this, newindex)
             implicit none
             class(t_io_txt_mligncomp), intent(inout):: this  !< dummy argument
             integer, dimension(:,:), intent(in) :: newindex !< indice des colonnes 
             integer nlig, ncol
             nlig = this%m_lig_nb
             ncol = this%m_col_nb
             
             allocate(this%m_ind2dreorder(1:ncol, nlig))
             this%m_ind2dreorder(:,:) = newindex(:,:)
             this%m_repeatfirstcomp = .false.
            end  subroutine io_txt_mligncomp_reorder

!***********************************************************************
! @brief ecriture des donnees dans le fichier lorsqu'aucun controle d'erreur n'est fait
!***********************************************************************
            subroutine io_txt_mligncomp_write_wo_mce(this, userdata,         &
     &                                             poswrite)
             implicit none
             class(t_buffer), intent(inout) :: this !< buffer de sortie
             class(t_io_txt_mligncomp), intent(inout) :: userdata    !< de type t_io_txt_mligncomp
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: j,k, p
             real(TREAL), allocatable, dimension(:) :: R
             integer nbcomp, nblig, nbinput
             character(len=40) fmt1001, fmt1002 
              
               nbcomp = userdata%m_col_nb
               nblig = userdata%m_lig_nb
               nbinput = userdata%m_comp_input
               write (fmt1001,'(A,I4,A,A,A)') "(",nbcomp,"(1X,",          &
     &    SFMT_TREAL,"))"
               write(fmt1002,'(A,I4,A,A,A)')"(A,",nbcomp,"(1X,",          &
     &    SFMT_TREAL,"))"

               allocate(R(nbinput))
               
               do j=1, poswrite
               call this%readdata(j, R)

               if (userdata%m_repeatfirstcomp.eqv..true.) then
               !!
               !!! repetition de la 1ere composante
               !!
                if (allocated(userdata%m_ind2dreorder)) then
                write(*,*) "io_txt_mligncomp_onbufferfull: repeat err."
                stop 1
                endif
              
               if (allocated(userdata%m_pref2dline)) then
               do p=0, nblig-1
               write(userdata%m_nf,fmt1002)                                 &
     &                 userdata%m_pref2dline(p+1)%m_prefixe,                &
     &                 R(1), (R(k+1+p*(nbcomp-1)),k=1,nbcomp-1) 
               enddo
               else
               do p=0, nblig-1
                write(userdata%m_nf,fmt1001)  R(1),                         &
     &                 (R(k+1+p*(nbcomp-1)),k=1,nbcomp-1) 
               enddo
               endif
               
               else
               !!
               !!! pas de repetition de la 1ere composante
               !!
               if (allocated(userdata%m_ind2dreorder)) then
               ! reordonnancement des colonnes
               if (allocated(userdata%m_pref2dline)) then
               do p=1, nblig
               write(userdata%m_nf,fmt1002)                                  &
     &                 userdata%m_pref2dline(p)%m_prefixe,                   &
     &                 (R(userdata%m_ind2dreorder(k,p)),k=1,nbcomp) 
               enddo
               else
               do p=1, nblig
               write(userdata%m_nf,fmt1001)                                  &
     &                 (R(userdata%m_ind2dreorder(k,p)),k=1,nbcomp) 
               enddo
               endif
               else
               !pas de reordonncement des colonnes
               if (allocated(userdata%m_pref2dline)) then
               do p=0, nblig-1
               write(userdata%m_nf,fmt1002)                                 &
     &                 userdata%m_pref2dline(p+1)%m_prefixe,                &
     &                 (R(k+p*nbcomp),k=1,nbcomp) 
               enddo
               else
               do p=0, nblig-1
               write(userdata%m_nf,fmt1001) (R(k+p*nbcomp),k=1,nbcomp) 
               enddo
               endif
               endif
               
               endif !!! pas de repetition de la 1ere composante
               
!1001  format(<nbcomp>(1X,FMT_TREAL))       !(etendu)
!1002  format(A,<nbcomp>(1X,FMT_TREAL))       !(etendu)
               enddo
               deallocate(R)
            
            end subroutine io_txt_mligncomp_write_wo_mce


!***********************************************************************
! @brief ecriture des donnees dans le fichier lorsqu'un controle d'erreur n'est fait par ligne
!***********************************************************************
            subroutine io_txt_mligncomp_write_w_mce(this, userdata,         &
     &                                             poswrite)
             implicit none
             class(t_buffer), intent(inout) :: this !< buffer de sortie
             class(t_io_txt_mligncomp), intent(inout) :: userdata    !< de type t_io_txt_mligncomp
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: j,k, p
             real(TREAL), allocatable, dimension(:) :: R
             integer nbcomp, nblig, nbinput
             character(len=40) fmt1001, fmt1002 
             
               nbcomp = userdata%m_col_nb
               nblig = userdata%m_lig_nb
               nbinput = userdata%m_comp_input
               write (fmt1001,'(A,I4,A,A,A)') "(",nbcomp,"(1X,",          &
     &    SFMT_TREAL,"))"
               write(fmt1002,'(A,I4,A,A,A)')"(A,",nbcomp,"(1X,",          &
     &    SFMT_TREAL,"))"

               allocate(R(nbinput))
               
               do j=1, poswrite
               call this%readdata(j, R)

               if (userdata%m_repeatfirstcomp.eqv..true.) then
               !!
               !!! repetition de la 1ere composante
               !!
                !if (allocated(userdata%m_ind2dreorder)) then
                !write(*,*) "io_txt_mligncomp_onbufferfull: repeat err."
                !stop 1
                !endif
              
               if (allocated(userdata%m_pref2dline)) then
               do p=0, nblig-1
                if (this%m_mce(p+1)%m_stop.eq.0) then
                 write(userdata%m_nf,fmt1002)                                 &
     &                 userdata%m_pref2dline(p+1)%m_prefixe,                  &
     &                 R(1), (R(k+1+p*(nbcomp-1)),k=1,nbcomp-1) 
                endif
               enddo
               else
               do p=0, nblig-1
                if (this%m_mce(p+1)%m_stop.eq.0) then
                 write(userdata%m_nf,fmt1001)  R(1),                          &
     &                 (R(k+1+p*(nbcomp-1)),k=1,nbcomp-1) 
                endif
               enddo
               endif
               
               else
               !!
               !!! pas de repetition de la 1ere composante
               !!
               if (allocated(userdata%m_ind2dreorder)) then
               ! reordonnancement des colonnes
               if (allocated(userdata%m_pref2dline)) then
               do p=1, nblig
               write(userdata%m_nf,fmt1002)                                  &
     &                 userdata%m_pref2dline(p)%m_prefixe,                   &
     &                 (R(userdata%m_ind2dreorder(k,p)),k=1,nbcomp) 
               enddo
               else
               do p=1, nblig
               write(userdata%m_nf,fmt1001)                                  &
     &                 (R(userdata%m_ind2dreorder(k,p)),k=1,nbcomp) 
               enddo
               endif
               else
               !pas de reordonncement des colonnes
               if (allocated(userdata%m_pref2dline)) then
               do p=0, nblig-1
               write(userdata%m_nf,fmt1002)                                 &
     &                 userdata%m_pref2dline(p+1)%m_prefixe,                &
     &                 (R(k+p*nbcomp),k=1,nbcomp) 
               enddo
               else
               do p=0, nblig-1
               write(userdata%m_nf,fmt1001) (R(k+p*nbcomp),k=1,nbcomp) 
               enddo
               endif
               endif
               
               endif !!! pas de repetition de la 1ere composante
               
!1001  format(<nbcomp>(1X,FMT_TREAL))       !(etendu)
!1002  format(A,<nbcomp>(1X,FMT_TREAL))       !(etendu)
               enddo
               deallocate(R)
            
            end subroutine io_txt_mligncomp_write_w_mce


!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
            subroutine io_txt_mligncomp_onbufferfull(this, userdata,         &
     &                                             poswrite)
             implicit none
             class(t_buffer), intent(inout) :: this !< buffer de sortie
             class(*), intent(inout) :: userdata    !< de type t_io_txt_mligncomp
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             
             !write(*,*) "t_io_txt_mligncomp_onbufferfull"
             select type(userdata)
              class is(t_io_txt_mligncomp)
               if (allocated(this%m_mce)) then
                call io_txt_mligncomp_write_w_mce(this, userdata,         &
     &          poswrite)
               else
                call io_txt_mligncomp_write_wo_mce(this, userdata,         &
     &          poswrite)
               endif
              
              class default
               stop 't_io_txt_mligncomp_onbufferfull : bad class'
             end select 
             
             call this%empty()
            
            end subroutine io_txt_mligncomp_onbufferfull

      end module mod_io_txt_mligncomp

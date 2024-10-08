!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file io_txt_ncomp.f 
!!  \brief ecriture en ascii sur fichier de n composantes sur 1 ligne
!!
!!     cette classe t_io_txt_ncomp gere la sortie de n composantes sous forme de n colonnes.\n
!!  Chaque ligne peut etre prefixee par un radical fourni.
!!     
!!
! history : creation 20/08/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture en ascii sur fichier de n composantes
!***********************************************************************
      module mod_io_txt_ncomp
       use mod_buffer
       
!***********************************************************************
!> @class t_io_txt_ncomp
!! classe d'ecriture en ascii sur fichier de n composantes
!!  
!***********************************************************************
      type :: t_io_txt_ncomp

          integer, private :: m_comp_nb !< nombre de composantes 
          integer, public :: m_nf !< numero fortran du fichier 
          character(len=512), private :: m_name !< nom du fichier
          type(t_buffer_consumer), pointer, public :: m_input_buffer    &
     &                                               =>NULL()  !< consommateur du buffer en entree
          character(len=40), public :: m_dotname = "I/O txt" !< nom du convertisseur pour le graphique dot
          
          character(len=:), allocatable, private :: m_prefixline !< prefixe au debut de chaque ligne
          integer, dimension(:), allocatable :: m_indexreorder !< indice pour le reordonnancement des colonnes

       contains
          
          
          procedure :: set_comp_nb =>io_txt_ncomp_set_comp_nb ! fixe le nombre de composantes
          procedure :: set_name =>io_txt_ncomp_set_name ! fixe le nom du fichier et le numero de fichier
          procedure :: fcreate => io_txt_ncomp_fcreate  ! cree le fichier et l'ouvre
          procedure :: fclose  => io_txt_ncomp_fclose   ! ferme le fichier
          
          procedure :: set_prefixline =>io_txt_ncomp_set_prefixline ! fixe le prefixe sur chaque ligne
        
          procedure :: getconsumer => io_txt_ncomp_getconsumer   ! retourne le consommateur du buffer  en entree
          
          final ::  io_txt_ncomp_destructor ! destructor       

      end type t_io_txt_ncomp  

      interface t_io_txt_ncomp
            module procedure t_io_txt_ncomp_dotname 
       end interface t_io_txt_ncomp 
          
      contains
         
!***********************************************************************
!> @brief constructeur qui specifie le nom pour le graphique DOT
!***********************************************************************
       function t_io_txt_ncomp_dotname (dotname) result (iotxt)
        implicit none
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        type ( t_io_txt_ncomp ) :: iotxt
        iotxt%m_dotname = dotname
       end function t_io_txt_ncomp_dotname
         
!***********************************************************************
!> @brief destructor
!***********************************************************************
       subroutine io_txt_ncomp_destructor (this)
        implicit none
        type ( t_io_txt_ncomp ) :: this  !< dummy argument
             if (associated(this%m_input_buffer)) then
              deallocate(this%m_input_buffer)
             endif
       end subroutine io_txt_ncomp_destructor

!***********************************************************************
!> @brief retourne le consommateur du buffer d'entree
!***********************************************************************
            function io_txt_ncomp_getconsumer(this) result(cs)
             implicit none
             class(t_io_txt_ncomp), intent(inout):: this  !< dummy argument
#if __INTEL_COMPILER <= 1500
             type(t_buffer_consumer), pointer :: cs 
#else
             class(t_buffer_consumer), pointer :: cs 
#endif
             procedure(t_procbufferfull), pointer :: pfcnfull
             procedure(t_procgraphdot), pointer :: pfcndot
             
             if (.not.associated(this%m_input_buffer)) then
              allocate(this%m_input_buffer)
              pfcnfull => io_txt_ncomp_onbufferfull
              pfcndot => io_txt_ncomp_ongraphdot
              call this%m_input_buffer%set_owner(this,                    &
     &         pfcnfull, pfcnfull, pfcndot)
             endif
             
             cs => this%m_input_buffer
            end  function io_txt_ncomp_getconsumer

!***********************************************************************
!> @brief fixe le nombre de composantes par ligne
!***********************************************************************
            subroutine io_txt_ncomp_set_comp_nb(this, nb)
             implicit none
             class(t_io_txt_ncomp), intent(inout):: this  !< dummy argument
             integer, intent(in) :: nb  !< nombre de composantes par ligne
             procedure(t_procbufferfull), pointer :: pfcnfull
             procedure(t_procgraphdot), pointer :: pfcndot
             this%m_comp_nb = nb
             if (.not.associated(this%m_input_buffer)) then
              allocate(this%m_input_buffer)
              pfcnfull => io_txt_ncomp_onbufferfull
              pfcndot => io_txt_ncomp_ongraphdot
              call this%m_input_buffer%set_owner(this,                   &
     &         pfcnfull, pfcnfull, pfcndot)
             endif
            end  subroutine io_txt_ncomp_set_comp_nb


!***********************************************************************
!> @brief fixe le nom du fichier et le numero fortran
!***********************************************************************
            subroutine io_txt_ncomp_set_name(this, name, nf)
             implicit none
             class(t_io_txt_ncomp), intent(inout):: this  !< dummy argument
             character(len=*), intent(in) :: name  !< nom du fichier
             integer, intent(in) :: nf  !< numero fortran du fichier
             this%m_name = name
             this%m_nf = nf
            end  subroutine io_txt_ncomp_set_name

!***********************************************************************
!> @brief cree le fichier et l'ouvre
!***********************************************************************
            subroutine io_txt_ncomp_fcreate(this)
             implicit none
             class(t_io_txt_ncomp), intent(inout):: this  !< dummy argument
             open(this%m_nf,file=this%m_name,status='unknown')
            end  subroutine io_txt_ncomp_fcreate

!***********************************************************************
!> @brief ferme le fichier
!***********************************************************************
            subroutine io_txt_ncomp_fclose(this)
             implicit none
             class(t_io_txt_ncomp), intent(inout):: this  !< dummy argument
             close(this%m_nf)
            end  subroutine io_txt_ncomp_fclose

!***********************************************************************
!> @brief indique le prefixe qui sera ajoute au debut de chaque ligne
!***********************************************************************
            subroutine io_txt_ncomp_set_prefixline(this, prefix)
             implicit none
             class(t_io_txt_ncomp), intent(inout):: this  !< dummy argument
             character(len=*), intent(in) :: prefix  !< prefix sur chaque ligne
             integer l
             l = len_trim(prefix)
             if (allocated(this%m_prefixline)) then
              deallocate(this%m_prefixline)
             endif
             allocate(character(len=l) :: this%m_prefixline)
             this%m_prefixline = trim(prefix)
            end  subroutine io_txt_ncomp_set_prefixline

!***********************************************************************
!> @brief fixe l'ordre des colonnes pour l'ecriture
!! newindex(i) indique que la colonne i contient les valeurs de l'entree newindex(i)
!***********************************************************************
            subroutine io_txt_ncomp_reorder(this, newindex)
             implicit none
             class(t_io_txt_ncomp), intent(inout):: this  !< dummy argument
             integer, dimension(:), intent(in) :: newindex !< indice des colonnes
             
             allocate(this%m_indexreorder(1:1+this%m_comp_nb))
             this%m_indexreorder(:) = newindex(:)
            end  subroutine io_txt_ncomp_reorder

!***********************************************************************
! @brief fonction appellee lorsque le buffer de sortie est plein
!***********************************************************************
            subroutine io_txt_ncomp_onbufferfull(this, userdata,         &
     &                                             poswrite)
             implicit none
             class(t_buffer), intent(inout) :: this !< buffer de sortie
             class(*), intent(inout) :: userdata    !< de type t_io_txt_ncomp
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) :: j,k
             real(TREAL), allocatable, dimension(:) :: R
             integer nbcomp
             character(len=40) fmt1001, fmt1002 
             
             !write(*,*) "t_io_txt_ncomp_onbufferfull"
             select type(userdata)
              class is(t_io_txt_ncomp)
               nbcomp = userdata%m_comp_nb
               write (fmt1001,'(A,I4,A,A,A)') "(",nbcomp,"(1X,",          &
     &    SFMT_TREAL,"))"
               write(fmt1002,'(A,I4,A,A,A)')"(A,",nbcomp,"(1X,",          &
     &    SFMT_TREAL,"))"
               allocate(R(nbcomp))
               do j=1, poswrite
               call this%readdata(j, R)
               if (allocated(userdata%m_indexreorder)) then
               ! reordonnancement des colonnes
               if (allocated(userdata%m_prefixline)) then
               write(userdata%m_nf,fmt1002) userdata%m_prefixline,             & 
     &                 (R(userdata%m_indexreorder(k)),k=1,nbcomp) 
               else
               write(userdata%m_nf,fmt1001)                                     & 
     &                 (R(userdata%m_indexreorder(k)),k=1,nbcomp) 
               endif
               else
               !pas de reordonncement des colonnes
               if (allocated(userdata%m_prefixline)) then
               write(userdata%m_nf,fmt1002) userdata%m_prefixline,             &
     &                        (R(k),k=1,nbcomp) 
               else
               write(userdata%m_nf,fmt1001) (R(k),k=1,nbcomp) 
               endif
               endif
!1001  format(<nbcomp>(1X,FMT_TREAL))       !(etendu)
!1002  format(A,<nbcomp>(1X,FMT_TREAL))       !(etendu)
               enddo
               deallocate(R)
              class default
               stop 't_io_txt_ncomp_onbufferfull : bad class'
             end select 
             
             call this%empty()
            
            end subroutine io_txt_ncomp_onbufferfull

!***********************************************************************
!> @brief  fonction appellee pour generer le graph au format dot
!***********************************************************************
       subroutine io_txt_ncomp_ongraphdot(this,userdata, nodesrc,              &
     &  nodecounter, num_file)
        implicit none
        class(t_buffer), intent(inout) :: this !< buffer ayant genere l'appel
        class(*), intent(inout) :: userdata !< donnee utilisateur de type t_io_txt_ncomp
        integer, intent(in) :: nodesrc !< numero du noeud  source 
        integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
        integer, intent(in) :: num_file !<numero du fichier de sortie
        
        integer curnode
        
        curnode = nodecounter+1

          select type(userdata)
           class is(t_io_txt_ncomp)
             call buffer_dot_writelabel(num_file,curnode,                 &
     &          trim(userdata%m_dotname))
             call buffer_dot_writenode(num_file, nodesrc,curnode)
             nodecounter = curnode
           class default
            stop 'io_txt_ncomp_ongraphdot : bad class'
          end select 

       end subroutine io_txt_ncomp_ongraphdot

      end module mod_io_txt_ncomp

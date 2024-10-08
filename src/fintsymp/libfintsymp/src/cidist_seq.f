!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file cidist_seq.f 
!!  \brief Gere la distribution de condition initiale et leur execution.
!!
!!  version sequentielle mono-coeur.
!!
!!  Pour chaque condition initiale lue par un lecteur de condition initiale, 
!!         la fonction "run" execute la fonction fournie.
!!          execution sequentielle de la fonction "run".
!!
!!
! history : creation 17/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la distribution d'un fichier texte de conditions initiales quelconques en sequentielle
!***********************************************************************
      module mod_cidist_seq
       use mod_ciread_base
       use mod_ci_base
!***********************************************************************
!> @class t_cidist_seq
!! classe de distribution et execution d'un fichier texte de conditions initiales quelconques
!!  
!***********************************************************************
      type :: t_cidist_seq
          
          logical :: m_ismaster = .true. !< = true si le processus est le maitre
          integer :: m_manyblocksize !< nombre de conditions initiales traitees simultanement par runmany

       contains
                    
          procedure :: set_manyblocksize =>cidist_seq_set_manyblocksize ! fixe le nombre de ci  traitees simultanement par runmany

          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          procedure :: run => cidist_seq_run ! execute toutes les conditions initiales du fichier (une par une)
          procedure :: runmany => cidist_seq_runmany ! execute toutes les conditions initiales du fichier (plusieurs en meme temps)
      end type t_cidist_seq  


      interface
      
       !> prototype des fonctions appellees par la fonction run
       subroutine t_cidist_userfunc(ci, userdata)
        import t_cidist_seq
        import t_ci_base
        class(t_ci_base), intent(inout) :: ci !< condition initiale
        class(*), intent(inout) :: userdata !< donnee utilisateur
       end subroutine t_cidist_userfunc
       
       !> prototype des fonctions appellees par la fonction runmany
       subroutine t_cidist_manyuserfunc(nbci, tabci, userdata)
        import t_cidist_seq
        import t_ci_base_ptr
        integer, intent(in) :: nbci  !< nombre d'elements valides dans ci
        type(t_ci_base_ptr), dimension(:), intent(inout) ::   tabci !< tableau de conditions initiales
        class(*), intent(inout) :: userdata !< donnee utilisateur
       end subroutine t_cidist_manyuserfunc
      end interface
    
      contains
      
!***********************************************************************
!> @brief fixe le nombre de conditions initiales qui seront execeutees simultanement par runmany
!***********************************************************************
            subroutine cidist_seq_set_manyblocksize(this, blocksize)
             implicit none
             class(t_cidist_seq), intent(inout):: this  !< distributeur de ci
             integer, intent(in) :: blocksize !< nombre de conditions initiales 
             this%m_manyblocksize = blocksize
            end subroutine cidist_seq_set_manyblocksize

!***********************************************************************
!> @brief lit chaque condition initiale via ciread 
!! puis appelle la fonction func avec la condition initiale lue et les donnees utilisateurs
!! Elle suppose que cireader est initialise correctement (fichier ouvert, ...)
!***********************************************************************
            subroutine cidist_seq_run(this, cireader, userdata, func)
             implicit none
             class(t_cidist_seq), intent(inout):: this  !< distributeur de ci
             class(t_ciread_base), intent(inout) :: cireader  !<lecteur du fichier de ci
             class(*), intent(inout) :: userdata !< donnee utilisateur
             procedure(t_cidist_userfunc), pointer, intent(in) :: func !< fonction a appeler pour chaque condition initiale
             
             class(t_ci_base), pointer :: ligneci ! condition initiale lue
             
             do while(cireader%freadline(ligneci))
               
               write(*,*) "cidist_seq_run : lecture d'une nouvelle ci"
               call ligneci%debug()
               call func(ligneci, userdata)
               deallocate(ligneci)
               
             enddo
             
            end  subroutine cidist_seq_run

!***********************************************************************
!> @brief lit au plus "m_manyblocksize" conditions initiales via ciread 
!! puis appelle la fonction func avec le tableau des conditions initiales lues 
!! et les donnees utilisateurs
!! Elle suppose que cireader est initialise correctement (fichier ouvert, ...)
!***********************************************************************
            subroutine cidist_seq_runmany(this,cireader,userdata,func)
             implicit none
             class(t_cidist_seq), intent(inout):: this  !< distributeur de ci
             class(t_ciread_base), intent(inout) :: cireader  !<lecteur du fichier de ci
             class(*), intent(inout) :: userdata !< donnee utilisateur
             procedure(t_cidist_manyuserfunc),pointer,intent(in):: func !< fonction a appeler pour chaque condition initiale
             
             type(t_ci_base_ptr), dimension(1:this%m_manyblocksize) ::    &
     &                    ligneci ! tableau de condition initiale lue
             logical bcontinue ! =false si la fin est atteinte
             integer nci ! nombre de ci lue
             integer maxnci
             integer k
             
             bcontinue=.true.
             maxnci = this%m_manyblocksize
             do while(bcontinue.eqv..true.)
               nci = 0
               do while((nci.lt.maxnci).and.(bcontinue.eqv..true.))
                  nci = nci+1
                  bcontinue = cireader%freadline(ligneci(nci)%m_ptr)
               enddo
               
               if (bcontinue.eqv..false.) then
                   nci = nci-1
               endif
               
               if (nci.gt.0) then
                  write(*,*) "cidist_seq_runmany : traitement de ",nci,  &
     &                       "ci"
                  do k=1, nci
                       call ligneci(k)%m_ptr%debug()
                  enddo 
                  call func(nci, ligneci, userdata)
                  do k=1, nci
                      deallocate(ligneci(k)%m_ptr)
                  enddo 
               endif

             enddo
             
            end  subroutine cidist_seq_runmany


      end module mod_cidist_seq

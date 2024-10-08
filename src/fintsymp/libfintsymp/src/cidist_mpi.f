!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file cidist_mpi.f 
!!  \brief Gere la distribution de condition initiale et leur execution.
!!
!!  version parallele mpi uniquement si USE_MPI=1
!!
!!  Pour chaque condition initiale lue par un lecteur de condition initiale, 
!!         la fonction "run" execute la fonction fournie.
!!          execution sequentielle de la fonction "run".
!!
!! 
!!
! history : creation 17/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la distribution d'un fichier texte de conditions initiales quelconques en mpi
!***********************************************************************
      module mod_cidist_mpi
#if USE_MPI
       use mpi
#endif
       use mod_cidist_seq
!***********************************************************************
!> @class t_cidist_mpi
!! classe de distribution via mpi des conditions initiales
!! et execution d'un fichier texte de conditions initiales quelconques
!! repartit la charge entre les noeuds
!!  
!***********************************************************************
      type, extends(t_cidist_seq) :: t_cidist_mpi
          private
          
          integer m_nb_ci !< nombre de conditions initiales que doit gerer le maitre MPI et distribuer aux esclaves
          integer m_sendblocksize !< nombre de conditions initaiaes envoyees a un seul esclave par requete
          integer m_rank !< rang mpi

       contains
          procedure :: set_nb_ci => cidist_mpi_set_nb_ci ! fixe le nombre de conditions initiales a traiter
          procedure :: set_sendblocksize =>cidist_mpi_set_sendblocksize ! nombre de ci a envoyer par tranche aux esclaves
          procedure, private :: run_master =>cidist_mpi_run_master ! fonction run du master
          procedure, private :: run_slave  =>cidist_mpi_run_slave ! fonction run des esclaves
         
          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          procedure :: run => cidist_mpi_run ! execute toutes les conditions initiales du fichier
          procedure :: runmany => cidist_mpi_runmany ! execute toutes les conditions initiales du fichier (plusieurs en meme temps)
          procedure :: init => cidist_mpi_init ! initialise le type
          procedure :: get_rank => cidist_mpi_get_rank ! recupere le rang du distributeur

          
      end type t_cidist_mpi  

      integer, private, parameter :: etiquette_request=111 !< etiquette mpi pour les requetes vers le serveur
      integer, private,  parameter :: etiquette_answer=110 !< etiquette mpi pur les requetes vers les esclaves
    
      contains



!***********************************************************************
!> @brief fixe le nombre de conditions initiales 
!! que doit distribuer le maitre MPI aux esclaves.
!***********************************************************************
            subroutine cidist_mpi_set_nb_ci(this, nbci)
             implicit none
             class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
             integer, intent(in) :: nbci !< nombre de conditions initiales 
             this%m_nb_ci = nbci
            end subroutine cidist_mpi_set_nb_ci

!***********************************************************************
!> @brief fixe le nombre de conditions initiales qui sera envoyee a un 
!! esclave en un seul envoi lorsqu'un esclave reclame des donnees a traiter
!***********************************************************************
            subroutine cidist_mpi_set_sendblocksize(this, blocksize)
             implicit none
             class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
             integer, intent(in) :: blocksize !< nombre de conditions initiales 
             this%m_sendblocksize = blocksize
            end subroutine cidist_mpi_set_sendblocksize
            
!***********************************************************************
!> @brief initialise le type: recupere son rang
!***********************************************************************
            subroutine cidist_mpi_init(this)
             implicit none
             class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
#if USE_MPI
             integer code

             call MPI_COMM_RANK ( MPI_COMM_WORLD ,this%m_rank,code)
             if (this%m_rank.ne.0) then 
             ! esclave
               this%m_ismaster = .false.
             endif 
#endif             
            end  subroutine cidist_mpi_init
      
!***********************************************************************
!> @brief retourne le rang parmi les distributeurs
!***********************************************************************
            function cidist_mpi_get_rank(this) result(r)
             implicit none
             class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
             integer r
             r = this%m_rank
            end  function cidist_mpi_get_rank

!***********************************************************************
!> @brief Fonction executee par le maitre MPI. 
!! Elle  envoie aux esclaves les tranches a traiter lorsque ceux-ci demande des donnees.
!! Le maitre envoie un couple(lig1, lig2) qui sont les indices de la 1ere ligne et de la derniere ligne que devra traiter l'esclave.
!! Les numeros de ligne commence a 1.
!! Si le maitre envoie (-1,-1) alors il n'y a plus de donnee a traiter
!***********************************************************************
            subroutine cidist_mpi_run_master(this, nbslavempi)
             implicit none
             class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
             integer, intent(in):: nbslavempi  !< nombre d'esclaves mpi
            
             integer lig1 ! numero de la prochaine a traiter
             integer lig2 
             integer j
             
#if USE_MPI
             integer code ! code retour mpi
             integer procslave !numero de l'esclave
             integer, dimension(1:2) :: line_send
             integer, dimension( MPI_STATUS_SIZE ) :: statut
#endif
             
             lig1 = 1
             ! boucle sur les donnees a traiter
             do while (lig1<=this%m_nb_ci)
              lig2 = min(lig1+this%m_sendblocksize-1, this%m_nb_ci)
              
#if USE_MPI
              line_send(1) = lig1
              line_send(2) = lig2
              write(*,*) "maitre se met en attente de requete  : ",         &
     &         lig1, "/",  this%m_nb_ci
              call MPI_RECV (procslave,1, MPI_INTEGER , MPI_ANY_SOURCE,     &
     &             etiquette_request, MPI_COMM_WORLD ,statut,code)
              ! Test du code de retour du sous-programme MPI_RECV
               if (code .ne. MPI_SUCCESS) call MPI_ABORT(MPI_COMM_WORLD,    &
     &             2,code)
              write(*,*) "maitre a recu la requete de  :", procslave

              write(*,*) "maitre envoie :", lig1, " a ", lig2
              call MPI_SEND (line_send(1),2, MPI_INTEGER , procslave ,      &
     &          etiquette_answer, MPI_COMM_WORLD, MPI_STATUS_IGNORE,        &
     &          code)
              ! Test du code de retour du sous-programme MPI_SEND
               if (code .ne. MPI_SUCCESS) call MPI_ABORT(MPI_COMM_WORLD,    &
     &             2,code)
#endif
              lig1 = lig2+1
             enddo
             
              write(*,*) "maitre a termine la distribution  :"
             ! il n'y a plus de ci a traiter
             ! le maitre arrete les esclaves en envoyant (-1,-1)
             do j=1, nbslavempi
#if USE_MPI
              line_send(1) = -1
              line_send(2) = -1
              call MPI_RECV (procslave,1, MPI_INTEGER , MPI_ANY_SOURCE,     &
     &             etiquette_request, MPI_COMM_WORLD ,statut,code)
              if (code .ne. MPI_SUCCESS) call MPI_ABORT(MPI_COMM_WORLD,     &
     &             2,code)
              write(*,*) "maitre a recu la requete de  :", procslave
              write(*,*) "maitre arrete l'esclave  :", procslave
              call MPI_SEND (line_send(1),2, MPI_INTEGER , procslave ,      &
     &          etiquette_answer, MPI_COMM_WORLD, MPI_STATUS_IGNORE,        &
     &          code)
              ! Test du code de retour du sous-programme MPI_SEND
              if (code .ne. MPI_SUCCESS) call MPI_ABORT(MPI_COMM_WORLD,     &
     &             2,code)
#endif
            
             enddo
             
         write(*,*) "maitre mpi arrete"
            end  subroutine cidist_mpi_run_master

!***********************************************************************
!> @brief Fonction executee par un esclave MPI. 
!! Elle demande au maitre un bloc de conditions initiales a traiter.
!! Elle recoit un couple(lig1, lig2) qui sont les indices de la 1ere ligne et de la derniere ligne que devra traiter cet esclave.
!! Les numeros de ligne commence a 1.
!! Si le maitre lui a envoye (-1,-1) alors l'esclave s'arrete.
!!
!!  lit chaque condition initiale de la ligne lig1 a lig2  via ciread 
!! puis appelle la fonction func avec la condition initiale lue et les donnees utilisateurs
!! Elle suppose que cireader est initialise correctement (fichier ouvert, ...)
!!
!***********************************************************************
        subroutine cidist_mpi_run_slave(this, cireader, userdata, func)
         implicit none
         class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
         class(t_ciread_base), intent(inout) :: cireader  !<lecteur du fichier de ci
         class(*), intent(inout) :: userdata !< donnee utilisateur
         procedure(t_cidist_userfunc), pointer, intent(in) :: func !< fonction a appeler pour chaque condition initiale

         class(t_ci_base), pointer :: ligneci ! condition initiale lue
         logical notfinish
         integer lineno, lig1, lig2
         integer  procslave, master
#if USE_MPI
         integer  code
         integer, dimension(1:2) :: line_recv
#endif
         
         lig2 = 0
         notfinish = .true.
         procslave = this%m_rank
         
         do while(notfinish)
         ! demande au maitre des donnees
         write(*,*) "esclave", this%m_rank, "demande au maitre "
         master = 0
#if USE_MPI
         call MPI_SENDRECV (procslave,1, MPI_INTEGER, master,                   &
     &     etiquette_request, line_recv(1), 2, MPI_INTEGER, master,             &
     &     etiquette_answer, MPI_COMM_WORLD, MPI_STATUS_IGNORE,code)
     
         if (code .ne. MPI_SUCCESS) call MPI_ABORT(MPI_COMM_WORLD,              &
     &             2,code)
         lig1 = line_recv(1) 
         lig2 = line_recv(2) 
         write(*,*) "esclave", this%m_rank, "a recu", lig1, lig2
#endif
         if (lig1.eq.-1) notfinish = .false.
         ! execute les donnees a traiter
         if (notfinish) then
          call cireader%fgoto_line(lig1) 
          do lineno=lig1, lig2
            if (cireader%freadline(ligneci)) then 
              write(*,*) "slave : lecture d'une nouvelle ci", lineno
              call ligneci%debug()
              call func(ligneci, userdata)
              deallocate(ligneci)
            endif
          enddo
         endif
         end do
         
         write(*,*) "esclave mpi", this%m_rank," arrete"
       end  subroutine cidist_mpi_run_slave

!***********************************************************************
!> @brief lit chaque condition initiale via ciread 
!! puis appelle la fonction func avec la condition initiale lue et les donnees utilisateurs
!! Elle suppose que cireader est initialise correctement (fichier ouvert, ...)
!! Elle distribue entre les noeuds mpi la liste des conditions initiales
!***********************************************************************
        subroutine cidist_mpi_run(this, cireader, userdata, func)
         implicit none
         class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
         class(t_ciread_base), intent(inout) :: cireader  !<lecteur du fichier de ci
         class(*), intent(inout) :: userdata !< donnee utilisateur
         procedure(t_cidist_userfunc), pointer, intent(in) :: func !< fonction a appeler pour chaque condition initiale
             
#if USE_MPI
         integer :: nb_procs,code
         call MPI_COMM_SIZE ( MPI_COMM_WORLD ,nb_procs,code)
         write(*,*) "cidist_mpi_run start on proc", this%m_rank,"/",          &
     &        nb_procs
        if (this%m_rank.ne.0) then 
         ! esclave
         call cidist_mpi_run_slave(this, cireader, userdata, func)
        
        else
         ! maitre
         call cidist_mpi_run_master(this, nb_procs-1)
         
        endif 
#else
        call cidist_seq_run(this, cireader, userdata, func)
#endif             
             
        end  subroutine cidist_mpi_run

!***********************************************************************
!> @brief Fonction executee par un esclave MPI. 
!! Elle demande au maitre un bloc de conditions initiales a traiter.
!! Elle recoit un couple(lig1, lig2) qui sont les indices de la 1ere ligne et de la derniere ligne que devra traiter cet esclave.
!! Les numeros de ligne commence a 1.
!! Si le maitre lui a envoye (-1,-1) alors l'esclave s'arrete.
!!
!!  lit chaque condition initiale de la ligne lig1 a lig2  via ciread 
!! puis appelle la fonction func avec la condition initiale lue et les donnees utilisateurs
!! Elle suppose que cireader est initialise correctement (fichier ouvert, ...)
!!
!***********************************************************************
        subroutine cidist_mpi_runmany_slave(this, cireader, userdata,      &
     &             func)
         implicit none
         class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
         class(t_ciread_base), intent(inout) :: cireader  !<lecteur du fichier de ci
         class(*), intent(inout) :: userdata !< donnee utilisateur
         procedure(t_cidist_manyuserfunc),pointer,intent(in):: func !< fonction a appeler pour chaque condition initiale

         logical notfinish
         integer lineno, lig1, lig2
         integer  procslave, master
         type(t_ci_base_ptr), dimension(1:this%m_manyblocksize) ::          &
     &                    ligneci ! tableau de condition initiale lue
         logical bcontinue ! =false si la fin est atteinte
         integer nci ! nombre de ci lue
         integer k
#if USE_MPI
         integer  code
         integer, dimension(1:2) :: line_recv
#endif
         
         notfinish = .true.
         procslave = this%m_rank
         
         do while(notfinish)
         ! demande au maitre des donnees
         write(*,*) "esclave", this%m_rank, "demande au maitre "
         master = 0
#if USE_MPI
         call MPI_SENDRECV (procslave,1, MPI_INTEGER, master,                   &
     &     etiquette_request, line_recv(1), 2, MPI_INTEGER, master,             &
     &     etiquette_answer, MPI_COMM_WORLD, MPI_STATUS_IGNORE,code)
     
         if (code .ne. MPI_SUCCESS) call MPI_ABORT(MPI_COMM_WORLD,              &
     &             2,code)
         lig1 = line_recv(1) 
         lig2 = line_recv(2) 
         write(*,*) "esclave", this%m_rank, "a recu", lig1, lig2
#endif
         if (lig1.eq.-1) notfinish = .false.
         ! execute les donnees a traiter
         if (notfinish) then
             call cireader%fgoto_line(lig1) 
             nci = 0
             do lineno=lig1, lig2
               write(*,*) "avcidist_mpi_runmany_slave : ",lig1,  lig2
               nci = nci+1
               bcontinue = cireader%freadline(ligneci(nci)%m_ptr)
               write(*,*) "slave : ",lineno, bcontinue
               write(*,*) "cidist_mpi_runmany_slave : ",lig1,  lig2
             enddo
             
             write(*,*) "cidist_mpi_runmany_slave : traitement de ",nci,  &
     &                       "ci"
          
            do k=1, nci
                 call ligneci(k)%m_ptr%debug()
            enddo 
            call func(nci, ligneci, userdata)
            do k=1, nci
                deallocate(ligneci(k)%m_ptr)
            enddo 
         endif
         end do
         
         write(*,*) "esclave mpi", this%m_rank," arrete"
       end  subroutine cidist_mpi_runmany_slave

!***********************************************************************
!> @brief lit au plus "m_manyblocksize" conditions initiales via ciread 
!! puis appelle la fonction func avec le tableau des conditions initiales lues 
!! et les donnees utilisateurs
!! Elle suppose que cireader est initialise correctement (fichier ouvert, ...)
!! Elle distribue entre les noeuds mpi la liste des conditions initiales
!***********************************************************************
         subroutine cidist_mpi_runmany(this,cireader,userdata,func)
         implicit none
         class(t_cidist_mpi), intent(inout):: this  !< distributeur de ci
         class(t_ciread_base), intent(inout) :: cireader  !<lecteur du fichier de ci
         class(*), intent(inout) :: userdata !< donnee utilisateur
         procedure(t_cidist_manyuserfunc),pointer,intent(in):: func !< fonction a appeler pour chaque condition initiale
             
#if USE_MPI
         integer :: nb_procs,code
         call MPI_COMM_SIZE ( MPI_COMM_WORLD ,nb_procs,code)
         write(*,*) "cidist_mpi_runmany start on proc", this%m_rank,         &
     &        "/", nb_procs
        if (this%m_rank.ne.0) then 
         ! esclave
         call cidist_mpi_runmany_slave(this, cireader, userdata, func)
        
        else
         ! maitre
         call cidist_mpi_run_master(this, nb_procs-1)
         
        endif 
#else
        call cidist_seq_runmany(this, cireader, userdata, func)
#endif             
             
        end  subroutine cidist_mpi_runmany
      
      end module mod_cidist_mpi

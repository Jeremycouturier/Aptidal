!-----------------------------------------------------------------------------------
!
! fichier : GETPARFILENAME.F
!
! auteur : M. Gastineau
!
! M. GASTINEAU 14/01/08 : version 1.0
!
!-----------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------
! retourne le fichier de parametres a utiliser en inspectant la ligne de commandes
! pour les appels du type : 
!
!       monapplication monfichierparam.par
!
!  si la ligne de commande comprend : 0 argument => parfilename = defaultfilename
!                                     1 argument => parfilename = cet argument
!  en cas d'erreur ou nombre d'arguments >1 : arret avec message d'erreur
!
!
!  @param defaultfilename (in) nom par defaut du fichier de parametre
!  @param parfilename (out) nom du fichier de parametre retourne
!
! 
!-----------------------------------------------------------------------------------
       subroutine GETPARFILENAME(defaultfilename, parfilename)
           implicit none
           character(len=*), intent(in)  :: defaultfilename
           character(len=*), intent(out)  :: parfilename
           
           integer status, length
           integer cmd_count
           logical fatalerror
           
           fatalerror = .false.
           
           ! recupere le nombre d'arguments
           cmd_count = COMMAND_ARGUMENT_COUNT()
           
           if (cmd_count.eq.1) then
             ! recupere le 1er argument
             call GET_COMMAND_ARGUMENT(1, parfilename, length, status)
             if (status.ne.0) then
               fatalerror = .true.
             endif
           else if (cmd_count.eq.0) then
               parfilename = defaultfilename
           else
               fatalerror = .true.
           endif
           
           if (fatalerror.eqv..true.) then
             write(*,*) 
             write(*,*) "Usage : <nom du programme> ",defaultfilename 
             stop
           endif
                      
       return
       end

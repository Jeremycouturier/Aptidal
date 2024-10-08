!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file io_dump.f 
!!  \brief dump binaire de l'etat du systeme
!!
!!     cette classe t_io_dump gere le dump et la restauration
!!
!!
! history : creation 01/08/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture des integrales premieres dans un fichier
!***********************************************************************
      module mod_io_dump
       
!***********************************************************************
!> @class t_io_dump
!! classe d'ecriture des "dump" dans un fichier binaire
!!  
!***********************************************************************
      type :: t_io_dump

          integer(8), public :: m_outputstep = 0 !< pas de sortie du dump 
          integer, public :: m_nf !< numero fortran du fichier 
          character(len=512), private :: m_name !< nom du fichier

       contains
          
          procedure :: set_name =>io_dump_set_name ! fixe le nom du fichier et le numero de fichier
          procedure :: fcreate => io_dump_fcreate  ! cree le fichier et l'ouvre
          procedure :: fopen => io_dump_fopen      ! ouvre le fichier en lecture
          procedure :: fclose  => io_dump_fclose   ! ferme le fichier
          procedure :: set_output => io_dump_set_output !< fixe le pas de sortie
          procedure :: isoutputstep => io_dump_isoutputstep !< determine si le pas fourni est un pas de sortie
                  
      end type t_io_dump  

     
      contains
         

!***********************************************************************
!> @brief fixe le nom du fichier et le numero fortran
!***********************************************************************
            subroutine io_dump_set_name(this, name, nf)
             implicit none
             class(t_io_dump), intent(inout):: this  !< dummy argument
             character(len=*), intent(in) :: name  !< nom du fichier
             integer, intent(in) :: nf  !< numero fortran du fichier
             
             this%m_name = name
             this%m_nf = nf
            end  subroutine io_dump_set_name

!***********************************************************************
!> @brief cree le fichier et l'ouvre
!***********************************************************************
            subroutine io_dump_fcreate(this)
             implicit none
             class(t_io_dump), intent(inout):: this  !< dummy argument
             open(this%m_nf,file=this%m_name,status='unknown',  form='unformatted')
            end  subroutine io_dump_fcreate

!***********************************************************************
!> @brief cree le fichier et l'ouvre
!***********************************************************************
            subroutine io_dump_fopen(this)
                  implicit none
                  class(t_io_dump), intent(inout):: this  !< dummy argument
                  open(this%m_nf,file=this%m_name,status='old',  form='unformatted')
            end  subroutine io_dump_fopen
     
!***********************************************************************
!> @brief ferme le fichier
!***********************************************************************
            subroutine io_dump_fclose(this)
             implicit none
             class(t_io_dump), intent(inout):: this  !< dummy argument
             close(this%m_nf)
            end  subroutine io_dump_fclose

!***********************************************************************
!> @brief fixe le pas du dump
!***********************************************************************
            subroutine io_dump_set_output(this, outputstep)
                  implicit none
                  class(t_io_dump), intent(inout):: this  !< dummy argument
                  integer(8), intent(in) :: outputstep !< pas de sortie
                  this%m_outputstep  = outputstep
            end  subroutine io_dump_set_output
     
!***********************************************************************
!> @brief retourne true si curstep est multiple de m_outputstep
!***********************************************************************
            function io_dump_isoutputstep(this, curstep) result(mult)
              implicit none
              class(t_io_dump), intent(in):: this  !< dummy argument
              integer(8), intent(in) :: curstep !< pas de sortie
              logical :: mult
              mult =.false.
              if (this%m_outputstep.ne.0) then
                  if (mod(curstep, this%m_outputstep).eq.0) then
                        mult = .true.
                  endif
              endif
            end  function io_dump_isoutputstep

      end module mod_io_dump

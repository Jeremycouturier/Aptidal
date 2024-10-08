!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file schema.f 
!!  \brief gestion des schema d'integration 
!!
!!  il faut appeller la fonction 'create_schema' pour creer le schema d'integration
!!
!! toutes les autres routines sont privees
!
! history : creation 25/07/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour les schemas d'integrations
!***********************************************************************
      module mod_schema
          use mod_syssympbase

          integer, private,  parameter :: coef_max = 100 !< nombre de max pour les coefficients pour les integrateurs
          private create_schemaABA_an
          private create_schemaABA_fill
          private create_schemaABAH
          ! private create_schemaABAJ
          private create_schemaBABH
          private create_schemaBABJ
         
!***********************************************************************
!> @class t_schema
!! classe de base decrivant 1 schema d'integration
!!
!***********************************************************************
      type t_schema
          
          integer, private  :: m_ordre    !< ordre de l'integrateur
          integer, private  :: m_if_cor    !<  =1 => correcteur present
          real(TREAL), private :: m_coef_depart !< coefficient pour le demarrage et pas de sortie
          real(TREAL), private :: m_cor         !< correcteur
          real(TREAL), private,allocatable,dimension(:) :: cb,db !< valeur des coefficients

          procedure(procrun), public, pointer :: m_procdepart !<fonction appellee au depart ou lors du pas de sortie
          procedure(procrun), public, pointer :: m_procrun  !< fonction appellee lors d'un pas normal

         logical, private  :: m_helio    !<  =true => schema helio, = false => schema jacobi

       contains
          private

          procedure :: t_schema_pasdepartA  !< fonction de depart pas A
          procedure :: t_schema_pasdepartB  !< fonction de depart pas B
          procedure :: t_schema_pasdepartC  !< fonction de depart pas C
          procedure :: t_schema_pasAB  !< fonction pas AB
          procedure :: t_schema_pasBA  !< fonction pas BA
          procedure :: t_schema_pasBAB  !< fonction pas BAB
          procedure :: t_schema_pasdepartvide  !< fonction de depart vide

      end type t_schema  

           abstract interface
            subroutine procrun(this, system, dt)
            import t_schema
            import t_syssymp
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in) :: dt ! pas d'inetgration
            end subroutine procrun

          end interface
      
      contains
      
!***********************************************************************
!> fonction qui appelle rien du tout de system au demarrage
!***********************************************************************
           subroutine t_schema_pasdepartvide(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in)  :: dt ! pas d'integration
            
            end subroutine t_schema_pasdepartvide

!***********************************************************************
!> fonction qui appelle pasA de system au demarrage
!***********************************************************************
           subroutine t_schema_pasdepartA(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in)  :: dt ! pas d'integration
            
            call system%pasA(this%m_coef_depart*dt) 
            end subroutine t_schema_pasdepartA

!***********************************************************************
!> fonction qui appelle pasB de system au demarrage
!***********************************************************************
           subroutine t_schema_pasdepartB(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in)  :: dt ! pas d'integration
            
            call system%pasB(this%m_coef_depart*dt) 
            end subroutine t_schema_pasdepartB


!***********************************************************************
!> fonction qui appelle pasC de system au demarrage
!***********************************************************************
           subroutine t_schema_pasdepartC(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in)  :: dt ! pas d'integration
            
            call system%pasC(this%m_coef_depart*dt) 
            end subroutine t_schema_pasdepartC

!***********************************************************************
!> fonction qui appelle AB de system selon l'ordre
!***********************************************************************
           subroutine t_schema_pasAB(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in) :: dt ! pas d'integration
            integer :: i
            
            do i=1,this%m_ordre
               call system%pasA(this%cb(i)*dt)
               call system%pasB(this%db(i)*dt)
            end do
            end subroutine t_schema_pasAB
            
!***********************************************************************
!> fonction qui appelle BA de system selon l'ordre
!***********************************************************************
           subroutine t_schema_pasBA(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in) :: dt ! pas d'integration
            integer :: i
            
            do i=1,this%m_ordre
               call system%pasB(this%db(i)*dt)
               call system%pasA(this%cb(i)*dt)
            end do
            end subroutine t_schema_pasBA

!***********************************************************************
!> fonction qui appelle BAB de system selon l'ordre 
!***********************************************************************
           subroutine t_schema_pasBAB(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in) :: dt ! pas d'integration
            integer :: i
            
            call system%pasB(this%db(1)*dt)
            do i=2,this%m_ordre+1
               call system%pasA(this%cb(i)*dt)
               call system%pasB(this%db(i)*dt)
            end do
            end subroutine t_schema_pasBAB

!***********************************************************************
!> fonction qui appelle ABC de system selon l'ordre 
!***********************************************************************
           subroutine t_schema_pasABC(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in) :: dt ! pas d'integration
            integer :: i
            
            do i=1,this%m_ordre
               call system%pasA(this%cb(i)*dt)
               call system%pasB(this%db(i)*dt)
            end do
            call system%pasA(this%cb(this%m_ordre+1)*dt)
            call system%pasC(-this%m_cor*dt)
            end subroutine t_schema_pasABC

!***********************************************************************
!> fonction qui appelle BAC de system selon l'ordre 
!***********************************************************************
           subroutine t_schema_pasBAC(this, system, dt)
            implicit none
            class(t_schema), intent(in) :: this
            class(t_syssymp), intent(inout) :: system ! systeme qui est integre (doit fournir: pasA, pasB et pas C)
            real(TREAL), intent(in) :: dt ! pas d'integration
            integer :: i
            
            do i=1,this%m_ordre
               call system%pasB(this%db(i)*dt)
               call system%pasA(this%cb(i)*dt)
            end do
            call system%pasB(this%db(this%m_ordre+1)*dt)
            call system%pasC(-this%m_cor*dt)
            end subroutine t_schema_pasBAC


!***********************************************************************
!> creation d'un schema d'integrateur
!!                 ('ABA' ou 'BAB') +['H'] + ['C'] + nombre a 1 ou 2 chiffres
!!                 Remarque:
!!                 'ABA' ou 'BAB' indique la sequence des pas.
!!                 'H' indique ce sont des variables heliocentriques
!!                 'C' indique la presence ou non du correcteur
!!                 le nombre indique l'ordre de l'integrateur
!!  valeurs possibles de schema :
!!  'ABA1' ou 'BABC10'
!!  'ABA2' ou 'BAB2'
!***********************************************************************
      subroutine create_schema(nom, sc)
      implicit none
      character(len=*), intent(in) :: nom !< nom du schema
      type(t_schema), intent(out) :: sc   !< schema cree
      integer :: l1 
      real(TREAL),dimension(coef_max) :: c,d
      real(TREAL) :: cor, coef_depart
      
      sc%m_if_cor = 0
      sc%m_ordre = 0
      
      c = 0.d0
      d = 0.d0
      cor = 0
      coef_depart = 0
      
      l1 = len_trim(nom)
      
      ! cas de l'integrateur ABA en helio
      if (nom(1:4).eq.'ABAH') then
        call create_schemaABAH(nom, sc)
      ! cas de l'integrateur BAB en helio
      else if (nom(1:4).eq.'BABH') then
        call create_schemaBABH(nom, sc)
      ! cas de l'integrateur ABA en jacobi
      else if (nom(1:3).eq.'ABA') then
         call create_schemaABAJ(nom, sc)
      ! cas de l'integrateur BAB en jacobi
      else if (nom(1:3).eq.'BAB') then
        call create_schemaBABJ(nom, sc)
      endif  

      end subroutine create_schema

!***********************************************************************
!> creation d'un schema d'integrateur de type ABA commun en jacobi ou helio
! a partir de numero: remplit c,d, cor et ordre
!***********************************************************************
      subroutine create_schemaABA_an(numero, c,d, cor, ordre)
      implicit none
      integer, intent(in) :: numero  !< numero de l'integrateur
      real(TREAL),dimension(coef_max), intent(out) :: c !< coefficient pour le pas A
      real(TREAL),dimension(coef_max), intent(out) :: d !< coefficient pour le pas B
      real(TREAL), intent(out) :: cor !< coefficient pour le pas C
      integer, intent(out) :: ordre !<ordre de l'integrateur (nombre d'etapes)
      real(TREAL) :: aux0, aux1

      c = 0.d0
      d = 0.d0
      cor = 0
      
      ordre = numero
      
      select case (numero)
      
      case (1)
         d(1) = 1d0
         c(1) = .5d0
         cor = 1.d0/12.d0
      
      case (2)
         c(2) = sqrt(3d0)/3d0
         c(1) = .5d0*(1d0 - c(2))
         d(1) = .5d0
         cor = 0.0111645496846301127696897357705886513774d0
      
      case (3)
         d(2) = 4d0/9d0
         d(1) = 5d0/18d0
         c(2) = sqrt(15d0)/10d0
         c(1) = .5d0 - c(2)
         cor = 0.0056345933631228094022678237697975386716d0
      
      case (4)
         aux0 = sqrt(525d0 + 70d0*sqrt(30d0))
         aux1 = sqrt(525d0 - 70d0*sqrt(30d0))
         d(2)  = 1d0/4d0 + sqrt(30d0)/72d0
         d(1)  = 1d0/4d0 - sqrt(30d0)/72d0
         c(3) = aux1/35d0
         c(2)  = aux0/70d0 - aux1/70d0
         c(1)  = .5d0 - aux0/70d0     
         cor = 0.0033967750482086013315321577834921437970d0
      
      case (5)
         aux0 = sqrt(490d0 - 42d0*sqrt(105d0))
         aux1 = sqrt(490d0 + 42d0*sqrt(105d0))
         d(3) = 64d0/225d0
         d(2) = 13d0/1800d0*sqrt(70d0) + 161d0/900d0
         d(1) = -13d0/1800d0*sqrt(70d0) + 161d0/900d0
         c(3) = -aux0/84d0 + aux1/84d0
         c(2) = aux0/42d0
         c(1) = .5d0 -aux0/84d0 - aux1/84d0
         cor = 0.0022705431214192648194349550500391298793d0
      
      case (6)
         c(1)=0.033765242898423986093849222753002695432617131143855d0
         c(2)=0.13563006386844375707545097973704463106415858665856d0
         c(3)=0.21129510019153380251544893666959670579391896712757d0
         c(4)=0.23861918608319690863050172168071193541861063014002d0
         d(1)=0.085662246189585172520148071086366446763411250742022d0
         d(2)=0.18038078652406930378491675691885805583076094637337d0
         d(3)=0.23395696728634552369493517199477549740582780288461d0
         cor = 0.0016244598416242825214522585124636070897d0
      
      case (7)
         c(1)=0.025446043828620737736905157976074368799614531164691d0
         c(2)=0.10378836337168204233116245538353142766331164526461d0
         c(3)=0.16784301711099863647862918060191347186338281652101d0
         c(4)=0.20292257568869858345330320603848073167369100704969d0
         d(1)=0.064742483084434846635305716339541009164293701129973d0
         d(2)=0.1398526957446383339507338857118897912434625326133d0
         d(3)=0.19091502525255947247518488774448756693918254176693d0
         d(4)=256.d0/1225.d0
         cor = 0.0012196439127604184725792118223316450935d0
      
      case (8)
         c(1)=0.019855071751231884158219565715263504785882382849274d0 
         c(2)=0.081811689541954746046003466046821276795531751342744d0 
         c(3)=0.13556703374864887688690744364329204389760374424779d0 
         c(4)=0.17104888371033959043913145341453118418760321510371d0 
         c(5)=0.18343464249564980493947614236018398066675781291297d0 
         d(1)=0.050614268145188129576265677154981095057697045525842d0 
         d(2)=0.11119051722668723527217799721312044221506543502562d0 
         d(3)=0.15685332293894364366898110099330065663016449950137d0 
         d(4)=0.18134189168918099148257522463859780609707301994717d0 
         cor = 0.0009493081777456022347921775035350542471d0
          
      case (9)
         c(1)=0.015919880246186955082211898548163564975297599754037d0 
         c(2)=0.066064566090495147768073207416968996752649041183625d0 
         c(3)=0.11132983731302269849536387436413034587919305659206d0 
         c(4)=0.1445590046483907341350820123490687881068821721018d0 
         c(5)=0.16212671170190446451926900732166830428597813036849d0 
         d(1)=0.040637194180787205985946079055261825337830860391205d0 
         d(2)=0.09032408034742870202923601562145640475716891086602d0 
         d(3)=0.13030534820146773115937143470931642488592010221865d0 
         d(4)=0.15617353852000142003431520329222183279937743063095d0 
         d(5)=16384.d0 /99225.d0 
         cor = 0.0007598460228604366463581966741768150274d0
          
      case (10)
          c(1)=0.013046735741414139961017993957773973285865026653809d0 
          c(2)=0.054421580914093604672933661830479502450363465863526d0 
          c(3)=0.09282689919498005224888466165430973637912415156526d0 
          c(4)=0.12300708708488860771753071097454470678461143708928d0 
          c(5)=0.14226052757380798995721997101803208879125348911778d0 
          c(6)=0.14887433898163121088482600112971998461756485942069d0 
          d(1)=0.033335672154344068796784404946665896428932417160079d0 
          d(2)=0.074725674575290296572888169828848666201278319834714d0 
          d(3)= 0.10954318125799102199776746711408159622938593526134d0 
          d(4)= 0.13463335965499817754561346078473467642987996923044d0 
          d(5)= 0.14776211235737643508694649732566916471052335851343d0
          cor =  0.0006219343314861664264970498453586457591d0

      case (82)
!*****************************************************************
!     Laskar SABA_4 and McLahan (8,2) integrating schemes
!     - 4 stages / symetric / NO corrector
!      int_type ABA82
!*****************************************************************

         aux0 = sqrt(REALC(525e0) + REALC(70e0)*sqrt(REALC(30e0)))
         aux1 = sqrt(REALC(525e0) - REALC(70e0)*sqrt(REALC(30e0)))
   
         d(1) = REALC(1e0)/REALC(4e0) - sqrt(REALC(30e0))/REALC(72e0)
         d(2) = REALC(1e0)/REALC(4e0) + sqrt(REALC(30e0))/REALC(72e0)
         d(3) = d(2)
         d(4) = d(1)
   
         c(1) = REALC(.5e0) - aux0/REALC(70e0)     
         c(2) = aux0/REALC(70e0) - aux1/REALC(70e0)
         c(3) = aux1/REALC(35e0)
         c(4) = c(2)
         c(5) = REALC(2.e0)*c(1)
   
         ordre = 4

      case (84)
!***********************************************************************
! initialize  structure for SABA McLa(8,4) McLahlan coefficients
!   coefficients from table 3 - from ariadna report 28/07/2011 : integNBP.pdf
!  McLahan (8,4) integrating schemes
! 		- 5 stages
!		- symetric 
! 		- NO corrector
!
!***********************************************************************
          d(1)  = REALC(0.19022593937367661925e0)
          d(2)  = REALC(0.84652407044352625706e0)
          d(3)  = REALC(-1.07350001963440575260e0)
          c(1)  = REALC(0.07534696026989288842e0)    
          c(2)  = REALC(0.51791685468825678230e0)
          c(3)  = REALC(-0.09326381495814967072e0)
          ordre = 5
          
      case (844)
!***********************************************************************
!  Blanes (8,4,4) integrating schemes
!  - 6 stages / symetric / extra stage for killig  tau^3
!   int_type ABAH844
!***********************************************************************

          d(1) = REALC(0.6408857951625127178e0);  
          d(2) = REALC(-0.8585754489567828567e0);
          d(3) = REALC(0.7176896537942701389e0);  
          d(4) = d(3)
          d(5) = d(2)
          d(6) = d(1)
    
     
          c(1) = REALC(0.2741402689434018762e0);
          c(2) = REALC(-0.10756843844016423066e0);
          c(3) = REALC(-0.048018502590601692667e0);
          c(4) = REALC(0.7628933441747280943e0);
          c(5) = c(3)
          c(6) = c(2)
          c(7) = REALC(2e0)*c(1)
          ordre = 6
          
!***********************************************************************
!  MacLachlan general integrating schemes
!  - 5 stages / symetric / order tau^5 independent of epsilon
!   int_type ABAH45
!***********************************************************************
       case (45)
          d(1) = 2.d0/5.d0
          d(2) = -1.d0/10.d0
          d(3) = d(1)
          aux0 = sqrt(19.d0)
          c(1) = (14.d0-aux0)/108.d0
          c(2) = (20.d0-7*aux0)/108.d0
          c(3) = (20.d0+8*aux0)/108.d0
          ordre =5

          
!***********************************************************************
!  MacLachlan general integrating schemes
!  - 7 stages / symetric / order tau^7 independent of epsilon
!   int_type ABAH67
!***********************************************************************
       case (67)
          d(1) = REALC( 7.84513610477557263820e-1)
          d(2) = REALC( 2.35573213359358133680e-1)
          d(3) = REALC(-1.17767998417887100695e0)
          d(4) = REALC( 1.31518632068391121890e0)

          c(1) = REALC( 3.92256805238778631901e-1)
          c(2) = REALC( 5.10043411918457698750e-1)
          c(3) = REALC(-4.71053385409756436635e-1)
          c(4) = REALC( 6.87531682525201059750e-2)
          ordre = 7

          
!***********************************************************************
!  MacLachlan general integrating schemes
!  - 9 stages / symetric / order tau^7 independent of epsilon
!   int_type ABAH69
!***********************************************************************   
       case (69)
          d(1) = REALC( 1.86700000000000000000e-1)
          d(2) = REALC( 5.55497023712478399160e-1)
          d(3) = REALC( 1.29466948913475358060e-1)
          d(4) = REALC(-8.43265623387734608550e-1)
          d(5) = REALC( 9.43203301523561702660e-1)

          c(1) = REALC( 9.33500000000000000000e-2)
          c(2) = REALC( 3.71098511856239199580e-1)
          c(3) = REALC( 3.42481986312976878610e-1)
          c(4) = REALC(-3.56899337237129625245e-1)
          c(5) = REALC( 4.99688390679135470550e-2)
          ordre = 9
          
!***********************************************************************
!  MacLachlan general integrating schemes
!  - 15 stages / symetric / order tau^9 independent of epsilon
!   int_type ABAH815
!***********************************************************************   
          
       case (815)
          d(1) = REALC( 7.41670364350612953450e-1)
          d(2) = REALC(-4.09100825800031594000e-1)
          d(3) = REALC( 1.90754710296238379950e-1)
          d(4) = REALC(-5.73862471116082266660e-1)
          d(5) = REALC( 2.99064181303655923840e-1)
          d(6) = REALC( 3.34624918245298183780e-1)
          d(7) = REALC( 3.15293092396766596630e-1)
          d(8) = REALC(-7.96887939352916353980e-1)

          c(1) = REALC( 3.708351821753064767250e-1)
          c(2) = REALC( 1.662847692752906797250e-1)
          c(3) = REALC(-1.091730577518966070250e-1)
          c(4) = REALC(-1.915538804099219433550e-1)
          c(5) = REALC(-1.373991449062131714100e-1)
          c(6) = REALC( 3.168445497744770538100e-1)
          c(7) = REALC( 3.249590053210323902050e-1)
          c(8) = REALC(-2.407974234780748786750e-1)

          ordre = 15
          
        case default
         write(*,*)"l'integrateur ABA ",numero,
     &              " n'est pas au point "
         stop
        end select ! fin de l'ordre
      
      end subroutine create_schemaABA_an

!***********************************************************************
!> remplissage de sc a partir de c,d, cor et ordre pour un schema ABA
!***********************************************************************
      subroutine create_schemaABA_fill(sc, c,d, cor, ordre)
      implicit none
      type(t_schema), intent(inout) :: sc   !< schema cree
      real(TREAL),dimension(coef_max), intent(in) :: c !< coefficient pour le pas A
      real(TREAL),dimension(coef_max), intent(in) :: d !< coefficient pour le pas B
      real(TREAL), intent(in) :: cor !< coefficient pour le pas C
      integer, intent(in) :: ordre !<ordre de l'integrateur (nombre d'etapes)
      integer :: i,j,n
      real(TREAL),dimension(coef_max) ::  cb,db
      real(TREAL) coef_depart
      
         coef_depart = 0

         sc%m_ordre = ordre
         
         write (*,*) "nombre etapes :", sc%m_ordre
        
!--------- Initalisation des coeff cb et db      

         cb = 0.d0
         db = 0.d0
      

       if (mod(ordre,2).eq.0) then 

         n = ordre/2
         if (sc%m_if_cor.eq.0) then 
            do j=1,n 
               db(j) = d(j)
               db(n+j) = d(n-j+1)
               cb(j) = c(j+1)
               cb(n+j) = c(n-j+1)
            end do
            cb(2*n) = REALC(2.E0)*c(1)
            coef_depart = c(1)
            write(*,*)'          i     cb                      db'
            do i=1,ordre
               write(*,*)i,cb(i),db(i)
            end do
            write(*,*)'coef_depart',  coef_depart
            write(*,*)sum(cb),sum(db)
 
         else 
           do  j=1,n 
             db(j) = d(j)
             db(n+j) = d(n-j+1)
             cb(j) = c(j)
             cb(n+j+1) = c(n-j+1)
           end do
           cb(n+1) = c(n+1)        
           coef_depart = -cor/REALC(2.E0) ! pour utiliser coef_depart systematiquement
            
           write(*,*)'          i     cb                      db'
           do i=1,ordre
              write(*,*)i,cb(i),db(i)
           end do
           write(*,*)ordre+1,cb(ordre+1)
           write(*,*)'cor ',  cor
           write(*,*)'coef_depart',  coef_depart
           write(*,*)sum(cb),sum(db)
         end if  
         
       else if (mod(ordre,2).eq.1) then
         n = (ordre-1)/2
         if (sc%m_if_cor.eq.0) then 
           do j=1,n 
               db(j) = d(j)
               db(n+j+1) = d(n-j+1)
               cb(j) = c(j+1)
               cb(n+j+1) = c(n-j+1)
           end do
           cb(n+1) = c(n+1)
           db(n+1) = d(n+1)
           cb(ordre) = REALC(2.E0)*c(1)
           coef_depart = c(1)
           write(*,*)'          i     cb                      db'
           do i=1,ordre
              write(*,*)i,cb(i),db(i)
           end do
           write(*,*)'coef_depart',  coef_depart
           write(*,*)sum(cb),sum(db)
         else 
           do j=1,n 
               db(j) = d(j)
               db(n+j+1) = d(n-j+1)
           end do
           db(n+1) = d(n+1)
           do j=1,n+1 
               cb(j) = c(j)
               cb(n+j+1) = c(n-j+2)
           end do
           coef_depart = -cor/REALC(2.E0) ! pour utiliser coef_depart systematiquement

           write(*,*)'          i     cb                      db'
           do i=1,ordre
              write(*,*)i,cb(i),db(i)
           end do
           write(*,*)ordre+1,cb(ordre+1)
           write(*,*)'coef_depart',  coef_depart
           write(*,*)sum(cb),sum(db)  
         end if
       endif

         allocate(sc%cb(ordre+1))
         allocate(sc%db(ordre))
         sc%cb = cb
         sc%db = db
      
      sc%m_cor = cor
      sc%m_coef_depart = coef_depart

      end subroutine create_schemaABA_fill

!***********************************************************************
!> creation d'un schema d'integrateur de type ABA
!!                 'ABA' + ['C'] + nombre a 1 ou 2 chiffres
!!                 Remarque:
!!                 'ABA'  indique la sequence des pas.
!!                 'C' indique la presence ou non du correcteur
!!                 le nombre indique l'ordre de l'integrateur
!!  valeurs possibles de schema :
!!  'ABA4' ou 'ABAC4'
!***********************************************************************
      subroutine create_schemaABAH(nom, sc)
      implicit none
      character(len=*), intent(in) :: nom !< nom du schema
      type(t_schema), intent(out) :: sc   !< schema cree
      integer :: l1, ordre,  numero
      real(TREAL),dimension(coef_max) :: c,d
      real(TREAL) :: cor
      
      sc%m_if_cor = 0
      sc%m_ordre = 0
      
      c = 0.d0
      d = 0.d0
      cor = 0
      
      l1 = len_trim(nom)
      
      if (nom(5:5).eq.'C') then
       sc%m_if_cor = 1
       read(nom(6:l1),*) ordre
       sc%m_procdepart => t_schema_pasdepartC
       sc%m_procrun => t_schema_pasABC
      else
       read(nom(5:l1),*) ordre
       sc%m_procdepart => t_schema_pasdepartA
       sc%m_procrun => t_schema_pasBA
      endif

      write(*,*)'integrateur ',nom
      write (*,*) "if_cor :", sc%m_if_cor
      write (*,*) "ordre :", ordre
      
      
      select case (ordre)
      
!***********************************************************************
! integrateur communaux helio et jacobi
!***********************************************************************
         case default
           numero = ordre
           call create_schemaABA_an(numero, c,d, cor, ordre)

!***********************************************************************
! les integrateurs suivants sont specifiques aux heliocentriques
!***********************************************************************


      case (864)
!***********************************************************************
!  Blanes (8,6,4) integrating schemes
!  - 8 stages / symetric / extra stage for killig  tau^3
!   int_type ABAH864
!***********************************************************************

          d(1) = REALC(0.1684432593618954534310382697756917558148)
          d(2) = REALC(0.4243177173742677224300351657407231801453)
          d(3) = REALC(-0.5858109694681756812309015355404036521923)
          d(4) = REALC(0.4930499927320125053698281000239887162321)
          d(5) = REALC(0.4930499927320125053698281000239887162321)
          d(6) = REALC(-0.5858109694681756812309015355404036521923)
          d(7) = REALC(0.4243177173742677224300351657407231801453)
          d(8) = REALC(0.1684432593618954534310382697756917558148)

          c(1) = REALC(0.06810235651658372084723976682061164571212)
          c(2) = REALC(0.2511360387221033233072829580455350680082)
          c(3) = REALC(-0.07507264957216562516006821767601620052338)
          c(4) = REALC(-0.009544719701745007811488218957217113269121)
          c(5) = REALC(0.5307579480704471776340674235341732001443)
          c(6) = REALC(-0.009544719701745007811488218957217113269121)
          c(7) = REALC(-0.07507264957216562516006821767601620052338)
          c(8) = REALC(0.2511360387221033233072829580455350680082)
          c(9) = REALC(2.e0)*c(1) !!REALC(0.06810235651658372084723976682061164571)

          ordre = 8

      case (1064)
************************************************************************
*  Blanes (10,6,4) integrating schemes
*  - 9 stages / symetric / extra stage for killig  tau^3
*   int_type ABAH1064
************************************************************************
 
          c(1) = REALC(0.04731908697653382270404371796320813250988e0)
          c(2) = REALC(0.2651105235748785159539480036185693201078e0)
          c(3) = REALC(-0.009976522883811240843267468164812380613143e0)
          c(4) = REALC(-0.05992919973494155126395247987729676004016e0)
          c(5) = REALC(0.2574761120673404534492282264603316880356e0)
          c(6) = c(5)
          c(7) = c(4)
          c(8) = c(3)
          c(9) = c(2)
          c(10) = REALC(2.e0)*c(1)

          d(1) = REALC(0.1196884624585322035312864297489892143852e0)
          d(2) = REALC(0.3752955855379374250420128537687503199451e0)
          d(3) = REALC(-0.4684593418325993783650820409805381740605e0)
          d(4) = REALC(0.3351397342755897010393098942949569049275e0)
          d(5) = REALC(0.2766711191210800975049457263356834696055e0)
          d(6) = d(4)
          d(7) = d(3)
          d(8) = d(2)
          d(9) = d(1)
        
          ordre = 9

         end select ! fin de l'ordre
      
      
       ! remplissage final des coefs et ordre  de sc
       call create_schemaABA_fill(sc, c,d, cor, ordre)
      
      end subroutine create_schemaABAH
      
!***********************************************************************
!> creation d'un schema d'integrateur de type BAB commun en jacobi ou helio
! a partir de numero: remplit c,d, cor et ordre
!***********************************************************************
      subroutine create_schemaBAB_an(numero, c,d, cor, ordre)
      implicit none
      integer, intent(in) :: numero  !< numero de l'integrateur
      real(TREAL),dimension(coef_max), intent(out) :: c !< coefficient pour le pas A
      real(TREAL),dimension(coef_max), intent(out) :: d !< coefficient pour le pas B
      real(TREAL), intent(out) :: cor !< coefficient pour le pas C
      integer, intent(out) :: ordre !<ordre de l'integrateur (nombre d'etapes)
      real(TREAL) :: aux0

      c = 0.d0
      d = 0.d0
      cor = 0
      
      ordre = numero
      
      select case (numero)
            
       case (1)
             c(2) = 1.d0
             d(1) = .5d0
             cor = -1.d0/24.d0
      
       case (2)
        d(2) = 2.d0/3.d0
        d(1) = 1.d0/6.d0
        c(2) = .5d0
        cor = 1.d0/72.d0
      
       case (3)
         c(3) = 1.d0/sqrt(5.d0)
         c(2) = .5d0*(1.d0 - c(3))
         d(2) = 5.d0/12.d0
         d(1) = 1.d0/12.d0
         cor = 0.0063182642795175399928962904734153431347d0
      
       case (4)
        d(3) = 16.d0/45.d0
         d(2) = 49.d0/180.d0
        d(1) = 1.d0/20.d0
         c(3) = sqrt(3.d0/7.d0)/2.d0
         c(2) = .5d0 - sqrt(3.d0/7.d0)/2.d0 
        cor = 0.0036447936001532493022971399654497729200d0
      
       case (5)
        aux0 = sqrt(7d0)
        d(3)  = (14.d0 + aux0)/60.d0
        d(2)  = (14.d0 - aux0)/60.d0
        d(1)  = 1.d0/30.d0
        c(4) = sqrt((1.d0 - 2.d0/aux0)/3.d0)
        c(3) = (sqrt(3.d0 + 6.d0/aux0) - sqrt(3.d0 - 6.d0/aux0))/6.d0
        c(2) = .5d0 - sqrt(3.d0 + 6.d0/aux0)/6.d0
        cor = 0.0023814866729536341874703862321814533098d0
      
       case (6)
        c(2)=.5d0-sqrt((15.d0+2.d0*sqrt(15.d0))/33.d0)/2.d0
        c(3)=sqrt(5.d0/22.d0-sqrt(5.d0/33.d0)/2.d0)
        c(4)= sqrt(5.d0/44.d0-sqrt(5.d0/3.d0)/22.d0)
        d(1)=1.d0/42.d0
        d(2)=31.d0/175.d0-sqrt(3.d0/5.d0)/20.d0
        d(3)=31.d0/175.d0+sqrt(3.d0/5.d0)/20.d0
        d(4)=128.d0/525.d0
        cor = 0.0016813465120919063265636932152964340419d0
      
       case (7)
        c(2) = 0.064129925745196692331277119389668280948109665161508d0
        c(3) = 0.14001998353823215659646751491135512407903984007983d0
        c(4) = 0.19120048176533171668792673552630096732507779221533d0 
        c(5) = 0.20929921790247886876865726034535125529554540508668d0
        d(1) = 1.d0/56.d0
        d(2) = 0.10535211357175301969149603288787816222767308308052d0
        d(3) = 0.17056134624175218238212033855387408588755548780279d0
        d(4) = 0.20622939732935194078352648570110489474191428625954d0
            cor = 0.0012517656160394000030725161002511908159d0        
      
       case (8)
        c(2) = 0.050121002294269921343827377790831020974259852216948d0
        c(3) = 0.1112858579503612019332299086634977536702253239761d0
        c(4) = 0.15703440784227979736756667919134161884519368500651d0
        c(5) = 0.18155873191308907935537603435432960651032113880044d0
        d(1) = 1.d0/72.d0
        d(2) = 0.08274768078040276252316986001460415291955361461258d0
        d(3) = 0.13726935625008086764035280928968636297062561049325d0
        d(4) = 0.17321425548652317255756576606985914397376635312546d0
        d(5) = 2048.d0/11025.d0
            cor = 0.0009687979680736885716546842084629822742d0
      
       case (9)
        c(2) =0.040233045916770593085533669588830932923228462276862d0
        c(3) =0.090380021530476869412913242981253704568678987034675d0
        c(4) =0.13042445764753028967096554106428636356470610402318d0
        c(5) =0.15632299607202873551747766338654223232781127489054d0
        c(6) =0.16527895766638702462621976595817353323115034354948d0
        d(1) = 1.d0/90.d0
        d(2) = 0.066652995425535055563113585377696449054812937498899d0
        d(3) = 0.11244467103156322605972891086552392137695645698088d0
        d(4) = 0.14602134183984187893779112868722194610374419484385d0
        d(5) = 0.16376988059194872832825526395844657235337529956526d0
            cor = 0.0007723490239999520782276868102603232890d0
      
       case (10)
        c(2) =0.032999284795970432833862931950308182730041334945019d0 
        c(3) =0.074758978372457357854928159995462765516359937507949d0
        c(4) =0.10962407333346970607572692331535321961207265476602d0
        c(5) =0.13473859570463280751952622695934707767102285560671d0
        c(6) =0.1478790677934696957159557577795287544705032171743d0
        d(1) = 1.d0/110.d0
        d(2) =0.054806136633497432230701724790175355024737071823638d0
        d(3) = 0.093584940890152602054070760949717459784597118812564d0
        d(4) = 0.12402405213201415702004243321093637667183658482386d0
        d(5) = 0.14343956238950404433961120166576761559182677372797d0
        d(6) = 32768.d0/218295.d0
            cor = 0.0006303200441631678407986387626651115288d0


        case default
         write(*,*)"l'integrateur BAB ",numero,
     &              " n'est pas au point "
         stop
        end select ! fin de l'ordre
      
      end subroutine create_schemaBAB_an


!***********************************************************************
!> creation d'un schema d'integrateur de type BAB en heliocentrique
!!                 'BAB' + ['C'] + nombre a 1 ou 2 chiffres
!!                 Remarque:
!!                 'BAB'  indique la sequence des pas.
!!                 'C' indique la presence ou non du correcteur
!!                 le nombre indique l'ordre de l'integrateur
!!  valeurs possibles de schema :
!!  'BAB4' ou 'BABC4'
!***********************************************************************
      subroutine create_schemaBABH(nom, sc)
      implicit none
      character(len=*), intent(in) :: nom !< nom du schema
      type(t_schema), intent(out) :: sc   !< schema cree
      integer :: l1, ordre,  numero
      integer, parameter :: coef_max = 100
      real(TREAL),dimension(coef_max) :: c,d
      real(TREAL) :: cor, coef_depart
      
      sc%m_if_cor = 0
      sc%m_ordre = 0
      
      c = 0.d0
      d = 0.d0
      cor = 0
      coef_depart = 0
      
      l1 = len_trim(nom)
      
      if (nom(5:5).eq.'C') then
       sc%m_if_cor = 1
       read(nom(6:l1),*) ordre
       sc%m_procdepart => t_schema_pasdepartC
       sc%m_procrun => t_schema_pasBAC
      else
       read(nom(5:l1),*) ordre
       sc%m_procdepart => t_schema_pasdepartvide
       sc%m_procrun => t_schema_pasBAB
      endif
      
      
      write(*,*)'integrateur ',nom
      write (*,*) "if_cor :", sc%m_if_cor
      write (*,*) "ordre :", ordre
      
      select case (ordre)

!***********************************************************************
! integrateur communaux helio et jacobi
!***********************************************************************
         case default
           numero = ordre
           call create_schemaBAB_an(numero, c,d, cor, ordre)

!***********************************************************************
! les integrateurs suivants sont specifiques aux heliocentriques
!***********************************************************************


       case (82)
***********************************************************************
*  Laskar SBAB_4 and McLahan (8,2) integrating schemes
* 	- 4 stages / symetric / NO corrector
*   int_type BAB82
************************************************************************

         d(1) = REALC(1.e0)/REALC(20.e0)
         d(2) = REALC(49.e0)/REALC(180.e0)
         d(3) = REALC(16.e0)/REALC(45.e0)
         d(4) = d(2) 
         d(5) = d(1)

         c(2) = REALC(.5e0) - sqrt(REALC(3.e0)/REALC(7.e0))/REALC(2.e0)
         c(3) = sqrt(REALC(3.e0)/REALC(7.e0))/REALC(2.e0)
         c(4) = c(3)
         c(5) = c(2)
         
         ordre = 4

       case (84)
************************************************************************
*  McLahan (8,4) integrating schemes
*	- 5 stages / symetric / NO corrector
*   int_type BAB84
************************************************************************
       d(1) = REALC(0.81186273854451628884e0); 
       d(2) = REALC(-0.67748039953216912289e0);
       d(3) = REALC(0.36561766098765283405e0); 
       d(4) = d(3)
       d(5) = d(2)
       d(6) = d(1)


       c(2) = REALC(-0.00758691311877447385e0); 
       c(3) = REALC(0.31721827797316981388e0); 
       c(4) = REALC(0.38073727029120931994e0); 
       c(5) = c(3)
       c(6) = c(2)

         ordre = 5

       case (844)
************************************************************************
*  Blanes (8,4,4) integrating schemes
*  - 6 stages / symetric / extra stage for killig  tau^3
*   int_type BABH844
************************************************************************

       d(1) = REALC(0.1308424104615589109e0);  
       d(2) = REALC(-0.010864481464054482540e0);
       d(3) = REALC(1.0281780095953900777e0);  
       d(4) = REALC(-1.2963118771857890123e0);
       d(5) = d(3);  
       d(6) = d(2);  
       d(7) = d(1);  

       c(2) = REALC(-0.1639587030679243705e0); 
       c(3) = REALC(0.7795825181082894712e0);
       c(4) = REALC(-0.11562381504036510071e0);
       c(5) = c(4);
       c(6) = c(3);
       c(7) = c(2);

         ordre = 6

       case (864)
************************************************************************
*  Blanes (8,6,4) integrating schemes
*  - 8 stages / symetric / extra stage for killig  tau^3
*   int_type ABAH864
************************************************************************

       c(1) = REALC(0e0);
       d(1) = REALC(0.4909632009912256537459526620221413892904e0);
       c(2) = REALC(-0.002723605386947930589739382599575762161318e0);
       d(2) = REALC(-0.4210537416491733768903025006178189256072e0);
       c(3) = REALC(0.1891751066286274399263963781088028608800e0);
       d(3) = REALC(0.2543676819230724635619821352391375985529e0);
       c(4) = REALC(0.3376654790456454914631693687072100451686e0);
       d(4) = REALC(0.5702018058693710290287561986575697238903e0);
       c(5) = REALC(-0.02411698028732500079982636421643714388728e0);
       d(5) = REALC(-0.7889578942689915388927769906020595722527e0);
       c(6) = REALC(-0.02411698028732500079982636421643714388728e0);
       d(6) = REALC(0.5702018058693710290287561986575697238903e0);
       c(7) = REALC(0.3376654790456454914631693687072100451686e0);
       d(7) = REALC(0.2543676819230724635619821352391375985529e0);
       c(8) = REALC(0.1891751066286274399263963781088028608800e0);
       d(8) = REALC(-0.4210537416491733768903025006178189256072e0);
       c(9) = REALC(-0.002723605386947930589739382599575762161318e0);
       d(9) = REALC(0.4909632009912256537459526620221413892904e0); 


         ordre = 8

         end select ! fin de l'ordre
         
         sc%m_ordre = ordre
         
         write (*,*) "nombre etapes :", sc%m_ordre
        
      
         allocate(sc%db(ordre+1))
         allocate(sc%cb(ordre+1))
         sc%cb = c(1:ordre+1)
         sc%db = d(1:ordre+1)
      
         sc%m_cor = cor
         if (sc%m_if_cor.eq.1) then 
            coef_depart = -cor/REALC(2.E0) ! pour utiliser coef_depart systematiquement
         else
            coef_depart = REALC(0.E0)
         endif
         sc%m_coef_depart = coef_depart
      
       
      end subroutine create_schemaBABH
      
!***********************************************************************
!> creation d'un schema d'integrateur de type BAB en jacobi
!!                 'BAB' + ['C'] + nombre a 1 ou 2 chiffres
!!                 Remarque:
!!                 'BAB'  indique la sequence des pas.
!!                 'C' indique la presence ou non du correcteur
!!                 le nombre indique l'ordre de l'integrateur
!!  valeurs possibles de schema :
!!  'BAB4' ou 'BABC4'
!***********************************************************************
      subroutine create_schemaBABJ(nom, sc)
      implicit none
      character(len=*), intent(in) :: nom !< nom du schema
      type(t_schema), intent(out) :: sc   !< schema cree
      integer :: l1, ordre, j,i,n, numero
      integer, parameter :: coef_max = 100
      real(TREAL),dimension(coef_max) :: c,d, cb,db
      real(TREAL) :: cor, coef_depart
      
      sc%m_if_cor = 0
      sc%m_ordre = 0
      
      c = 0.d0
      d = 0.d0
      cor = 0
      coef_depart = 0
      
      l1 = len_trim(nom)
      
      if (nom(4:4).eq.'C') then
       sc%m_if_cor = 1
       read(nom(5:l1),*) ordre
       sc%m_procdepart => t_schema_pasdepartC
       sc%m_procrun => t_schema_pasBAC
      else
       read(nom(4:l1),*) ordre
       sc%m_procdepart => t_schema_pasdepartB
       sc%m_procrun => t_schema_pasAB
      endif
      
      
      write(*,*)'integrateur ',nom
      write (*,*) "if_cor :", sc%m_if_cor
      write (*,*) "ordre :", ordre
      
      select case (ordre)

!***********************************************************************
! integrateur communaux helio et jacobi
!***********************************************************************
         case default
           numero = ordre
           call create_schemaBAB_an(numero, c,d, cor, ordre)

!***********************************************************************
! les integrateurs suivants sont specifiques aux jacobis
!***********************************************************************
       case (82)
************************************************************************
*  Laskar SBAB_4 and McLahan (8,2) integrating schemes
* 	- 4 stages / symetric / NO corrector
*   int_type BAB82
************************************************************************

         d(1) = REALC(1.e0)/REALC(20.e0)
         d(2) = REALC(49.e0)/REALC(180.e0)
         d(3) = REALC(16.e0)/REALC(45.e0)
         d(4) = d(2) 
         d(5) = REALC(2e0)*d(1)

         c(2) = REALC(.5e0) - sqrt(REALC(3.e0)/REALC(7.e0))/REALC(2.e0)
         c(3) = sqrt(REALC(3.e0)/REALC(7.e0))/REALC(2.e0)
         c(4) = c(3)
         c(5) = c(2)
         
         ordre = 4

       case (84)
************************************************************************
*  McLahan (8,4) integrating schemes
*	- 5 stages / symetric / NO corrector
*   int_type BAB84
************************************************************************
         d(1) = REALC(0.81186273854451628884e0)
         d(2) = REALC(-0.67748039953216912289e0)
         d(3) = REALC(0.36561766098765283405e0)
         d(4) = d(3)
         d(5) = d(2)
         d(6) = REALC(2.e0)*d(1)

         c(2) = REALC(-0.00758691311877447385e0)
         c(3) = REALC(0.31721827797316981388e0)
         c(4) = REALC(0.38073727029120931994e0)
         c(5) = c(3)
         c(6) = c(2)

         ordre = 5
         
       case (864)
************************************************************************
*  Blanes (8,6,4) integrating schemes eo (10, 6, 4)
* - 9 stages / symetric / NO corrector
*   int_type BAB864
************************************************************************
         d(1) = REALC(0.0492150656634390294003903097390e0)
         d(2) = REALC(0.404400755945652352844233877907e0)
         d(3) = REALC(0.294714729388018745931498720662e0)
         d(4) = REALC(0.229051044384905913483491183170e0)
         d(5) = REALC(-0.477381595382016041659614091478e0)        
         !d(5) = REALC(0.5e0) - (d(1) + d(2) + d(3) + d(4))         
         d(6) = d(5)
         d(7) = d(4)
         d(8) = d(3)
         d(9) = d(2)
         d(10) = REALC(2e0)*d(1)
       
         c(2) = REALC(0.0492150656634390294003903097390e0)
         c(3) = REALC(0.375906575650518113208289535856e0)
         c(4) = REALC(0.220249195979751384369990554054e0)
         c(5) = REALC(-0.0892572733579359490419485657182e0)
         c(6) = REALC(-0.482448929390071663627890699658e0) 
         !c(6) = REALC(1e0) - REALC(2e0)*(c(1) + c(2) + c(3) + c(4))
         c(7) = c(5)
         c(8) = c(4)
         c(9) = c(3)
         c(10) = c(2)

         ordre = 9

         end select ! fin de l'ordre
         
         sc%m_ordre = ordre
         
         write (*,*) "nombre etapes :", sc%m_ordre
        
      
!--------- Initalisation des coeff cb et db      
         cb = 0.d0
         db = 0.d0

       if (mod(ordre,2).eq.0) then 

         n = ordre/2
         if (sc%m_if_cor.eq.0) then 
            do j=1,n 
               cb(j) = c(j+1)
               cb(n+j) = c(n-j+2)
               db(j) = d(j+1)
               db(n+j) = d(n-j+1)
            end do
            db(2*n) = 2*d(1)
            coef_depart = d(1)
            write(*,*)'          i     cb                      db'
            do i=1,ordre
               write(*,*)i,cb(i),db(i)
            end do
            write(*,*)'coef_depart',  coef_depart
            write(*,*)sum(cb),sum(db) 
         else 
            do j=1,n 
               cb(j) = c(j+1)
               cb(n+j) = c(n-j+2)
               db(j) = d(j)
               db(n+j+1) = d(n-j+1)
            end do
            db(n+1) = d(n+1)
            coef_depart = -cor/REALC(2.E0) ! pour utiliser coef_depart systematiquement
            write(*,*)'          i     cb                      db'
            do i=1,ordre
               write(*,*)i,cb(i),db(i)
            end do
            write(*,*)ordre+1,0.d0,db(ordre+1)
            write(*,*)'cor',  cor
            write(*,*)sum(cb),sum(db)
         end if  
         
       else if (mod(ordre,2).eq.1) then
         n = (ordre-1)/2
         if (sc%m_if_cor.eq.0) then 
            do j=1,n
              cb(j) = c(j+1)
              cb(n+j+1) = c(n-j+2)
              db(j) = d(j+1)
              db(n+j+1) = d(n-j+1)
            end do
            db(n+1) = d(n+1)
            cb(n+1) = c(n+2)
            db(ordre) = 2*d(1)
            coef_depart = d(1)
            write(*,*)'          i     cb                      db'
              do i=1,ordre
              write(*,*)i,cb(i),db(i)
            end do
            write(*,*)'coef_depart',  coef_depart
            write(*,*)sum(cb),sum(db)
         else 
            do j=1,n
              cb(j) = c(j+1)
              cb(n+j+1) = c(n-j+2)
            end do
            cb(n+1) = c(n+2)
            do j=1,n+1
              db(j) = d(j)
              db(n+j+1) = d(n-j+2)
            end do
           coef_depart = -cor/REALC(2.E0) ! pour utiliser coef_depart systematiquement

            write(*,*)'          i     cb                      db'
            do i=1,ordre
              write(*,*)i,cb(i),db(i)
            end do
            write(*,*)ordre+1,0.d0,db(ordre+1)
            write(*,*)sum(cb),sum(db)
            write(*,*)'cor',  cor
            write(*,*)sum(cb),sum(db)  
         end if
       endif

         allocate(sc%cb(ordre+1))
         allocate(sc%db(ordre))
         sc%cb = cb
         sc%db = db
      
      sc%m_cor = cor
      sc%m_coef_depart = coef_depart
            
      end subroutine create_schemaBABJ

!***********************************************************************
!> creation d'un schema d'integrateur de type ABA
!!                 'ABA' + ['C'] + nombre a 1 ou 2 chiffres
!!                 Remarque:
!!                 'ABA'  indique la sequence des pas.
!!                 'C' indique la presence ou non du correcteur
!!                 le nombre indique l'ordre de l'integrateur
!!  valeurs possibles de schema :
!!  'ABA4' ou 'ABAC4'
!***********************************************************************
      subroutine create_schemaABAJ(nom, sc)
      implicit none
      character(len=*), intent(in) :: nom !< nom du schema
      type(t_schema), intent(out) :: sc   !< schema cree
      integer :: l1, ordre, numero
      real(TREAL),dimension(coef_max) :: c,d
      real(TREAL) :: cor, aux0, aux1
      
      sc%m_if_cor = 0
      sc%m_ordre = 0
      
      c = 0.d0
      d = 0.d0
      cor = 0
      
      l1 = len_trim(nom)
      
      if (nom(4:4).eq.'C') then
       sc%m_if_cor = 1
       read(nom(5:l1),*) ordre
       sc%m_procdepart => t_schema_pasdepartC
       sc%m_procrun => t_schema_pasABC
      else
       read(nom(4:l1),*) ordre
       sc%m_procdepart => t_schema_pasdepartA
       sc%m_procrun => t_schema_pasBA
      endif

      write(*,*)  'integrateur ',nom
      write (*,*) 'if_cor :', sc%m_if_cor
      write (*,*) 'ordre :', ordre
      
      
      select case (ordre)
      
!***********************************************************************
! integrateur communaux helio et jacobi
!***********************************************************************
         case default
           numero = ordre
           call create_schemaABA_an(numero, c,d, cor, ordre)

!***********************************************************************
! les integrateurs suivants sont specifiques aux jacobis
!***********************************************************************


       case (82)
******************************************************************
*     Laskar SABA_4 and McLahan (8,2) integrating schemes
* 		- 4 stages / symetric / NO corrector
*      int_type ABA82
******************************************************************

         aux0 = sqrt(REALC(525e0) + REALC(70e0)*sqrt(REALC(30e0)))
         aux1 = sqrt(REALC(525e0) - REALC(70e0)*sqrt(REALC(30e0)))
   
         d(1) = REALC(1e0)/REALC(4e0) - sqrt(REALC(30e0))/REALC(72e0)
         d(2) = REALC(1e0)/REALC(4e0) + sqrt(REALC(30e0))/REALC(72e0)
         d(3) = d(2)
         d(4) = d(1)
   
         c(1) = REALC(.5e0) - aux0/REALC(70e0)     
         c(2) = aux0/REALC(70e0) - aux1/REALC(70e0)
         c(3) = aux1/REALC(35e0)
         c(4) = c(2)
         c(5) = REALC(2.e0)*c(1)
   
         ordre = 4 

       case (864)
************************************************************************
*  Blanes (8,6,4) integrating schemes eo(10,8,6)
*  - 9 stages / symetric / NO corrector
*   int_type ABA864
************************************************************************
 
         c(1) = REALC(0.0453712130326967590033317595992e0)  
         c(2) = REALC(0.266355488928810576888813408539e0)   
         c(3) = REALC(0.470996475404286449136399805662e0)   
         c(4) = REALC(-0.0426935662057334025124411638611e0) 
         c(5) = REALC(0.5e0) - (c(1)+c(2)+c(3)+c(4)) 
         c(6) = c(5)
         c(7) = c(4)
         c(8) = c(3)
         c(9) = c(2)
         c(10) = REALC(2e0)*c(1)

         d(1) = REALC(0.110697092141418037890462980579e0)
         d(2) = REALC(0.456621746800863158212216102725e0)
         d(3) = REALC(0.447019291364693625053008799194e0)
         d(4) = REALC(-0.575034109315983726141544334602e0)
         d(5) = REALC(1e0) - REALC(2e0)*(d(1)+d(2)+d(3)+d(4)) 
         d(6) = d(4)
         d(7) = d(3)
         d(8) = d(2)
         d(9) = d(1)
        
         ordre = 9

       case (1064)
************************************************************************
*  Blanes (10,6,4) integrating schemes 
*  - 8 stages / symetric / NO corrector
*   int_type ABA864
************************************************************************
         c(1) = REALC(0.03809449742241219545697532230863756534060)
         c(2) = REALC(0.1452987161169137492940200726606637497442)
         c(3) = REALC(0.2076276957255412507162056113249882065158)
         c(4) = REALC(0.4359097036515261592231548624010651844006)
         c(5) = REALC(-0.6538612258327867093807117373907094120024)
         c(6) = c(4) 
         c(7) = c(3) 
         c(8) = c(2) 
         c(9) = REALC(2e0)*c(1) 

         d(1) = REALC(0.09585888083707521061077150377145884776921)
         d(2) = REALC(0.2044461531429987806805077839164344779763)
         d(3) = REALC(0.2170703479789911017143385924306336714532)
         d(4) = REALC(-0.01737538195906509300561788011852699719871)
         d(5) = d(4)
         d(6) = d(3)
         d(7) = d(2)
         d(8) = d(1)

         ordre = 8

       case (104)
************************************************************************
*  Blanes (10,4) integrating schemes 
*  - 7 stages / symetric / NO corrector
*   int_type ABA104 (com el McLahan pero amb una etapa mes)
************************************************************************
         c(1) = REALC(0.04706710064597250612947887637243678556564e0)
         c(2) = REALC(0.1847569354170881069247376193702560968574e0)
         c(3) = REALC(0.2827060056798362053243616565541452479160e0)
         c(4) = REALC(-0.01453004174289681837857815229683813033908e0)
         c(5) = c(4) 
         c(6) = c(3) 
         c(7) = c(2) 
         c(8) = REALC(2e0)*c(1) 

         d(1) = REALC(0.1188819173681970199453503950853885936957e0)
         d(2) = REALC(0.2410504605515015657441667865901651105675e0)
         d(3) = REALC(-0.2732866667053238060543113981664559460630e0)
         d(4) = REALC(0.8267085775712504407295884329818044835997e0)
         d(5) = d(3)
         d(6) = d(2)
         d(7) = d(1)

         ordre = 7

        end select ! fin de l'ordre
      
      
       ! remplissage final des coefs et ordre  de sc
       call create_schemaABA_fill(sc, c,d, cor, ordre)
      
      end subroutine create_schemaABAJ

      end module mod_schema

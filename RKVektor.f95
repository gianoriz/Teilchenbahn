module Saturndata
  implicit none
  save

  ! Konstanten *****************************
  double precision  :: gamma = 6.67428e-11                            ![m³/kg*s²]                 (Wiki)
  double precision  :: d =  527040000.0                               ![m]   Entf. Rhea zu Saturn (Wiki) 
  double precision  :: MSaturn = 5.685e26                             ![kg]  Masse Saturn         (Wiki)
  double precision  :: MRhea = 0.000023166e26                         ![kg]  Masse Rhea 2.3166e21 (berechnet) 
  double precision  :: L1                                             !   Lagrangepunkt        (berechnet)
  double precision  :: RRhea = 764000.0                               ![m]   Rhearadius           (Wiki) 
  double precision  :: RSaturn = 57316000.0                           ![m]   Saturnradius         (Wiki)
  double precision  :: omegaz                                         ![1/s] Winkelgeschw.        (berechnet) 

  ! Globale Variablen **********************
  double precision  :: dt                          ! Zeitintervall
  double precision , dimension(3) :: VecX          ! Ortsvektor   
  double precision , dimension(3) :: VecV          ! Geschwindigkeitsvektor   

end module Saturndata

module vectors
  implicit none

  public :: vector, cross, scalar

  type vector
     double precision :: x,y,z 
     character(len=20) :: SIunit                        ! Einheit des Vektors in SI Einheiten  (z.B. m/s)
  end type vector

contains
  function cross(a, b)
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: cross

    cross%x = a%y * b%z - a%z * b%y
    cross%y = a%z * b%x - a%x * b%z
    cross%z = a%x * b%y - a%y * b%x
  end function cross

  function scalar(a, b)
    TYPE (vector), INTENT (in) :: a, b
    double precision :: scalar
    scalar = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)
  end function scalar

  function vabs(a)
    TYPE (vector), INTENT (in) :: a
    double precision :: vabs
    vabs = sqrt(a%x **2 + a%y **2 + a%z **2)
  end function vabs

  function vadd(a,b)
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: vadd
    vadd%x = a%x + b%x
    vadd%y = a%y + b%y
    vadd%z = a%z + b%z
  end function vadd

  function vsub(a,b)
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: vsub
    vsub%x = a%x - b%x
    vsub%y = a%y - b%y
    vsub%z = a%z - b%z
  end function vsub

  function vsmul(a,s)
    TYPE (vector), INTENT (in) :: a
    double precision :: s
    TYPE (vector) :: vsmul
    vsmul%x = a%x * s
    vsmul%y = a%y * s
    vsmul%z = a%z * s
  end function vsmul


end module vectors




! **********************************************************************
! P R O G R A M
! **********************************************************************

program TeilchenbahnNeu
  use Saturndata
  use vectors
  implicit none   

  ! Funktionen *****************************
  double precision :: Betrag3
  double precision :: Wurzel
  double precision :: Fc
  double precision :: Fz
  double precision :: Fg

  ! Variablen ******************************
  TYPE (vector) :: Nullvektor

  ! Berechnete Variablen *******************
  L1 = d - (d/( sqrt(MRhea/MSaturn)+1.0)) 
  omegaz = sqrt((gamma * MSaturn)/(d**3))

  dt = 42
  call Werteausgeben



end program TeilchenbahnNeu

! **********************************************************************
! P R O G R A M     E N D
! **********************************************************************

double precision  function Betrag3(x,y,z)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: x,y,z
  write(*,*) x,y,z
  Betrag3 = sqrt(x**2 + y**2 + z**2)
  return
! TESTFUNKTION:
! write(*,*) Betrag3(dble(1), dble(1), dble(0))
end function Betrag3



!---------------------------------------
double precision  function Fg(x,y,z)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: x,y,z


  return
end function Fg
!---------------------------------------



!---------------------------------------
double precision  function Fc(x,y,z)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: x,y,z


  return
end function Fc
!---------------------------------------



!---------------------------------------
double precision  function Fz(x,y,z)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: x,y,z


  return
end function Fz
!---------------------------------------



subroutine Werteausgeben
  use Saturndata
  use vectors
  implicit none
  write(*,*) "Hallo Du Eimer", dt
  return 
end subroutine Werteausgeben



! Wissen-Konserve:
! COMMON http://www.obliquity.com/computer/fortran/common.html
! MODULE https://srv.rz.uni-bayreuth.de/lehre/fortran90/vorlesung/V08/V08.html

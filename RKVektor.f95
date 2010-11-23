
module Saturndata
  implicit none
  save




  ! Konstanten *****************************
  double precision  :: gamma = 6.67428e-11         ![m³/kg*s²]                      (Wiki)
  double precision  :: d =  527040000.0            ![m]        Entf. Rhea zu Saturn (Wiki) 
  double precision  :: MSaturn = 5.685e26          ![kg]       Masse Saturn         (Wiki)
  double precision  :: MRhea = 0.000023166e26      ![kg]       Masse Rhea 2.3166e21 (berechnet) 
  double precision  :: MStein = 1.d0                  ![kg]       Masse Rhea 2.3166e21 (berechnet) 
  double precision  :: L1                          ![m]        Lagrangepunkt        (berechnet)
  double precision  :: RRhea = 764000.0            ![m]        Rhearadius           (Wiki) 
  double precision  :: RSaturn = 57316000.0        ![m]        Saturnradius         (Wiki)
  double precision  :: omegaz                      ![1/s]      Winkelgeschw.        (berechnet) 
  ! Globale Variablen **********************
  double precision  :: dt                          !Zeitintervall
  double precision , dimension(3) :: VecX          !Ortsvektor   
  double precision , dimension(3) :: VecV          !Geschwindigkeitsvektor   
end module Saturndata






! Vektorrechnungen ************************
module vectors
  implicit none

  public :: vector, cross, scalar, vsmul, svmul

  type vector
     double precision :: x,y,z                     !Komponenten des Vektors sind reelle Zahlen
     character(len=20) :: SIunit                   !Einheit des Vektors in SI Einheiten  (z.B. m/s)
  end type vector

  interface operator (-)
     module procedure vsub
  end interface

  interface operator (+)
     module procedure vadd
  end interface

  interface operator (*)
     module procedure scalar
  end interface

  interface operator (.kreuz.)
     module procedure cross
  end interface

contains
  function cross(a, b)                             !Kreuzprodukt von Vektoren
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: cross

    cross%x = a%y * b%z - a%z * b%y
    cross%y = a%z * b%x - a%x * b%z
    cross%z = a%x * b%y - a%y * b%x
  end function cross

  function scalar(a, b)                            !Skalarprodukt von Vektoren
    TYPE (vector), INTENT (in) :: a, b
    double precision :: scalar
    scalar = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)
    ! FIXME    vadd%SIunit = a%SIunit & "*" & b%SIunit  !! DIESE ZEILE MUSS ZEICHENKETTE ANEINANDERHAENGEN!
  end function scalar

  function vabs(a)                                 !Betrag eines Vektors
    TYPE (vector), INTENT (in) :: a
    double precision :: vabs
    vabs = sqrt(a%x **2 + a%y **2 + a%z **2)
  end function vabs

  function vadd(a,b)                               !Addition von Vektoren
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: vadd
    !   if (a%SIunit == b%SIunit) then
    vadd%x = a%x + b%x
    vadd%y = a%y + b%y
    vadd%z = a%z + b%z
    vadd%SIunit = a%SIunit
    !   else 
    !      write(*,*) "EE: Vektoren inkompatibel!"
    !      stop
    !   end if
  end function vadd

  function vsub(a,b)                               !Subtraktion von Vektoren
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: vsub
    !    if (a%SIunit == b%SIunit) then
    vsub%x = a%x - b%x
    vsub%y = a%y - b%y
    vsub%z = a%z - b%z
    vsub%SIunit = a%SIunit
    !   else 
    !       write(*,*) "EE: Vektoren inkompatibel!"
    !       stop
    !    end if
  end function vsub


  function vsmul(a,s)                              !Skalarmultiplikation von rechts
    TYPE (vector), INTENT (in) :: a
    double precision :: s
    TYPE (vector) :: vsmul
    vsmul%x = a%x * s
    vsmul%y = a%y * s
    vsmul%z = a%z * s
    vsmul%SIunit = a%SIunit
  end function vsmul



 function svmul(s,a)                              !Skalarmultiplikation von links
    TYPE (vector), INTENT (in) :: a
    double precision :: s
    TYPE (vector) :: svmul
    svmul%x = s * a%x 
    svmul%y = s * a%y 
    svmul%z = s * a%z 
    svmul%SIunit = a%SIunit
  end function svmul






  ! Benoetigen hier eigentlich noch Constructor und Methoden zum Datenaustausch
  ! Das MUSS nachgebessert werden wir haben jetzt aber Hunger!




  function newvec(a, b, c, SIunit)
    TYPE (vector) :: newvec
    double precision, intent(IN) :: a, b, c
    character(len=20) :: SIunit 
    newvec%x = a
    newvec%y = b
    newvec%z = c
    newvec%SIunit = SIunit     
  end function newvec


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
  type (vector) :: Fc
  type (vector) :: Fz
  type (vector) :: Fg 
  type (vector) :: Fges


  ! Variablen ******************************
  TYPE (vector) :: Nullvektor
  TYPE (vector) :: vektora 
  TYPE (vector) :: vektorb
  TYPE (vector) :: vektorc
  TYPE (vector) :: vektord
  double precision :: dzahl 


  ! Berechnete Variablen *******************
  L1 = d - (d/( sqrt(MRhea/MSaturn)+1.0)) 
  omegaz = sqrt((gamma * MSaturn)/(d**3))


  vektora = newvec(4.d0, 2.d0, 8.d0, "m/s                 ")
  vektorb = newvec(1.d0, 3.d0, 9.d0, "m/s                 ")

  !Test:
  !vektorc = cross(vektora, vektorb) !OK
  !dzahl = scalar(vektora, vektorb) !OK
  !vektord = vsmul(vektora, 3.d0) !OK
  !vektord = svmul(3.d0, vektora) !OK


  !Fges = Fg + Fc + Fz


  write(*,*) Fges

  !write(*,*) MStein, vektora, vektorb
  !write(*,*) cross(vektora,  vektorb)



  !write(*,*) Fc(MStein, vektorb, vektora)
  !write(*,*) Fz(MStein, vektorb, vektora)



end program TeilchenbahnNeu

! **********************************************************************
! P R O G R A M     E N D
! **********************************************************************

double precision  function Betrag(x,y,z)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: x,y,z
  write(*,*) x,y,z
  Betrag = sqrt(x**2 + y**2 + z**2)
  return
  ! TESTFUNKTION:
  ! write(*,*) Betrag(dble(1), dble(1), dble(0))
end function Betrag



type(vector)  function Fg(vecr) !Wir definieren hier die Gravitationskraft als Summe der Gravikraefte von Saturn und Rhea
  use Saturndata
  use vectors
  implicit none
  Type (ec), intent(in) :: vecr!Hier weiter machen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Fg = -gamma * MSaturn * MTeilchen * r/r**3    !-gamma * MRhea * MTeilchen
  return
end function Fg


!-------------------------------------------------
type(vector)  function Fc(m, w, v)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: m
  type(vector) :: w, v

  write(*,*) w .kreuz. v

  Fc = vsmul((w .kreuz. v), (2 * m))
  return
end function Fc


type(vector)  function Fz(m, w, r)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: m
  type(vector) :: w, r

  Fz = vsmul((w .kreuz. (w .kreuz. r)), m)

  return
end function Fz



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



! Wie es weiter geht:
! DONE Constructor und Interface fuer Vektor erzeugen
! >>>>> Datentyp Vektor testen Beispielaufgaben fuer jede Rechnung!
! >>>>>>>>> Alle Befehle im programm verstehen
! DONE Alle Funktionen zur Vektorrechnung von double precision auf vector umstellen und mit Rechnung fuellen
! 






!Protokoll:
!Vektor anlegen funktioniert

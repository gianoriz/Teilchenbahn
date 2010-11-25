



!KONSTANTEN######################################################################
module Saturndata
  implicit none
  save
  double precision :: gamma = 6.67428e-11    ![m³/kg*s²]                (Wiki)
  double precision :: MSaturn = 5.685e26     ![kg] Masse Saturn         (Wiki)
  double precision :: MRhea = 0.000023166e26 ![kg] Masse Rhea 2.3166e21 (berechnet)
  double precision :: MTeilchen = 1.d0       ![kg] Masse Rhea 2.3166e21 (berechnet)
  double precision :: L1                     ![m] Lagrangepunkt         (berechnet)
  double precision :: RRhea = 764000.0       ![m] Rhearadius            (Wiki)
  double precision :: RSaturn = 57316000.0   ![m] Saturnradius          (Wiki)

  !Globale Variablen **********************
!  double precision :: dt !Zeitintervall
!  double precision , dimension(3) :: VecX    !Ortsvektor
!  double precision , dimension(3) :: VecV    !Geschwindigkeitsvektor
 
end module Saturndata
!##################################################################################





!VEKTORRECHNUNGEN##################################################################
module vectors
  implicit none

  public :: vector, cross, scalar, vsmul, svmul

  type vector
     double precision :: x,y,z   !Komponenten des Vektors sind reelle Zahlen
     character(len=5) :: SIunit  !Einheit des Vektors in SI Einheiten (z.B. m/s)
  end type vector


contains

  function vadd(a,b)                                   !Addition von Vektoren
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: vadd
    vadd%x = a%x + b%x
    vadd%y = a%y + b%y
    vadd%z = a%z + b%z
    vadd%SIunit = a%SIunit
  end function vadd


  function vsub(a,b)                                  !Subtraktion von Vektoren
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: vsub
    vsub%x = a%x - b%x
    vsub%y = a%y - b%y
    vsub%z = a%z - b%z
    vsub%SIunit = a%SIunit
  end function vsub


  function vsmul(a,s)                                 !Skalarmultiplikation von rechts
    TYPE (vector), INTENT (in) :: a   
    double precision, intent(in) ::s
    TYPE (vector) :: vsmul
    vsmul%x = a%x * s
    vsmul%y = a%y * s
    vsmul%z = a%z * s
    vsmul%SIunit = a%SIunit
  end function vsmul


  function svmul(s,a)                                 !Skalarmultiplikation von links
    TYPE (vector), INTENT (in) :: a
    double precision, intent(in) ::s
    TYPE (vector) :: svmul
    svmul%x = s * a%x
    svmul%y = s * a%y
    svmul%z = s * a%z
    svmul%SIunit = a%SIunit
  end function svmul


 function scalar(a, b)                                !Skalarprodukt von Vektoren
    TYPE (vector), INTENT (in) :: a, b
    double precision :: scalar
    scalar = (a%x * b%x) + (a%y * b%y) + (a%z * b%z)
  end function scalar


  function cross(a, b)                                !Kreuzprodukt von Vektoren
    TYPE (vector), INTENT (in) :: a, b
    TYPE (vector) :: cross
    cross%x = a%y * b%z - a%z * b%y
    cross%y = a%z * b%x - a%x * b%z
    cross%z = a%x * b%y - a%y * b%x
  end function cross


  function vabs(a)                                    !Betrag eines Vektors
    TYPE (vector), INTENT (in) :: a
    double precision :: vabs
    vabs = sqrt(a%x **2.d0 + a%y **2.d0 + a%z **2.d0)
  end function vabs


  function newvec(a, b, c, SIunit)                    !Erstellt neue Vektoren mit SI-Einheit
    TYPE (vector) :: newvec
    double precision, intent(IN) :: a, b, c
    character(len=5) :: SIunit
    newvec%x = a
    newvec%y = b
    newvec%z = c
    newvec%SIunit = SIunit
  end function newvec

end module vectors
!##################################################################################









!##################################################################################
! P R O G R A M
!##################################################################################



program TeilchenbahnNeu

  use Saturndata
  use vectors
  implicit none

  ! Vektor Funktionen *****************************
  double precision :: Betrag

  TYPE (vector) :: Fgges,   A        !Ergebnis der Berechnung der Gesamtgravikraefte wird in der Variable A gespeichert  
  TYPE (vector) :: Fgsat,   B        !Ergebnis der Saturn-Gravikraft wird in der Variablen B gespeichert 
  TYPE (vector) :: Fgrhe,   C        !Ergebnis der Rhea-Gravikraft wird in der Variablen C gespeichert 
  TYPE (vector) :: Fc,      D        !Ergebnis der Cori.-Kraft wird in der Variablen D gespeichert 
  TYPE (vector) :: Fz,      E        !Ergebnis der Zentrip.-Kraft wird in der Variablen E gespeichert 
  TYPE (vector) :: Fgesamt, F        !Gesamtkraft, die auf das Teilchen wirkt wird in der Variablen F gespeichert 

  TYPE (vector) :: K1r               !Runge-Kutta-Konstante fuer den Ort 
  TYPE (vector) :: K2r               !Runge-Kutta-Konstante fuer den Ort 
  TYPE (vector) :: K3r               !Runge-Kutta-Konstante fuer den Ort 
  TYPE (vector) :: K4r               !Runge-Kutta-Konstante fuer den Ort 

  TYPE (vector) :: K1v               !Runge-Kutta-Konstante fuer die Geschw. 
  TYPE (vector) :: K2v               !Runge-Kutta-Konstante fuer die Geschw. 
  TYPE (vector) :: K3v               !Runge-Kutta-Konstante fuer die Geschw. 
  TYPE (vector) :: K4v               !Runge-Kutta-Konstante fuer die Geschw. 

  ! Vektor Variablen ******************************
  TYPE (vector) :: vecr              !Positionsvektor des Teilchens bezgl. Rhea 
  TYPE (vector) :: vecs              !Positionsvektor des Teilchens bezgl. Saturn
  TYPE (vector) :: vecv              !Geschwindigkeitsvektor 
  TYPE (vector) :: veca              !Beschleunigung, die auf das Teilchen wirkt
  TYPE (vector) :: AbstRheaSaturnVec !Abstandsvektor zwischen Saturn und Rhea
  TYPE (vector) :: omega             ![1/s] Winkelgeschw. Vektor


  !Parameter**************************************
  double precision :: Abstand        !Abstand Saturn-Rhea
  double precision :: wz             ![1/s] z-Komponente der Winkelgeschw. 
  double precision :: dzahl          !Testzahl 

  double precision t, dt             !Flugzeit & Zeitintervall
  integer(kind=8) :: i               !Laufvariable
  integer(kind=8) :: n               !Endwert





  !Anfangswerte des Teilchens:
  AbstRheaSaturnVec = newvec(0.d0, 527040000.d0, 0.d0, "m    ")!Verbindungsvektor zwischen Saturn und Rhea
  vecr = newvec(4.d0, 500000.d0, 8.d0, "m    ")                !Positionsvektor des Teilchens bezgl. Rhea   
  vecs = vsub(vecr, AbstRheaSaturnVec)                         !Positionvektor des Teilchens im Ruhesystem von Saturn 
  vecv = newvec(1.d0, 1.d0, 1.d0, "m/s  ")                     !Geschwindigkeitsvektor des Teilchen


  !Berechne den Vektor der Winkelgeschwindigkeit: OK
  Abstand = vabs(AbstRheaSaturnVec)
  wz = sqrt((gamma * MSaturn)/Abstand**3)
  omega = newvec(0.d0, 0.d0, wz, "1/s  ") !Vektor Omega
  write(*,*) "omega   =", omega


  !Berechne die Gesamt-Gravikraft: OK
  B = Fgsat(vecr, AbstRheaSaturnVec) 
  C = Fgrhe(vecr)
  A = vadd(B, C)
  write(*,*) "Fgges   =", A !Gesamtgravikraft


  !Berechne die Corioliskraft: OK
  D = Fc(MTeilchen, omega, vecv)
  write(*,*) "Fc      =", D 


  !Berechne die Zentripetalkraft: 
  E = Fz(MTeilchen, omega, vecs)
  write(*,*) "Fz      =", E 


  !Berechne die Gesamtkraft, die auf das Teichen wirkt: OK
  Fgesamt = vadd(vadd(A,D), E)  
  write(*,*) "Fgesamt =", Fgesamt

  !Berechne die Beschleunigung, die auf das Teichen wirkt: OK
  veca = vsmul(Fgesamt, 1.d0/MTeilchen)
  write(*,*) "veca    =", veca





  !RUNGE-KUTTA-SOLVER:

  !K1r = dt * v fuer die Ortskoordinaten


  n = 5
  dt = 1.0

  do i = 1, n
     K1r = svmul(dt*i, vecv)                   
     write(*,*) K1r
     !end do

     ! do i = 1, n
     K2r = svmul(dt*i, vadd(vecv, vsmul(K1r,0.5))) !Ergebnisse werden nicht in der zweiten Schleife Uebernommen wieso???
     write(*,*) K2r
  end do

  !  do i = 1, n
  !     K3r = svmul(dt*i, vecv)
  !     write(*,*) K3r
  !  end do

  !  do i = 1, n
  !     K4r = svmul(dt*i, vecv)
  !     write(*,*) K4r
  !  end do


  !K1v = dt * v fuer die Geschw-Koordinaten
  !.
  !.
  !.
  !kommt noch.....
  !.
  !.
  !.


end program TeilchenbahnNeu

! **********************************************************************
! P R O G R A M E N D E
! **********************************************************************
!(Wichtig: Implecit None wird immer nach use definiert)

!DEFINIERE DIE FUNKTIONEN FUER DIE KRAEFTE:
type(vector) function Fgsat(vecs, AbstRheaSaturnVec) !Wichtig: vecs soll den Vektor: (r - d) darstellen
  use Saturndata
  use vectors
  implicit none
  Type (vector), intent(in) :: vecs
  Type (vector), intent(in) :: AbstRheaSaturnVec
  Fgsat = svmul(-gamma * MSaturn * MTeilchen *(1.d0/vabs(vsub(vecs, AbstRheaSaturnVec))**3.d0), vsub(vecs, AbstRheaSaturnVec))    
  return
end function Fgsat


type(vector) function Fgrhe(vecr) 
  use Saturndata
  use vectors
  implicit none
  Type (vector), intent(in) :: vecr
  Fgrhe = svmul(-gamma * MRhea * MTeilchen *(1.d0/vabs(vecr)**3.d0), vecr) 
  return
end function Fgrhe


type(vector) function Fc(m, w, v)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: m!Skalar
  type(vector), intent(in) :: w, v  !2 Vektoren
  Fc = svmul(2.d0 * m, cross(v, w))
  return
end function Fc


type(vector) function Fz(m, w, vecs)
  use Saturndata
  use vectors
  implicit none
  double precision , intent(in) :: m
  type(vector), intent(in) :: w, vecs
  Fz = vsmul(cross(w,cross(w, vecs)), m)
  return
end function Fz


!DEFINIERE DIE FUNKTIONEN FUER DIE BERECHNUNG DER RUNGE-KUTTA-KONSTANTEN:

!type(vector) function Fz(m, w, vecs)
!  use Saturndata
!  use vectors
!  implicit none
!  double precision , intent(in) :: m
!  type(vector), intent(in) :: w, vecs
!  Fz = vsmul(cross(w,cross(w, vecs)), m)
!  return
!end function Fz









!type(vector) function trafo(vec, RheaSaturnVec) !Trafo-Matrix (Spaeter)
!  use Saturndata
!  use vectors
!  implicit none
!  Type (vector), intent(in) :: vec, RheaSaturnVec  
!  trafo = vsub(vec, RheaSaturnVec)!Momentan ist trafo noch auf vec gesetzt so, dass die Umrechnung spaeter erfolgt 
!  return
!end function trafo



























!subroutine Werteausgeben
!  use Saturndata
!  use vectors
!  implicit none
!  write(*,*) "Hallo Du Eimer", dt
!  return
!end subroutine Werteausgeben



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

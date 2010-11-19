module Saturndata
  implicit none
  save

  ! Konstanten *****************************
  integer, parameter :: mem = selected_real_kind(10,12)            ! Variablengroesse
  real(kind=mem) :: gamma = 6.67428e-11                            ![m³/kg*s²]                 (Wiki)
  real(kind=mem) :: d =  527040000.0                               ![m]   Entf. Rhea zu Saturn (Wiki) 
  real(kind=mem) :: MSaturn = 5.685e26                             ![kg]  Masse Saturn         (Wiki)
  real(kind=mem) :: MRhea = 0.000023166e26                         ![kg]  Masse Rhea 2.3166e21 (berechnet) 
  real(kind=mem) :: L1                                             !   Lagrangepunkt        (berechnet)
  real(kind=mem) :: RRhea = 764000.0                               ![m]   Rhearadius           (Wiki) 
  real(kind=mem) :: RSaturn = 57316000.0                           ![m]   Saturnradius         (Wiki)
  real(kind=mem) :: omegaz                                         ![1/s] Winkelgeschw.        (berechnet) 

  ! Globale Variablen **********************
  real(kind=mem) :: dt                          ! Zeitintervall

end module Saturndata



! **********************************************************************
! P R O G R A M
! **********************************************************************

program TeilchenbahnNeu
  use Saturndata
  implicit none   

  ! Variablen ******************************

  real(kind=mem), dimension(3) :: VecX          ! Ortsvektor   
  real(kind=mem), dimension(3) :: VecV          ! Geschwindigkeitsvektor   


  ! Berechnete Variablen *******************
  L1 = d - (d/( sqrt(MRhea/MSaturn)+1.0)) 
  omegaz = sqrt((gamma * MSaturn)/(d**3))

  dt = 42
  call Werteausgeben

end program TeilchenbahnNeu


subroutine Werteausgeben
  use Saturndata
  implicit none


  write(*,*) "Hallo Du Eimer", dt

  return 
end subroutine Werteausgeben
! real(kind=mem) function ()



! Wissen-Konserve:
! COMMON http://www.obliquity.com/computer/fortran/common.html
! MODULE https://srv.rz.uni-bayreuth.de/lehre/fortran90/vorlesung/V08/V08.html

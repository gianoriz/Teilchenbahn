program TeilchenbahnNeu
  implicit none   

  !PROGRAMM: RUNGE-KUTTA-4TER-ORDNUNG (MIT DEN RICHTIGEN BEWEGUNGSGLEICHUNGEN)
  ! Der Letzte Anlauf


  !DEFINITION DER VARIABLEN:
  integer, parameter :: mem = selected_real_kind(10,12)    ! Variablengroesse





 !SETZE KONSTANTEN & PARAMETER:

   real(kind=mem) :: gamma = 6.67428e-11                            ![m³/kg*s²]                 (Wiki)
   real(kind=mem) :: d =  527040000.0                                ![m]   Entf. Rhea zu Saturn (Wiki) 
   real(kind=mem) :: Ms = 5.685e26                                  ![kg]  Masse Saturn         (Wiki)
   real(kind=mem) :: Mr = 0.000023166e26                            ![kg]  Masse Rhea 2.3166e21 (berechnet) 
   real(kind=mem) :: L1 = d - (d/( sqrt(Mr/Ms)+1.0))                ![m]   Lagrangepunkt        (berechnet)
   real(kind=mem) :: Rr = 764000.0                                  ![m]   Rhearadius           (Wiki) 
   real(kind=mem) :: Rs = (120536000.0 + 108728000.0)/4             ![m]   Saturnradius         (Wiki)
   real(kind=mem) :: omegaz = sqrt((gamma * Ms)/(d**3))             ![1/s] Winkelgeschw.        (berechnet) 





!!$  real(kind=mem) ::  dt                            !Flugzeit & Zeitintervall
!!$  real(kind=mem) ::  omegaz                        !Winkelgeschw. in z-Richtung 
!!$  real(kind=mem) ::  gamma                         !Gravitationskonstante       
!!$  real(kind=mem) ::  Ms, Mr                        !Masse Saturn, Masse Rhea    
!!$  real(kind=mem) :: , dimension(9) :: W            !W(i)=(W(1)=x,...,W(4)=v_x,...,W(7)=a_x) Vektor   
!!$  real(kind=mem) :: , dimension(4,6) :: K          !RungeKuttaKonst.
!!$  real(kind=mem) ::  d, L1                         !Abstand Saturn-Rhea;L1,L2=Lagrangepunkte
!!$  real(kind=mem) ::  Rs                            !Saturnradius
!!$  real(kind=mem) ::  Rr                            !Rhearadius              
!!$  integer(kind=8) :: i, t                        !Laufvariablen
!!$  real(kind=mem) ::  Zeit

  !#####################################
  !Koordinatenursprung liegt bei Rhea
  !#####################################


 







end program TeilchenbahnNeu





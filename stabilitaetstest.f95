program TeilchenTrajektorie

  !PROGRAMM: RUNGE-KUTTA-4TER-ORDNUNG (MIT DEN RICHTIGEN BEWEGUNGSGLEICHUNGEN)
  !Jonas schau dir bitte meinen Rungekuttasolver an und vergleiche dies mit dem online numerical recepies!!!
  implicit none   

  !DEFINITION DER VARIABLEN:
  integer, parameter :: rk = SELECTED_REAL_KIND(10,30)
  double precision xr, v_xr, a_xr       !x-Koord. & Geschw. im Rhea-Ruhesystem  
  double precision yr, v_yr, a_yr       !y-Koord. & Geschw. im Rhea-Ruhesystem
  double precision zr, v_zr, a_zr       !z-Koord. & Geschw. im Rhea-Ruhesystem 
  double precision t, dt                !Flugzeit & Zeitintervall
  double precision v_0                  !Anfangsgeschw.              
!  double precision alpha                !Abwurfwinkel                    
  double precision omegaz               !Winkelgeschw. in z-Richtung 
  double precision gamma                !Gravitationskonstante       
  double precision Ms, Mr               !Masse Saturn, Masse Rhea    
  double precision v_r                  !Geschw. von Rhea um Saturn  
  double precision Mj, Mk               !Hilfsgroessen: Mj = 3 bzw. 6 & Mk = 2 bzw. 1
  double precision A                    !Hilfsgr.: eliminiert nichtvorhandene Matrixelemente 

  double precision, dimension(6) :: k1 
  double precision, dimension(6) :: k2
  double precision, dimension(6) :: k3
  double precision, dimension(6) :: k4


  double precision, dimension(9) :: W   !W(i)=(W(1)=x,...,W(4)=v_x,...,W(7)=a_x) Vektor   
  real :: Pi = 3.14159265358979323846   !Ludolphsche Zahl 
  double precision d, L1, L2            !Abstand Saturn-Rhea;L1,L2=Lagrangepunkte
  double precision S                    !Gemeinsamer Schwerpunkt
  double precision Rs                   !Saturnradius
  double precision R                    !Rhearadius              
  integer :: q, n, i, j, schrittzaehler !Laufvariablen


  !#####################################
  !Koordinatenursprung liegt bei Rhea
  !#####################################


  !SETZE KONSTANTEN & PARAMETER:
  gamma = 6.67428e-11                               ![m³/kg*s²]                 (Wiki)
  d = 527040000.0                                   ![m]   Entf. Rhea zu Saturn (Wiki) 
  Ms = 5.685e26                                     ![kg]  Masse Saturn         (Wiki)
  Mr = 0.000023166e26                               ![kg]  Masse Rhea 2.3166e21 (berechnet)
  omegaz = 0.0!sqrt((gamma * Ms)/(d**3))            ![1/s] Winkelgeschw.        (berechnet) 
  L1 = d - (d/( sqrt(Mr/Ms)+1))                     ![m]                        (berechnet)
  R = 764000.0                                      ![m] Rhearadius             (Wiki) 
  Rs = (120536000.0 + 108728000.0)/4                ![m] Saturnradius



  dt   = 1.0                                       ![s]
  n    = 1000000                                   ![Anzahl der Iterationen]

!  Schrittweite: 
do schrittzaehler = 1, 10

  dt = schrittzaehler * 1.0 
  write(*,*) dt



write(*,*) "#######################################################################"



     !write(*,*) L1

     !STARTWERTE DES TEILCHENS IM RUHESYSTEM VON RHEA: 
     W(1) = 0.0  !x  Vertikalachse                                      ![m]    Startpunkt
     W(2) = L1   !y  Horizontalachse                                    ![m]    Startpunkt             
     W(3) = 0.0  !z                                                     ![m]    Startpunkt                                  
     W(4) = 0.0  !omegaz * d * 1.9                                      ![m/s]  Startgeschwindigkeit
     W(5) = 0.0                                                         ![m/s]  Startgeschwindigkeit
     W(6) = 0.0                                                         ![m/s]  Startgeschwindigkeit 
     W(7) = -(gamma * Mr * W(1))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
          -(gamma * Ms * W(1))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) & !Teilchen in Rhea Test
          + 2 * omegaz * W(5) + omegaz**2 * W(1)
     W(8) = -(gamma * Mr * W(2))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
          -(gamma * Ms * (W(2) - d))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) &
          - 2 * omegaz * W(4) + omegaz**2 * (W(2) - d)
     W(9) = -(gamma * Mr * W(3))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) &
          -(gamma * Ms * W(3))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5)




     !Runge-Kutta-Solver:


     do q = 0, n !Grosse do-Zeitschleife!Muss ueberhaupt eine do Zeitschleife rein?


        do  i = 1,6
           k1(i) = dt * W(i+3)
        end do
        W(7) = -(gamma * Mr * W(1))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
             -(gamma * Ms * W(1))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) & !Teilchen in Rhea Test
             + 2 * omegaz * W(5) + omegaz**2 * W(1)
        W(8) = -(gamma * Mr * W(2))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
             -(gamma * Ms * (W(2) - d))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) &
             - 2 * omegaz * W(4) + omegaz**2 * (W(2) - d)
        W(9) = -(gamma * Mr * W(3))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) &
             -(gamma * Ms * W(3))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5)





        do  i = 1,6
           if ((i > 3) .or. (j == 1)) then            !  Beseitigt nichtvorhandene Matrixelem.                   
              A = 0.0                                 !  K(1,7),K(1,8),...,K(1,10)                              
           else                                       
              A = 1.0
           end if
           k2(i) = dt * (W(i+3) + 0.5 * k1(i+3) * A)
        end do
        W(7) = -(gamma * Mr * W(1))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
             -(gamma * Ms * W(1))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) & !Teilchen in Rhea Test
             + 2 * omegaz * W(5) + omegaz**2 * W(1)
        W(8) = -(gamma * Mr * W(2))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
             -(gamma * Ms * (W(2) - d))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) &
             - 2 * omegaz * W(4) + omegaz**2 * (W(2) - d)
        W(9) = -(gamma * Mr * W(3))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) &
             -(gamma * Ms * W(3))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5)





        do i = 1,6
           if (i > 3) then                            !  Beseitigt nichtvorhandene Matrixelem.                   
              A = 0.0                                 !  K(1,7),K(1,8),...,K(1,10)                                 
           else                                       
              A = 1.0
           end if
           k3(i) = dt * (W(i+3) + 0.5 * k2(i+3) * A)
        end do
        W(7) = -(gamma * Mr * W(1))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
             -(gamma * Ms * W(1))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) & !Teilchen in Rhea Test
             - 2 * omegaz * W(5) - omegaz**2 * W(1)
        W(8) = -(gamma * Mr * W(2))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
             -(gamma * Ms * (W(2) - d))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) &
             + 2 * omegaz * W(4) - omegaz**2 * (W(2) - d)
        W(9) = -(gamma * Mr * W(3))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) &
             -(gamma * Ms * W(3))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) 


        do i = 1,6
           if (i > 3)  then                           !  Beseitigt nichtvorhandene Matrixelem.                   
              A = 0.0                                 !  K(1,7),K(1,8),...,K(1,10)                               
           else                                       
              A = 1.0
           end if
           k4(i) = dt * (W(i+3) + k3(i+3) * A)
           W(i) = W(i) + k1(i)/6. + k2(i)/3. + k3(i)/3. + k4(i)/6. !Algorithmus
        end do




        t = q * dt          

        if (W(2) > L1+R) then
           exit 
        end if

        if (W(2) < L1-R) then
           exit
        end if



        !write(*,*) W(2)/R, W(1)/R!W(2)=y-Koordinate in Rhearadien, W(1)=x-Koordinate in Rhearadien
        !write(*,*) W(2)/R, W(5) * (3600/1000) 
         write(*,format(3A, 3X, F3.0, F3.2, 3X, F5.3 )) "dt:", schrittzaehler, t/3600, W(2)/1000
        !write(*,*) t, W(2)-L1!/R


     end do !Ende der grossen do-Zeitschleife 

end do !Schrittweite


end program TeilchenTrajektorie

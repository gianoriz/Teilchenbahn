program TeilchenUmSaturnTest

  !PROGRAMM: RUNGE-KUTTA-4TER-ORDNUNG (MIT DEN RICHTIGEN BEWEGUNGSGLEICHUNGEN)
  implicit none   

  !DEFINITION DER VARIABLEN:
  integer, parameter :: mem = selected_real_kind(10,12)    ! Variablengroesse

  real(kind=mem) :: x, v_x, a_x       !x-Koord. & Geschw. im Rhea-Ruhesystem  
  real(kind=mem) :: y, v_y, a_y       !y-Koord. & Geschw. im Rhea-Ruhesystem
  real(kind=mem) :: z, v_z, a_z       !z-Koord. & Geschw. im Rhea-Ruhesystem 
  real(kind=mem) :: t, dt                !Flugzeit & Zeitintervall
  real(kind=mem) :: omegaz               !Winkelgeschw. in z-Richtung 
  real(kind=mem) :: gamma                !Gravitationskonstante       
  real(kind=mem) :: Ms, Mr               !Masse Saturn, Masse Rhea    



  real(kind=mem) :: k11, k12, k13, k14, k15, k16
  real(kind=mem) :: k21, k22, k23, k24, k25, k26
  real(kind=mem) :: k31, k32, k33, k34, k35, k36
  real(kind=mem) :: k41, k42, k43, k44, k45, k46


  real(kind=mem) :: Pi = 3.14159265358979323846   !Ludolphsche Zahl 
  real(kind=mem) :: d, L1                         !Abstand Saturn-Rhea;L1
  real(kind=mem) :: Saturnradius                  !Saturnradius
  real(kind=mem) :: Rhearadius                    !Rhearadius              

  integer(kind=8) :: counter

  !#####################################
  !Koordinatenursprung liegt bei Rhea
  !#####################################


  !SETZE KONSTANTEN & PARAMETER:
  gamma = 6.67428e-11                               ![m³/kg*s²]                 (Wiki)
  d = 527040000.0                                   ![m]   Entf. Rhea zu Saturn (Wiki) 
  Ms = 5.685e26                                     ![kg]  Masse Saturn         (Wiki)
  Mr = 0.000023166e26                               ![kg]  Masse Rhea 2.3166e21 (berechnet) 
  L1 = d - (d/( sqrt(Mr/Ms)+1.0))                   ![m]                        (berechnet)
  Rhearadius = 764000.0                                      ![m]   Rhearadius           (Wiki) 
  Saturnradius = (120536000.0 + 108728000.0)/4                ![m]   Saturnradius         (Wiki)
  omegaz = sqrt((gamma * Ms)/(d**3))               ![1/s] Winkelgeschw.        (berechnet) 





  !STARTWERTE DES TEILCHENS IM RUHESYSTEM VON RHEA: 
  x   = 0.0           !x  Vertikalachse                                      ![m]    Startpunkt
  y   = Rhearadius    !y  Horizontalachse                                    ![m]    Startpunkt             
  z   = 0.0           !z                                                     ![m]    Startpunkt                                  
  v_x = 0.0                                                         ![m/s]  Startgeschwindigkeit
  v_y = 0.0                                                         ![m/s]  Startgeschwindigkeit
  v_z = 0.0                                                         ![m/s]  Startgeschwindigkeit 
!  a_x = -(gamma * Ms * x)/(x**2 + (y - d)**2 + z**2)**(1.5) &
!       + 2 * omegaz * v_y + omegaz**2 * x
!  a_y = -(gamma * Ms * (y - d))/(x**2 + (y - d)**2 + z**2)**(1.5) &
!       - 2 * omegaz * v_x + omegaz**2 * (y - d)
!  a_z = -(gamma * Ms * z)/(x**2 + (y - d)**2 + z**2)**(1.5)


  dt   = 1000.0                                       ![s]

  !Runge-Kutta-Solver:

  Zeitschleife : do counter = 0, 10000                                    ![Anzahl der Iterationen]

     a_x = -(gamma * Ms * x)/(x**2 + (y - d)**2 + z**2)**(1.5) &
          + 2 * omegaz * v_y + omegaz**2 * x
     a_y = -(gamma * Ms * (y - d))/(x**2 + (y - d)**2 + z**2)**(1.5) &
          - 2 * omegaz * v_x + omegaz**2 * (y - d)
     a_z = -(gamma * Ms * z)/(x**2 + (y - d)**2 + z**2)**(1.5)



     k11 = dt * v_x
     k12 = dt * v_y
 !    k13 = dt * v_z 
     k14 = dt * a_x  
     k15 = dt * a_y
 !    k16 = dt * a_z

     k21 = dt * (v_x + k14/2.)  
     k22 = dt * (v_y + k15/2.)  
  !   k23 = dt * (v_z + k16/2.) 
     k24 = dt * (a_x)
     k25 = dt * (a_y)
!     k26 = dt * (a_z) 

     k31 = dt * (v_x + k24/2.)
     k32 = dt * (v_y + k25/2.) 
 !    k33 = dt * (v_z + k26/2.)
     k34 = dt * (a_x)
     k35 = dt * (a_y)
  !   k36 = dt * (a_z)

     k41 = dt * (v_x + k34)
     k42 = dt * (v_y + k35)
  !   k43 = dt * (v_z + k36)
     k44 = dt * (a_x) 
     k45 = dt * (a_y)
  !   k46 = dt * (a_z)


     !###################################################################
     x   = x   + k11/6. + k21/3. + k31/3. + k41/6.  
     y   = y   + k12/6. + k22/3. + k32/3. + k42/6.  
!     z   = z   + k13/6. + k23/3. + k33/3. + k43/6.  
     v_x = v_x + k14/6. + k24/3. + k34/3. + k44/6.  
     v_y = v_y + k15/6. + k25/3. + k35/3. + k45/6.    
!     v_z = v_z + k16/6. + k26/3. + k36/3. + k46/6.           !Algorithmus 
     !###################################################################


     t = counter * dt          


     write(*,*) y/Rhearadius, x/Rhearadius

  end do Zeitschleife

end program TeilchenUmSaturnTest

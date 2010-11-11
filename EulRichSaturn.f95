program SaturnRheaTeilchen

  implicit none   

  !Definition der verwendeten Variablen:
  double precision x, y, z                  !Ortskoordinaten
  double precision v_x, v_y, v_z            !Geschwindigkeitskoordinaten
  double precision a_x, a_y, a_z            !Beschleunigungskoordinaten
  double precision x_mid, y_mid, z_mid      !Mittelwerte der Ortskoordinaten
  double precision v_xmid, v_ymid, v_zmid   !Mittelwerte der Geschwindigkeitskoord.
  double precision a_xmid, a_ymid, a_zmid   !Mittelwerte der Beschleunigungskoord.
  double precision x_neu, y_neu, z_neu      !Neue Koordinaten
  double precision v_xneu, v_yneu, v_zneu   !Neue Koordinaten
  double precision dt                       !Zeitintervall
  double precision t                        !Zeit
  double precision gamma                    !Gravitationskonstante
  double precision Mr, Ms                   !Masse Rhea, Masse Saturn
  double precision Rr, Rs                   !Radius Rhea, Radius Saturn
  double precision d                        !Abstand Saturn-Rhea  
  double precision L1                       !Lagrangepunkt
  double precision omegaz                   !Winkelgeschw. in z-Richt.
  integer*8 :: i, n                         !Laufvariablen
  integer :: schrittzaehler                 ! Fuer die Hauptschleife mit dt

  !Werte einsetzten: 
  gamma = 6.67428e-11                       ![m³/kg*s²]                 (Wiki)  
  Ms = 5.685e26                             ![kg]  Masse Saturn         (Wiki)
  Mr = 0.000023166e26                       ![kg]  Masse Rhea           (Wiki)
  Rr = 764000.0                             ![m] Rhearadius             (Wiki) 
  Rs = (120536000.0 + 108728000.0)/4        ![m] Saturnradius           (Wiki) 
  d  = 527040000.0                          ![m] Entf. Rhea zu Saturn   (Wiki) 
  L1 = d - (d/( sqrt(Mr/Ms)+1))             ![m] Lagrangepunkt




  ! Startwerte: (sind angegeben in s, m, m/s, m/s²)
  !t = 0.0 
  x = 0.0 
  y = L1  
  z = 0.0
  v_x = 0.0
  v_y = 0.0
  v_z = 0.0 
  a_x = -(gamma * Mr * x)/(x**2 + y**2 + z**2)**(1.5) & 
       -(gamma * Ms * x)/(x**2 + (y - d)**2 + z**2)**(1.5) & 
       - 2 * omegaz * v_y - omegaz**2 * x
  a_y = -(gamma * Mr * y)/(x**2 + y**2 + z**2)**(1.5) & 
       -(gamma * Ms * (y - d))/(x**2 + (y - d)**2 + z**2)**(1.5) &
       + 2 * omegaz * v_x - omegaz**2 * (y - d)
  a_z = -(gamma * Mr * z)/(x**2 + y**2 + z**2)**(1.5) &
       -(gamma * Ms * z)/(x**2 + (y - d)**2 + z**2)**(1.5)




  dt = 1.0      ![s]       
  n = 100000    !Anzahl der Iterationen
  !  Schrittweite: 
  do schrittzaehler = 1, 10

     dt = schrittzaehler * 1.0 
     write(*,*) dt



     write(*,*) "#######################################################################"



     Zeitschleife : do i = 0, n                    


        a_x = -(gamma * Mr * x)/(x**2 + y**2 + z**2)**(1.5) & 
             -(gamma * Ms * x)/(x**2 + (y - d)**2 + z**2)**(1.5) & 
             + 2 * omegaz * v_y + omegaz**2 * x
        a_y = -(gamma * Mr * y)/(x**2 + y**2 + z**2)**(1.5) & 
             -(gamma * Ms * (y - d))/(x**2 + (y - d)**2 + z**2)**(1.5) &
             - 2 * omegaz * v_x + omegaz**2 * (y - d)
        a_z = -(gamma * Mr * z)/(x**2 + y**2 + z**2)**(1.5) &
             -(gamma * Ms * z)/(x**2 + (y - d)**2 + z**2)**(1.5)



        v_xmid =  v_x + a_x * dt/2 !mittlere Geschw. 
        v_ymid =  v_y + a_y * dt/2
        v_zmid =  v_z + a_z * dt/2


        x_mid = x +  v_x * dt/2    !mittlere Ort
        y_mid = y +  v_y * dt/2
        z_mid = z +  v_z * dt/2


        a_xmid = -(gamma * Mr * x_mid)/(x_mid**2 + y_mid**2 + z_mid**2)**(1.5) &        !mittlere Beschl. 
             -(gamma * Ms * x_mid)/(x_mid**2 + (y_mid - d)**2 + z_mid**2)**(1.5) & 
             + 2 * omegaz * v_ymid + omegaz**2 * x_mid

        a_ymid = -(gamma * Mr * y_mid)/(x_mid**2 + y_mid**2 + z_mid**2)**(1.5) & 
             -(gamma * Ms * (y_mid - d))/(x_mid**2 + (y_mid - d)**2 + z_mid**2)**(1.5) &
             - 2 * omegaz * v_xmid + omegaz**2 * (y_mid - d)

        a_zmid = -(gamma * Mr * z_mid)/(x_mid**2 + y_mid**2 + z_mid**2)**(1.5) &
             -(gamma * Ms * z_mid)/(x_mid**2 + (y_mid - d)**2 + z_mid**2)**(1.5)


        v_xneu = v_x + a_xmid * dt !neue Geschw.
        v_yneu = v_y + a_ymid * dt
        v_zneu = v_z + a_zmid * dt


        x_neu = x + v_xmid * dt    !neue Koordinaten
        y_neu = y + v_ymid * dt
        z_neu = z + v_zmid * dt






        x = x_neu                  !Puffer, aktualisieren des Ortsvektors
        y = y_neu
        z = z_neu



        v_x = v_xneu               !Puffer, aktualisieren des Geschwindigkeitsvektors
        v_y = v_yneu 
        v_z = v_zneu


        t = i * dt                         



        !write(*,*) r(1)/Radius, r(2)/Radius           !Zuerst x-Achse, dann die y-Achse 
        !write(*,*) r(1)/Radius, v(1)* (3600/1000)     !Geschw. vs. Ort
        !write(*,*) t/3600, (y-L1)/Rr                  !Zeit[h], Abweichung [Rhearadien]


        if (y < L1+Rr) then
           exit 
        end if

        if (y > L1-Rr) then
           exit
        end if




        write(*,523) "dt:", schrittzaehler, t/3600, y/1000
523     format(A3, 1X, I3, 3X, F10.3, 3X, F10.3)




     end do Zeitschleife

  end do !Schrittweite

end program SaturnRheaTeilchen

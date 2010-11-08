program TeilchenTrajektorie

  !PROGRAMM: RUNGE-KUTTA-4TER-ORDNUNG (MIT DEN RICHTIGEN BEWEGUNGSGLEICHUNGEN)

  implicit none   

  !DEFINITION DER VARIABLEN:
  integer, parameter :: rk = SELECTED_REAL_KIND(10,30)
  double precision xr, v_xr, a_xr       !x-Koord. & Geschw. im Rhea-Ruhesystem  

  !***********************************
  ! wie gross ist double precision? Am besten zur Erinnerung in Kommentar schreiben. KEINE AHNUNG, WO SEHE ICH NACH???
  !***********************************

  double precision yr, v_yr, a_yr       !y-Koord. & Geschw. im Rhea-Ruhesystem
  double precision zr, v_zr, a_zr       !z-Koord. & Geschw. im Rhea-Ruhesystem 
  double precision t, dt                !Flugzeit & Zeitintervall
  double precision v_0                  !Anfangsgeschw.              
  double precision alpha                !Abwurfwinkel wird spaeter benutzt                    
  double precision omegaz               !Winkelgeschw. in z-Richtung 
  double precision v_r                  !Geschw. von Rhea um Saturn  
  double precision L1 
  double precision A                    !Hilfsgr.: eliminiert nichtvorhandene Matrixelemente 
  !***********************************
  ! verstehe ich nicht. sieht mir sehr sehr seltsam aus. Flicken wir spaeter. OK, DAS IST WIRKLICH EIN PROBLEM ERKLAERE ICH SPAETER 
  !***********************************


  double precision, dimension(6) :: k1 
  double precision, dimension(6) :: k2
  double precision, dimension(6) :: k3
  double precision, dimension(6) :: k4
  double precision, dimension(9) :: W   !W(i)=(W(1)=x,...,W(4)=v_x,...,W(7)=a_x) Vektor   


  !Konstanten:
  real, parameter:: PI = 3.14159265358979323846 !Ludolphsche Zahl
  real, parameter:: gamma = 6.67428*10**(-11.0) ![m^3/(kg*s^2)]   Gravitationskonstante        (Wiki)   
  real, parameter:: Ms = 5.685*10**(26.0)       ![kg] Masse Saturn                             (Wiki)
  real, parameter:: Mr = 0.000023166*10**(26.0) ![kg] Masse Rhea                               (Wiki)
  real, parameter:: R = 764000.0                ![m]  Rhearadius                               (Wiki) 
  real, parameter:: Rs = 57316000               ![m]  Saturnradius=(120536000.0+108728000.0)/4 (Wiki)  
  real, parameter:: d = 527040000.0             ![m]  Entf. Rhea zu Saturn                     (Wiki) 

  integer*8 :: q, n, i, j                       !Laufvariablen


  !***********************************
  ! wie gross ist Datentyp REAL? Genau nachpruefen. Ist Pi so OK??? WEISS ICH NICHT!
  !***********************************


 



  !#####################################
  !Koordinatenursprung liegt in Rhea
  !#####################################

                              
  omegaz = 0.0!sqrt((gamma * Ms)/(d**3.0))            ![1/s] Winkelgeschw.      
  L1 = d - (d/( sqrt(Mr/Ms)+1.0))                     ![m]                       

  !***********************************
  ! Ausrechnen, einsetzen! NEIN WIRD ZU UNGENAU
  !***********************************



  dt   = 10.0                                        ![s]
  n    = 1000.0                                       ![Anzahl der Iterationen]



  !STARTWERTE DES TEILCHENS IM RUHESYSTEM VON RHEA: 
  W(1) = 0.0  !x  Vertikalachse                                      ![m]    Startpunkt
  W(2) = L1   !y  Horizontalachse                                    ![m]    Startpunkt             
  W(3) = 0.0  !z                                                     ![m]    Startpunkt                                  
  W(4) = 0.0  !omegaz * d * 1.9                                      ![m/s]  Startgeschwindigkeit
  W(5) = 0.0                                                         ![m/s]  Startgeschwindigkeit
  W(6) = 0.0                                                         ![m/s]  Startgeschwindigkeit 
  W(7) = -(gamma * Mr * W(1))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
       -(gamma * Ms * W(1))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) & !Teilchen in Rhea Test
       - 2 * omegaz * W(5) - omegaz**2 * W(1)
  W(8) = -(gamma * Mr * W(2))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
       -(gamma * Ms * (W(2) - d))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) &
       + 2 * omegaz * W(4) - omegaz**2 * (W(2) - d)
  W(9) = -(gamma * Mr * W(3))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) &
       -(gamma * Ms * W(3))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5)




  !Runge-Kutta-Solver:


  do q = 0, n !Grosse do-Zeitschleife!Muss ueberhaupt eine do Zeitschleife rein?

     !***********************************
     ! Stufe 3
     ! Keine Ahnung... Das sollte der Programmierer wissen... Unbedingt erst ein Konzept auf Papier machen.
     ! Aus dem Bauch raus vermute ich aber mal dass Du diese Schleife brauchst. Aendere Dein Programm aber 
     ! besser bald mal um, so dass Abbruch, wenn Delta(Ergebnis) < Zielgenauigkeit  
     !***********************************


     do  i = 1,6
        k1(i) = dt * W(i+3)
     end do

     !***********************************
     ! Den W's widmen wir uns, in Stufe 2 WAS MEINST DU MIT STUFE 2?
     !***********************************

     W(7) = -(gamma * Mr * W(1))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
          -(gamma * Ms * W(1))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) & !Teilchen in Rhea Test
          - 2 * omegaz * W(5) - omegaz**2 * W(1)
     W(8) = -(gamma * Mr * W(2))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) & 
          -(gamma * Ms * (W(2) - d))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5) &
          + 2 * omegaz * W(4) - omegaz**2 * (W(2) - d)
     W(9) = -(gamma * Mr * W(3))/(W(1)**2 + W(2)**2 + W(3)**2)**(1.5) &
          -(gamma * Ms * W(3))/(W(1)**2 + (W(2) - d)**2 + W(3)**2)**(1.5)



     do  i = 1,6 
        if (i > 3) then                            !  Beseitigt nichtvorhandene Matrixelem.                   
           A = 0.0                                 !  K(1,7),K(1,8),...,K(1,10)                              
        else                                       
           A = 1.0
        end if
        k2(i) = dt * (W(i+3) + 0.5 * k1(i+3) * A)
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


     !write(*,*) W(2)/R, W(1)/R!W(2)=y-Koordinate in Rhearadien, W(1)=x-Koordinate in Rhearadien
     !write(*,*) W(2)/R, W(5) * (3600/1000) 
     write(*,*) t/3600, W(2)-L1!/R

     !***********************************
     ! was soll da Ausgegeben werden? Spalte 2 ist hier immer 0 DIE SPALTE SOLL AUCH NULL SEIN ICH WILL NUR DIE AENDERUNG SEHEN WILL!
     ! lass t erstmal in Sekunden. Umrechnen erst wenn alles laeuft
     !***********************************


  end do !Ende der grossen do-Zeitschleife 

end program TeilchenTrajektorie

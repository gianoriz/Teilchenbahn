program Pendel
  implicit none


  double precision t, dt             !Flugzeit & Zeitintervall
  integer(kind=8) :: i               !Laufvariable
 


double precision, dimension(1000000) :: x
double precision :: y
double precision :: K1, K2, K3, K4
double precision :: k = -9.81 ! 1/s^2
double precision :: v
double precision :: a


x(1) = 5.d0
x(2) = 5.d0

dt = 0.0001

i = 1


do while (i .lt. 1000000)   

i = i+1
t = dt * i

a = k * x(i)
v = (x(i-1)- x(i-2) )/ dt

! K Ort

!!$K1 = x(i)  
!!$K2 = x(i) + dt * K1/2.d0
!!$K3 = x(i) + dt * K2/2.d0
!!$K4 = x(i) + dt * K3 

x(i+1) = x(i) + v * dt + a * dt**2

! x(i+1) = x(i) + dt * (K1 + 2.d0* K2 + 2.d0 *  K3 + K4)/6.d0

!t = i * dt
write(*,*) x(i), t  

end do 

! write(*,*) x


 

!dt = 1.d0

!do i = 1, 10000

!t = dt * i
!write(*,*) y
!end do 











end program Pendel

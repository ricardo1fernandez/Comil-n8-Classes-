! water_rocket_real.f90
! Caso real: botella de 1L, 1/3 agua, presión inicial 3 bar gauge (4 bar absolutos).

program water_rocket
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: g = 9.80665_dp
  real(dp), parameter :: rho_air = 1.225_dp
  real(dp), parameter :: rho_water = 997.0_dp
  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: Pa = 101325.0_dp

  ! ---- Parámetros físicos del cohete ----
  real(dp) :: mass_dry = 0.15_dp            ! kg (botella + nose + aletas)
  real(dp) :: mass_water0 = 0.33_dp         ! 1/3 L = 0.33 kg
  real(dp) :: P0 = 4.0e5_dp                 ! 4 bar absolutos (3 bar manométricos)
  real(dp) :: V_bottle = 1.0e-3_dp          ! 1 L = 0.001 m^3
  real(dp) :: An = 5.0e-5_dp                ! área tobera (8 mm diámetro)
  real(dp) :: Cd = 0.8_dp
  real(dp) :: Cd_drag = 0.5_dp
  real(dp) :: A_front = 0.005_dp            ! área frontal ~ 5 cm^2
  real(dp) :: launch_angle_deg = 60.0_dp
  real(dp) :: launch_height = 0.0_dp

  real(dp) :: dt = 0.001_dp, tmax = 8.0_dp

  ! ---- Variables ----
  real(dp) :: t,x,y,vx,vy,theta
  real(dp) :: mass_water,mass_total
  real(dp) :: V_air0,V_air,P,const_PV
  real(dp) :: mdot,Ve,thrust,Pe
  real(dp) :: v,drag,Fx,Fy
  integer :: i,nsteps

  integer :: ios

  ! RK4 temp vars
  real(dp) :: k1x,k1y,k1vx,k1vy
  real(dp) :: k2x,k2y,k2vx,k2vy
  real(dp) :: k3x,k3y,k3vx,k3vy
  real(dp) :: k4x,k4y,k4vx,k4vy

  theta = launch_angle_deg * acos(-1.0_dp) / 180.0_dp

  ! Iniciales
  x=0; y=launch_height; vx=0; vy=0
  mass_water = mass_water0
  mass_total = mass_dry + mass_water

  V_air0 = V_bottle - mass_water0/rho_water
  const_PV = P0 * V_air0**gamma
  P = P0
  Pe = Pa

  open(10,file="trajectory.csv",status="replace",action="write")

  write(10,'(A)') "t,x,y,vx,vy,mass_total,P"

  nsteps = int(tmax/dt)
  t = 0.0_dp

  do i=1,nsteps

    write(10,'(F8.4,",",F10.6,",",F10.6,",",F9.4,",",F9.4,",",F10.6,",",F10.2)') &
      t,x,y,vx,vy,mass_total,P

    if (y <= 0.0_dp .and. t > 0.1_dp) exit

    ! ::::: EMPUJE ::::::
    if (mass_water > 1e-8_dp .and. P > Pa) then
      V_air = V_bottle - mass_water/rho_water
      Ve = Cd * sqrt( max(0.0_dp, 2.0_dp*(P-Pa)/rho_water) )
      mdot = rho_water * An * Ve
      if (mdot*dt > mass_water) mdot = mass_water / dt
      thrust = mdot * Ve
    else
      thrust = 0.0_dp
      mdot = 0.0_dp
    end if

    ! ::::: AERODINÁMICA ::::::
    v = sqrt(vx*vx + vy*vy)
    if (v > 1e-12_dp) then
      drag = 0.5_dp * rho_air * Cd_drag * A_front * v*v
      Fx = -drag * vx/v
      Fy = -drag * vy/v
    else
      Fx = 0; Fy = 0
    end if

    ! ::::: INTEGRACIÓN (RK4) ::::::
    ! k1
    k1x = vx
    k1y = vy
    k1vx = (thrust*cos(theta) + Fx)/mass_total
    k1vy = (thrust*sin(theta) + Fy)/mass_total - g

    ! k2
    k2x = vx + 0.5_dp*dt*k1vx
    k2y = vy + 0.5_dp*dt*k1vy
    k2vx = k1vx
    k2vy = k1vy

    ! k3
    k3x = vx + 0.5_dp*dt*k2vx
    k3y = vy + 0.5_dp*dt*k2vy
    k3vx = k2vx
    k3vy = k2vy

    ! k4
    k4x = vx + dt*k3vx
    k4y = vy + dt*k3vy
    k4vx = k3vx
    k4vy = k3vy

    ! Actualización
    x  = x  + dt*(k1x + 2*k2x + 2*k3x + k4x)/6
    y  = y  + dt*(k1y + 2*k2y + 2*k3y + k4y)/6
    vx = vx + dt*(k1vx + 2*k2vx + 2*k3vx + k4vx)/6
    vy = vy + dt*(k1vy + 2*k2vy + 2*k3vy + k4vy)/6

    ! ::::: MASA Y PRESIÓN ::::::
    mass_water = mass_water - mdot*dt
    if (mass_water < 0) mass_water = 0
    V_air = V_bottle - mass_water/rho_water
    if (V_air < 1e-8_dp) V_air = 1e-8_dp
    P = const_PV / V_air**gamma

    mass_total = mass_dry + mass_water
    t = t + dt

  end do

  close(10)
  print *, "Simulación finalizada. Datos guardados en trajectory.csv"

end program water_rocket


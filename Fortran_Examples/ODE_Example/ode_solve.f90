program ode_solve
    use runge_kutta_2, only: rk2
    use real_type_mod, only: dp
    implicit none
    integer :: num_dims, beginning, end
    real(dp), allocatable :: y0(:)
    real(dp) :: t0, tend, dt, rate
    character(len=30) :: filename_solution, filename_times

    call system_clock(beginning, rate)

    ! Set up files for outputting data
    filename_solution = "fortran_soln.txt"
    filename_times = "fortran_soln_times.txt"

    ! Number of equations to solve
    num_dims = 1

    ! Set simulation parameters
    allocate(y0(num_dims)) ! Set the dimension of the initial data vector
    y0 = 0.01_dp
    t0 = 0.0_dp
    tend = 5.0_dp
    dt = 0.01_dp

    
    ! Save initial data values
    open(unit=10, file=filename_solution, status="replace", action="write")
    write(10, *) y0

    open(unit=11, file=filename_times, status="replace", action="write")
    write(11, *) t0

    
    ! Main rk2 loop
    do while (t0 < tend)
        call rk2(t0,y0,dt,num_dims,ODE)
        t0 = t0 + dt

        ! write solution data to txt files
        write(10, *) y0
        write(11, *) t0
    
    end do
    close(10)
    close(11)

    call system_clock(end)
    print *, "elapsed time: ", real(end - beginning) / real(rate)

contains

! Exponential growth/decay ODE
! Solution is exp(2*t)
subroutine ODE(t,y,neqn,dydt)
    implicit none
    integer, intent(in) :: neqn
    real(dp), intent(in) :: t, y(neqn)
    real(dp), intent(out) :: dydt(neqn)

    dydt = 2*y

end subroutine ODE

end program ode_solve
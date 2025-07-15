program evolve_advection
    use real_type_mod, only: dp
    use time_stepping, only: rk4
    use boundary_conditions, only: periodic_bc
    use finite_difference, only: FD_C4
    implicit none

    integer :: Nx, Ngz, i, k, beginning, end, steps
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: x_min, x_max, dx, CFL, t0, tend, dt, rate
    real(dp), allocatable ::  grid(:)
    real(dp), allocatable :: y0(:,:)
    character(len=30) :: filename_solution, filename_times

    call system_clock(beginning, rate)

    !#########################################
    ! Set simulation parameters
    !########################################
    Nx = 501 ! Grid points
    Ngz = 2 ! Ghost points

    ! Boundaries of grid [x_min,x_max]
    x_min = 0.0_dp 
    x_max = 2.0_dp*pi

    ! Spatial step size, CFL number, timestep size
    dx = (x_max - x_min)/(real(Nx,dp)-1.0_dp)
    CFL = 0.5_dp
    dt = CFL*dx ! Set timestep with CFL condition

    ! Compute grid point values
    allocate(grid(1-Ngz:Nx+Ngz))
    do i = 1-Ngz, Nx+Ngz
        grid(i) = x_min + (real(i,dp) -1.0_dp)*dx
    end do

    !#########################################
    ! Set initial data
    !########################################
    allocate(y0(1-Ngz:Nx+Ngz,1))

    y0(:,1) = sin(grid)

    call periodic_bc(Nx,Ngz,y0)

    !#########################################
    ! Time interval 
    !########################################
    t0 = 0.0_dp
    tend = 20.0_dp

    !#########################################
    ! Set up files for saving data
    !########################################

    ! Name files
    filename_solution = "advection_soln.dat"
    filename_times = "advection_soln_times.dat"

    ! Save initial data
    open(unit=10, file=filename_solution, status="replace", action="write")
    do k = 1, Nx
        write(10, '(es21.12E3,a1)',advance='NO') y0(k,1)
    end do
    write(10,*)

    open(unit=11, file=filename_times, status="replace", action="write")
    write(11,'(es21.12E3,a1)') t0

    ! close(10)
    ! close(11)

    !#########################################
    ! Evolution Routine
    !########################################
    steps = 0
    do while (t0<tend)

        steps = steps + 1
        call rk4(Nx,Ngz,dx,t0,y0,dt)
        t0 = t0 + dt

        if (steps == 100) then
            do k = 1, Nx
                write(10, '(es21.12E3,a1)',advance='NO') y0(k,1)
            end do
            write(10, *)
            ! write(10) y0
            write(11,'(es21.12E3,a1)') t0
            steps = 0
        end if

    end do

    close(10)
    close(11)

    call system_clock(end)
    print *, "elapsed time: ", real(end - beginning) / real(rate)

    

end program evolve_advection
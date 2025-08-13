program evolve_advection
    use real_type_mod, only: dp
    use time_stepping_mod, only: rk4
    use boundary_conditions_mod, only: periodic_bc
    use burgers_mod, only: compute_CS_tstep
    use read_parameters_mod
    implicit none

    integer :: beginning, steps, end
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: CFL, dt, CS, rate
    real(dp), allocatable :: y0(:,:)
    character(len=30) :: filename_solution, filename_times


    ! Buffer for efficient I/O
    integer, parameter :: buffer_size = 100
    real(dp), allocatable :: solution_buffer(:,:)
    real(dp) :: time_buffer(buffer_size)
    integer :: buffer_index = 0

    call system_clock(beginning, rate)


    !-------------------------------------
    ! Load Simulation Parameters
    !-------------------------------------
    call read_params()

    ! Initial timestep size
    CFL = 0.5
    dt = CFL*dx

    ! Allocate I/O Buffer
    allocate(solution_buffer(Nx, buffer_size))


    !#########################################
    ! Set initial data
    !########################################
    allocate(y0(1-Ngz:Nx+Ngz,Nvar))

    ! y0(:,Nvar) = 2.0_dp*sin(grid)

    y0(:,Nvar) = 10.0_dp*exp(-20*(grid-pi/2.0_dp)**2)

    call periodic_bc(Nx,Ngz,Nvar,y0)

    print *, t0, tend


    !#########################################
    ! Set up files for saving data
    !########################################

    ! Name files
    filename_solution = "Burgers_soln.dat"
    filename_times = "Burgers_soln_times.dat"


    ! Create binary files
    open(unit=10, file=filename_solution, status="replace", action="write", &
         access="stream", form="unformatted")
    open(unit=11, file=filename_times, status="replace", action="write", &
         access="stream", form="unformatted")

    ! Save initial data to buffer arrays
    buffer_index = 1
    solution_buffer(:, buffer_index) = y0(1:Nx, 1)
    time_buffer(buffer_index) = t0

    

    !#########################################
    ! Evolution Routine
    !########################################
    steps = 0
    do while (t0<tend)

        steps = steps + 1

        ! Adjust timestep
        call compute_CS_tstep(y0,CS)
        dt = min(CFL*dx/CS,CFL*dx)

        call rk4(t0,y0,dt)
        t0 = t0 + dt


        ! Save solution data every 100 timesteps
        if (steps == 10) then 
            buffer_index = buffer_index + 1
            solution_buffer(:, buffer_index) = y0(1:Nx, 1)
            time_buffer(buffer_index) = t0
            
            ! Save buffer array to file when
            ! it is full
            if (buffer_index == buffer_size) then
                call write_buffer_to_file(10, 11, solution_buffer, time_buffer, &
                                        Nx, buffer_size, buffer_index)
                buffer_index = 0
            end if
            
            steps = 0
        end if

    end do

    ! Save any remaining data in buffer array
    if (buffer_index > 0) then
        call write_buffer_to_file(10, 11, solution_buffer, time_buffer, &
                                Nx, buffer_size, buffer_index)
    end if

    close(10)
    close(11)

    call system_clock(end)
    print *, "elapsed time: ", real(end - beginning) / real(rate)




contains 

    subroutine write_buffer_to_file(unit_soln, unit_time, soln_buffer, time_buffer, &
                                   Nx, buffer_size, N_entries)
        integer, intent(in) :: unit_soln, unit_time, Nx, buffer_size, N_entries
        real(dp), intent(in) :: soln_buffer(nx, buffer_size)
        real(dp), intent(in) :: time_buffer(buffer_size)
        integer :: j
        
        ! Write solution data
        do j = 1, n_entries
            write(unit_soln) soln_buffer(:, j)
        end do
        
        ! Write time data
        do j = 1, n_entries
            write(unit_time) time_buffer(j)
        end do
        
        ! Flush buffers to ensure data is written
        flush(unit_soln)
        flush(unit_time)

    end subroutine write_buffer_to_file

    

end program evolve_advection
module read_parameters_mod
    use real_type_mod, only: dp
    implicit none
    real(dp) :: x_min, x_max, dx, t0, tend
    integer :: Nx, Ngz, Nvar, i
    real(dp), allocatable :: grid(:)

    public x_min, x_max, dx, t0, tend, Nx, Ngz, Nvar, grid, read_params
    private i

    contains

    subroutine read_params 

        ! Load in parameters from namelist
        namelist /Grid_Params/ x_min, x_max, Nx, Ngz, Nvar
        namelist /Sim_Params/ t0, tend

        open(file='params.nml',status='old',unit=10,action='read')
        read(10,nml=Grid_Params)
        read(10,nml=Sim_Params)
        close(10)

        ! Created derived parameters
        dx = (x_max - x_min)/(real(Nx,dp))

        allocate(grid(1-Ngz:Nx+Ngz))

        do i = 1-Ngz, Nx+Ngz
            grid(i) = x_min + (real(i,dp) -0.5_dp)*dx
        end do

    
    end subroutine read_params



end module read_parameters_mod
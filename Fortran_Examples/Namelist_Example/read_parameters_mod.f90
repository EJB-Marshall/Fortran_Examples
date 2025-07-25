module read_parameters
    use real_type_mod, only: dp
    implicit none
    real(dp) :: x_min, x_max, dx
    integer :: Nx, Ngz, i
    real(dp), allocatable :: grid(:)

    contains

    subroutine read_params 

        ! Load in parameters from namelist
        namelist /Grid_Params/ x_min, x_max, Nx, Ngz

        open(file='params.nml',status='old',unit=10,action='read')
        read(10,nml=Grid_Params)
        close(10)

        ! Created derived parameters
        dx = (x_max - x_min)/(real(Nx,dp))

        allocate(grid(1-Ngz:Nx+Ngz))

        do i = 1-Ngz, Nx+Ngz
            grid(i) = x_min + (real(i,dp) -0.5_dp)*dx
        end do

    
    end subroutine read_params



end module read_parameters
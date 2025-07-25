program namelist_read
    use read_parameters
    implicit none

    call read_params

    ! Print out the parameters:
    print *, "The grid domain is: ", x_min, x_max
    print *, "The number of cell centres is: ", Nx
    print *, "The number of ghost points at each boundary is: ", Ngz
    print *, "The step size is: ", dx
    print *, "The cell centres are: ", grid

end program namelist_read
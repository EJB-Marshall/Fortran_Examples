program io
    use real_type_mod, only: dp
    implicit none
    character(len=20) :: filename
    real(dp) ::  x,y

    x = 1.0_dp
    y = 12.0_dp

    filename = 'result.txt'

    ! Write to a text file
    open(unit=1,file=filename,status='replace',form='formatted')
    write(1,*) x,y
    close(1)

    ! Reset the variables
    x = 0.0_dp
    y = 0.0_dp

    ! Read in from text file
    open(unit=2,file=filename,status='old')
    read(2,*) x,y
    close(2)

    ! Print the variables
    print *, 'x = ', x, 'y = ', y

end program io
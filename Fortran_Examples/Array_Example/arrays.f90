module arrays
    implicit none

    contains

    subroutine array_operations()
        use real_type_mod, only: dp
        real(dp) :: x(3)
        real(dp) :: z(2,2), y(2,2)

        x = 1.0_dp ! Set all entries of x to 1.0
        print *,"x is", x

        x(1) = 0.62_dp ! Set the first entry to 0.62

        print *,"The updated x is", x

        z = 1.0_dp ! Set all entries of z to 1.0

        print *, "z is ", z

        ! Set the first row of y to be -0.5 and the second to 0.5
        y(1:2,1) = -0.5_dp
        y(1:2,2) = 0.5_dp 

        print *, "y is ", y

        ! Add y and z elementwise
        y = y + z 

        print *, "y+z is ", y

    
    end subroutine array_operations
    
end module arrays
! Module containing the boundary conditions
module boundary_conditions
    use real_type_mod, only: dp
    implicit none

    contains 

    subroutine periodic_bc(Nx,Ngz,u)
        ! Inputs:            
        ! dx - The spatial step size, real(dp)
        ! Nx - The number of grid points, integer
        ! u - The function we apply the BCs to, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        ! Output:
        ! u - The function with updated boundaries, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        implicit none
        integer, intent(in) :: Nx, Ngz
        real(dp), intent(inout) :: u(1-Ngz:Nx+Ngz,1)
        ! integer :: i 

        u(Nx,1) = u(1,1) ! First and last grid points are identified

        ! Right Boundary
        u(Nx+1:Nx+Ngz,1) = u(1+1:1+Ngz,1)

        ! Left Boundary
        u(1-Ngz:0,1) = u(Nx-Ngz:Nx-1,1)

        ! do i = 1,Ngz

        !     ! Right Boundary
        !     u(1,Nx+i) = u(1,1+i)

        !     ! Left Boundary
        !     u(1,1-i) = u(1,Nx-i)
        
        ! end do

    end subroutine periodic_bc

end module boundary_conditions
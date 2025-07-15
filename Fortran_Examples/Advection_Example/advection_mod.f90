! Module containing the advection equation subroutine
module advection
    use real_type_mod, only: dp
    use finite_difference, only: FD_C4
    implicit none

    contains 

    subroutine compute_rhs(Nx,Ngz,dx,t,u,dtu)
        implicit none
        integer, intent(in) :: Nx, Ngz
        real(dp), intent(in) :: t, dx
        real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,1)
        real(dp) :: du(1-Ngz:Nx+Ngz,1)
        real(dp), intent(out) :: dtu(1-Ngz:Nx+Ngz,1)
        

       call  FD_C4(dx,Nx,Ngz,u,du) ! Compute the derivative

       dtu = -1.0_dp*du

        

    end subroutine compute_rhs

end module advection
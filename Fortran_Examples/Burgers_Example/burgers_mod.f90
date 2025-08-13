! Module containing the advection equation subroutine
module burgers_mod
    use real_type_mod, only: dp
    use read_parameters_mod
    use diff_op_mod
    implicit none

    contains 

    subroutine compute_rhs(t,u,dtu)
        implicit none
        real(dp), intent(in) :: t 
        real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,Nvar)
        real(dp) :: u_plus(1-Ngz:Nx+Ngz,Nvar), u_minus(1-Ngz:Nx+Ngz,Nvar)
        real(dp), intent(out) :: dtu(1-Ngz:Nx+Ngz,Nvar)
        real(dp) :: flux(1-Ngz:Nx+Ngz,Nvar)
        real(dp) :: char_speed(1-Ngz:Nx+Ngz,1)

        call linear_reconstruct(u,u_plus,u_minus)
        
        ! Set local max characteristic speed
        where (abs(u_plus) > abs(u_minus))
            char_speed = abs(u_plus)
        elsewhere (abs(u_plus) < abs(u_minus))
            char_speed = abs(u_minus)
        end where


        call local_lax_friedrichs(0.5_dp*u_plus**2,0.5_dp*u_minus**2,&
                                    u_plus,u_minus,char_speed,flux)


        dtu(1:Nx,:) = (-1/dx)*(flux(1:Nx,:) - flux(0:Nx-1,:))


    end subroutine compute_rhs


    subroutine compute_CS_tstep(u,CS)
        real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,Nvar)
        real(dp), intent(out) :: CS

        CS = maxval(abs(u(1:Nx,:)))

        CS = max(CS, 1.0e-12_dp) ! Make sure CS doesn't return zero

    end subroutine compute_CS_tstep

end module burgers_mod
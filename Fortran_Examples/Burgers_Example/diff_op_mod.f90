module diff_op_mod
    use real_type_mod, only: dp
    use read_parameters_mod
    implicit none


    contains 


    subroutine linear_reconstruct(u,u_plus,u_minus)
        real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,Nvar)
        real(dp), intent(out) :: u_plus(1-Ngz:Nx+Ngz,Nvar), u_minus(1-Ngz:Nx+Ngz,Nvar)
        real(dp) :: ratio(1-Ngz:Nx+Ngz,Nvar), phi(1-Ngz:Nx+Ngz,Nvar)


        ! Initialise arrays to zero
        ratio = 0.0_dp
        phi = 0.0_dp
        u_plus = 0.0_dp
        u_minus = 0.0_dp

        ! Compute ratio of slopes for limiter
        ! NB: We add small number to the denominator to 
        ! stop NaN issues when neighbouring points are close to equal
        ratio(0:Nx+1,:) = (u(0:Nx+1,:) - u(-1:Nx,:))/(u(1:Nx+2,:) - u(0:Nx+1,:) + 1e-16_dp) 

        ! Minmod limiter
        phi = max(0.0_dp,min(1.0_dp,ratio))

        ! Compute the slope-limited reconstruction
        u_minus(-1:Nx+1,:) = u(-1:Nx+1,:) + 0.5_dp*phi(-1:Nx+1,:)*(u(0:Nx+2,:) - u(-1:Nx+1,:))

        u_plus(-1:Nx,:) = u(0:Nx+1,:) - 0.5_dp*phi(0:Nx+1,:)*(u(1:Nx+2,:) - u(0:Nx+1,:))

    end subroutine linear_reconstruct

    subroutine local_lax_friedrichs(flux_plus,flux_minus,cons_plus,cons_minus,char_speed,flux)
        real(dp), intent(in) :: flux_plus(1-Ngz:Nx+Ngz,Nvar), flux_minus(1-Ngz:Nx+Ngz,Nvar)
        real(dp), intent(in) :: cons_plus(1-Ngz:Nx+Ngz,Nvar), cons_minus(1-Ngz:Nx+Ngz,Nvar)
        real(dp), intent(in) :: char_speed(1-Ngz:Nx+Ngz,1)
        real(dp), intent(out) :: flux(1-Ngz:Nx+Ngz,Nvar)


        flux = 0.5*( (flux_plus + flux_minus)  - char_speed*(cons_plus  - cons_minus)) 

    end subroutine local_lax_friedrichs


end module diff_op_mod
! Module containing the time stepping routines
module time_stepping_mod
    use real_type_mod, only: dp
    use burgers_mod, only: compute_rhs
    use boundary_conditions_mod, only: periodic_bc
    use read_parameters_mod, only: Nx, Ngz, Nvar


    contains
    
    subroutine rk4(t0,y0,dt)
        real(dp), intent(in) :: t0,dt ! Current time, timestep 
        real(dp), intent(inout) :: y0(1-Ngz:Nx+Ngz,Nvar) ! Current values
        real(dp) :: k1(1-Ngz:Nx+Ngz,Nvar), k2(1-Ngz:Nx+Ngz,Nvar), k3(1-Ngz:Nx+Ngz,Nvar), k4(1-Ngz:Nx+Ngz,Nvar) ! Runge-Kutta Substeps
        real(dp) :: y1(1-Ngz:Nx+Ngz,Nvar), y2(1-Ngz:Nx+Ngz,Nvar), y3(1-Ngz:Nx+Ngz,Nvar) 

        ! Step 1
        call compute_rhs(t0,y0,k1)
        y1 = y0 + 0.5*dt*k1
        call periodic_bc(Nx,Ngz,Nvar,y1) ! Update Boundary Conditions
        
        ! Step 2
        call compute_rhs(t0+0.5_dp*dt,y1,k2)
        y2 = y0 + 0.5*dt*k2
        call periodic_bc(Nx,Ngz,Nvar,y2)

        ! Step 3
        call compute_rhs(t0+0.5_dp*dt,y2,k3)
        y3 = y0 + dt*k3
        call periodic_bc(Nx,Ngz,Nvar,y3)

        ! Step 4
        call compute_rhs(t0+dt,y3,k4)

        
        ! Update y
        y0 = y0 + dt*((1.0_dp/6.0_dp)*k1 + (1.0_dp/3.0_dp)*k2 + (1.0_dp/3.0_dp)*k3  + (1.0_dp/6.0_dp)*k4)
        call periodic_bc(Nx,Ngz,Nvar,y0)


    end subroutine rk4

end module time_stepping_mod
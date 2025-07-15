! Module containing the time stepping routines
module time_stepping
    use real_type_mod, only: dp
    use advection, only: compute_rhs
    use boundary_conditions, only: periodic_bc


    contains
    
    subroutine rk4(Nx,Ngz,dx,t0,y0,dt)
        integer, intent(in) :: Nx,Ngz ! Num of eqns to solve
        real(dp), intent(in) :: t0,dt, dx ! Current time, timestep, spatial step size
        real(dp), intent(inout) :: y0(1-Ngz:Nx+Ngz,1) ! Current values
        real(dp) :: k1(1-Ngz:Nx+Ngz,1), k2(1-Ngz:Nx+Ngz,1), k3(1-Ngz:Nx+Ngz,1), k4(1-Ngz:Nx+Ngz,1) ! Runge-Kutta Substeps
        real(dp) :: y1(1-Ngz:Nx+Ngz,1), y2(1-Ngz:Nx+Ngz,1), y3(1-Ngz:Nx+Ngz,1) 

        ! Step 1
        call compute_rhs(Nx,Ngz,dx,t0,y0,k1)
        y1 = y0 + 0.5*dt*k1
        call periodic_bc(Nx,Ngz,y1) ! Update Boundary Conditions
        
        ! Step 2
        call compute_rhs(Nx,Ngz,dx,t0+0.5_dp*dt,y1,k2)
        y2 = y0 + 0.5*dt*k2
        call periodic_bc(Nx,Ngz,y2)

        ! Step 3
        call compute_rhs(Nx,Ngz,dx,t0+0.5_dp*dt,y2,k3)
        y3 = y0 + dt*k3
        call periodic_bc(Nx,Ngz,y3)

        ! Step 4
        call compute_rhs(Nx,Ngz,dx,t0+dt,y3,k4)
        
        ! Update y0
        y0 = y0 + dt*((1.0_dp/6.0_dp)*k1 + (1.0_dp/3.0_dp)*k2 + (1.0_dp/3.0_dp)*k3  + (1.0_dp/6.0_dp)*k4)
        call periodic_bc(Nx,Ngz,y0)


    end subroutine rk4

end module time_stepping
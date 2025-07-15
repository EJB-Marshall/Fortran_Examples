module runge_kutta_2
    use real_type_mod, only: dp
    implicit none

    ! First we construct an abstract interface
    ! This allows us to pass arbitrary ODEs to 
    ! the Runge-Kutta method
    abstract interface
        subroutine rhs_interface(t, y, neqn, dydt)
            import :: dp
            integer, intent(in) :: neqn ! The number of eqns to solve
            real(dp), intent(in) :: t, y(neqn)
            real(dp), intent(out):: dydt(neqn)
        end subroutine rhs_interface
    end interface
    
    contains

    ! 2nd order accurate Runge-Kutta method (Heun's method)
    subroutine rk2(t0,y0,dt,neqn,rhs)
        integer :: neqn ! Num of eqns to solve
        real(dp), intent(in) :: t0,dt ! Current time, timestep
        real(dp), intent(inout) :: y0(neqn) ! Current values
        real(dp) :: k1(neqn), k2(neqn) ! Runge-Kutta Substeps
        procedure(rhs_interface) :: rhs

        call rhs(t0,y0,neqn,k1)
        call rhs(t0 + 0.5_dp*dt, y0 + dt*k1, neqn, k2)

        y0 = y0 + dt*(0.5_dp*k1 + 0.5_dp*k2)

    end subroutine rk2

end module runge_kutta_2
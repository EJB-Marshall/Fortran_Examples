module compute_area
    implicit none
    real (kind=8), parameter :: pi = acos(-1.d0) 
    public pi, area
    private

    contains
    
    real (kind = 8) function area(r)
        real (kind=8), intent(in) :: r
        
        area = pi*r**2
        
    end function area

end module compute_area
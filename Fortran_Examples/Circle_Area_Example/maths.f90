program maths
    use compute_area, only: area, pi
    implicit none
    real (kind = 8) :: r

    r = 2.0

    print *, "Pi = ", pi
    print *, "The area of the circle is ", area(r)

end program maths
real function degree_to_radian(theta)
    implicit none 
    real, intent(in) :: theta
    real, parameter :: PI = 3.14159

    degree_to_radian = theta * PI / 180.
end function

subroutine radian_to_degree(radian, degree)
    implicit none 
    real, intent(in) :: radian
    real, intent(out) :: degree 
    real, parameter :: PI = 3.14159

    degree = radian * 180. / PI
end subroutine radian_to_degree

program ex11
    implicit none 
    integer, parameter :: N = 6
    real degree_to_radian
    integer :: i
    real, dimension(N+1) ::             &
        angles = [ (i * 15., i=0,N)],   &
        radian_ouput,                   &
        degree_ouput = 0.            

    do i=1,N+1
        radian_ouput(i) = degree_to_radian(angles(i))
    enddo

    print *, 'Radians'
    write(*, '(f8.2)') radian_ouput

    do i=1,N+1
        call radian_to_degree(radian_ouput(i), degree_ouput(i))
    enddo

    print *, 'Back to degrees'
    write(*, '(f8.2)') degree_ouput

    do i=1,N+1
        call radian_to_degree(degree_to_radian(angles(i)), degree_ouput(i))
    enddo

    print '(A)', 'Radians and then back to degrees'
    print '(A)', 'Is it the same?'
    write(*, '(f8.2)') degree_ouput

end program


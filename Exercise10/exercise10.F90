program ex10
    implicit none
    integer, parameter :: NUM_DRAWS=1e7
    integer :: circle_draw_counts=0, i 
    real( kind=4 ):: random_position(2), PI=3.1415926, error = 1., delta=1e-7
    integer :: seed_array(16)

    call date_and_time(values=seed_array(1:8))
    call date_and_time(values=seed_array(9:16))
    ! print *, seed_array

    call random_seed(put=seed_array)

    mc_pi_loop: do
        if (( i > NUM_DRAWS ) .or. (abs(error) <= delta)) then 
            exit mc_pi_loop
        else
            call random_number(random_position)
            if ( (random_position(1)**2 + random_position(2)**2) < 1.0) then
                circle_draw_counts = circle_draw_counts + 1
            end if
            i = i + 1
            error = (4.0*(real(circle_draw_counts) / real(NUM_DRAWS))) - PI
        end if
    end do mc_pi_loop

    print *, "estimated pi =", &
        4.0*(real(circle_draw_counts) / real(NUM_DRAWS))
    
end program ex10

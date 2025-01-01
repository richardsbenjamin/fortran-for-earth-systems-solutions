program timestep_filename_construction
    implicit none
    character(40) :: aux_string, width_string
    integer :: i, num_time_steps = 10, speed_file_id, width

    width = aint(log10(real(num_time_steps))) + 1
    do i=1, num_time_steps
        write(width_string, 9) width
9        format (I4)
        write(aux_string, '(i0.'//width_string//')') i

        print *, aux_string

        open(newunit=speed_file_id, &
                file="speed_"// trim(adjustl(aux_string)) //".dat", &
                action="write")

        close(speed_file_id)

    end do
end program timestep_filename_construction
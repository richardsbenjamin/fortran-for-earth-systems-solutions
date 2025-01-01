program ex7
    implicit none
    integer, parameter :: int_exp_range = 45, real_exp_range = 5500, prec_limit = 60
    integer :: large_int
    real :: large_real
    integer :: i, j, k
    integer :: int_exp_file_id, real_exp_file_id

    open(newunit=int_exp_file_id, &
            file="int_exp_ranges.txt", &
            action="write")
    open(newunit=real_exp_file_id, &
            file="real_exp_ranges.txt", &
            action="write")

    do i = 0, int_exp_range
        large_int = selected_int_kind(i)
        write(int_exp_file_id, '(i0, a i0)') i, ", ", large_int
    end do

    do j = 0, real_exp_range
        do k = 0, prec_limit
            large_real = selected_real_kind(k, j)
            write(real_exp_file_id, '(i0, a, i0, a, f0.2)') &
                    k, ", ", j, ", ", large_real
        end do
    end do

end program ex7

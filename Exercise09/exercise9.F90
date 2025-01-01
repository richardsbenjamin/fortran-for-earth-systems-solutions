program ex9
    implicit none
    integer, parameter :: N = 813, N_rep = 30
    integer :: iter, i, j, k, stat_code
    real( kind=4 ), dimension(:,:,:), allocatable :: a, b

    allocate( a(N,N,N), b(N,N,N), stat=stat_code)

    a = 1
    b = 1

    do iter=1,N_rep
        do k=1,N
            do j=1,N
                do i=1,N
                    a(i,j,k) = a(i,j,k) + b(i,j,k)
                enddo
            enddo
        enddo
    enddo

end program ex9

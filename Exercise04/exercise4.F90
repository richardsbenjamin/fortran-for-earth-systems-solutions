program ex4
    implicit none

    integer :: i, numTimesteps = 80, speedFileID
    real :: Si, pi = 3.1459, radianN, radianS
    character(len=256) :: buffer

    open(newunit=speedFileID, file="coordinates", action="write")

    do i = 0, numTimesteps, 10
        radianS = i * pi / 180.
        radianN = (i + 10.) * pi / 180
        Si = 10. * (sin(radianN) - sin(radianS))
        write (speedFileID, '(f10.6)') Si
    end do

    close(speedFileID)

    open(newunit=speedFileID, file="coordinates", action="read")

    do i = 0, numTimesteps, 10
        read(speedFileID, *) buffer
        print *, buffer
    end do    

end program ex4

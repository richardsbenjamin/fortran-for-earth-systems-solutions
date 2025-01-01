program ex5
    implicit none

    real :: sigma, G, arc
    integer, parameter :: Nx = 100, Ny = 1000, Nz = 4000
    integer, parameter :: X = 100, L = 1000, H = 4
    real, parameter :: dx = 1., dy = 1., dz = 0.001
    integer :: ix, iy, iz
    real :: y, z, total_volume, prop_vol

    prop_vol = 0

    do ix = 0, Nx
        do iy = 0, Ny
            y = iy
            do iz = 0, Nz
                z = iz * dz

                G = (2.0 - y / H)**2 + (0.1 - z / H)**2

                sigma = 0.9184 * (-sqrt(G) + 1.0) ** 2 &
                         + (0.9184 * acos((1.0 / sqrt(G)) * (2.0 - y / H)) ** 2) &
                         + 26.57

                total_volume = total_volume + (dx * dy * dz)

                if (27.68 <= sigma .and. sigma <= 27.74) then
                    prop_vol = prop_vol + (dx * dy * dz)

                end if
            end do
        end do
    end do

    print *, prop_vol, total_volume, prop_vol / total_volume

end program ex5


! command line : mpifort -I sprng5/include/ main.F90 -L sprng5/lib/ -lsprng -lstdc++

module constants
    implicit none
    doubleprecision, parameter :: pi = dacos(-1d0), epsilon = 1d-6
    doubleprecision, parameter :: d = 1d0, d_half = 0.5d0
    doubleprecision, parameter :: two_pi = 2d0 * pi, drho = 0.01d0
end module constants

module particle
    implicit none
    doubleprecision :: x, y, z, r_wos
    integer :: hit
end module particle

module container
    implicit none
    integer, parameter :: n_bin = 300
    integer up(n_bin), down(n_bin), total(n_bin), ind
end module container

!****************************************************************
program main
    use constants
    use particle
    use container
    implicit none

#define SIMPLE_SPRNG
#include "sprng_f.h"
    SPRNG_POINTER stream

    character(len=20) :: name_file
    integer :: i, n_total, n_file
    integer :: t1, t2, clock_rate, clock_max
    doubleprecision :: rho_middle, area
    write(*, *) "Give the name of csv file"
    read(*, *) n_file
    write(name_file, "('csv/IPPWOS',i0,'.csv')") n_file

    write(*, *) "Enter the number of particle"
    read(*, *) n_total

    stream = init_sprng(12345, SPRNG_DEFAULT, 5)

    call system_clock(t1, clock_rate, clock_max)
    do i = 1, n_total
        call firstpassage
    end do
    do ind = 1, n_bin
        total(ind) = up(ind) + down(ind)
    end do
    call system_clock(t2, clock_rate, clock_max)
    write(*, *) "Running time: ", real(t2-t1)/real(clock_rate)

    open(unit = 1, file= name_file, status="new")
    write(1, *) "Total number: ", n_total
    write(1, *) "Running time: ", real(t2-t1)/real(clock_rate)

    do ind = 1, n_bin
        write(1, *) total(ind)
    end do
    close(1)

end program main

!*********************************
subroutine wos
    use particle
    use constants
    implicit none

#define SIMPLE_SPRNG
#include "sprng_f.h"
    SPRNG_POINTER stream

    doubleprecision xxx, yyy
    doubleprecision theta, phi
    xxx = sprng(stream)
    yyy = sprng(stream)
    theta = dacos(2d0 * xxx - 1d0)
    phi = two_pi * yyy
    x = x + r_wos * dsin(theta) * dcos(phi)
    y = y + r_wos * dsin(theta) * dsin(phi)
    z = z + r_wos * dcos(theta)
end subroutine wos

subroutine initialize
    use particle
    use constants
    implicit none

    x = 0d0
    y = 0d0
    z = d_half
    r_wos = d_half
    hit = 0
end subroutine initialize

subroutine firstpassage
    use particle
    use constants
    use container
    implicit none

    call initialize
    do while (hit == 0)
        call wos
        r_wos = min(z, d-z)
        if (r_wos < epsilon) then
            hit = 1
            ind = dsqrt(x**2+y**2)/drho + 1
            if (ind < n_bin) then
                if (z > d_half) then
                    up(ind) = up(ind) + 1
                else
                    down(ind) = down(ind) + 1
                end if
            end if
        end if
    end do
end subroutine firstpassage

! command line : mpifort -I sprng5/include/ IPPWOH.F90 -o IPPWOH -L sprng5/lib/ -lsprng -lstdc++

module constants
    implicit none
    doubleprecision, parameter :: pi = dacos(-1d0)
    doubleprecision, parameter :: d = 1d0, d_half = 0.5d0
    doubleprecision, parameter :: two_pi = 2d0 * dacos(-1d0), drho = 0.01d0
end module constants

module particle
    implicit none
    doubleprecision :: x, y, z
    integer :: hit
end module particle

module container
    implicit none
    integer, parameter :: n_bin = 300
    integer :: up(n_bin), down(n_bin), total(n_bin), ind
end module container

module woh_variables
    implicit none
    doubleprecision :: alpha, beta, gamma, sigma, theta
end module woh_variables

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

    write(*, *) "Give the name of csv file"
    read(*, *) n_file
    write(name_file, "('csv/IPPWOH',i0,'.csv')") n_file

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

    open(unit = 1, file=name_file, status="new")
    write(1, *) "Total number: ", n_total
    write(1, *) "Running time: ", real(t2-t1)/real(clock_rate)

    do ind = 1, n_bin
        write(1, *) total(ind)
    end do
    close(1)

end program main

!*********************************
subroutine woh
    use particle
    use constants
    use woh_variables
    implicit none

#define SIMPLE_SPRNG
#include "sprng_f.h"
    SPRNG_POINTER stream

    doubleprecision :: rnd, rn1, rn2, phi, temp, tempr

    rnd = sprng(stream)
    tempr = nint(rnd)

    z = (d-z) * tempr + z * (1d0 - tempr)
    rn1= sprng(stream)
    phi = two_pi * rn1

    gamma = z/d
    alpha = gamma ** 2d0 + 1d0
    beta =  gamma ** 2d0 - 1d0
    sigma = (beta + dsqrt(alpha))/(gamma * dsqrt(alpha))

    rn2= sprng(stream)

    if (rn2 < sigma) then
        call hemisphere_angle
        x = x + d * sin(theta) * cos(phi)
        y = y + d * sin(theta) * sin(phi)
        z = d * cos(theta)
        z = z * tempr +(d-z) * (1d0 - tempr)
    else
        gamma = (1d0-gamma)/(1d0+gamma)
        alpha = gamma ** 2d0 + 1d0
        beta =  gamma ** 2d0 - 1d0
        sigma = (beta + dsqrt(alpha))/(gamma * dsqrt(alpha))
        call disk_angle
        x = x + tan(theta*0.5d0) * dcos(phi)
        y = y + tan(theta*0.5d0) * dsin(phi)
        z = 0d0
        z = z * tempr +(d-z) * (1d0 - tempr)
    end if
end subroutine woh

subroutine hemisphere_angle
    use constants
    use woh_variables
    implicit none

#define SIMPLE_SPRNG
#include "sprng_f.h"
    SPRNG_POINTER stream

    doubleprecision :: rnd, a, b, costh
    rnd = sprng(stream)
    a = (2d0*(-1d0+gamma*sigma*rnd)/beta)**2d0
    b = -2d0*alpha/a +alpha**2d0 -2d0*(1d0+dsqrt(1d0+2d0*a*alpha))/a**2d0
    costh = dsqrt(b/(4d0*gamma**2))
    theta = dacos(costh)
end subroutine hemisphere_angle

subroutine disk_angle
    use particle
    use constants
    use woh_variables
    implicit none

#define SIMPLE_SPRNG
#include "sprng_f.h"
    SPRNG_POINTER stream
    integer m
    doubleprecision :: rnd, a, b, c, rej, test, costh
    do while (.true.)
        rnd = sprng(stream)
        a = (2d0*(-1d0+gamma*sigma*rnd)/beta)**2d0
        b = -2d0*alpha/a +alpha**2d0 -2d0*(1d0+dsqrt(1d0+2d0*a*alpha))/a**2d0
        costh = dsqrt(b/(4d0*gamma**2))
        theta = dacos(costh)
        rej = dsqrt(0.5d0)/dcos(theta*0.5d0)
        test = sprng(stream)
        if (test < rej) then
            exit
        end if
    end do
end subroutine disk_angle



subroutine initialize
    use particle
    use constants
    implicit none

    x = 0d0
    y = 0d0
    z = d_half
    hit = 0
end subroutine initialize

subroutine firstpassage
    use particle
    use constants
    use container
    implicit none

    call initialize
    do while (hit == 0)
        call woh
        if (z==0d0 .or. z==d) then
            hit = 1
            ind = dsqrt(x**2d0+y**2d0)/drho + 1
            if (ind < n_bin) then
                if (z == 0) then
                    down(ind) = down(ind) + 1
                else
                    up(ind) = up(ind) + 1
                end if
            end if
        end if
    end do
end subroutine firstpassage

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

    call system_clock(t1, clock_rate, clock_max)
    do i = 1, n_total
        call firstpassage
    end do
!
!    do ind = 1, n_bin
!        total(ind) = up(ind) + down(ind)
!    end do

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

    doubleprecision :: rn1, rn2, phi

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
    else
        gamma = (1d0-gamma)/(1d0+gamma)
        alpha = gamma ** 2d0 + 1d0
        beta =  gamma ** 2d0 - 1d0
        sigma = (beta + dsqrt(alpha))/(gamma * dsqrt(alpha))

        call disk_angle
        x = x + tan(theta*0.5d0) * dcos(phi)
        y = y + tan(theta*0.5d0) * dsin(phi)
        z = 0d0
    end if
end subroutine woh

subroutine hemisphere_angle
    use constants
    use woh_variables
    implicit none

#define SIMPLE_SPRNG
#include "sprng_f.h"
    SPRNG_POINTER stream

    doubleprecision :: rnd, xxx, yyy
    rnd = sprng(stream)
    xxx = 4d0*(-1d0+gamma*sigma*rnd)**2d0/beta**2d0
    yyy = -2d0*(1d0+xxx*alpha-0.5d0*(alpha*xxx)**2d0+dsqrt(1d0+2d0*xxx*alpha))
    theta = dacos(0.5d0/gamma/xxx*dsqrt(yyy))
end subroutine hemisphere_angle

subroutine disk_angle
    use constants
    use woh_variables
    implicit none

#define SIMPLE_SPRNG
#include "sprng_f.h"
    SPRNG_POINTER stream

    doubleprecision :: xxx, yyy, rnd, rej
    do while (.true.)
        rnd = sprng(stream)
        xxx = 4d0*(-1d0+gamma*sigma*rnd)**2d0/beta**2d0
        yyy = -2d0*(1d0+xxx*alpha-0.5d0*(alpha*xxx)**2d0+dsqrt(1d0+2d0*xxx*alpha))
        theta = dacos(0.5d0/gamma/xxx*dsqrt(yyy))
        rej = sprng(stream)
        if (rej < sqrt(0.5d0)/cos(theta*0.5d0)) then
            exit
        end if
    end do
end subroutine disk_angle

subroutine initialize_z
    use particle, only:z
    use constants, only:d
#define SIMPLE_SPRNG
#include "sprng_f.h"
    SPRNG_POINTER stream
    test = sprng(stream)
    if (test > 0.5d0) then
        z = d - z
    end if
endsubroutine initialize_z


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
    doubleprecision rnd

    call initialize
    do while (hit == 0)
        z = min(z, d-z)
        call woh
        if (z == 0d0) then
            hit = 1
            ind = dsqrt(x**2d0+y**2d0)/drho + 1
            total(ind) = total(ind) + 1
        end if
    end do
end subroutine firstpassage

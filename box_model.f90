program box_model
use micro
use thermo
use solver
use util
implicit none

    integer, parameter :: ntime = 3600, nbin = 51
    double precision, parameter :: dlm = 0.1d0, dt = 1.d0, rhoi = 920.d0, pi = 3.14159265d0
    integer :: k, i
    double precision, dimension(nbin) :: mass, ch_mass, ncon, reqv, dmass, a, c, rhoeff, vol_fac
    double precision, dimension(nbin+1) :: mass_e
    double precision :: alp, bet, ntot, ncon_fac
    double precision :: tmpk, pres, si, vapi, vapl

    ! set intial psd of all spherical ice particles
    alp = 18.d0
    bet = 9.d4
    ntot = 6.d3
    mass_e(1) = 10.d0**(-13.d0)
    !ncon = 1.d0
    !ncon(1) = ntot

    do i = 1, nbin
        ! set mass bins
        mass_e(i+1) = 10.d0**(i*dlm-12.d0)
        mass(i) = 10.d0**(0.5d0*(dlog10(mass_e(i))+dlog10(mass_e(i+1))))
        dmass(i) = mass_e(i+1)-mass_e(i)

        reqv(i) = (3.d0*mass(i)/(4.d0*rhoi*pi))**(1.d0/3.d0)
        !ncon(i) = ntot*bet**alp*reqv(i)**(alp-3.d0)/dgamma(alp)*dexp(-bet*reqv(i))/(4.d0*pi*rhoi)
        ncon(i) = dexp(-bet*reqv(i))

        ! set volumes and axis lengths (all spherical initially)
        a(i) = reqv(i)
        c(i) = reqv(i)
        print *, i, a(i), mass(i), ncon(i)
    end do

    !ncon(1) = 100.d0

    ncon = ntot*ncon/(0.5d0*sum(ncon(2:nbin)*dmass(2:nbin)+ncon(1:nbin-1)*dmass(1:nbin-1)))

    ! set environment
    tmpk = 258.15d0
    pres = 80000.d0
    vapi = svp_ice(tmpk)
    vapl = svp_liq(tmpk)
    si = vapl/vapi-1.d0
    !si = 0.15d0
    print *, si

    ! time loop
    open(unit=11, file="ndist.txt")
    open(unit=21, file="a.txt")
    open(unit=31, file="c.txt")
    open(unit=41, file="mass.txt")
    write(11,*) ncon
    write(21,*) a
    write(31,*) c
    write(41,*) mass

    do k = 1, ntime
        ! call vapor growth subroutine
        call vapor_vec(tmpk, si, pres, dt, nbin, a, c, ch_mass)
        rhoeff = (mass+ch_mass)/(4.d0/3.d0*pi*a**2.d0*c)

        ! interpolate particle properties to old mass grid
        call interp_linear_pw(mass, mass+ch_mass, a, a, nbin, nbin, mass_e(1), mass_e(nbin+1),&
                              (3.d0*mass_e(1)/(4.d0*rhoi*pi))**(1.d0/3.d0), a(nbin))
        call interp_linear_pw(mass, mass+ch_mass, c, c, nbin, nbin, mass_e(1), mass_e(nbin+1),&
                              (3.d0*mass_e(1)/(4.d0*rhoi*pi))**(1.d0/3.d0), c(nbin))
        call interp_linear_pw(mass, mass+ch_mass, ncon, ncon, nbin, nbin, mass_e(1), mass_e(nbin+1),&
                              1.d-3, ncon(nbin))
        call interp_linear_pw(mass, mass+ch_mass, rhoeff, rhoeff, nbin, nbin, mass_e(1), mass_e(nbin+1),&
                              920.d0, rhoeff(nbin))

        ! adjust by effective density
        vol_fac = rhoeff/(mass/(4.d0/3.d0*pi*a**2.d0*c))
        a = a/(vol_fac)**(1.d0/3.d0)
        c = c/(vol_fac)**(1.d0/3.d0)

        ! adjust total number
        !ncon = ntot*ncon/(0.5d0*sum(ncon(2:nbin)*dmass(2:nbin)+ncon(1:nbin-1)*dmass(1:nbin-1)))

        ! advect properties between mass bins
        !call tridiagonal(-ch_mass(1:nbin-1)/dmass(1:nbin-1),&
        !                 &1.d0+dabs(ch_mass)/dmass, ch_mass(2:nbin)/dmass(2:nbin),&
        !                 &ncon, nbin, ncon_new)

        !ncon_new(2:nbin) = (ncon(2:nbin)+ch_mass(1:nbin-1)/dmass(1:nbin-1)*ncon(1:nbin-1))/&
        !                   &(1.d0+ch_mass(1:nbin-1)/dmass(1:nbin-1))
        !ncon_new(1) = ncon(1)
        !va_new(2:nbin) = (va(2:nbin)+ch_mass(1:nbin-1)/dmass(1:nbin-1)*va(1:nbin-1))/&
        !                   &(1.d0+ch_mass(1:nbin-1)/dmass(1:nbin-1))
        !va_new(1) = a(1)
        !vc_new(2:nbin) = (vc(2:nbin)+ch_mass(1:nbin-1)/dmass(1:nbin-1)*vc(1:nbin-1))/&
        !                   &(1.d0+ch_mass(1:nbin-1)/dmass(1:nbin-1))
        !vc_new(1) = c(1)


        !ncon_fac = 0.5d0*sum(ncon_new(2:nbin)*dmass(2:nbin)+ncon_new(1:nbin-1)*dmass(1:nbin-1))/&
        !       (0.5d0*sum(ncon(2:nbin)*dmass(2:nbin)+ncon(1:nbin-1)*dmass(1:nbin-1)))
        !ncon = ncon_new/ncon_fac

        !va_fac = 0.5d0*sum(va_new(2:nbin)*dmass(2:nbin)+va_new(1:nbin-1)*dmass(1:nbin-1))/&
        !       (0.5d0*sum(va(2:nbin)*dmass(2:nbin)+va(1:nbin-1)*dmass(1:nbin-1)))
        !vc_fac = 0.5d0*sum(vc_new(2:nbin)*dmass(2:nbin)+vc_new(1:nbin-1)*dmass(1:nbin-1))/&
        !       (0.5d0*sum(vc(2:nbin)*dmass(2:nbin)+vc(1:nbin-1)*dmass(1:nbin-1)))
        !va = va_new
        !vc = vc_new

        !print *, (ncon(2)+ch_mass(1)/dmass(1)*ncon(1))/(1.d0+ch_mass(1)/dmass(1)), ncon(2), ncon(1)
        !ncon(2:nbin) = (ncon(2:nbin)+ch_mass(1:nbin-1)/dmass(1:nbin-1)*ncon(1:nbin-1))/&
        !               &(1.d0+ch_mass(1:nbin-1)/dmass(1:nbin-1))
        !ncon(1) = ncon(2)
        !call tridiagonal(-ch_mass(1:nbin-1)/dmass(1:nbin-1),&
        !                 &1.d0+dabs(ch_mass)/dmass, ch_mass(2:nbin)/dmass(2:nbin),&
        !                 &a, nbin, a)
        !call tridiagonal(-ch_mass(1:nbin-1)/dmass(1:nbin-1),&
        !                 &1.d0+dabs(ch_mass)/dmass, ch_mass(2:nbin)/dmass(2:nbin),&
        !                 &c, nbin, c)

        ! recalculate a and c in each bin
        !a = (va/ncon)**(1.d0/3.d0)
        !c = (vc/ncon)**(1.d0/3.d0)

        !do i = 1, nbin
            !print *, i, mass(i), 4.d0/3.d0*rhoeff(i)*pi*a(i)**2.d0*c(i), rhoeff(i)
        !    print *, i, a(i), a(i)/c(i), rhoeff(i)
        !    print *, i, a(i), mass(i), ch_mass(i)/dmass(i), ncon(i)
            !print *, ch_mass(i)
        !end do
        write(11,*) ncon
        write(21,*) a
        write(31,*) c
        !write(41,*) mass_new
    end do

    close(11)
    close(21)
    close(31)
    close(41)

end program

module micro
use thermo
implicit none
contains


    ! growth of a spheroid over one time step
    subroutine spheroid_vap(tmpk, si, pres, dt, a, c)

        double precision :: tmpk, si, pres, dt, a, c
        double precision :: dmdt, rhodep, dvdt, dphidt, vol, phi, gam, reqv, pi
        double precision :: k1v, k2v, k3v, k4v, k1p, k2p, k3p, k4p
        double precision :: v1, v2, v3, v4, p1, p2, p3, p4

        ! common parameters
        pi = 3.14159265d0
        gam = 0.4d0
        rhodep = rhodepbranched(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0)
        print *, rhodep
        !rhodep = 920.d0
        vol = 4.d0/3.d0*pi*a**2.d0*c
        phi = c/a
        mass_i = vol*rhodep

        ! mass growth
        !dmdt = 4.d0*pi*capac(a,c)*si*g_diff(tmpk, pres, si, (3.d0*vol/(4.d0*pi))**(1.d0/3.d0))
        !dvdt = dmdt/rhodep
        !vol = vol+dvdt*dt

        ! change in shape
        !dphidt = phi/vol*(gam-1.d0)/(gam+2.d0)*dvdt
        !dphidt = 1.d0/vol*(gam-phi)/(gam+2.d0*phi)*dvdt
        !phi = phi+dphidt*dt
        !vol = vol+dvdt*dt

        ! use runge-kutta method
        !-------------------------------
        ! step 1
        k1v = 4.d0*pi*capac(a,c)*si*g_diff(tmpk, pres, si, (3.d0*vol/(4.d0*pi))**(1.d0/3.d0))/&
              &rhodep
        k1p = phi/vol*(gam-1.d0)/(gam+2.d0)*k1v

        ! step 2
        v1 = vol+dt*k1v/2.d0
        p1 = phi+dt*k1p/2.d0
        a = (3.d0*v1/(4.d0*pi*p1))**(1.d0/3.d0)
        c = a*p1
        rhodep = rhodepbranched(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0)

        k2v = 4.d0*pi*capac(a,c)*si*g_diff(tmpk, pres, si, (3.d0*v1/(4.d0*pi))**(1.d0/3.d0))/&
              &rhodep
        k2p = p1/v1*(gam-1.d0)/(gam+2.d0)*k2v

        ! step 3
        v2 = vol+dt*k2v/2.d0
        p2 = phi+dt*k2p/2.d0
        a = (3.d0*v2/(4.d0*pi*p2))**(1.d0/3.d0)
        c = a*p2
        rhodep = rhodepbranched(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0)

        k3v = 4.d0*pi*capac(a,c)*si*g_diff(tmpk, pres, si, (3.d0*v2/(4.d0*pi))**(1.d0/3.d0))/&
              &rhodep
        k3p = p2/v2*(gam-1.d0)/(gam+2.d0)*k2v

        ! step 4
        v3 = vol+dt*k3v
        p3 = phi+dt*k3p
        a = (3.d0*v3/(4.d0*pi*p3))**(1.d0/3.d0)
        c = a*p3
        rhodep = rhodepbranched(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0)

        k4v = 4.d0*pi*capac(a,c)*si*g_diff(tmpk, pres, si, (3.d0*v3/(4.d0*pi))**(1.d0/3.d0))/&
              &rhodep
        k4p = p3/v3*(gam-1.d0)/(gam+2.d0)*k3v


        vol = vol+1.d0/6.d0*dt*(k1v+2.d0*k2v+2.d0*k3v+k4v)
        phi = phi+1.d0/6.d0*dt*(k1p+2.d0*k2p+2.d0*k3p+k4p)

        ! update sizes
        a = (3.d0*vol/(4.d0*pi*phi))**(1.d0/3.d0)
        c = a*phi

    end subroutine

    function capac(a, c)

        double precision :: capac, a, c, ecc

        ! oblate spheroid
        if (a.gt.c) then
            ecc = dsqrt(1.d0-(c/a)**2.d0)
            capac = a*ecc/(dasin(ecc))

        ! prolate spheroid
        else if (a.lt.c) then
            ecc = dsqrt(1.d0-(a/c)**2.d0)
            capac = c*ecc/(dlog(c/a*(1.d0+ecc)))

        ! sphere
        else
            capac = a
        end if
    
    end function

    function rhodepbranched(a,ac,fb,ft,fi,rhoi)

        double precision :: a,ac,asc,fb,rhoi,rhodepbranched
        double precision :: fi,ft,amax
        double precision :: wt,wsb,wmb,fmb,ai
        integer nsb

        ! set density parameters
        amax = 3.d0
        nsb = 5

        ! calculate fmb and ai
        wt = ft/2.d0*amax
        wsb = 1.d0/(dble(nsb-1)/fb+1.d0)*(amax-ac)
        wmb = min(max(wsb/2.d0, ac/2.d0), min(wsb/2.d0, ac/2.d0))
        fmb = wmb/(ac/2.d0)
        ai = fi*amax+(1.d0-fi)*ac

        ! calculate deposition density (scaled linear fraction x rhoi)
        if (a.le.ac) then
            rhodepbranched = rhoi
        elseif ((a.gt.ac).and.(a.le.ai)) then
            rhodepbranched = rhoi*(fb+ac/a*fmb*(1.d0-fb))
        else
            rhodepbranched = rhoi*(fb/(amax-ai)*(ft*amax-ai+ai*amax*(1.d0-ft)/a)+&
                                  (1.d0-fb)*fmb*ac/a)
        endif

      end function

end module

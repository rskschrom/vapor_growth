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

    ! vectorized vapor growth of an array of spheroids
    subroutine vapor_vec(tmpk, si, pres, dt, npar, a, c, dmass)

        integer :: npar
        double precision :: tmpk, si, pres, dt
        double precision, dimension(npar) :: a, c, dmass
        double precision, dimension(npar) :: dmdt, rhodep, dvdt, dphidt, vol, phi, gam, reqv, pi
        double precision, dimension(npar) :: k1v, k2v, k3v, k4v, k1p, k2p, k3p, k4p
        double precision, dimension(npar) :: v1, v2, v3, v4, p1, p2, p3, p4, voli

        ! common parameters
        pi = 3.14159265d0
        gam = 0.28d0
        rhodep = rhodep_vec(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0,npar)
        !rhodep = 920.d0
        vol = 4.d0/3.d0*pi*a**2.d0*c
        voli = vol
        phi = c/a

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
        k1v = 4.d0*pi*capac_vec(a,c,npar)*si*g_vec(tmpk,pres,si,(3.d0*vol/(4.d0*pi))**(1.d0/3.d0),npar)/&
              &rhodep
        k1p = phi/vol*(gam-1.d0)/(gam+2.d0)*k1v

        ! step 2
        v1 = vol+dt*k1v/2.d0
        p1 = phi+dt*k1p/2.d0
        a = (3.d0*v1/(4.d0*pi*p1))**(1.d0/3.d0)
        c = a*p1
        rhodep = rhodep_vec(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0,npar)

        k2v = 4.d0*pi*capac_vec(a,c,npar)*si*g_vec(tmpk,pres,si,(3.d0*v1/(4.d0*pi))**(1.d0/3.d0),npar)/&
              &rhodep
        k2p = p1/v1*(gam-1.d0)/(gam+2.d0)*k2v

        ! step 3
        v2 = vol+dt*k2v/2.d0
        p2 = phi+dt*k2p/2.d0
        a = (3.d0*v2/(4.d0*pi*p2))**(1.d0/3.d0)
        c = a*p2
        rhodep = rhodep_vec(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0,npar)

        k3v = 4.d0*pi*capac_vec(a,c,npar)*si*g_vec(tmpk,pres,si,(3.d0*v2/(4.d0*pi))**(1.d0/3.d0),npar)/&
              &rhodep
        k3p = p2/v2*(gam-1.d0)/(gam+2.d0)*k2v

        ! step 4
        v3 = vol+dt*k3v
        p3 = phi+dt*k3p
        a = (3.d0*v3/(4.d0*pi*p3))**(1.d0/3.d0)
        c = a*p3
        rhodep = rhodep_vec(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0,npar)

        k4v = 4.d0*pi*capac_vec(a,c,npar)*si*g_vec(tmpk,pres,si,(3.d0*v3/(4.d0*pi))**(1.d0/3.d0),npar)/&
              &rhodep
        k4p = p3/v3*(gam-1.d0)/(gam+2.d0)*k3v


        vol = vol+1.d0/6.d0*dt*(k1v+2.d0*k2v+2.d0*k3v+k4v)
        phi = phi+1.d0/6.d0*dt*(k1p+2.d0*k2p+2.d0*k3p+k4p)

        ! update sizes
        a = (3.d0*vol/(4.d0*pi*phi))**(1.d0/3.d0)
        c = a*phi
        rhodep = rhodep_vec(a*1.d3,0.1d0,0.3d0,0.4d0,0.7d0,920.d0,npar)
        dmass = (vol-voli)*rhodep

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

    function capac_vec(a, c, npar)

        integer :: npar
        double precision, dimension(npar) :: capac_vec, a, c, ecc

        capac_vec = a
        where (a.gt.c) ecc = dsqrt(1.d0-(c/a)**2.d0)
        where (a.gt.c) capac_vec = a*ecc/(dasin(ecc))

        where (a.lt.c) ecc = dsqrt(1.d0-(a/c)**2.d0)
        where (a.lt.c) capac_vec = c*ecc/(dlog(c/a*(1.d0+ecc)))
    
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

    ! vectorized deposition density function
    function rhodep_vec(a,ac,fb,ft,fi,rhoi,npar)

        integer :: npar
        double precision, dimension(npar) :: a,rhodep_vec
        double precision :: ac,fb,rhoi,fi,ft,amax,wt,wsb,wmb,fmb,ai
        integer :: nsb

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
        rhodep_vec = rhoi*(fb/(amax-ai)*(ft*amax-ai+ai*amax*(1.d0-ft)/a)+&
                     (1.d0-fb)*fmb*ac/a)
        where (a.le.ac) rhodep_vec = rhoi
        where ((a.gt.ac).and.(a.le.ai)) rhodep_vec = rhoi*(fb+ac/a*fmb*(1.d0-fb))

      end function

end module

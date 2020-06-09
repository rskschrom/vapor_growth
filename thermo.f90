module thermo
implicit none
contains

    ! Goff-Gatch down to -100C for ice (Pa)
    function polysvp_ice(tmpk)

        double precision :: polysvp_ice, tmpk

        polysvp_ice = 10.d0**(-9.09718d0*(273.16d0/tmpk-1.d0)-3.56654d0*&
                      &dlog10(273.16d0/tmpk)+0.876793d0*(1.d0-tmpk/273.16d0)+&
                      &dlog10(6.1071d0))*100.d0
    
    end function

    ! Murphy + Koop (2005) RMET (Pa; T > 100 K)
    function svp_ice(tmpk)

        double precision :: svp_ice, tmpk

        svp_ice = dexp(9.550426d0-5723.265d0/tmpk+3.53068*dlog(tmpk)-0.00728332*tmpk)
    
    end function

    ! Murphy + Koop (2005) RMET (Pa; 332 K > T > 123 K)
    function svp_liq(tmpk)

        double precision :: lsvp_liq, svp_liq, tmpk

        lsvp_liq = 54.842763d0-6763.22d0/tmpk-4.210d0*dlog(tmpk)+0.000367d0*tmpk+&
                   &dtanh(0.0415d0*(tmpk-218.8d0))*&
                   &(53.878d0-1331.22d0/tmpk-9.44523d0*dlog(tmpk)+0.014025d0*tmpk)
        svp_liq = dexp(lsvp_liq)
    
    end function

    ! Murphy + Koop (2005) RMET (J mol-1; 273 K > T > 236 K)
    function lheat_vap(tmpk)

        double precision :: lheat_vap, tmpk

        lheat_vap = 56579.d0-42.212d0*tmpk+dexp(0.1149d0*(281.6d0-tmpk))
    
    end function

    ! Murphy + Koop (2005) RMET (J mol-1; T > 30 K)
    function lheat_sub(tmpk)

        double precision :: lheat_sub, tmpk

        lheat_sub = 46782.5d0+35.8925*tmpk-0.07414*tmpk**2.d0+&
                    541.5d0*dexp(-(tmpk/123.75d0)**2.d0)
    
    end function

    ! modified diffusivity P+K97 (eqn. 13-14)
    function dv_mod(tmpk, pres, rad)

        double precision :: dv_mod, dv, tmpk, pres, rad
        double precision ::ta, delv, alpc, pi

        pi = 3.14159265d0
        ta = tmpk
        delv = 1.3d0*8.d-8
        alpc = 0.036d0

        ! regular diffusivity (m2 s-1; PK97 eqn. 13-3)
        dv = 0.211d-4*(tmpk/273.15)**1.94d0*(101325.d0/pres)

        ! modified diffusivity
        dv_mod = dv/(rad/(rad+delv)+dv/(rad*alpc)*dsqrt(2.d0*pi*18.d-3/(8.314d0*ta)))
    
    end function

    ! vectorized modified diffusivity P+K97 (eqn. 13-14)
    function dv_vec(tmpk, pres, rad, npar)

        integer :: npar
        double precision, dimension(npar) :: dv_vec, rad
        double precision :: dv, tmpk, pres
        double precision ::ta, delv, alpc, pi

        pi = 3.14159265d0
        ta = tmpk
        delv = 1.3d0*8.d-8
        alpc = 0.036d0

        ! regular diffusivity (m2 s-1; PK97 eqn. 13-3)
        dv = 0.211d-4*(tmpk/273.15)**1.94d0*(101325.d0/pres)

        ! modified diffusivity
        dv_vec = dv/(rad/(rad+delv)+dv/(rad*alpc)*dsqrt(2.d0*pi*18.d-3/(8.314d0*ta)))
    
    end function

    ! modified conductivity P+K97 (eqn. 13-20)
    function ka_mod(tmpk, pres, rad)

        double precision :: ka_mod, ka, tmpk, pres, rad
        double precision ::ta, delt, alpt, cpa, rho, pi

        pi = 3.14159265d0
        ta = tmpk
        delt = 2.16d-7
        alpt = 0.7d0
        cpa = 1004.d0
        rho = pres/(287.d0*tmpk)

        ! regular conductivity (J m-1 s-1 k-1; PK97 eqn. 13-18a)
        ka = (5.69d0+0.017d0*(tmpk-273.15))*1.d-3

        ! modified diffusivity
        ka_mod = 4.184d0*ka/(rad/(rad+delt)+4.184d0*ka/(rad*alpt*rho*cpa)*dsqrt(2.d0*pi*29.d-3/(8.314d0*ta)))
    
    end function

    ! vectorized modified conductivity P+K97 (eqn. 13-20)
    function ka_vec(tmpk, pres, rad, npar)

        integer :: npar
        double precision, dimension(npar) :: ka_vec, rad
        double precision :: ka, tmpk, pres
        double precision :: ta, delt, alpt, cpa, rho, pi

        pi = 3.14159265d0
        ta = tmpk
        delt = 2.16d-7
        alpt = 0.7d0
        cpa = 1004.d0
        rho = pres/(287.d0*tmpk)

        ! regular conductivity (J m-1 s-1 k-1; PK97 eqn. 13-18a)
        ka = (5.69d0+0.017d0*(tmpk-273.15))*1.d-3

        ! modified diffusivity
        ka_vec = 4.184d0*ka/(rad/(rad+delt)+4.184d0*ka/(rad*alpt*rho*cpa)*dsqrt(2.d0*pi*29.d-3/(8.314d0*ta)))
    
    end function

    ! combined diffusivity
    function g_diff(tmpk, pres, si, rad)

        double precision :: g_diff, ka, dv, rv, vapi, tmpk, pres, si, rad, ls

        ka = ka_mod(tmpk, pres, rad)
        dv = dv_mod(tmpk, pres, rad)
        vapi = (si+1.d0)*svp_ice(tmpk)
        rv = 8.314d0/(18.d-3)
        ls = lheat_sub(tmpk)/(18.d-3)

        g_diff = 1.d0/(rv*tmpk/(vapi*dv)+ls/(tmpk*ka)*(ls/(tmpk*rv)-1.d0))        

    end function

    ! combined diffusivity vectorized
    function g_vec(tmpk, pres, si, rad, npar)

        integer :: npar
        double precision, dimension(npar) :: g_vec, ka, dv, rad
        double precision :: rv, vapi, tmpk, pres, si, ls

        ka = ka_vec(tmpk, pres, rad, npar)
        dv = dv_vec(tmpk, pres, rad, npar)
        vapi = (si+1.d0)*svp_ice(tmpk)
        rv = 8.314d0/(18.d-3)
        ls = lheat_sub(tmpk)/(18.d-3)

        g_vec = 1.d0/(rv*tmpk/(vapi*dv)+ls/(tmpk*ka)*(ls/(tmpk*rv)-1.d0))        

    end function

end module

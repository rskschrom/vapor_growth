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

        double precision :: dv_mod, tmpk, pres, rad
    
    end function


end module

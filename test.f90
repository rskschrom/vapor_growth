program test
use thermo

    double precision :: tk, vapi, vapl, lv, ls

    ! test thermo
    tk = 240.d0
    vapi = svp_ice(tk)
    vapl = svp_liq(tk)
    lv = lheat_vap(tk)
    ls = lheat_sub(tk)

    print *, vapi, vapl, vapl/vapi*100.d0, vapl-vapi
    print *, lv, ls

end program

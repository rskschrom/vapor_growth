program test
use thermo
use micro
implicit none

    double precision :: tk, vapi, vapl, lv, ls, a, c, si, pres, dt
    integer :: ntime, k

    ! test thermo
    tk = 258.15d0
    pres = 80000.d0
    vapi = svp_ice(tk)
    vapl = svp_liq(tk)
    lv = lheat_vap(tk)
    ls = lheat_sub(tk)
    !si = vapl/vapi-1.d0
    si = 0.15d0

    !print *, vapi, vapl, vapl/vapi*100.d0, vapl-vapi
    !print *, lv, ls

    !print *, dv_mod(tk, 80000.d0, 1.d-6)
    !print *, ka_mod(tk, 80000.d0, 1.d-6)
    !print *, g_diff(tk, 80000.d0, 0.04d0, 1.d-6)

    ! test micro
    a = 1.d-6
    c = 1.d-6
    dt = 1.d0
    ntime = 5400
    !print *, capac(a, c)/c

    do k = 1, ntime
        call spheroid_vap(tk, si, pres, dt, a, c)
        print *, k, a, c, c/a, si, tk
    end do

end program

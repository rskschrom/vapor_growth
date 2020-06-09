module util
implicit none
contains

    ! piecewise linear interpolation 
    subroutine interp_linear_pw(xgrid, xpoint, ygrid, ypoint, ngrid, npoint, xmin, xmax, ymin, ymax)

        integer :: ngrid, npoint, i
        double precision :: xmin, xmax, ymin, ymax
        double precision, dimension(ngrid) :: xgrid, ygrid
        double precision, dimension(npoint) :: xpoint, ypoint
        integer, dimension(1) :: close_ind
        double precision :: xpl, xph, wgh, wgl, ypl, yph

        ! loop over grid locs and find closest point bounds
        do i = 1, ngrid
            close_ind = minloc(dabs(xgrid(i)-xpoint))
            xpl = xpoint(close_ind(1))
            ypl = ypoint(close_ind(1))

            if (xpl.le.xgrid(i)) then
                if (close_ind(1).eq.npoint) then
                    xph = xmax
                    yph = ymax
                else
                    xph = xpoint(close_ind(1)+1)
                    yph = ypoint(close_ind(1)+1)
                end if
            else
                xph = xpl
                yph = ypoint(close_ind(1))
                if (close_ind(1).eq.1) then
                    xph = xmin
                    yph = ymin
                else
                    xpl = xpoint(close_ind(1)-1)
                    ypl = ypoint(close_ind(1)-1)
                end if
            end if

            ! calculate point weights
            wgh = (xgrid(i)-xpl)/(xph-xpl)
            wgl = 1.d0-wgh
            ygrid(i) = wgl*ypl+wgh*yph
        end do

        ! output to text for testing
        !open(unit=11, file="xg.txt")
        !open(unit=21, file="yg.txt")
        !open(unit=101, file="xp.txt")
        !open(unit=201, file="yp.txt")
        !write(11,*) xgrid
        !write(21,*) ygrid
        !write(101,*) xpoint
        !write(201,*) ypoint
        !close(11)
        !close(21)
        !close(101)
        !close(201)

    end subroutine
end module

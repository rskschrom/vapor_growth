module solver
implicit none
contains

    subroutine tridiagonal(avec, bvec, cvec, dvec, n, xsol)

        integer :: n, i
        double precision, dimension(n) :: bvec, dvec, xsol, dprime
        double precision, dimension(n-1) :: avec, cvec, cprime

        ! forward sweep
        cprime(1) = cvec(1)/bvec(1)
        dprime(1) = dvec(1)/bvec(1)
        do i = 2, n-1
            cprime(i) = cvec(i)/(bvec(i)-avec(i-1)*cprime(i-1))
            dprime(i) = (dvec(i)-avec(i-1)*dprime(i-1))/(bvec(i)-avec(i-1)*cprime(i-1))
        end do

        dprime(n) = (dvec(n)-avec(n-1)*dprime(n-1))/(bvec(n)-avec(n-1)*cprime(n-1))

        ! backward sweep
        xsol(n) = dprime(n)
        do i = 1, n-1
            xsol(n-i) = dprime(n-i)-cprime(n-i)*xsol(n-i+1)
        end do

    
    end subroutine
end module

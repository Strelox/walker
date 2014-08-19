! * HEigensystem.F
! * diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
! * code adapted from the "Handbook" routines for complex A
! * (Wilkinson, realinsch: Handbook for Automatic Computation, p. 202)
! * this file is part of the Diag library
! * last modified 9 Aug 11 th

! #include "diag-f.h"


! ************************************************************************
! ** HEigensystem diagonalizes a Hermitian n-by-n matrix.
! ** Input: n, A = n-by-n matrix, Hermitian
! ** (only the upper triangle of A needs to be filled).
! ** Output: d = vector of eigenvalues, U = transformation matrix
! ** these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.
module diag
    use accuracy

contains

    real(dp) function Sq(c)
        complex(dp) c
        Sq = real(c*conjg(c))
    end function Sq
    
    subroutine HEigensystem(n, A,ldA, d, U,ldU, sort)
        implicit none
        integer n, ldA, ldU, sort
        complex(dp) A(ldA,*), U(ldU,*)
        real(dp) d(*)

        integer p, q, j
        integer :: MAXDIM = 16
        real(dp) red, off, thresh
        real(dp) t, delta, invc, s,SYM_EPS
        complex(dp) x, y, Apq
        real(dp) :: ev(2,n)
        
        integer sweep
        common /nsweeps/ sweep

        SYM_EPS = err


        if( n .gt. MAXDIM ) then
        print *, "Dimension too large"
        d(1) = -999
        return
        endif

        do p = 1, n
        ev(1,p) = 0
        ev(2,p) = real(A(p,p))
        d(p) = ev(2,p)
        enddo

        do p = 1, n
        do q = 1, n
            U(q,p) = 0
        enddo
        U(p,p) = 1
        enddo

        red = .04D0/n**4

        do sweep = 1, 50
        off = 0
        do q = 2, n
            do p = 1, q - 1
            off = off + Sq(A(p,q))
            enddo
        enddo
        if( .not. off .gt. SYM_EPS ) goto 1

        thresh = 0
        if( sweep .lt. 4 ) thresh = off*red

        do q = 2, n
            do p = 1, q - 1
            Apq = A(p,q)
            off = Sq(Apq)
            if( sweep .gt. 4 .and. off .lt.&
            &              SYM_EPS*(ev(2,p)**2 + ev(2,q)**2) ) then
                A(p,q) = 0
            else if( off .gt. thresh ) then
                t = .5D0*(ev(2,p) - ev(2,q))
                t = 1/(t + sign(sqrt(t**2 + off), t))

                delta = t*off
                ev(1,p) = ev(1,p) + delta
                ev(2,p) = d(p) + ev(1,p)
                ev(1,q) = ev(1,q) - delta
                ev(2,q) = d(q) + ev(1,q)

                invc = sqrt(delta*t + 1)
                s = t/invc
                t = delta/(invc + 1)

                do j = 1, p - 1
                x = A(j,p)
                y = A(j,q)
                A(j,p) = x + s*(conjg(Apq)*y - t*x)
                A(j,q) = y - s*(Apq*x + t*y)
                enddo

                do j = p + 1, q - 1
                x = A(p,j)
                y = A(j,q)
                A(p,j) = x + s*(Apq*conjg(y) - t*x)
                A(j,q) = y - s*(Apq*conjg(x) + t*y)
                enddo

                do j = q + 1, n
                x = A(p,j)
                y = A(q,j)
                A(p,j) = x + s*(Apq*y - t*x)
                A(q,j) = y - s*(conjg(Apq)*x + t*y)
                enddo

                A(p,q) = 0

                do j = 1, n
                x = U(p,j)
                y = U(q,j)
                U(p,j) = x + s*(Apq*y - t*x)
                U(q,j) = y - s*(conjg(Apq)*x + t*y)
                enddo
            endif
            enddo
        enddo

        do p = 1, n
            ev(1,p) = 0
            d(p) = ev(2,p)
        enddo
        enddo

        print *, "Bad convergence in HEigensystem"

        1 if( sort .eq. 0 ) return

!         * sort the eigenvalues

        do p = 1, n - 1
        j = p
        t = d(p)
        do q = p + 1, n
            if( sort*(t - d(q)) .gt. 0 ) then
            j = q
            t = d(q)
            endif
        enddo

        if( j .ne. p ) then
            d(j) = d(p)
            d(p) = t
            do q = 1, n
            x = U(p,q)
            U(p,q) = U(j,q)
            U(j,q) = x
            enddo
        endif
        enddo
    end subroutine HEigensystem
    
!     * ApplyEigen.F
!     * transform the n-by-n matrix A by a similarity transformation U,
!     * i.e. form U A U^+.
!     * The original matrix A is overwritten with the result.
!     * this file is part of the Diag library
!     * last modified 9 Aug 11 th
    subroutine ApplyEigen(n, A,ldA, U,ldU)
    implicit none
    integer n, ldA, ldU
    complex(dp) A(ldA,*), U(ldU,*)

    integer i, j, k
    complex(dp) AUdag(n,n), s

    do j = 1, n
      do i = 1, n
        s = 0
        do k = 1, n
          s = s + A(i,k)*conjg(U(j,k))
        enddo
        AUdag(i,j) = s
      enddo
    enddo

    do j = 1, n
      do i = 1, n
        s = 0
        do k = 1, n
          s = s + U(i,k)*AUdag(k,j)
        enddo
        A(i,j) = s
      enddo
    enddo
    end subroutine ApplyEigen
    
!     * ReverseHEigen.F
!     * recover the original matrix A from the eigenvalues and
!     * the transformation matrix U by forming A = U^+ diag(d) U.
!     * this file is part of the Diag library
!     * last modified 9 Aug 11 th
    subroutine ReverseHEigen(n, A,ldA, d, U,ldU)
    implicit none
    integer n, ldA, ldU
    complex(dp) A(ldA,*), U(ldU,*)
    real(dp) d(*)

    integer i, j, k
    complex(dp) s

    do j = 1, n
      do i = 1, n
        s = 0
        do k = 1, n
          s = s + d(k)*conjg(U(k,i))*U(k,j)
        enddo
        A(i,j) = s
      enddo
    enddo
    end subroutine ReverseHEigen
    
end module diag
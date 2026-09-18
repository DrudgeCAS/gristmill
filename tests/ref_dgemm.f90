! Reference implementation of the BLAS DGEMM interface for tests.
!
! Only the 'N' and 'T' transposition flags are supported.  The routine is kept
! external with an implicit interface, exactly like the reference BLAS.

subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  implicit none
  character, intent(in) :: transa, transb
  integer, intent(in) :: m, n, k, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta
  double precision, intent(in) :: a(lda, *), b(ldb, *)
  double precision, intent(inout) :: c(ldc, *)
  integer :: i, j, l
  double precision :: tmp, aval, bval

  do j = 1, n
    do i = 1, m
      tmp = 0.0d0
      do l = 1, k
        if (transa == 'N') then
          aval = a(i, l)
        else
          aval = a(l, i)
        end if
        if (transb == 'N') then
          bval = b(l, j)
        else
          bval = b(j, l)
        end if
        tmp = tmp + aval * bval
      end do
      if (beta == 0.0d0) then
        c(i, j) = alpha * tmp
      else
        c(i, j) = alpha * tmp + beta * c(i, j)
      end if
    end do
  end do
end subroutine dgemm

program check_pis
   use regime
   use scattering_calculation
   use utils
   use constants
   use contexts
   use spheroidal_scatterer
   use spheroidal_indicatrix
   implicit none

   integer :: m, n, l, f, md, r, lnum, n0, l0, ssize, i
   complex(knd), allocatable :: p(:,:), appp(:,:)
   type(SpheroidalCalculation) :: inside, outside
   type(SpheroidalScatterer) :: scat
   complex(knd) :: ri, ddd, fl, sl, fc, sc
   complex(knd) :: ksi1, ksi2, approx, approx_pi
   real(knd) :: arg(1), alpha
   real(knd), allocatable :: p_coef(:)
   complex(knd), allocatable :: pi1(:,:), pi3(:,:)

   m = 1
   lnum = 300
   allocate(p(lnum, lnum), appp(lnum, lnum))
   n0 = 150
   l0 = 100
   f = 1
   ri = 1q0 !cmplx(1.7, 0.0, knd)
   alpha = 0.5q0
   arg(1) = cos(alpha)
   ! call scat%set(f, [3.0q0, 1.0q0], [2.0q0, 1.01q0], alpha, 2)
   call scat%set(f, [1.5q0, 1.2q0]/2.06209445549790393957521879532901848q0*1.5q0, [2.0q0, 2.0q0], alpha, 2, 2q0*PI)
   ! write(*,*) 'd = ', scat%d
   ! write(*,*) 'c = ', scat%c0
   ! c1^2/4 = 22
   ! write(*,*) 'sigma = ', (scat%d(1) + sqrt(scat%d(1)**2 - scat%d(2)**2)) / 2q0

   ! call inside%calculate(m, lnum, scat%c0(2) * ri, scat%ksi(2), 1, arg, f)
   ! call outside%calculate(m, lnum, scat%c0(1) * ri, scat%ksi(1), 1, arg, f)

   ! md = min(inside%maxd, outside%maxd)
   ! allocate(p_coef(0:md))
   ! p_coef = 1.0_knd / calculate_legendre_coef(m, md + 1)

   !  call compare_bb(m, n0, -inside%c ** 2 / outside%c ** 2)
   ! call compare_bb_direct(m, l0, -inside%c ** 2 / outside%c ** 2)
   ! call single_dep_integral(p, m, outside%legendre, inside%legendre, &
   ! p_coef, delta_coef)
   ! call log_matrix('pinit', p)

   ! p = 0
   ! do l = 1, lnum
   !     do n = 2 - mod(l, 2), lnum, 2
   !         do r = mod(n + 1, 2), md, 2
   !             p(n,l) = p(n,l) + abs(outside%legendre(r,n) * inside%legendre(r,l) * p_coef(r))
   !             if (n < l) then
   !                 appp(n,l) = approx_pi(m, l + m - 1, n + m - 1, abs(outside%c), abs(inside%c))
   !             else
   !                 appp(n,l) = approx_pi(m, n + m - 1, l + m - 1, abs(outside%c), abs(inside%c))
   !             endif
   !         enddo
   !     enddo
   ! enddo

   ! ssize = lnum
   ! call log_matrix('pcalc', p(:ssize,:ssize))
   ! call log_matrix('papprox', appp(:ssize,:ssize))
   ! call log_matrix('approx / p', appp(:ssize,:ssize) / p(:ssize,:ssize))
   ! n = n0 - m + 1
   ! l = l0 - m + 1
   ! write(*,*) 'm = ', m, 'n0 = ', n0, 'n = ', n
   ! do r = 1,5
   ! write(*,*) 'leg = ', outside%legendre(:20,r)
   ! enddo
   ! write(*,*) 'leg = ', outside%legendre(:20,n)
   ! do r = 0, md, 2
   !     ! approx = ddd(m, n0, r + m, outside%c)
   !     ! write(*,*) r, outside%legendre(r, n), approx, approx / outside%legendre(r,n)
   !     write(*,*) r, outside%legendre(r,n) * inside%legendre(r,l) * p_coef(r)
   ! enddo

   call compare_ratios(scat, alpha)
   ! call compare_convergeance(scat, alpha)

end program check_pis

subroutine compare_ratios(scat, alpha)
   use regime
   use scattering_calculation
   use utils
   use constants
   use contexts
   use spheroidal_scatterer
   implicit none

   type(SpheroidalScatterer) :: scat
   type(SpheroidalCalculation) :: inside, outside
   integer :: m, lnum, f, i, riind, xvind
   integer :: md
   real(knd), allocatable :: p_coef(:)
   real(knd) :: alpha, k, d(2), sigma2

   complex(knd) :: fl, sl, fc, sc, ri, c(2), ris(3)
   complex(knd), allocatable :: pi1(:,:), pi3(:,:)
   real(knd) :: etf, ets, xv(2), ab(2), xvmain, xvs(3)
   m = 1
   lnum = 60
   f = 1
   ab = [2q0, 2q0]
   ris = [cmplx(1.3, 0.0, knd), cmplx(1.7, 0.0, knd), cmplx(2.5, 0.0, knd)]
   xvs = [0.01, 0.03, 0.09]
   do riind=1,3
      do xvind = 1,3
   ri = ris(riind)
   ! write(*,*) 'd = ', scat%d
   ! write(*,*) 'sigma^2 = ', (scat%d(1) + sqrt(scat%d(1)**2 - scat%d(2)**2))**2 *0.25
   xvmain = xvs(xvind)
   ! core = 1/2 of volume
   xv = [xvmain, xvmain / 2**(1.0q0/3.0q0)]
   ab = [4q0, 4q0]
   k = 1q0
   d = [1.5q-2, 1.2q-2]
   ! c = d / 2 * k
   ! c = k * d / 2q0
   c = xv * sqrt(ab**2q0 - 1q0) * (1q0 / ab)**(1q0 / 3q0) * ri
   ! c(2) = c(1) / 1.5q0
   d = 2q0 * c / k
   write(*,*) 'xv = ', xv
   write(*,*) 'ab = ', ab
   write(*,*) 'ri = ', ri
   write(*,*) 'c = ', c
   write(*,*) 'k = ', k
   write(*,*) 'd = ', d

   sigma2 = (d(1) + sqrt(d(1)**2 - d(2)**2))**2 *0.25
   write(*,*) 'sigma^2 = ', sigma2

   call inside%calculate(m, lnum, c(2), 1.1_knd, 1, [cos(alpha)], f)
   call outside%calculate(m, lnum, c(1), 1.1_knd, 1, [cos(alpha)], f)
   ! write(*,*) 'c = ', inside%c, outside%c

   ! call inside%calculate(m, lnum, 1.5 scat%ksi(2), 1, [cos(alpha)], f)
   ! call outside%calculate(m, lnum, 1.2, scat%ksi(1), 1, [cos(alpha)], f)

   md = min(inside%maxd, outside%maxd)
   allocate(p_coef(0:md))
   p_coef = 1.0_knd / calculate_legendre_coef(m, md + 1)
   allocate(pi1(lnum, lnum), pi3(lnum, lnum))
   call single_dep_integral(pi1(:,:), m, outside%legendre, inside%legendre, p_coef, delta_coef)
   call mult_matr_by_i(pi1(:,:), lnum)
   ! call log_matrix("pi1", pi1)
   ! call log_matrix("pi1*pi^T", matmul(pi1, transpose(pi1)))
   pi3 = pi1
   ! call multiply_by_diag_left(pi1(:,:), lnum, inside%r1(1:lnum))
   ! call multiply_by_diag_right(pi1(:,:), lnum, 1.0_knd / outside%r1(1:lnum))
   ! call multiply_by_diag_left(pi3(:,:), lnum, inside%r3(1:lnum))
   ! call multiply_by_diag_right(pi3(:,:), lnum, 1.0_knd / outside%r3(1:lnum))
   ! call log_matrix('pi1', pi1)
   write(*,*) 'element ratios pi1'
   100 format('#',10A24)
   101 format(' ',4E24.15)
   write(*, 100) 'first_row', 'second_row', 'first_column', 'second_column'
   do i = 1, lnum - 3, 2
      fl = pi1(1, i + 2) / pi1(1, i)! * i**2! * (i + 2q0) / i
      sl = pi1(2, i + 3) / pi1(2, i + 1)! * i**2! * (i + 3q0) / (i + 1q0)
      fc = pi1(i + 2, 1) / pi1(i, 1)! * i**2 !* (i + 2q0) / i
      sc = pi1(i + 3, 2) / pi1(i + 1, 2)!* i**2 !* (i + 3q0) / (i + 1q0)
      ! etf = ((i + 2q0) / i)**(m-1) * sigma2
      ! ets = ((i + 3q0) / (i + 1q0))**(m-1) * sigma2
      write(*, 101) real(fl), real(sl), real(fc), real(sc) !, etf, ets
   enddo
   ! write(*,*) 'element ratios pi3'
   ! write(*, 100) 'fl_re', 'fl_im', 'sl_re', 'sl_im', 'fc_re', 'fc_im', 'sc_re', 'sc_im'
   ! do i = 1, lnum - 3, 2
   !    fl = pi3(1, i + 2) / pi3(1, i)!  * i**2
   !    sl = pi3(2, i + 3) / pi3(2, i + 1)!  * i**2
   !    fc = pi3(i + 2, 1) / pi3(i, 1)!  * i**2
   !    sc = pi3(i + 3, 2) / pi3(i + 1, 2)!  * i**2
   !    etf = ((i + 2q0) / i)**(m-1) * sigma2
   !    ets = ((i + 3q0) / (i + 1q0))**(m-1) * sigma2
   !    write(*, 101) real(fl), imag(fl), real(sl), imag(sl), real(fc), imag(fc), real(sc), imag(sc), etf, ets
   ! enddo
   write(*,*)
   deallocate(p_coef, pi1, pi3)
enddo
enddo
end subroutine compare_ratios

subroutine compare_spheroidal_functions()
end subroutine compare_spheroidal_functions

real(knd) function rel_diff1(n, a, b)
use regime
implicit none
real(knd), parameter :: eps = 1q-256
complex(knd) :: a(n), b(n)
integer :: i, n
rel_diff1 = 0
do i = 1, size(a)
   if (abs(a(i)) < eps) cycle
   rel_diff1 = max(rel_diff1, abs(a(i) - b(i)) / (abs(a(i)) + abs(b(i))))
enddo
end function rel_diff1

subroutine compare_convergeance(scat, alpha)
   use regime
   use scattering_calculation
   use utils
   use constants
   use contexts
   use spheroidal_scatterer
   implicit none

   type(SpheroidalScatterer) :: scat
   type(SpheroidalCalculation) :: inside, outside
   integer :: m, lnum, f, i
   integer :: md
   real(knd), allocatable :: p_coef(:)
   real(knd) :: alpha, rel_diff1
   logical :: found
   complex(knd) :: ri

   complex(knd), allocatable :: pi1(:,:), pi3(:,:), fl(:)
   m = 1
   lnum = 150
   f = 1
   ri = cmplx(1.7, 0.0, knd)
   fl = [cmplx(1q0, 0q0, knd)]
   found = .false.
   ! do lnum = 4, 200, 2

      call inside%calculate(m, lnum, scat%c0(2) * ri, scat%ksi(2), 1, [cos(alpha)], f)
      call outside%calculate(m, lnum, scat%c0(1) * ri, scat%ksi(1), 1, [cos(alpha)], f)

      md = min(inside%maxd, outside%maxd)
      ! allocate(p_coef(0:md))
      p_coef = 1.0_knd / calculate_legendre_coef(m, md + 1)
      allocate(pi1(lnum, lnum), pi3(lnum, lnum))
      call single_dep_integral(pi1(:,:), m, inside%legendre, outside%legendre, p_coef, delta_coef)
      call mult_matr_by_i(pi1(:,:), lnum)
      pi3 = pi1
      call multiply_by_diag_left(pi1(:,:), lnum, inside%r1(1:lnum))
      call multiply_by_diag_right(pi1(:,:), lnum, 1.0_knd / outside%r1(1:lnum))
      call multiply_by_diag_left(pi3(:,:), lnum, inside%r3(1:lnum))
      call multiply_by_diag_right(pi3(:,:), lnum, 1.0_knd / outside%r3(1:lnum))
      ! call log_matrix('pi1', pi1)
   !    if (rel_diff1(size(fl), fl, pi1(:size(fl), 1)) < 1q-32) then
   !       found = .true.
   !       exit
   !    endif
   !    fl = pi1(:,1)
   !    deallocate(pi1, pi3)
   ! enddo
   ! if (found) then
      ! write(*,*) 'found convergeance: lnum=', lnum
      call log_matrix('pi1', pi1(:,1:1))
   ! else
      ! write(*,*) 'does not converge'
   ! endif
end subroutine compare_convergeance

complex(knd) function ddd(m, n, r, c)
   use regime
   implicit none

   integer :: m, n, r, i
   complex(knd) :: c, deg
   real(knd) :: coef

   coef = 2.0q0 * n + 1.0q0
   do i = n - m + 1, m + n
      coef = coef / i
   enddo
   coef = sqrt(coef)

   deg = 1q0
   do i = 1, max(n - m, r) - min(n - m, r)
      deg = deg * c / 4q0
   enddo

   ! write(*,*) m, n, r, c, coef, deg
   if (r <= n - m) then
      ddd = coef * deg * gamma((n + m + r + 1q0) / 2q0) / (gamma((n - m - r + 2q0) / 2q0) * gamma(n + 0.5q0))
   else
      ddd = coef * deg * gamma(n + 1.5q0) / (gamma((r - n + m + 2q0) / 2q0) * gamma((n + m + r + 3q0) / 2q0))
   endif
end function ddd

real(knd) function nnorm(m, n)
   use regime
   implicit none
   integer, intent(in) :: m, n
   integer :: i

   nnorm = 2q0 / (2q0 * n + 1q0)
   do i = n - m + 1, n + m
      nnorm = nnorm * i
   enddo
   nnorm = sqrt(nnorm)

end function nnorm

complex(knd) function approx_pi(m, n, l, c1, c2)
   use regime
   implicit none
   real(knd) :: nlfac, nfac, lfac, nnorm
   integer :: i, m, n, l
   real(knd) :: c1, c2, c

   nlfac = 1q0
   do i = n - l, 1, -2
      nlfac = nlfac * i
   enddo

   nfac = 1q0
   do i = 2 * n - 1, 1, -4
      nfac = nfac * i
   enddo

   lfac = 1q0
   do i = 2 * l - 3, 1, -4
      lfac = lfac * i
   enddo

   c = c2 ** 2 / c1 ** 2

   ! write(*,*) n, l, nfac, lfac, nlfac, c1, c2, c, nnorm(m, l), nnorm(m, n)
   approx_pi = nnorm(m, l) / nnorm(m, n) * (c1 / 2q0) ** (n - l) * lfac / nfac / nlfac * &
      (1q0 + ((2q0 * c) ** ((n - l) / 2 + 1) - 1q0) / (2q0 * c - 1q0))
end function approx_pi

real(knd) function alp(n,l)
   use regime
   implicit none
   integer :: n,l

   alp = (n + l - 1q0) * (n - l + 2q0) / (2q0 * l + 1q0) / (2q0 * l - 1)

end function alp

real(knd) function bet(n,l)
   use regime
   implicit none
   integer :: n, l

   bet = (n - l - 2q0) * (n + l + 3q0) / (2q0 * l + 1q0) / (2q0 * l + 3q0)

end function bet

real(knd) function phi(n,l)
   use regime
   implicit none
   integer :: n,l

   phi = (n + l - 1q0) * (n - l - 2q0) / (2q0 * n + 1q0) / (2q0 * n - 1)

end function phi

real(knd) function psi(n,l)
   use regime
   implicit none
   integer :: n, l

   psi = (n - l + 2q0) * (n + l + 3q0) / (2q0 * n + 1q0) / (2q0 * n + 3q0)

end function psi

real(knd) function fac(n, k)
   use regime
   implicit none
   integer :: n, k, i
   fac = 1q0
   do i = n, 1, -k
      fac = fac * i
   enddo
end function fac

complex(knd) function yy(m,n,l,s,c)
   use regime
   implicit none
   integer :: m,n,l,s, i
   complex(knd) :: c
   real(knd) :: fac

   yy = 1q0
   do i = s - m + 1, s + m
      yy = yy * i
   enddo
   do i = l + s + 3, n + s - 1, 2
      yy = yy * i
   enddo
   do i = n + s + 1, l + s + 1, 2
      yy = yy / i
   enddo
   yy = yy / fac((n - s) / 2, 1) / fac((s - l) / 2, 1) &
      / (2q0 * s + 1q0) * c ** (s / 2)
end function yy

complex(knd) function BB(m,n,l,c)
   use regime
   implicit none
   integer :: m,n,l,s
   complex(knd) :: c, yy
   BB = 0
   do s = l, n, 4
    ! write(*,*) 'add', yy(m,n,l,s,c)
      BB = BB + yy(m,n,l,s,c)
   enddo
   do s = l + 2, n, 4
    ! write(*,*) 'add', yy(m,n,l,s,c)
      BB = BB + yy(m,n,l,s,c)
   enddo
end function BB

complex(knd) function SS(m,n,l,c, a)
   use regime
   implicit none
   integer :: m,n,l,s, a
   complex(knd) :: c, yy
   SS = 0
   do s = l, n, 4
      SS = SS + yy(m,n,l,s,c) / (s - a)
   enddo
   do s = l + 2, n, 4
      SS = SS + yy(m,n,l,s,c) / (s - a)
   enddo
end function SS

! complex(knd) function BB(m,n,l,c)
!    use regime
!    implicit none
!    integer :: m,n,l, i, s
!    complex(knd) :: c
!    real(knd) :: smfac, nfac, lfac, nsfac, lsfac

!    BB = 0

!    smfac = 1q0
!    do i = l - m + 1, l + m
!       smfac = smfac * i
!    enddo
!    nfac = 1q0
!    do i = n + l - 1, 1, -2
!       nfac = nfac * i
!    enddo
!    lfac = 1
!    do i = 2 * l + 1, 1, -2
!       lfac = lfac * i
!    enddo
!    nsfac = 1
!    do i = 1, (n - l) / 2
!       nsfac = nsfac * i
!    enddo
!    lsfac = 1

!    do s = l, n, 2
!       BB = BB + smfac * nfac / lfac / nsfac / lsfac * c ** (s / 2q0) / (2q0 * s + 1)

!       smfac = smfac * (s + m + 1q0) * (s + m + 2q0) / (s - m + 1q0) / (s - m + 2q0)
!       nfac = nfac * (n + s + 1q0)
!       lfac = lfac * (l + s + 3q0)
!       nsfac = nsfac / (n - s) * 2q0
!       lsfac = lsfac * ((s - l) / 2q0 + 1q0)

!    enddo
! end function BB

subroutine compare_bb(m, n, c)
   use regime
   implicit none
   real(knd) :: alp, bet
   complex(knd) :: c, BB, prev, prevprev, calc, etalon, alpha, beta, SS, yy, left, right
   integer :: m, n, l, s
   write(*,*) 'm = ', m, 'n = ', n
   !  write(*,*) 'compare frac'
   !  do l = n - 2, m, -2
   !     do s = l, n, 2
   !        left = (n + s + 1q0) * (n - s) / (s + l + 3q0) / (s - l + 2q0) * (2q0 * l + 1q0)
   !        right = (-2q0 * l - 1q0 - (n - l - 2q0) * (n + l + 3q0) / (s + l + 3q0) + (n + l - 1q0) * (n - l + 2q0) / (s - l + 2q0))

   !        write(*,*) l,s,left,right,left / right
   !     enddo
   !  enddo
   !  write(*,*)
   !  write(*,*) 'compare frac'
   !  do l = n - 2, m, -2
   !     do s = l, n, 2
   !        left = (n + s + 1q0) * (n - s) / (s + l + 3q0) / (s - l + 2q0)
   !        right = (n + l - 1q0) * (n - l + 2q0) / (2q0 * l + 1q0) / (s - l + 2q0) &
   !        - (n - l - 2q0) * (n + l + 3q0) / (2q0 * l + 1q0) / (s + l + 3q0) - 1

   !        write(*,*) l,s,left,right,left / right
   !     enddo
   !  enddo
   !  write(*,*)
   !  write(*,*) 'compare y s, s + 2 before expansion'
   !  do l = n - 2, m, -2
   !     do s = l, n - 2, 2
   !        left = yy(m,n,l,s + 2,c)
   !        right = c * yy(m,n,l,s,c) * (n + s + 1q0) * (n - s) / (l + s + 3q0) / (s - l + 2)!&
   !       !  * (2q0 * s + 1q0) / (2q0 * s + 5q0) * (s + m + 2q0) * (s + m + 1q0) / (s - m + 2q0) / (s - m + 1q0)
   !       !  * (1q0 - 4q0 / (2q0 * s + 5) + 2q0 * m / (s + 2q0 - m) + 2q0 * m / (s + 1q0 - m))
   !        write(*,*) l,s,left,right,left / right
   !     enddo
   !  enddo
   !  write(*,*) 'compare y s, s + 2'
   !  do l = n - 2, m, -2
   !     do s = l, n - 2, 2
   !        left = yy(m,n,l,s + 2,c)
   !        right = c * yy(m,n,l,s,c) * &
   !           (-(n - l - 2q0) * (n + l + 3q0) / (s + l + 3q0) / (2q0 * l + 1)&
   !            + (n + l - 1q0) * (n - l + 2q0) / (s - l + 2q0)/ (2q0 * l + 1) - 1)

   !        write(*,*) l,s,left,right,left / right
   !     enddo
   !  enddo
   !  write(*,*)
   !  write(*,*) 'compare SS l - 2'
   !  do l = n, m, -2
   !     left = SS(m,n,l,c,l-2) * 2q0 * (2q0 * l - 1)
   !     right = BB(m,n,l-2,c) - yy(m,n,l-2,l-2,c) - 2q0 * BB(m,n,l,c)
   !     write(*,*) l, left, right, left / right
   !  enddo
   !  write(*,*)
   !  write(*,*) 'compare SS -l - 3'
   !  do l = n, m, -2
   !     left = SS(m,n,l,c,-l-3) / 2q0 * (2q0 * l + 3q0)
   !     right = -BB(m,n,l+2,c) + 0.5q0 * BB(m,n,l,c)
   !     write(*,*) l, left, right, left / right
   !  enddo
   !  write(*,*)
   !  write(*,*) 'compare recurrent 1'
   !  do l = n - 2, m + 2, -2
   !     left = BB(m,n,l,c) - yy(m,n,l,l,c)
   !     right = c * ((n + l - 1q0) * (n - l + 2q0) / (2q0 * l + 1) * (SS(m,n,l,c,l - 2) - yy(m,n,l,n,c) / (l - n + 2q0)) - &
   !     (n - l - 2q0) * (n + l + 3q0) / (2q0 * l + 1) * (SS(m,n,l,c,-l-3) - yy(m,n,l,n,c) / (n + l + 3q0)) - &
   !     (BB(m,n,l,c) - yy(m,n,l,n,c)))

   !     write(*,*) l, left, right, left / right
   !  enddo

   !  write(*,*) 'compare recurrent 2'
   !  do l = n - 2, m + 2, -2
   !     left = BB(m,n,l,c) - yy(m,n,l,l,c)
   !     right = c * ((n + l - 1q0) * (n - l + 2q0) / (2q0 * l + 1) * (SS(m,n,l,c,l - 2)) - &
   !     (n - l - 2q0) * (n + l + 3q0) / (2q0 * l + 1) * (SS(m,n,l,c,-l-3)) - &
   !     (BB(m,n,l,c)))

   !     write(*,*) l, left, right, left / right
   !  enddo

   !  write(*,*) 'compare recurrent 3'
   !  do l = n - 2, m + 2, -2
   !     left = BB(m,n,l,c) - yy(m,n,l,l,c)
   !     right = c * ((n + l - 1q0) * (n - l + 2q0) / (2q0 * l + 1) / (2q0 * l - 1q0) / 2q0 * &
   !     (BB(m,n,l-2,c) - yy(m,n,l-2,l-2,c) - 2q0 * BB(m,n,l,c)) - &
   !     2q0 * (n - l - 2q0) * (n + l + 3q0) / (2q0 * l + 1) / (2q0 * l + 3q0) * (0.5q0 * BB(m,n,l,c) - BB(m,n,l+2,c)) - &
   !     (BB(m,n,l,c)))

   !     write(*,*) l, left, right, left / right
   !  enddo
   !  write(*,*) 'compare recurrent 4'
   !  do l = n - 2, m + 2, -2
   !     left = BB(m,n,l,c)
   !     right = c * ((n + l - 1q0) * (n - l + 2q0) / (2q0 * l + 1) / (2q0 * l - 1q0) / 2q0 * &
   !     (BB(m,n,l-2,c) - 2q0 * BB(m,n,l,c)) - &
   !     2q0 * (n - l - 2q0) * (n + l + 3q0) / (2q0 * l + 1) / (2q0 * l + 3q0) * (0.5q0 * BB(m,n,l,c) - BB(m,n,l+2,c)) - &
   !     (BB(m,n,l,c)))

   !     write(*,*) l, left, right, left / right
   !  enddo

   !  write(*,*) 'compare recurrent 5'
   !  do l = n - 2, m + 2, -2
   !     left = BB(m,n,l,c) * (1q0 + c * (1q0 + (n + l - 1q0) * (n - l + 2q0) / (2q0 * l + 1) / (2q0 * l - 1q0) + &
   !     (n - l - 2q0) * (n + l + 3q0) / (2q0 * l + 1) / (2q0 * l + 3q0)))
   !     right = c * ((n + l - 1q0) * (n - l + 2q0) / (2q0 * l + 1) / (2q0 * l - 1q0) / 2q0 * &
   !     (BB(m,n,l-2,c)) + &
   !     2q0 * (n - l - 2q0) * (n + l + 3q0) / (2q0 * l + 1) / (2q0 * l + 3q0) * ( BB(m,n,l+2,c)))

   !     write(*,*) l, left, right, left / right
   !  enddo

   write(*,*) 'compare recurrent 6'
   do l = n - 2, m + 2, -2
      alpha = alp(n,l)
      beta = bet(n,l)
      left = BB(m,n,l,c) * (1q0 + c * (1q0 + alpha + beta))
      right = c * (alpha / 2q0 * (BB(m,n,l-2,c)) + 2q0 * beta * ( BB(m,n,l+2,c)))

      write(*,*) l, left, right, left / right
   enddo
   write(*,*) 'compare recurrent'
   prevprev = BB(m,n,n,c)
   prev = BB(m,n,n-2,c)
   do l = n - 2, n / 2, -2
      etalon = BB(m,n,l-2,c)
      alpha = alp(n,l)
      beta = bet(n,l)

      calc = -4q0 * beta / alpha * prevprev + 2q0 / alpha  * (1q0 / c + 1q0 + alpha + beta ) * prev

      ! write(*,*) l, calc, etalon, calc / etalon
      left = 2q0 / alpha  * (1q0 / c + 1q0 + alpha + beta ) * prev
      right = 4q0 * beta / alpha * prevprev
      write(*,*) l, left, right, left - right, (left - right) / left
      prevprev = prev
      prev = calc
   enddo
end subroutine compare_bb

subroutine compare_bb_direct(m, l, c)
   use regime
   implicit none
   real(knd) :: phi, psi
   complex(knd) :: c, BB, prev, prevprev, calc, etalon, alpha, beta, SS, yy, left, right
   integer :: m, n, l, s, maxl
   maxl = 4 * l
   write(*,*) 'm = ', m, 'l = ', l
   !  write(*,*) 'compare frac'
   !  do n = l, 4 * l, 2
   !     do s = l, n, 2
   !        left = (l + s + 1q0) * (s-l) / (n + s - 1q0) / (n - s + 2q0)
   !        right = ( (n - l - 2q0) * (n + l - 1q0) / (s + n - 1q0) / (2q0 * n + 1q0) + &
   !        (n + l + 3q0) * (n - l + 2q0) / (n - s + 2q0) / (2q0 * n + 1q0) - 1q0)

   !        write(*,*) n,s,left,right,left / right
   !     enddo
   !  enddo
   !  write(*,*)

   !  write(*,*) 'compare y s, s - 2 before expansion'
   !  do n = l + 2, maxl, 2
   !     do s = l + 2, n, 2
   !        left = yy(m,n,l,s - 2,c)
   !        right = 1q0 / c * yy(m,n,l,s,c) * (l + s + 1q0) * (s - l) / (n + s - 1q0) / (n - s + 2)!&
   !       !  * (2q0 * s + 1q0) / (2q0 * s + 5q0) * (s + m + 2q0) * (s + m + 1q0) / (s - m + 2q0) / (s - m + 1q0)
   !       !  * (1q0 - 4q0 / (2q0 * s + 5) + 2q0 * m / (s + 2q0 - m) + 2q0 * m / (s + 1q0 - m))
   !        write(*,*) n,s,left,right,left / right
   !     enddo
   !  enddo
   !  write(*,*) 'compare y s, s - 2'
   !  do n = l + 2, maxl, 2
   !     do s = l+2, n, 2
   !        left = yy(m,n,l,s - 2,c)
   !        right = 1q0 / c * yy(m,n,l,s,c) * &
   !           ((n - l - 2q0) * (n + l - 1q0) / (s + n - 1q0) / (2q0 * n + 1q0) + &
   !                  (n + l + 3q0) * (n - l + 2q0) / (n - s + 2q0) / (2q0 * n + 1q0) - 1q0)

   !        write(*,*) n,s,left,right,left / right
   !     enddo
   !  enddo
   !  write(*,*)
  !  write(*,*) 'compare SS n + 2'
  !  do n = l, maxl, 2
  !     left = -SS(m,n,l,c,n + 2) * 2q0 * (2q0 * n + 3q0)
  !     right = BB(m,n+2,l,c) - yy(m,n+2,l,n+2,c) + 2q0 * BB(m,n,l,c)
  !     ! write(*,*) l, BB(m,n+2,l,c), yy(m,n+2,l,n+2,c), 2q0 * BB(m,n,l,c)
  !     write(*,*) n, left, right, left / right
  !  enddo
  !  write(*,*)
    ! write(*,*) 'compare SS -n + 1'
    ! do n = l, maxl, 2
    !    left = SS(m,n,l,c,1- n) * (2q0 * n - 1q0)
    !    right = 2q0 * BB(m,n-2,l,c) + BB(m,n,l,c)
    !    write(*,*) n, left, right, left / right
    ! enddo
    ! write(*,*)
    ! write(*,*) 'compare recurrent 1'
    ! do n = l + 2, maxl, 2
    !    left = BB(m,n,l,c) - yy(m,n,l,n,c)
    !    right = 1q0 / c * ((n + l - 1q0) * (n - l - 2q0) / (2q0 * n + 1) * (SS(m,n,l,c,1-n) - yy(m,n,l,l,c) / (n + l - 1q0)) + &
    !    (n - l + 2q0) * (n + l + 3q0) / (2q0 * n + 1) * (-SS(m,n,l,c,n+2) - yy(m,n,l,l,c) / (n - l + 2q0)) - &
    !    (BB(m,n,l,c) - yy(m,n,l,l,c)))

    !    write(*,*) l, left, right, left / right
    ! enddo

    ! write(*,*) 'compare recurrent 2'
    ! do n = l + 2, maxl, 2
    !    left = BB(m,n,l,c) - yy(m,n,l,n,c)
    !    right = 1q0 / c * ((n + l - 1q0) * (n - l - 2q0) / (2q0 * n + 1) * (SS(m,n,l,c,1-n)) + &
    !       (n - l + 2q0) * (n + l + 3q0) / (2q0 * n + 1) * (-SS(m,n,l,c,n+2)) - &
    !       (BB(m,n,l,c)))

    !    write(*,*) l, left, right, left / right
    ! enddo

    ! write(*,*) 'compare recurrent 3'
    ! do n = l + 2, maxl, 2
    !    left = BB(m,n,l,c) - yy(m,n,l,n,c)
    !    right = 1q0 / c * ((n + l - 1q0) * (n - l - 2q0) / (2q0 * n + 1)  / (2q0 * n - 1) * (2q0 * BB(m,n-2,l,c) + BB(m,n,l,c)) + &
    !       (n - l + 2q0) * (n + l + 3q0) / (2q0 * n + 1) / (2q0 * n + 3q0) / 2q0 * (BB(m,n+2,l,c) + 2q0 * BB(m,n,l,c) &
    !       - yy(m,n+2,l,n+2,c)) - &
    !       (BB(m,n,l,c)))
    !    write(*,*) l, left, right, left / right
    ! enddo
    ! write(*,*) 'compare recurrent 4'
    ! do n = l + 2, maxl, 2
    !    left = BB(m,n,l,c)
    !    right = 1q0 / c * ((n + l - 1q0) * (n - l - 2q0) / (2q0 * n + 1)  / (2q0 * n - 1) * (2q0 * BB(m,n-2,l,c) + BB(m,n,l,c)) + &
    !       (n - l + 2q0) * (n + l + 3q0) / (2q0 * n + 1) / (2q0 * n + 3q0) / 2q0 * (BB(m,n+2,l,c) + 2q0 * BB(m,n,l,c)) - &
    !       (BB(m,n,l,c)))

    !    write(*,*) l, left, right, left / right
    ! enddo

    ! write(*,*) 'compare recurrent 5'
    ! do n = l + 2, maxl, 2
    !    left = BB(m,n,l,c) * (1q0 - 1q0 / c * (n + l - 1q0) * (n - l - 2q0) / (2q0 * n + 1)  / (2q0 * n - 1) - &
    !    1q0 / c * (n - l + 2q0) * (n + l + 3q0) / (2q0 * n + 1) / (2q0 * n + 3q0) + 1q0 / c)
    !    right = 1q0 / c * ((n + l - 1q0) * (n - l - 2q0) / (2q0 * n + 1)  / (2q0 * n - 1) * (2q0 * BB(m,n-2,l,c)) + &
    !       (n - l + 2q0) * (n + l + 3q0) / (2q0 * n + 1) / (2q0 * n + 3q0) / 2q0 * (BB(m,n+2,l,c)) )

    !    write(*,*) l, left, right, left / right
    ! enddo
    ! write(*,*) 'compare recurrent 6'
    ! do n = l + 2, maxl, 2
    !         alpha = phi(n,l)
    !   beta = psi(n,l)
    !    left = BB(m,n,l,c) * (1q0 - 1q0 / c * alpha - &
    !    1q0 / c * beta + 1q0 / c)
    !    right = 1q0 / c * (alpha * (2q0 * BB(m,n-2,l,c)) + &
    !       beta / 2q0 * (BB(m,n+2,l,c)) )

    !    write(*,*) l, left, right, left / right
    ! enddo

  !  write(*,*) 'compare recurrent'
  !  prevprev = BB(m,l,l,c)
  !  prev = BB(m,l + 2,l,c)
  !  do n = l + 2, 4 * l, 2
  !     etalon = BB(m,n + 2,l,c)
  !     alpha = phi(n,l)
  !     beta = psi(n,l)

  !     calc = -4q0 * alpha / beta * prevprev + 2q0 / beta  * (c + 1q0 - alpha - beta ) * prev

  !     ! write(*,*) l, calc, etalon, calc / etalon

  !     write(*,*) l, etalon, calc, calc / etalon
  !     prevprev = prev
  !     prev = calc
  !  enddo
   write(*,*) 'compare last element'
   do n = l, 4 * l, 2
    write(*,*) 'sum = ', BB(m,n,l,c), 'elem = ', yy(m,n,l,n,c)
   enddo
end subroutine compare_bb_direct

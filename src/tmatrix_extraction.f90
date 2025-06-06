module tmatrix_extraction
    use regime
    use constants
    use contexts
    use scattering_calculation
    use iso_c_binding, only: C_FLOAT128, C_FLOAT128_COMPLEX, C_INT32_T
    implicit none
contains

    subroutine GetSpheroidalTmatrix(nol_in, ab_in, rv_in, lambda_in, lnum_in, ri_in, m_in, res)
        real(C_FLOAT128) :: ab_in(nol_in), rv_in(nol_in), lambda_in
        integer(C_INT32_T) :: lnum_in, m_in, nol_in
        complex(C_FLOAT128_COMPLEX) :: ri_in(nol_in), res(2*lnum_in, 2*lnum_in)

        type(ScatteringContext) :: global_context
        integer :: nol
        integer :: f, lnum, spherical_lnum, m, result_num, i
        real(knd), parameter :: alpha = 0.0_knd !PI / 4.0_knd
        real(knd) :: xv(nol_in), ab(nol_in), lambda
        complex(knd) :: ri(0:nol_in)
        type(Node), allocatable:: queue(:)
        type(ModeCalculationResult), allocatable :: mode_res(:)

        f = -1
        ab = ab_in
        if (ab_in(1) < 1.0q0) then
            f = 1
            ab = 1q0 / ab
        end if
        nol = nol_in
        xv = rv_in * 2.0_knd * PI / lambda_in
        lambda = lambda_in
        ri(0) = cmplx(1.0_knd, 0.0_knd, knd)
        ri(1:nol) = [(cmplx(real(ri_in(i), knd), imag(ri_in(i)), knd), i=1,nol)]
        lnum = lnum_in
        spherical_lnum = lnum_in
        if (mod(lnum_in, 2) == 1) then
            lnum = lnum_in + 1
            spherical_lnum = lnum_in + 1
        endif
        if (lnum < 8) then
            lnum = 8
            spherical_lnum = 8;
        endif
        m = m_in

        write(*,*) 'fortran f = ', f, ' ab = ', ab, 'xv = ', xv, 'lambda = ', lambda, 'lnum = ', lnum, &
        'lnum_in = ', lnum_in, 'ri = ', ri, 'm = ', m

        call global_context%initialize(f, nol, xv, ab, alpha, lambda, ri, lnum, spherical_lnum, m, m, &
        0, 0.0_knd, 0.0_knd, 0, 0.0_knd, 0.0_knd)
        if (m == 0) then
            queue = [&
                Node(MODE_SPH_TE_PQ, 1, lnum, [integer::], .false., .false.), &
                Node(MODE_SPH_TM_PQ, 1, lnum, [integer::], .false., .false.), &
                Node(MODE_FAR_TE_PQ, 1, spherical_lnum, [1], .false., .false.), &
                Node(MODE_FAR_TM_PQ, 1, spherical_lnum, [2], .false., .false.), &
                Node(MODE_FAR_TETM, 1, spherical_lnum, [3, 4], .false., .false.), &
                Node(MODE_BARBER, m, spherical_lnum, [5], .false., .false.), &
                Node(MODE_MISHCH, m, spherical_lnum, [6], .false., .false.) &
                ! Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .false.), &
                ! Node(MODE_FAR_TM_UV, m, spherical_lnum, [8], .false., .false.), &
                ! Node(MODE_BARBER, m, spherical_lnum, [9], .false., .false.), &
                ! Node(MODE_MISHCH, m, spherical_lnum, [10], .false., .false.) &
            ]
            result_num = 7
            m = 1
        else
            queue = [&
                Node(MODE_SPH_TM_UV, m, lnum, [integer::], .true., .false.), &
                Node(MODE_FAR_TM_UV, m, spherical_lnum, [1], .false., .false.), &
                Node(MODE_BARBER, m, spherical_lnum, [2], .false., .false.), &
                Node(MODE_MISHCH, m, spherical_lnum, [3], .false., .false.) &
            ]
            result_num = 4
        endif
        
        mode_res = calculate_m(global_context, queue, m, lnum, spherical_lnum)

        ! call log_mode_factors('Q SPH_TM', mode_res(1)%factors)
        res(1:lnum_in, 1:lnum_in) = transpose(mode_res(result_num)%tmatrix(1:lnum_in, 1:lnum_in))
        res(lnum_in + 1:2 *lnum_in, lnum_in + 1:2 *lnum_in) = &
            transpose(mode_res(result_num)%tmatrix(lnum + 1:lnum + lnum_in, lnum + 1:lnum + lnum_in))
        res(1:lnum_in, lnum_in + 1:2 *lnum_in) = -transpose(mode_res(result_num)%tmatrix(lnum + 1:lnum + lnum_in, 1:lnum_in))
        res(lnum_in + 1:2 *lnum_in, 1:lnum_in) = -transpose(mode_res(result_num)%tmatrix(1:lnum_in, lnum + 1:lnum + lnum_in))
        ! res = transpose(mode_res(result_num)%tmatrix(1:2*lnum_in, 1:2*lnum_in))
        ! write(*,*) 'inside te result = ', mode_res(1)%tmatrix(:,1)
        ! write(*,*) 'inside result = ', res(:,1)
        do i = 1, size(queue)
            if (allocated(queue(i)%previous)) deallocate(queue(i)%previous)
            call delete_mode_calculation_result(mode_res(i))
        enddo
        deallocate(queue)
        deallocate(mode_res)
        ! deallocate(queue, mode_res)
    end subroutine GetSpheroidalTmatrix
end module tmatrix_extraction
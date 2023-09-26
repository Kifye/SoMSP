! Created by drakosha on 13.07.2021.

module spheroidal_tmatrix
    use regime
    use wavelength_point
    use spheroidal_scatterer
    use integrals
    use spheroidal
    use matrix
    use contexts
    use logging
!    use spherical
    implicit none

contains
    ! formulas (23)-(26) AppliedOptics 2012
    ! A_ik = W(left_R Delta - mu_j/mu_{j+1}Delta right_R - (mu_j/mu_{j+1}-1)ksi/(ksi^20f)Delta)P
    function get_part(f, ksi, mu, left_R, right_R, W, P, &
            Delta, matrix_size) result(res)
        complex(knd) :: mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f, i
        complex(knd) :: left_R(:), right_R(:), W(:), val
        type(NoZeroesMatrix) :: P, Delta, tmp, x, res

        tmp = Delta
        ! call log_nzmatrix('tmp', tmp)
        ! call nzmultiply_by_diag_left(tmp, left_R)
        x = left_R .dm. tmp

        tmp = Delta
        tmp = tmp .dm. right_R
        ! call nzmultiply_by_diag_right(tmp, right_R)
        x = x - mu(0) / mu(1) * tmp - (mu(0) / mu(1) - 1.0_knd)*ksi / (ksi**2 - f) * Delta
        x = W .dm. x
        ! call nzmultiply_by_diag_left(x, W)

        res = x * P
        ! call log_nzmatrix('res', res)
    end function get_part

    subroutine calculate_tmatrix_sph_pq(scatterering_context, computation_context, &
            is_te_mode, tmatrix)

        type(ScatteringContext), intent(in) :: scatterering_context
        type(SpheroidalContext), intent(in) :: computation_context
        logical, intent(in) :: is_te_mode
        type(NoZeroesMatrix), intent(out) :: tmatrix(1,1)

        integer :: lnum, nol
        integer :: n, i, l, j
        ! by layer
        complex(knd) :: mu(0:computation_context%nol)
        ! by expansion
        complex(knd), dimension(computation_context%lnum * 2) :: R11, R31, R12, W1, R32
        type(NoZeroesMatrix) :: adder, A11, A31, A31inv, P
        type(NoZeroesMatrix) :: big_matr(2, 1)
        type(NoZeroesMatrix) :: tmp(2, 2)

        real(knd), allocatable :: mult_coef(:)

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate tmatrix sph pq'

        nol = computation_context%nol
        lnum = computation_context%lnum

        mu = scatterering_context%calculation_point%mu
        if (.not. is_te_mode) then
            mu = scatterering_context%calculation_point%eps
        endif

        call assert(computation_context%layers(0, nol)%lnum >= 2 * lnum, 'too small calc to pq')
        ! setting up core
        R11 = computation_context%layers(0, nol)%r1d(1:2*lnum) / computation_context%layers(0, nol)%r1(1:2*lnum)
        R31 = computation_context%layers(0, nol)%r3d(1:2*lnum) / computation_context%layers(0, nol)%r3(1:2*lnum)
        R12 = computation_context%layers(1, nol)%r1d(1:2*lnum) / computation_context%layers(1, nol)%r1(1:2*lnum)
        W1 = -1.0_knd / (R31 - R11)
        call get_identity_nzmatrix(P, lnum)

        ! big_matr = 0
        big_matr(2,1) = get_part( &
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(nol), &
                mu(nol - 1:nol), &
                R11, R12, W1, &
                computation_context%layer_contexts(nol)%Pi1, &
                computation_context%layer_contexts(nol)%Delta, &
                lnum)
        big_matr(1,1) = -get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(nol), &
                mu(nol - 1:nol), &
                R31, R12, W1, &
                computation_context%layer_contexts(nol)%Pi1, &
                computation_context%layer_contexts(nol)%Delta, &
                lnum)

        do j = nol - 1, 1, -1
            R11 = computation_context%layers(0, j)%r1d(1:2*lnum) / computation_context%layers(0, j)%r1(1:2*lnum)
            R31 = computation_context%layers(0, j)%r3d(1:2*lnum) / computation_context%layers(0, j)%r3(1:2*lnum)
            R12 = computation_context%layers(1, j)%r1d(1:2*lnum) / computation_context%layers(1, j)%r1(1:2*lnum)
            R32 = computation_context%layers(1, j)%r3d(1:2*lnum) / computation_context%layers(1, j)%r3(1:2*lnum)
            W1 = -1.0_knd / (R31 - R11)
            ! tmp = 0

            tmp(2, 1) = get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(j), &
                mu(j - 1:j), &
                R11, R12, W1, &
                computation_context%layer_contexts(j)%Pi1, &
                computation_context%layer_contexts(j)%Delta, &
                lnum)
            tmp(1, 1) = -get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(j), &
                mu(j - 1:j), &
                R31, R12, W1, &
                computation_context%layer_contexts(j)%Pi1, &
                computation_context%layer_contexts(j)%Delta, &
                lnum)

            tmp(2, 2) = get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(j), &
                mu(j - 1:j), &
                R11, R32, W1, &
                computation_context%layer_contexts(j)%Pi3, &
                computation_context%layer_contexts(j)%Delta, &
                lnum)
            tmp(1, 2) = -get_part(&
                scatterering_context%scatterer%spheroidal_type, &
                scatterering_context%scatterer%ksi(j), &
                mu(j - 1:j), &
                R31, R32, W1, &
                computation_context%layer_contexts(j)%Pi3, &
                computation_context%layer_contexts(j)%Delta, &
                lnum)
            call log_nzmatrix2('big_matr', big_matr)
            call log_nzmatrix2('tmp', tmp)
            big_matr = tmp * big_matr
            call log_nzmatrix2('big_matr2', big_matr)
        enddo

        ! A11 = 0
        ! A31 = 0
        ! A31inv = 0
        A31 = big_matr(1,1)
        A11 = (1.0_knd / computation_context%layers(0,1)%r3(1:2 * lnum)) .dm. big_matr(2,1)
        call log_nzmatrix('a11', a11)
        call log_nzmatrix('a31', a31)
        ! call multiply_by_diag_left(A11, lnum, 1.0_knd / computation_context%layers(0,1)%r3(1:lnum))
        A31inv = nzinverse(A31) .dm. computation_context%layers(0,1)%r1(1:2*lnum)
        ! call quick_inverse_matrix(A31, lnum, A31inv)
        ! call multiply_by_diag_right(A31inv, lnum, computation_context%layers(0,1)%r1(1:lnum))

        tmatrix = -(A11* A31inv)

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate tmatrix sph pq'

    end subroutine calculate_tmatrix_sph_pq

    !  Nonaxisymmetric part

    function get_part_11(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            lc, matrix_size) result(res)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(:), right_multiplier(:)
        type(LayerContext), intent(in) :: lc
        type(NoZeroesMatrix) :: adder, res, identity

        call get_identity_nzmatrix(identity, matrix_size)

        ! tmp2 = -nzmatmul(Q01, Epsilon)
        ! ee = ((eps(1) / eps(0) - mu(0) / mu(1)) * f * ksi / (ksi ** 2 - f))
        ! tmp1 = ee * tmp2
        ! tmp2 = (eps(1) / eps(0) - 1q0) * ksi * (Q01 - 2q0 * ksi**2 * Q01Q11)
        ! res = tmp1 - tmp2
        res = -(eps(1) / eps(0) - mu(0) / mu(1)) * f * ksi / (ksi ** 2 - f) * nzmatmul(lc%Q01, lc%Epsilon) - &
                (eps(1) / eps(0) - 1q0) * ksi * (lc%Q01 - 2q0 * ksi**2 * lc%Q01Q11)
        ! call log_nzmatrix('res', res)
        adder = (mu(0) / mu(1) - 1q0) * ksi**2 * lc%Q01 - mu(0) / mu(1) * lc%Delta
        adder = (adder .dm. right_multiplier)
        ! call nzmultiply_by_diag_right(adder, right_multiplier)
        res = res + adder

        adder = lc%Delta + (eps(1) / eps(0) - 1q0) * ksi**2 * lc%Q01
        ! call nzmultiply_by_diag_left(adder, left_multiplier)
        res = res + (left_multiplier .dm. adder)
        ! call log_nzmatrix('res', res)
    end function get_part_11

    function get_part_12(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            lc, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(:), right_multiplier(:)
        type(LayerContext), intent(in) :: lc
        type(NoZeroesMatrix) :: adder, result, identity

        call get_identity_nzmatrix(identity, matrix_size)

        result = -(eps(1) / eps(0) - mu(0) / mu(1)) * f / (ksi ** 2 - f) * &
                (nzmatmul((ksi**2 * lc%Q01 - lc%Delta), lc%Kappa) + nzmatmul(lc%Delta, lc%Gamma11)) + &
                (eps(1) / eps(0) - 1q0) * f * ksi**2 * 2q0 * nzmatmul(lc%Q01Q11, lc%Gamma11)
                ! call log_nzmatrix('res', result)
        adder = (mu(0) / mu(1) - 1q0) * f * ksi * nzmatmul(lc%Q01, lc%Gamma11)
        ! call log_nzmatrix('res1', adder)
        ! call nzmultiply_by_diag_right(adder, right_multiplier)
        adder = (adder .dm. right_multiplier)
        ! call log_nzmatrix('res1.5', adder)
        result = result + adder
        ! call log_nzmatrix('res2', result)
        adder = (eps(1) / eps(0) - 1q0) * f * ksi * nzmatmul(lc%Q01, lc%Gamma11)
        ! call log_nzmatrix('res3', adder)
        ! call nzmultiply_by_diag_left(adder, left_multiplier)
        result = result + (left_multiplier .dm. adder)
        ! call log_nzmatrix('res', result)
    end function get_part_12

    function get_part_21(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            lc, matrix_size) result(result)

        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(:), right_multiplier(:)
        type(LayerContext), intent(in) :: lc
        type(NoZeroesMatrix) :: result, tmp1, tmp2, adder

        tmp1 = nzmatmul(lc%Q01, lc%Kappa)
        tmp1 = ((eps(1) / eps(0) - mu(0) / mu(1)) * ksi**2 / (ksi ** 2 - f)) * tmp1
        tmp2 = nzmatmul(lc%Q01Q11, lc%Gamma11)
        tmp2 = ((eps(1) / eps(0) - 1q0) * ksi**2 * 2q0) * tmp2
        result = tmp1 - tmp2
        ! result = ((eps(1) / eps(0) - mu(0) / mu(1)) * ksi**2 / (ksi ** 2 - f)) * nzmatmul(Q01, Kappa) - &
        !         ((eps(1) / eps(0) - 1q0) * ksi**2 * 2q0) * nzmatmul(Q01Q11, Gamma11)
        adder = -(mu(0) / mu(1) - 1q0) * ksi * nzmatmul(lc%Q01, lc%Gamma11)
        ! call nzmultiply_by_diag_right(adder, right_multiplier)
        result = result + (adder .dm. right_multiplier)

        adder = -(eps(1) / eps(0) - 1q0) * ksi * nzmatmul(lc%Q01, lc%Gamma11)
        ! call nzmultiply_by_diag_left(adder, left_multiplier)
        result = result + (left_multiplier .dm. adder)

    end function get_part_21

    function get_part_22(f, ksi, eps, mu, left_multiplier, right_multiplier, &
            lc, matrix_size) result(result)
        complex(knd) :: eps(0:1), mu(0:1)
        real(knd) :: ksi
        integer :: matrix_size, f
        complex(knd) :: left_multiplier(:), right_multiplier(:)
        type(LayerContext), intent(in) :: lc
        type(NoZeroesMatrix) :: adder, identity, &
                result

        call get_identity_nzmatrix(identity, matrix_size)

        result = (eps(1) / eps(0) - mu(0) / mu(1)) * ksi / (ksi ** 2 - f) * (f * nzmatmul(lc%Q01, lc%Epsilon) + lc%Delta) + &
                (eps(1) / eps(0) - 1q0) * ksi * (lc%Q01 - 2q0 * ksi**2 * lc%Q01Q11)

        adder = -(mu(0) / mu(1) - 1q0) * ksi**2 * lc%Q01 - lc%Delta
        ! call nzmultiply_by_diag_right(adder, right_multiplier)
        result = result + (adder .dm. right_multiplier)

        adder = eps(1) / eps(0) * lc%Delta - (eps(1) / eps(0) - 1q0) * ksi**2 * lc%Q01
        ! call nzmultiply_by_diag_left(adder, left_multiplier)
        result = result + (left_multiplier .dm. adder)

    end function get_part_22

    subroutine set_full_matrix(f, ksi, eps, mu, first_multiplier, left_multiplier, right_multiplier, last_multiplier, &
            layer_context, matrix_size, result, ddd)

        complex(knd), intent(in) :: eps(0:1), mu(0:1)
        real(knd), intent(in) :: ksi
        integer, intent(in) :: matrix_size, f
        real(knd), intent(in), optional :: ddd
        complex(knd), intent(in) :: first_multiplier(:), left_multiplier(:), right_multiplier(:)
        type(LayerContext), intent(in) :: layer_context
        type(NoZeroesMatrix), intent(in) :: last_multiplier
        type(NoZeroesMatrix), intent(out) :: result(2, 2)
        
        integer :: i, j


        ! result = 0
        result(1,1) = get_part_11(f, ksi, eps, mu, left_multiplier, right_multiplier, layer_context, matrix_size)
        result(1,2) = get_part_12(f, ksi, eps, mu, left_multiplier, right_multiplier, layer_context, matrix_size)
        result(2,1) = get_part_21(f, ksi, eps, mu, left_multiplier, right_multiplier, layer_context, matrix_size)
        result(2,2) = get_part_22(f, ksi, eps, mu, left_multiplier, right_multiplier, layer_context, matrix_size)
        do i = 1,2
            do j = 1,2
                ! if (i == 1 .and. j == 2) then
                !     call log_nzmatrix('tmp1', (first_multiplier .dm. result(i,j)))
                !     call log_nzmatrix('lm', last_multiplier)
                ! endif
                result(i,j) = (first_multiplier .dm. result(i,j)) * last_multiplier
            end do
            if (present(ddd)) then
                result(i, 2) = ddd * result(i, 2)
            end if
        end do
        ! if (present(ddd)) then
        !     result(:, 2) = ddd * result(:, 2)
        ! end if
        ! call log_nzmatrix2('result', result)

    end subroutine set_full_matrix

    subroutine calculate_tmatrix_spheroidal_uv(scattering_context, computation_context, is_te_mode, tmatrix)

        type(ScatteringContext), intent(in) :: scattering_context
        type(SpheroidalContext), intent(in) :: computation_context
        logical, intent(in) :: is_te_mode
        type(NoZeroesMatrix), intent(out) :: tmatrix(2, 2)

        integer :: n, i, l, nol, j, lnum
        real(knd) :: k1, ddd, start, finish, global_start
        complex(knd) :: c1
        complex(knd), dimension(0:scattering_context%scatterer%number_of_layers)  :: mu, eps
        complex(knd), dimension(2*computation_context%lnum) :: R11, R31, R12, R32, W1
        complex(knd), dimension(4 * computation_context%lnum) :: initial_corrector, solution_corrector
        type(NoZeroesMatrix), dimension(2,2) :: A11, A31, A31inv
        type(NoZeroesMatrix) :: big_matr(4,2), tmp(4,4), pip, stmp(2,2)

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate tmatrix sph uv'

        call cpu_time(global_start)
        mu = scattering_context%calculation_point%mu
        eps = scattering_context%calculation_point%eps
        if (.not. is_te_mode) then
            mu = scattering_context%calculation_point%eps
            eps = scattering_context%calculation_point%mu
        endif

        k1 = scattering_context%calculation_point%k
        c1 = scattering_context%scatterer%c0(1)
        nol = scattering_context%scatterer%number_of_layers
        lnum = computation_context%lnum

        ! setting up core
        call assert(computation_context%layers(0, nol)%lnum >= 2 * lnum, 'too small spheroidal calc')
        R11 = computation_context%layers(0, nol)%r1d(1:2*lnum) / computation_context%layers(0, nol)%r1(1:2*lnum)
        R31 = computation_context%layers(0, nol)%r3d(1:2*lnum) / computation_context%layers(0, nol)%r3(1:2*lnum)
        R12 = computation_context%layers(1, nol)%r1d(1:2*lnum) / computation_context%layers(1, nol)%r1(1:2*lnum)
        W1 = -1.0_knd / (R31 - R11)

        ! big_matr = 0
        call cpu_time(start)
        ! A31
        call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(nol), &
                eps(nol - 1: nol), &
                mu(nol - 1: nol), &
                W1, &
                R31, &
                R12, &
                computation_context%layer_contexts(nol)%Pi1, &
                computation_context%layer_contexts(nol),&
                lnum, &
                stmp)
        big_matr(1 : 2,1:2) = stmp
        ! A11
        call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(nol), &
                eps(nol - 1: nol), &
                mu(nol - 1: nol), &
                W1, &
                R11, &
                R12, &
                computation_context%layer_contexts(nol)%Pi1, &
                computation_context%layer_contexts(nol),&
                lnum, &
                stmp)
                big_matr(3: 4 , :) = stmp
        big_matr(1 : 2,:) = -big_matr(1 : 2,:)
        call cpu_time(finish)
        call log_time('tmatrix core', finish - start)    

        call cpu_time(start)
        do j = nol - 1, 1, -1
            R11 = computation_context%layers(0, j)%r1d(1:2*lnum) / computation_context%layers(0, j)%r1(1:2*lnum)
            R31 = computation_context%layers(0, j)%r3d(1:2*lnum) / computation_context%layers(0, j)%r3(1:2*lnum)
            R12 = computation_context%layers(1, j)%r1d(1:2*lnum) / computation_context%layers(1, j)%r1(1:2*lnum)
            R32 = computation_context%layers(1, j)%r3d(1:2*lnum) / computation_context%layers(1, j)%r3(1:2*lnum)
            W1 = -1.0_knd / (R31 - R11)
            ddd = scattering_context%scatterer%d(j) / scattering_context%scatterer%d(j + 1)

            ! tmp = 0

            ! A31
            call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(j), &
                eps(j - 1:j), &
                mu(j - 1:j), &
                W1, &
                R31, &
                R12, &
                computation_context%layer_contexts(j)%Pi1, &
                computation_context%layer_contexts(j),&
                lnum, &
                tmp(1:2, 1:2), &
                ddd)

            ! A11
            call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(j), &
                eps(j - 1:j), &
                mu(j - 1:j), &
                W1, &
                R11, &
                R12, &
                computation_context%layer_contexts(j)%Pi1, &
                computation_context%layer_contexts(j),&
                lnum, &
                tmp(3: 4, 1:2), &
                ddd)

            ! A33
            call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(j), &
                eps(j - 1:j), &
                mu(j - 1:j), &
                W1, &
                R31, &
                R32, &
                computation_context%layer_contexts(j)%Pi3, &
                computation_context%layer_contexts(j),&
                lnum, &
                tmp(1:2 , 3: 4), &
                ddd)
            ! A13
            call set_full_matrix(&
                scattering_context%scatterer%spheroidal_type, &
                scattering_context%scatterer%ksi(j), &
                eps(j - 1:j), &
                mu(j - 1:j), &
                W1, &
                R11, &
                R32, &
                computation_context%layer_contexts(j)%Pi3, &
                computation_context%layer_contexts(j),&
                lnum, &
                tmp(3:4,3:4), &
                ddd)
            tmp(1 : 2,:) = -tmp(1 : 2,:)
            big_matr = tmp * big_matr

        enddo
        call cpu_time(finish)
        call log_time('tmatrix layers', finish - start)

        call cpu_time(start)
        ! A11 = 0
        ! A31 = 0
        ! A31inv = 0
        A31 = big_matr(1:2, :)
        A11 = big_matr(3:4, :)
        ! call log_nzmatrix2('a11', a11)
        ! call log_nzmatrix2('a31', a31)
        initial_corrector(1:2*lnum) = 1q0 / (computation_context%layers(0, 1)%r1(1:2*lnum) * k1)
        initial_corrector((2*lnum + 1):(4 * lnum)) = 1q0 / (computation_context%layers(0, 1)%r1(1:2*lnum) * c1)
        solution_corrector(1:2*lnum) = 1q0 / (computation_context%layers(0, 1)%r3(1:2*lnum) * k1)
        solution_corrector((2*lnum + 1):(4 * lnum)) = 1q0 / (computation_context%layers(0, 1)%r3(1:2*lnum) * c1)

        A31 = initial_corrector .dm. A31
        A11 = solution_corrector .dm. A11
        ! call multiply_by_diag_left(A31, 2 * lnum, initial_corrector)
        ! call multiply_by_diag_left(A11, 2 * lnum, solution_corrector)
        call cpu_time(finish)
        call log_time('tmatrix basis correction', finish - start)
        
        ! call log_nzmatrix2('a11', a11)
        ! call log_nzmatrix2('a31', a31)
        call cpu_time(start)
        A31inv = nzinverse2(A31)
        ! call quick_inverse_matrix(A31, 2 * lnum, A31inv)
        call cpu_time(finish)
        call log_time('tmatrix inversion', finish - start)

        call cpu_time(start)
        tmatrix = A11 * A31inv
        call cpu_time(finish)
        call log_time('tmatrix multiplication', finish - start)
        call log_time('tmatrix total', finish - global_start)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate tmatrix sph uv'

    end subroutine calculate_tmatrix_spheroidal_uv


end module spheroidal_tmatrix
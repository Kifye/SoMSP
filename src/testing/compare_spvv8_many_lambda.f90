program compare_spvv8_many_lambda

    use regime
    use scattering_calculation
    use utils
    use constants
    use contexts
    use spheroidal_scatterer
    use spheroidal_indicatrix
    implicit none

    integer :: f, nol, matrix_size, spherical_lnum, minm, maxm, ntheta, nphi, rvind, alpind, riind
    real(knd), allocatable :: rv(:), xv(:), ab(:)
    real(knd) :: alpha, lambda, theta0, theta1, phi0, phi1, prev_ext, prev_sca, start, finish
    complex(knd), allocatable :: ri(:), ris(:)
    ! type(ScatteringQuery) :: query
    type(ScatteringResult):: result, res1, res2
    type(ScatteringContext) :: global_context
    type(SpheroidalShape) :: shape
    character(32) :: model
    character(1024) :: input_file, scatmatr_file, arg, outfile
    integer :: i,j,k
    logical :: need_far
    type(ModeFactors) :: fact, avgfact, diffact
    real(knd), allocatable :: rvs(:)
    real(knd), allocatable :: alphas(:), lambdas(:)
    call cpu_time(start)
    if (LOG_INFO) open(LOG_FD, FILE=LOG_FILENAME, status='replace')
    input_file = 'input.txt'
    scatmatr_file = 'scattering_matrix.txt'
    outfile = 'total_out.txt'
    open(115, FILE=outfile, status='replace')
    do i = 1, command_argument_count(), 2
        call assert(i < command_argument_count(), 'Number of arguments must be even')
        call get_command_argument(i, arg)
        if (trim(arg) == '--input') then
            call get_command_argument(i + 1, input_file) 
        elseif (trim(arg) == '--scat-matr') then
            call get_command_argument(i + 1, scatmatr_file) 
        else
            call assert(.false., 'unknown parameter '//trim(arg))
        endif
    enddo

    call read_input(input_file, f, nol, rv, xv, ab, alpha, lambda, ri, matrix_size, spherical_lnum, minm, maxm, model, &
    ntheta, theta0, theta1, nphi, phi0, phi1)

    ! open(120, file='compout')
    ! prev_ext = 0
    ! prev_sca = 0
    ! do matrix_size = 8, 220, 4
        spherical_lnum = matrix_size * 1.3
        call cpu_time(finish)
    call log_time('prefix', finish - start)
    rvs = [  0.16000q0]
    lambdas = [ 1.002q+02,&
    2.992q+01 ,&
    1.001q+01,&
    3.002q+00   ,&
    1.000q+00   ,&
    2.999q-01   ,&
    9.481q-02   ,&
    4.990q-02  ,&
    4.005q-02   ,&
    2.996q-02 ]
    ris = [cmplx( 3.14758e+00 , 5.05666e+00, knd), &
    cmplx( 1.33952e+00, 1.86188e+00, knd), &
    cmplx(9.91696e-01, 8.17848e-01, knd ), &
    cmplx(7.93805e-01, 3.12897e-01, knd), &
    cmplx(7.73674e-01, 1.96323e-01, knd), &
    cmplx(5.98916e-01, 4.89562e-01, knd), &
    cmplx(2.89623e-02 , 4.35010e-01, knd), &
    cmplx(-9.42650e-02, 1.59667e-01, knd), &
    cmplx( -7.24923e-02, 6.80800e-02, knd), &
    cmplx(-2.73124e-02, 5.36397e-02, knd)]
    alphas = [89.999999q0, 45.0q0, 0.000001q0]
    do rvind = 1, size(rvs)
        do riind = 1, size(ris)
        do alpind = 1, size(alphas)
    lambda = lambdas(riind)
    ri = [cmplx(1.0q0, 0.0q0, knd), ris(riind) + cmplx(1.0q0, 0.0q0, knd)]
    rv(1) = rvs(rvind)
    xv = 2 * PI * rv / lambda
    alpha = alphas(alpind) * PI / 180.0q0
    matrix_size = size_of_matrices(f, xv(1), ab(1), ri(1))
    spherical_lnum = matrix_size * 1.3
    write(*, *) 'Read input:'
    write(*, *) 'f = ', f
    write(*, *) 'rv = ', rv
    write(*, *) 'xv = ', xv
    write(*, *) 'ab = ', ab
    write(*, *) 'alpha = ', alpha, 'radian = ', alpha / PI * 180_knd, 'degrees'
    write(*, *) 'lambda = ', lambda
    write(*, *) 'ri = ', ri
    write(*, *) 'lnum = ', matrix_size
    write(*, *) 'spherical_lnum = ', spherical_lnum
    write(*, *) 'm = ', minm, ':', maxm
    write(*,*) 'model = ', model
    call cpu_time(start)
    call global_context%initialize(f, nol, xv, ab, alpha, lambda, ri, matrix_size, spherical_lnum, minm, maxm, &
    ntheta, theta0, theta1, nphi, phi0, phi1)
    call cpu_time(finish)
    call log_time('initialization', finish - start)
    call cpu_time(start)
    res1 = calculate_indicatrix(global_context, 1,1, model, matrix_size, spherical_lnum, scatmatr_file)
    res2 = calculate_indicatrix(global_context, 1,1, model, matrix_size + 4, spherical_lnum + 4, scatmatr_file)
    if (abs(res1%sph_te%Qext - res2%sph_te%Qext) / abs(res1%sph_te%Qext + res2%sph_te%Qext) > 1q-6) then
        write(*,*) 'inaccurate! diff is', abs(res1%sph_te%Qext - res2%sph_te%Qext) / abs(res1%sph_te%Qext + res2%sph_te%Qext)
    endif
    result = calculate_indicatrix(global_context, minm, maxm, model, matrix_size, spherical_lnum, scatmatr_file)
    call cpu_time(finish)
    call log_time('calculation', finish - start)
    call cpu_time(start)
    call shape%set(f, rv(1), ab(1), alpha)
    ! fact = get_normalized_c_factors_from_q(result%sph_tm, shape)
    ! write(120, '(I5,1x,5E24.12)') matrix_size, fact%Qext, fact%Qsca, abs(fact%Qext - fact%Qsca) / abs(fact%Qext + fact%Qsca), &
    ! abs(fact%Qext - prev_ext) / abs(fact%Qext + prev_ext), abs(fact%Qsca - prev_sca) / abs(fact%Qsca + prev_sca)
    ! write(*, '(I5,1x,5E24.12)') matrix_size, fact%Qext, fact%Qsca, abs(fact%Qext - fact%Qsca) / abs(fact%Qext + fact%Qsca), &
    ! abs(fact%Qext - prev_ext) / abs(fact%Qext + prev_ext), abs(fact%Qsca - prev_sca) / abs(fact%Qsca + prev_sca)
    ! prev_ext = fact%Qext
    ! prev_sca = fact%Qsca
    need_far = (model == 'uv_pq_te_from_tm_with_far')
    call log_mode_factors('Q SPH_TM', result%sph_tm)
    if (need_far) call log_mode_factors('Q FAR_TM', result%far_tm)
    call log_mode_factors('Q SPH_TE', result%sph_te)
    if (need_far) call log_mode_factors('Q FAR_TE', result%far_te)
    call log_mode_factors('C_SPH_TM', get_c_factors_from_q(result%sph_tm, shape))
    if (need_far) call log_mode_factors('C FAR_TM', get_c_factors_from_q(result%far_tm, shape))
    call log_mode_factors('C  SPH_TE', get_c_factors_from_q(result%sph_te, shape))
    if (need_far) call log_mode_factors('C FAR_TE', get_c_factors_from_q(result%far_te, shape))
    call log_mode_factors('C_norm SPH_TM', get_normalized_c_factors_from_q(result%sph_tm, shape))
    if (need_far) call log_mode_factors('C_norm FAR_TM', get_normalized_c_factors_from_q(result%far_tm, shape))
    call log_mode_factors('C norm SPH_TE', get_normalized_c_factors_from_q(result%sph_te, shape))
    if (need_far) call log_mode_factors('C norm FAR_TE', get_normalized_c_factors_from_q(result%far_te, shape))
    avgfact%Qext = (result%sph_te%Qext + result%sph_tm%Qext) / 2.0q0
    avgfact%Qsca = (result%sph_te%Qsca + result%sph_tm%Qsca) / 2.0q0
    avgfact = get_normalized_c_factors_from_q(avgfact, shape)
    call log_mode_factors('C_norm SPH_AVG', avgfact)
    diffact%Qext = (-result%sph_te%Qext + result%sph_tm%Qext) / 2.0q0
    diffact%Qsca = (-result%sph_te%Qsca + result%sph_tm%Qsca) / 2.0q0
    diffact = get_normalized_c_factors_from_q(diffact, shape)
    call log_mode_factors('C_norm SPH_POL', diffact)
    write(*,*) 
    write(115,'(3F10.6, 4E25.15)') rv(1), xv(1), alphas(alpind), avgfact%Qext, diffact%Qext, avgfact%Qsca, avgfact%qabs()
    flush(115)
        enddo
    enddo
    enddo
    ! enddo
    deallocate(rv, xv, ab, ri)
    ! close(120)
    call cpu_time(finish)
    call log_time('postfix', finish - start)
    if (LOG_INFO) close(LOG_FD)
    close(115)
end program compare_spvv8_many_lambda

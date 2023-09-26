module spheroidal_context_module
    use regime
    use constants
    use integrals
    use matrix
    use scattering_context_module
    use logging
    use legendre_functions
    use no_zeroes

    implicit none

    type :: LayerContext
        type(NoZeroesMatrix) :: Delta, Q01, Q11, Q01Q11, Kappa, Gamma11, Epsilon, Pi1, Pi3
    end type LayerContext

    type :: SpheroidalContext
        ! general
        ! fnum - number of function values to calculate
        ! lnum - number of field/potential expansion 
        ! must be fnum >= lnum
        integer :: m, fnum, lnum, nol
        ! functions
        ! layers(0,j) = outside_layer for layer j (with ref_ind related to the outside and c, ksi to this layer's size
        ! layers(1,j) - inside_layer
        ! [0..1][1..nol]
        type(SpheroidalCalculation), allocatable :: layers(:,:)
        integer, private :: maxd
        ! integrals [1..lnum, 1..lnum, 1..nol]
        type(LayerContext), allocatable, dimension(:) :: layer_contexts
        ! state: 0 - not initialized, 1 - only functions are calculated, 2 - for PQ, 3 - for UV
        integer :: state
    contains
        procedure, private :: check_spheroidal_functions, check_base_matrices
        procedure :: initialize => initialize_spheroidal_context
        procedure :: reset => reset_spheroidal_context
        final :: delete_spheroidal_context
    end type SpheroidalContext

contains
    subroutine reset_spheroidal_context(this)
        class(SpheroidalContext), intent(out) :: this

        this%state = 0
    end subroutine reset_spheroidal_context

    subroutine initialize_spheroidal_context(this, m, fnum, lnum, scattering_context, target_state)
        class(SpheroidalContext), intent(inout) :: this
        ! expansion parameters
        integer, intent(in) :: m, fnum, lnum
        ! physical scattering state for functions arguments
        type(ScatteringContext), intent(in) :: scattering_context
        ! what needs to be calculated: 
        !   1 - only functions, 
        !   2 - functions and delta integral matrix (for pq solution)
        !   3 - functions and all integral matrices (for uv solution)
        integer, intent(in) :: target_state

        if (target_state < 1 .or. target_state > 3) then
            write(LOG_FD, *) '{ERROR} invalid target_state = ', target_state
            call exit(1)
        endif

        if (max(this%fnum, fnum) < 2 * lnum) then
            write(LOG_FD, *) '{ERROR} fnum = ', fnum, ' is less than lnum = ', lnum
            call exit(1)
        endif

        ! allocated check is used to see whether the state variable is valid, if it is not, it is equivalent to it being 0
        if (.not. allocated(this%layers) .or. m /= this%m .or. &
            any([fnum, lnum, scattering_context%scatterer%number_of_layers] > [this%fnum, this%lnum, this%nol])) then
            this%state = 0
        endif

        if (this%state >= target_state) return

        ! if (LOG_INFO) write(LOG_FD,*) 'initialize spheroidal context with target state =', target_state

        if (this%state == 0) then
            this%m = m
            this%fnum = fnum
            this%lnum = lnum
            this%nol = scattering_context%scatterer%number_of_layers

            call check_spheroidal_functions(this, scattering_context)

            this%state = 1
        endif

        if (target_state > 1) then
            call check_base_matrices(this, target_state)
        endif

    end subroutine initialize_spheroidal_context

    subroutine check_spheroidal_functions(this, scattering_context)
        class(SpheroidalContext), intent(inout) :: this
        type(ScatteringContext), intent(in) :: scattering_context

        integer :: j
        real(knd) :: start, finish

        if (this%state /= 0) return

        if (allocated(this%layers) .and. (size(this%layers(0,:)) /= this%nol)) then
            deallocate(this%layers)
        end if
        if (.not. allocated(this%layers)) then
            allocate(this%layers(0:1, 1:this%nol))
        end if

        call cpu_time(start)

        ! layer(0, j) outside layer of the j layer of the particle
        ! layer(1, j) - inside layer of j layer
        ! j = 1 - mantle
        ! j = 2 - core
        do j = 1, this%nol
            call this%layers(0,j)%calculate( &
                this%m, &
                this%fnum, &
                scattering_context%scatterer%c0(j) * scattering_context%calculation_point%get_refractive_index(j - 1), &
                scattering_context%scatterer%ksi(j), &
                1, &
                (/ scattering_context%directions%alpha%angle_cos /), &
                scattering_context%scatterer%spheroidal_type &
            )
            call this%layers(1,j)%calculate( &
                this%m, &
                this%fnum, &
                scattering_context%scatterer%c0(j) * scattering_context%calculation_point%get_refractive_index(j), &
                scattering_context%scatterer%ksi(j), &
                1, &
                (/ scattering_context%directions%alpha%angle_cos /), &
                scattering_context%scatterer%spheroidal_type &
            )
        end do

        this%state = 1

        call cpu_time(finish)

        call log_time('context spheroidal functions', finish - start)

    end subroutine check_spheroidal_functions

    subroutine calculate_base_matrices(layer0, layer1, matrix_size, &
        layer_context, mult_coef)
        type(SpheroidalCalculation), intent(in) :: layer0, layer1
        integer, intent(in) :: matrix_size
        type(LayerContext), intent(out) :: layer_context

        type(NoZeroesMatrix) :: tmp, result, identity
        real(knd) :: ksi
        integer :: m, full_size, i, j, k
        real(knd), intent(in) :: mult_coef(:)
        real(knd) :: start, finish

        if (layer0%m /= layer1%m) then
            write(*,*) 'different m in layers!'
            return
        endif

        m = layer1%m
        ksi = layer0%ksi

        call layer_context%Delta%setzeroes(1, matrix_size)
        call layer_context%Q01%setzeroes(1, matrix_size)
        call layer_context%Q11%setzeroes(1, matrix_size)
        call layer_context%Q01Q11%setzeroes(1, matrix_size)
        call layer_context%Kappa%setzeroes(-1, matrix_size)
        call layer_context%Gamma11%setzeroes(-1, matrix_size)
        call layer_context%Epsilon%setzeroes(1, matrix_size)

    !        call calculate_kappa(layer1, layer1, Kappa, matrix_size)
        call double_dep_integral(layer_context%Kappa, m, layer1%legendre, layer1%legendre, &
        mult_coef, kappa_c_lower, kappa_c_upper)
        call log_nzmatrix('kappa', layer_context%Kappa)
    !        call calculate_gamma(layer1, layer1, Gamma11, matrix_size)
        call double_dep_integral(layer_context%Gamma11, m, layer1%legendre, layer1%legendre, &
        mult_coef, gamma_c_lower, gamma_c_upper)
    !        call calculate_epsilon(layer1, layer1, Epsilon, matrix_size)
        call triple_dep_integral(layer_context%Epsilon, m, layer1%legendre, layer1%legendre, &
        mult_coef, epsilon_c_lower, epsilon_c_middle, &
                epsilon_c_upper)

        full_size = min(get_full_matrix_size(matrix_size), min(layer0%lnum, layer1%lnum))

        call tmp%setzeroes(1, full_size)
        call identity%setzeroes(1, full_size)
        call result%setzeroes(1, full_size)

        call get_identity_nzmatrix(identity, full_size)
    !        call calculate_omega(layer1, layer1, tmp, full_size)
        call triple_dep_integral(tmp, m, layer1%legendre, layer1%legendre, mult_coef, omega_c_lower, omega_c_middle, &
                omega_c_upper)
                
        tmp = (ksi**2 - layer0%spheroidal_type) * identity + layer0%spheroidal_type * tmp
        call cpu_time(start)
        result = nzinverse(tmp)
        call cpu_time(finish)
        ! write(*,*) 'inverse time = ', finish - start
        layer_context%Q11 = cut(result, matrix_size)

        call identity%setzeroes(-1, matrix_size)
    !        call calculate_delta(layer0, layer1, identity, full_size)
        call single_dep_integral(identity, m, layer0%legendre, layer1%legendre, mult_coef, delta_coef)
    !        Delta = 0
        layer_context%Delta = cut(identity, matrix_size)
        call cpu_time(start)
        tmp = nzmatmul(identity, result)
        ! do j = 1, full_size
        !     do k = 1, full_size
        !         do i = 1, matrix_size
        !             tmp(i,j) = tmp(i,j) + identity(i,k) * result(k,j)
        !         end do
        !     end do
        ! end do
        call cpu_time(finish)
        ! write(*,*) 'first matmul time = ', finish - start
        layer_context%Q01 = cut(tmp, matrix_size)
        call cpu_time(start)
        result = nzmatmul(tmp, result)
        ! do j = 1, matrix_size
        !     do k = 1, full_size
        !         do i = 1, matrix_size
        !             Q01Q11(i,j) = Q01Q11(i,j) + tmp(i,k) * result(k,j)
        !         end do
        !     end do
        ! end do
        call cpu_time(finish)
        ! write(*,*) 'second matmul time = ', finish - start
        layer_context%Q01Q11 = cut(result, matrix_size)

    end subroutine calculate_base_matrices

    subroutine mult_matr_by_i(a, lnum)
        complex(knd), intent(inout) :: a(lnum, lnum)
        integer, intent(in) :: lnum

        integer :: i, j

        do i = 1, lnum
            do j = 2 - mod(i,2), lnum, 2
                a(i,j) = a(i,j) * IDEG(mod(j-i + 4 * lnum, 4))
            enddo
        enddo
    end subroutine mult_matr_by_i

    subroutine check_base_matrices(this, target_state)
        class(SpheroidalContext), intent(inout) :: this
        integer, intent(in) :: target_state

        real(knd), allocatable :: mult_coef(:)
        integer :: j, nol, lnum, md, i, k
        real(knd) :: start, finish
        real(knd), allocatable :: p_coef(:)

        if (this%state >= target_state) return

        nol = this%nol
        lnum = this%lnum

        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{BEGIN} calculate integrals'
        if (LOG_INFO) write(LOG_FD, *) '{INFO} calculate base matrices of sizes ', lnum, 'x', lnum, 'x', nol

        if (allocated(this%layer_contexts)) then
            if (size(this%layer_contexts) /= nol) deallocate(this%layer_contexts)
        end if

        if (.not. allocated(this%layer_contexts)) then
            allocate(this%layer_contexts(nol))
        end if

        call cpu_time(start)

        do j = 1, this%nol
            md = max(this%layers(0, j)%maxd, this%layers(1, j)%maxd)
            p_coef = 1.0_knd / calculate_legendre_coef(this%m, md + 1)
            call fill_common_multiplier(this%m, md, mult_coef)

            if (target_state == 3) then
                call calculate_base_matrices(this%layers(0, j), this%layers(1, j), lnum, &
                        this%layer_contexts(j), mult_coef)
            endif

            if (this%state < 2) then
                call this%layer_contexts(j)%Pi1%setzeroes(1, lnum)
                call this%layer_contexts(j)%Pi3%setzeroes(1, lnum)

                if (j < nol) then
                    call single_dep_integral(this%layer_contexts(j)%Pi1, this%m, this%layers(1,j)%legendre, &
                    this%layers(0,j + 1)%legendre, &
                            p_coef, delta_coef)
                    call mult_nzmatr_by_i(this%layer_contexts(j)%Pi1)
                    this%layer_contexts(j)%Pi3 = this%layer_contexts(j)%Pi1
                    this%layer_contexts(j)%Pi1 = (this%layers(1,j)%r1(1:2*lnum) .dm. &
                        this%layer_contexts(j)%Pi1) .dm. &
                         (1.0_knd / this%layers(0,j + 1)%r1(1:2*lnum))
                    ! call nzmultiply_by_diag_left(this%layer_contexts(j)%Pi1, this%layers(1,j)%r1(1:2*lnum))
                    ! call nzmultiply_by_diag_right(this%layer_contexts(j)%Pi1, 1.0_knd / this%layers(0,j + 1)%r1(1:2*lnum))
                    this%layer_contexts(j)%Pi3 = this%layers(1,j)%r3(1:2*lnum) .dm. this%layer_contexts(j)%Pi3 .dm. &
                        (1.0_knd / this%layers(0,j + 1)%r3(1:2*lnum))
                    call log_nzmatrix('pi1', this%layer_contexts(j)%Pi1)
                    call log_nzmatrix('pi3', this%layer_contexts(j)%Pi3)
                    ! call nzmultiply_by_diag_left(this%layer_contexts(j)%Pi3, this%layers(1,j)%r3(1:2*lnum))
                    ! call nzmultiply_by_diag_right(this%layer_contexts(j)%Pi3, 1.0_knd / this%layers(0,j + 1)%r3(1:2*lnum))
                else
                    call get_identity_nzmatrix(this%layer_contexts(j)%Pi1, lnum)
                    this%layer_contexts(j)%Pi3 = this%layer_contexts(j)%Pi1
                end if

                if (target_state == 2) then
                    call this%layer_contexts(j)%Delta%setzeroes(1, lnum)
                    call single_dep_integral(this%layer_contexts(j)%Delta, this%m, &
                    this%layers(0, j)%legendre, this%layers(1, j)%legendre, &
                        mult_coef, delta_coef)
                endif
            endif
        end do

        deallocate(mult_coef, p_coef)

        this%state = max(this%state, target_state) 

        call cpu_time(finish)

        call log_time('context base matrices', finish - start)
        if (LOG_BLOCKS) write(LOG_FD, *) '{BLOCK}{END} calculate integrals'

    end subroutine check_base_matrices


    subroutine delete_spheroidal_context(this)
        type(SpheroidalContext), intent(inout) :: this

        if (allocated(this%layers)) then
            deallocate(this%layers)
        endif

        if (allocated(this%layer_contexts)) then
            deallocate(this%layer_contexts)
        endif

        ! if (allocated(this%Q01)) then
        !     deallocate(this%Q01, this%Q11, this%Q01Q11, this%Kappa, this%Gamma11, this%Epsilon, this%Pi1, this%Pi3)
        ! endif

    end subroutine delete_spheroidal_context

end module spheroidal_context_module
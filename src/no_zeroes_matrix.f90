module no_zeroes
    use regime
    ! use utils
    use constants
    use matrix
    implicit none

    type :: NoZeroesMatrix
        integer :: ms
        complex(knd), allocatable, dimension(:,:), private :: first, second
        integer :: half ! 1/-1
    contains
        procedure :: assign => assign_no_zeroes_matrix
        procedure :: setzeroes
        procedure :: log => log_nz
        ! procedure :: to_normal
        ! procedure :: set_part => set_quater_part
    end type NoZeroesMatrix

    interface operator(*)
        procedure :: multiply_by_integer_left
        procedure :: multiply_by_integer_right
        procedure :: multiply_by_real_left
        procedure :: multiply_by_real_right
        procedure :: multiply_by_complex_left
        procedure :: multiply_by_complex_right
        procedure :: nzmultiply_by_linear_left
        procedure :: nzmultiply_by_linear_right
        procedure :: nzmultiply2_by_linear_right
        procedure :: nzmatmul
        procedure :: nzmatmul2
    end interface

    interface operator(.dm.)
        procedure :: nzmultiply_by_diag_left
        procedure :: nzmultiply_by_diag_right
        procedure :: nzmultiply2_by_diag_left
        procedure :: nzmultiply2_by_diag_right
    end interface

    interface to_normal
        module procedure to_normal_matr, to_normal_sing
    end interface

    interface operator(+)
        procedure :: sum_no_zeroes
    end interface

    interface assignment(=)
        procedure :: assign_nz_matrix
    end interface

    interface operator(-)
        procedure :: minus_no_zeroes
        procedure :: opposite_no_zeroes
        procedure :: opposite_nz_matr
    end interface
    private :: log_matrix_here, assert_here
contains
subroutine assert_here(expression, message)
    logical, intent(in) :: expression
    character(*), intent(in) :: message

    if (.not. expression) then
        write(LOG_FD,*) message
        write(*,*) 'assert failed: ', message, ' see detailes in scattering.log'
        close(LOG_FD)
        call backtrace()
        call exit(1)
    endif

end subroutine assert_here
subroutine log_matrix_here(name, matrix, border)
    character(*), intent(in) :: name
    integer, optional, intent(in) :: border
    complex(knd), dimension(:,:), intent(in) :: matrix

    logical :: by_columns

    integer :: border_row, border_column, i

    if (.not. LOG_MATRICES) return

    border_row = size(matrix, 1)
    border_column = size(matrix, 2)
    if (present(border)) then
        if (border > 0) then
        border_row = border
        border_column = border

        end if
    end if

    write(LOG_FD, *) '{MATRIX} ', name, ' [', border_row, 'x', border_column, ']:'
        do i = 1, border_column
            write(LOG_FD, *) matrix(1:border_row, i)
        end do
end subroutine log_matrix_here
    subroutine log_nz(this, name, border)
        class(NoZeroesMatrix) :: this
        character(*), intent(in) :: name
        integer, optional, intent(in) :: border

        call log_matrix_here(name//'_first', this%first, border)
        call log_matrix_here(name//'_second', this%second, border)
    end subroutine log_nz

    subroutine setzeroes(this, half, ms)
        class(NoZeroesMatrix), intent(out) :: this
        integer, intent(in) :: half, ms

        if (allocated(this%first) .and. this%ms /= ms) then
            deallocate(this%first, this%second)
        endif 
        if (.not. allocated(this%first)) then
            allocate(this%first(ms, ms), this%second(ms, ms))
        endif
        this%first = 0
        this%second = 0
        this%half = half
        this%ms = ms
        
    end subroutine setzeroes

    subroutine assign_nz_matrix(this, other)
        type(NoZeroesMatrix), intent(out) :: this
        type(NoZeroesMatrix), intent(in) :: other

        this%ms = other%ms
        this%half = other%half
        this%first = other%first
        this%second = other%second

    end subroutine assign_nz_matrix

    function to_normal_sing(this) result(res)
        type(NoZeroesMatrix), intent(in) :: this
        complex(knd) :: res(2 * this%ms, 2 * this%ms)
        integer :: i, j

        res = 0
        if (this%half == 1) then
            do i = 1, 2 * this%ms, 2
                do j = 1, 2 * this%ms, 2
                    res(i,j) = this%first((i + 1) / 2, (j + 1) / 2)
                enddo
            enddo
            do i = 2, 2 * this%ms, 2
                do j = 2, 2 * this%ms, 2
                    res(i,j) = this%second(i/ 2, j / 2)
                enddo
            enddo
        else
            do i = 1, 2 * this%ms, 2
                do j = 2, 2 * this%ms, 2
                    res(i,j) = this%first((i + 1) / 2, j / 2)
                enddo
            enddo
            do i = 2, 2 * this%ms, 2
                do j = 1, 2 * this%ms, 2
                    res(i,j) = this%second(i/ 2, (j + 1) / 2)
                enddo
            enddo
        endif
    end function to_normal_sing

    function to_normal_matr(this) result(res)
        type(NoZeroesMatrix), intent(in) :: this(:,:)
        complex(knd), allocatable :: res(:,:)
        integer :: i, j, rows, columns, ms

        rows = size(this, 1)
        columns = size(this, 2)
        ms = this(1,1)%ms
        allocate(res(rows * 2 * ms, columns * 2 * ms))
        do i = 1, rows
            do j = 1, columns
                res((i - 1) * 2*ms + 1:i*2*ms, (j-1)*2*ms + 1:j*2*ms) = to_normal(this(i,j))
            enddo
        enddo
    end function to_normal_matr

    ! subroutine set_quater_part(this, other, left, right)
    !     class(NoZeroesMatrix), intent(inout) :: this
    !     type(NoZeroesMatrix) :: other
    !     integer :: left, right

    !     call assert(this%ms == 2 * other%ms, 'wrong other ms')

    !     if (left == 1 .and. right == 1) then
    !     elseif (left == 2 .and. right == 1) then
    !     elseif (left == 1 .and. right == 2) then
    !     elseif (left == 2 .and. right == 2) then
    !     else
    !         call assert(.false., 'wrong left and right')
    !     endif
        
    ! end subroutine set_quater_part

    subroutine assign_no_zeroes_matrix(this, ms, matrix, half)
        class(NoZeroesMatrix), intent(inout) :: this
        integer, intent(in) :: ms
        complex(knd), intent(in) :: matrix(:,:)
        integer, intent(in) :: half
        integer :: n, l, l0

        this%half = half
        this%ms = ms
        do n = 1, ms
            if (half == 1) then
                l0 = 1
            else
                l0 = 0
            endif
            do l = 1, ms
                this%first(n, l) = matrix(2 * n - 1, 2 * l - l0)
            enddo
        enddo
        do n = 1, ms
            if (half == 1) then
                l0 = 0
            else
                l0 = 1
            endif
            do l = 1, ms
                this%second(n, l) = matrix(2 * n, 2 * l - l0)
            enddo
        enddo
        
    end subroutine assign_no_zeroes_matrix

    function nzmatmul(left, right) result(res)
        type(NoZeroesMatrix), intent(in) :: left, right
        type(NoZeroesMatrix) :: res

        if (left%ms /= right%ms) then
            call assert_here(.false., 'different sizes if nz matrices')
        endif
        res%ms = left%ms
        res%half = left%half * right%half
        if (left%half == 1) then
            res%first = matmul(left%first, right%first)
            res%second = matmul(left%second, right%second)
        else
            res%first = matmul(left%first, right%second)
            res%second = matmul(left%second, right%first)
        endif
    end function nzmatmul

    function nzinverse(matrix) result(res)
        type(NoZeroesMatrix), intent(in) :: matrix
        type(NoZeroesMatrix) :: res

        if (matrix%half /= 1) then
            call assert_here(.false., 'can only inverse nz type 1')
        endif

        res%ms = matrix%ms
        res%half = 1
        
        allocate(res%first(res%ms, res%ms), res%second(res%ms, res%ms))
        call quick_inverse_matrix(matrix%first, matrix%ms, res%first)
        call quick_inverse_matrix(matrix%second, matrix%ms, res%second)
    end function nzinverse

    function nzinverse2(a) result(res)
        type(NoZeroesMatrix), intent(in) :: a(2,2)
        type(NoZeroesMatrix) :: res(2,2)
        type(NoZeroesMatrix) :: a11inv, a22inv, r12, r21

        a11inv = nzinverse(a(1,1))
        a22inv = nzinverse(a(2,2))
        r12 = -a11inv * a(1,2)
        r21 = -a22inv * a(2,1)
        res(1,1) = nzinverse(a(1,1) + a(1,2) * r21)
        res(2,2) = nzinverse(a(2,2) + a(2,1) * r12)
        res(1,2) = r12 * res(2,2)
        res(2,1) = r21 * res(1,1)

    end function nzinverse2

    function nzmatmul2(a, b) result(res)
        type(NoZeroesMatrix), intent(in) :: a(:,:), b(:,:)
        type(NoZeroesMatrix), allocatable :: res(:,:)

        integer :: i, j,k

        call assert_here(size(a,2) == size(b,1), 'bad sizes in nzmult2')
        allocate(res(size(a, 1), size(b,2)))
        do i = 1, size(a,1)
            do j = 1, size(b,2)
                res(i,j) = a(i,1) * b(1,j)
                do k = 2, size(a,2)
                    res(i,j) = res(i,j) + a(i,k) * b(k,j)
                enddo
            enddo
        enddo

    end function nzmatmul2

    function cut(matrix, border) result(res)
        type(NoZeroesMatrix), intent(in) :: matrix
        integer, intent(in) :: border
        type(NoZeroesMatrix) :: res

        res%ms = border
        res%half = matrix%half
        res%first = matrix%first(:border, :border)
        res%second = matrix%second(:border, :border)
    end function cut

    subroutine get_identity_nzmatrix(this, ms)
        type(NoZeroesMatrix), intent(out) :: this
        integer, intent(in) :: ms
        complex(knd) :: first(ms, ms)

        this%ms = ms
        
        call get_identity_matrix(first, ms)
        this%first = first
        this%second = first
        this%half = 1

    end subroutine get_identity_nzmatrix

    function multiply_by_integer_left(a, b) result(c)
        integer, intent(in) :: a
        type(NoZeroesMatrix), intent(in) :: b
        type(NoZeroesMatrix) :: c

        c = b
        c%first = c%first * a
        c%second = c%second * a
    end function multiply_by_integer_left

    function multiply_by_integer_right(b, a) result(c)
        integer, intent(in) :: a
        type(NoZeroesMatrix), intent(in) :: b
        type(NoZeroesMatrix) :: c

        c = multiply_by_integer_left(a, b)
    end function multiply_by_integer_right

    function multiply_by_real_left(a, b) result(c)
        real(knd), intent(in) :: a
        type(NoZeroesMatrix), intent(in) :: b
        type(NoZeroesMatrix) :: c

        c = b
        c%first = c%first * a
        c%second = c%second * a
    end function multiply_by_real_left

    function multiply_by_real_right(b, a) result(c)
        real(knd), intent(in) :: a
        type(NoZeroesMatrix), intent(in) :: b
        type(NoZeroesMatrix) :: c

        c = multiply_by_real_left(a, b)
    end function multiply_by_real_right

    function multiply_by_complex_left(a, b) result(c)
        complex(knd), intent(in) :: a
        type(NoZeroesMatrix), intent(in) :: b
        type(NoZeroesMatrix) :: c

        c = b
        c%first = c%first * a
        c%second = c%second * a
    end function multiply_by_complex_left

    function multiply_by_complex_right(b, a) result(c)
        complex(knd), intent(in) :: a
        type(NoZeroesMatrix), intent(in) :: b
        type(NoZeroesMatrix) :: c

        c = multiply_by_complex_left(a, b)
    end function multiply_by_complex_right

    function sum_no_zeroes(a, b) result(c)
        type(NoZeroesMatrix), intent(in) :: a, b
        type(NoZeroesMatrix) :: c

        call assert_here(a%half == b%half, 'different half in sum')

        c = a
        c%first = c%first + b%first
        c%second = c%second + b%second

    end function sum_no_zeroes

    function minus_no_zeroes(a, b) result(c)
        type(NoZeroesMatrix), intent(in) :: a, b
        type(NoZeroesMatrix) :: c

        call assert_here(a%half == b%half, 'different half in sum')

        c = a
        c%first = c%first - b%first
        c%second = c%second - b%second

    end function minus_no_zeroes

    function opposite_no_zeroes(a) result(c)
        type(NoZeroesMatrix), intent(in) :: a
        type(NoZeroesMatrix) :: c

        c%first = -a%first
        c%second = -a%second
        c%ms = a%ms
        c%half = a%half

    end function opposite_no_zeroes

    function opposite_nz_matr(a) result(c)
        type(NoZeroesMatrix), intent(in) :: a(:,:)
        type(NoZeroesMatrix), allocatable :: c(:,:)
        integer :: i,j

        c = a
        do i = 1, size(c(:,1))
            do j = 1, size(c(1,:))
                c(i,j) = -c(i,j)
            enddo
        enddo

    end function opposite_nz_matr

    function nzmultiply_by_diag_left(b, a) result(c)
        type(NoZeroesMatrix), intent(in) :: a
        integer :: ms
        complex(knd), intent(in) :: b(:)
        type(NoZeroesMatrix) :: c
        complex(knd) :: b1(a%ms), b2(a%ms)
        integer :: i

        ms = size(b) / 2
        call assert_here(a%ms == ms, 'wrong ms in multiply_by_diag_left')

        b1 = [(b(i), i = 1, 2 * ms, 2)]
        b2 = [(b(i), i = 2, 2 * ms, 2)]
        c = a
        if (a%half == 1) then
            call multiply_by_diag_left(c%first, ms, b1)
            call multiply_by_diag_left(c%second, ms, b2)
        else
            call multiply_by_diag_left(c%first, ms, b1)
            call multiply_by_diag_left(c%second, ms, b2)
        endif
    end function nzmultiply_by_diag_left

    function nzmultiply_by_diag_right(a,b) result(c)
        type(NoZeroesMatrix), intent(in) :: a
        integer :: ms
        complex(knd), intent(in) :: b(:)
        type(NoZeroesMatrix) :: c
        complex(knd) :: b1(a%ms), b2(a%ms)
        integer :: i

        ms = size(b) / 2
        call assert_here(a%ms == ms, 'wrong ms in multiply_by_diag_right')

        b1 = [(b(i), i = 1, 2 * ms, 2)]
        b2 = [(b(i), i = 2, 2 * ms, 2)]
        c = a
        if (a%half == 1) then
            call multiply_by_diag_right(c%first, ms, b1)
            call multiply_by_diag_right(c%second, ms, b2)
        else
            call multiply_by_diag_right(c%first, ms, b2)
            call multiply_by_diag_right(c%second, ms, b1)
        endif
    end function nzmultiply_by_diag_right

    function nzmultiply_by_linear_left(b, a) result(c)
        type(NoZeroesMatrix), intent(in) :: a
        integer :: ms
        complex(knd), intent(in) :: b(:)
        complex(knd) :: c(2 * a%ms)
        complex(knd) :: tmp(2, a%ms)
        integer :: i

        ms = size(b) / 2
        call assert_here(a%ms == ms, 'wrong ms in multiply_by_diag_left')

        tmp(1,:) = matmul([(b(i), i = 1, 2 * a%ms, 2)], a%first)
        tmp(2,:)  = matmul([(b(i), i = 2, 2 * a%ms, 2)], a%second)
        c = reshape(tmp, (/2 * a%ms/))
    end function nzmultiply_by_linear_left

    function nzmultiply_by_linear_right(a, b) result(c)
        type(NoZeroesMatrix), intent(in) :: a
        integer :: ms
        complex(knd), intent(in) :: b(:)
        complex(knd) :: c(2 * a%ms)
        complex(knd) :: tmp(2, a%ms)
        integer :: i

        ms = size(b) / 2
        call assert_here(a%ms == ms, 'wrong ms in multiply_by_diag_left')

        if (a%half == 1) then
            tmp(1,:) = matmul(a%first, [(b(i), i = 1, 2 * a%ms, 2)])
            tmp(2,:)  = matmul(a%second, [(b(i), i = 2, 2 * a%ms, 2)])
        else
            tmp(1,:) = matmul(a%first, [(b(i), i = 2, 2 * a%ms, 2)])
            tmp(2,:)  = matmul(a%second, [(b(i), i = 1, 2 * a%ms, 2)])
        endif
        c = reshape(tmp, (/2 * a%ms/))
    end function nzmultiply_by_linear_right

    function nzmultiply2_by_linear_right(a,b) result(c)
        type(NoZeroesMatrix), intent(in) :: a(:,:)
        integer :: ms
        complex(knd), intent(in) :: b(:)
        complex(knd), allocatable :: c(:)
        complex(knd), allocatable :: tmp(:,:)

        integer :: i, rows, columns, j
        rows = size(a,1)
        columns = size(a,2)
        ms = a(1,1)%ms
        call assert_here(size(b) >= 2 * ms * columns, 'small b in nzmult2_by_linear')
        allocate(tmp(2 * ms, rows))
        do i = 1, rows
            tmp(:,i) = 0
            do j = 1, columns
                tmp(:,i) = tmp(:,i) + a(i,j) * b(2 * ms * (j - 1) + 1 : 2 * ms * j)
            enddo
        enddo

        c = reshape(tmp, (/2 * ms * rows/))
    end function nzmultiply2_by_linear_right

    function nzmultiply2_by_diag_left(b, a) result(c)
        type(NoZeroesMatrix), intent(in) :: a(2,2)
        integer :: ms
        complex(knd), intent(in) :: b(:)
        type(NoZeroesMatrix) :: c(2,2)
        integer :: i,j

        call assert_here(a(1,1)%ms * 4 == size(b), 'wrong ms in multiply_by_diag_left')
        ms = a(1,1)%ms
        do i = 1, 2
            do j = 1, 2
                c(i,j) = b((i - 1) * 2 * ms + 1:i * 2 * ms) .dm. a(i,j)
            enddo
        enddo
        
    end function nzmultiply2_by_diag_left
    function nzmultiply2_by_diag_right(a,b) result(c)
        type(NoZeroesMatrix), intent(in) :: a(2,2)
        integer :: ms
        complex(knd), intent(in) :: b(:)
        type(NoZeroesMatrix) :: c(2,2)
        integer :: i,j,k

        call assert_here(a(1,1)%ms * 4 == size(b), 'wrong ms in multiply_by_diag_left')
        ms = a(1,1)%ms
        do i = 1, 2
            do j = 1, 2
                c(i,j) = a(i,j) .dm. b((j - 1) * 2 * ms + 1:j * 2 * ms)
            enddo
        enddo
        
    end function nzmultiply2_by_diag_right

    subroutine mult_nzmatr_by_i(a)
        type(NoZeroesMatrix), intent(inout) :: a

        integer :: i, j, left, right, lnum

        lnum = 2 * a%ms
        
        if (a%half == 1) then
            do i = 1, 2 * a%ms
                do j = 2 - mod(i,2), lnum, 2
                    if (mod(i, 2) == 1) then
                        left = (i + 1) / 2
                        right = (j + 1) / 2
                        a%first(left,right) = a%first(left,right) * IDEG(mod(j-i + 4 * lnum, 4))
                    else
                        left = i / 2
                        right = j / 2
                        a%second(left,right) = a%second(left,right) * IDEG(mod(j-i + 4 * lnum, 4))
                    endif
                enddo
            enddo
        else
            do i = 1, lnum
                do j = 2 - mod(i,2), lnum, 2
                    if (mod(i, 2) == 1) then
                        left = (i + 1) / 2
                        right = j / 2
                        a%first(left,right) = a%first(left,right) * IDEG(mod(j-i + 4 * lnum, 4))
                    else
                        left = i / 2
                        right = (j + 1) / 2
                        a%second(left,right) = a%second(left,right) * IDEG(mod(j-i + 4 * lnum, 4))
                    endif
                enddo
            enddo
        endif

    end subroutine mult_nzmatr_by_i
end module no_zeroes
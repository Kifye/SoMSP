module logging
    use regime
    use utils
    use constants

    implicit none
contains
    subroutine log_matrix(name, matrix, bc, border)
        character(*), intent(in) :: name
        integer, optional, intent(in) :: border
        complex(knd), dimension(:,:), intent(in) :: matrix
        logical, optional, intent(in) :: bc

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

        by_columns = .false.
        if (present(bc)) then
            by_columns = bc
        endif

        write(LOG_FD, *) '{MATRIX} ', name, ' [', border_row, 'x', border_column, ']:'
        if (by_columns) then
            do i = 1, border_column
                write(LOG_FD, *) matrix(1:border_row, i)
            end do
        else
            do i = 1, border_row
                write(LOG_FD, *) matrix(i, 1:border_column)
            end do
        end if
    end subroutine log_matrix

    subroutine log_nzmatrix(name, matrix, aznz, border)
        character(*), intent(in) :: name
        integer, optional, intent(in) :: border
        type(NoZeroesMatrix), intent(in) :: matrix
        logical, optional, intent(in) :: aznz
        complex(knd), allocatable :: normal(:,:)

        logical :: az_nz

        integer :: border_row, border_column, i

        if (.not. LOG_MATRICES) return

        border_row = matrix%ms
        border_column = matrix%ms
        if (present(border)) then
            if (border > 0) then
            border_row = border
            border_column = border

            end if
        end if

        az_nz = .false.
        if (present(aznz)) then
            az_nz = aznz
        endif

        if (az_nz) then
            write(LOG_FD, *) '{MATRIX} as nz: ', name, ' [', border_row, 'x', border_column, '] x 2:'
            call matrix%log(name, border)
        else
            normal = to_normal(matrix)
            call log_matrix(name, normal, border=border)
        endif

    end subroutine log_nzmatrix

    subroutine log_nzmatrix2(name, matrix, aznz, border)
        character(*), intent(in) :: name
        integer, optional, intent(in) :: border
        type(NoZeroesMatrix), intent(in) :: matrix(:,:)
        logical, optional, intent(in) :: aznz
        complex(knd), allocatable :: normal(:,:)

        logical :: az_nz

        integer :: border_row, border_column, i
        integer :: j

        if (.not. LOG_MATRICES) return

        border_row = matrix(1,1)%ms
        border_column = matrix(1,1)%ms
        if (present(border)) then
            if (border > 0) then
            border_row = border
            border_column = border

            end if
        end if

        az_nz = .false.
        if (present(aznz)) then
            az_nz = aznz
        endif
        write(LOG_FD, *) '{MATRIX} nz: ', name, ' of shape [', size(matrix, 1), 'x', size(matrix, 2), ']:'
        do i = 1, size(matrix, 1)
            do j = 1, size(matrix, 2)
                call log_nzmatrix(name//'_element_'//to_string(i)//','//to_string(j), matrix(i,j), az_nz, border)    
            enddo
        enddo        

    end subroutine log_nzmatrix2

    subroutine log_time(name, duration)
        character(*), intent(in) :: name
        real(knd), intent(in) :: duration

        if (LOG_TIMES) then
            write(LOG_FD, *) '{TIME} ', name, ': ', duration
        endif
    end subroutine log_time

    subroutine log_node_queue(queue)
        type(Node), intent(in) :: queue(:)

        integer :: qlen, i

        if (.not. LOG_INFO) return
        qlen = size(queue)
        write(LOG_FD, *) 'Queue: ['
    901 format(4x,I3, ' ', A5,' ',A20, ' ', A3, ' ', I4, ' ', A6, ' ', I4)
    902 format(4x,I3, ' ', A5,' ',A20, ' ', A3, ' ', I4, ' ', A6, ' ', I4, ' ', A9, ' ', I3)
    903 format(4x,I3, ' ', A5,' ',A20, ' ', A3, ' ', I4, ' ', A6, ' ', I4, ' ', A9, ' ', I3, ' ', A3, ' ', I3)
        do i = 1, qlen
            if (size(queue(i)%previous) == 0) then
            write(LOG_FD,901) i, 'mode:', queue(i)%info%to_string(), 'm =', queue(i)%item%m, 'lnum =', queue(i)%item%lnum
            elseif (size(queue(i)%previous) == 1) then
                write(LOG_FD,902) i, 'mode:', queue(i)%info%to_string(), 'm =', queue(i)%item%m, 'lnum =', queue(i)%item%lnum, &
                'from node', queue(i)%previous(1)
            elseif (size(queue(i)%previous) == 2) then
                write(LOG_FD,903) 'mode:', queue(i)%info%to_string(), 'm =', queue(i)%item%m, 'lnum =', queue(i)%item%lnum, &
                'from nodes', queue(i)%previous(1), 'and', queue(i)%previous(2)
            else
                call assert(.false., 'too many previous nodes')
            endif
        enddo
        write(LOG_FD, *) ']'
    end subroutine log_node_queue
end module logging
! Created by drakosha on 18.02.2021.

module bessel_functions

    use regime
    use constants

    use complex_prolate_swf
!    use complex_oblate_swf

    implicit none

    type, public :: BesselCalculation
        integer :: lnum
        complex(knd) :: c
        real(knd) :: x
        ! functions and their derivatives for n = 0,.., lnum - 1, arrays [1:lnum]
        real(knd), allocatable, dimension(:) :: sr, sdr, sn, sdn

        logical :: calculated
    contains
        !  initializes, allocates arrays, does not calculate
        procedure :: set
        !  call to calculate all the values, can be long
        procedure :: calculate => calculate_bessel_functions !, calculate_q
        !  destructor, deallocate arrays
        final :: delete_calculation
    end type BesselCalculation

contains

    subroutine set(this, n, c, x)
        class(BesselCalculation) :: this
        integer, intent(in) :: n
        complex(knd), intent(in) :: c
        real(knd), intent(in) :: x

        this%lnum = n
        this%c = c
        this%x = x

        this%calculated = .false.
    end subroutine set

    subroutine calculate_bessel_functions(this)
        class(BesselCalculation), intent(inout) :: this

        integer :: i, j, maxlp, maxj
        complex(vb_knd) :: sbesdf(this%lnum), sbesdr(this%lnum), sbesf(this%lnum), sbesn(this%lnum)
        integer ibese(this%lnum)

        maxlp = this%lnum
        maxj = this%lnum
        if (allocated(this%sr)) then
            deallocate(this%sr, this%sdr, this%sn, this%sdn)
        endif
        allocate(this%sr(0:maxlp), this%sdr(0:maxlp), this%sn(0:maxlp), this%sdn(0:maxlp))

        call sphbes (this%c, this%x, maxj, maxj, maxlp, precision(sbesn(1)), range(sbesn(1)) - 1, sbesf, sbesdf, &
        sbesn, ibese, sbesdr)

        this%sr(0:this%lnum-1) = sbesn * 10q0**ibese
        this%sdr(0:this%lnum-1) = this%sr(0:this%lnum-1)*sbesdr

        call sphneu (this%c, this%x, maxj-2, maxj, maxlp, precision(sbesn(1)), range(sbesn(1)) - 1, maxlp, sbesf, &
        sbesn, ibese, sbesdf, sbesdr)

        this%sn(0:this%lnum-1) = sbesn * 10q0**ibese
        this%sdn(0:this%lnum-1) = this%sn(0:this%lnum-1)*sbesdr

        this%calculated = .true.

    end subroutine calculate_bessel_functions

    subroutine delete_calculation(this)
        type(BesselCalculation), intent(inout) :: this

        if (allocated(this%sr)) then
            deallocate(this%sr, this%sdr)
        endif

    end subroutine delete_calculation

end module bessel_functions
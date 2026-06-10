! Small formatted table/matrix printers.

module matrix_print_mod
    use kind_mod, only: dp
    implicit none
    private

    public :: print_real_matrix
    public :: print_vol_table

contains

    ! Print a labeled real matrix. Labels can refer to assets, models, or any row/column names.
    subroutine print_real_matrix(title, labels, mat, value_format, label_width, value_width)
        character(len=*), intent(in) :: title
        character(len=*), intent(in) :: labels(:)
        real(dp), intent(in) :: mat(:, :)
        character(len=*), intent(in), optional :: value_format
        integer, intent(in), optional :: label_width, value_width
        integer :: i, j, lw, vw
        logical :: scientific
        character(len=32) :: label_fmt, col_fmt, val_fmt
        character (len=*), parameter :: fmt_title_ = "(/,a)"

        if (size(mat, 1) /= size(labels) .or. size(mat, 2) /= size(labels)) then
            error stop "print_real_matrix: label count differs from matrix shape"
        end if
        lw = 16
        vw = 12
        if (present(label_width)) lw = label_width
        if (present(value_width)) vw = value_width
        write(label_fmt, '("(A",I0,")")') lw
        write(col_fmt, '("(1X,A",I0,")")') vw
        scientific = maxval(abs(mat)) < 1.0e-3_dp
        if (present(value_format)) then
            val_fmt = value_format
        else if (scientific) then
            write(val_fmt, '("(1X,ES",I0,".4)")') vw
        else
            write(val_fmt, '("(1X,F",I0,".6)")') vw
        end if
        print fmt_title_, trim(title)
        write(*, label_fmt, advance='no') ""
        do j = 1, size(labels)
            write(*, col_fmt, advance='no') labels(j)(1:min(vw, len_trim(labels(j))))
        end do
        write(*, *)
        do i = 1, size(labels)
            write(*, label_fmt, advance='no') labels(i)
            do j = 1, size(labels)
                write(*, val_fmt, advance='no') mat(i, j)
            end do
            write(*, *)
        end do
    end subroutine print_real_matrix

    ! Print realized variance and annualized volatility by label.
    subroutine print_vol_table(labels, rv, vol_ann)
        character(len=*), intent(in) :: labels(:)
        real(dp), intent(in) :: rv(:), vol_ann(:)
        integer :: i

        if (size(rv) /= size(labels) .or. size(vol_ann) /= size(labels)) then
            error stop "print_vol_table: label and value sizes differ"
        end if
        print '(A)', "Realized volatility by asset"
        print '(A16,1X,A14,1X,A12)', "Asset", "sum_ret2", "vol_ann%"
        do i = 1, size(labels)
            print '(A16,1X,ES14.5,1X,F12.4)', labels(i), rv(i), vol_ann(i)
        end do
    end subroutine print_vol_table

end module matrix_print_mod

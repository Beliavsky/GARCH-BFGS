module model_selection_mod
    use kind_mod, only: dp
    implicit none
    private

    public :: model_summary_order, print_model_selection_counts, print_model_fit_times

contains

    subroutine model_summary_order(aic_wins, bic_wins, order)
        integer, intent(in)  :: aic_wins(:), bic_wins(:)
        integer, intent(out) :: order(:)

        integer :: i, j, candidate

        if (size(order) /= size(aic_wins) .or. size(bic_wins) /= size(aic_wins)) then
            error stop "model_summary_order: inconsistent array sizes"
        end if

        do i = 1, size(order)
            order(i) = i
        end do

        do i = 2, size(order)
            candidate = order(i)
            j = i - 1
            do
                if (j < 1) exit
                if (.not. comes_before(candidate, order(j), aic_wins, bic_wins)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = candidate
        end do
    end subroutine model_summary_order

    subroutine print_model_selection_counts(model_names, param_count, aic_wins, bic_wins, aic_symbols, &
                                            aic_rank_sum, bic_rank_sum, rank_count, sec_per_asset)
        character(len=16), intent(in) :: model_names(:)
        integer, intent(in) :: param_count(:), aic_wins(:), bic_wins(:)
        character(len=256), intent(in) :: aic_symbols(:)
        integer, intent(in), optional :: aic_rank_sum(:), bic_rank_sum(:), rank_count(:)
        real(dp), intent(in), optional :: sec_per_asset(:)

        integer, allocatable :: order(:)
        integer :: iord, imod
        logical :: have_rank, have_sec
        real(dp) :: aic_avg, bic_avg

        if (size(param_count) /= size(model_names) .or. size(aic_wins) /= size(model_names) .or. &
            size(bic_wins) /= size(model_names) .or. size(aic_symbols) /= size(model_names)) then
            error stop "print_model_selection_counts: inconsistent array sizes"
        end if

        have_rank = present(aic_rank_sum) .and. present(bic_rank_sum) .and. present(rank_count)
        have_sec = present(sec_per_asset)
        if (have_rank) then
            if (size(aic_rank_sum) /= size(model_names) .or. size(bic_rank_sum) /= size(model_names) .or. &
                size(rank_count) /= size(model_names)) then
                error stop "print_model_selection_counts: inconsistent rank array sizes"
            end if
        end if
        if (have_sec) then
            if (size(sec_per_asset) /= size(model_names)) then
                error stop "print_model_selection_counts: inconsistent timing array sizes"
            end if
        end if

        allocate(order(size(model_names)))
        call model_summary_order(aic_wins, bic_wins, order)

        write(*,'(/,A)') "Model selection counts:"
        if (have_sec .and. have_rank) then
            write(*,'(A)') "Model            #param sec/asset  AIC_wins  BIC_wins AIC_avg_rank BIC_avg_rank  AIC_symbols"
            write(*,'(A)') repeat("-", 111)
        else if (have_rank) then
            write(*,'(A)') "Model            #param  AIC_wins  BIC_wins AIC_avg_rank BIC_avg_rank  AIC_symbols"
            write(*,'(A)') repeat("-", 100)
        else
            write(*,'(A)') "Model            #param  AIC_wins  BIC_wins  AIC_symbols"
            write(*,'(A)') repeat("-", 72)
        end if

        do iord = 1, size(order)
            imod = order(iord)
            if (have_rank) then
                if (rank_count(imod) > 0) then
                    aic_avg = real(aic_rank_sum(imod), dp) / real(rank_count(imod), dp)
                    bic_avg = real(bic_rank_sum(imod), dp) / real(rank_count(imod), dp)
                else
                    aic_avg = 0.0_dp
                    bic_avg = 0.0_dp
                end if
            else
                aic_avg = 0.0_dp
                bic_avg = 0.0_dp
            end if

            if (have_sec .and. have_rank) then
                write(*,'(A16,I8,F10.4,2I10,2F13.2,2X,A)') trim(model_names(imod)), param_count(imod), &
                    sec_per_asset(imod), aic_wins(imod), bic_wins(imod), aic_avg, bic_avg, trim(aic_symbols(imod))
            else if (have_rank) then
                write(*,'(A16,I8,2I10,2F13.2,2X,A)') trim(model_names(imod)), param_count(imod), &
                    aic_wins(imod), bic_wins(imod), aic_avg, bic_avg, trim(aic_symbols(imod))
            else
                write(*,'(A16,I8,2I10,2X,A)') trim(model_names(imod)), param_count(imod), &
                    aic_wins(imod), bic_wins(imod), trim(aic_symbols(imod))
            end if
        end do
    end subroutine print_model_selection_counts

    subroutine print_model_fit_times(model_names, param_count, fit_seconds)
        character(len=16), intent(in) :: model_names(:)
        integer, intent(in) :: param_count(:)
        real(dp), intent(in) :: fit_seconds(:)

        integer, allocatable :: order(:)
        integer :: i, j, candidate, imod
        real(dp) :: total_seconds, cumul_seconds, cumul_frac

        if (size(param_count) /= size(model_names) .or. size(fit_seconds) /= size(model_names)) then
            error stop "print_model_fit_times: inconsistent array sizes"
        end if

        allocate(order(size(model_names)))
        do i = 1, size(order)
            order(i) = i
        end do

        do i = 2, size(order)
            candidate = order(i)
            j = i - 1
            do
                if (j < 1) exit
                if (fit_seconds(candidate) <= fit_seconds(order(j))) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = candidate
        end do

        write(*,'(/,A)') "Total fitting time by model:"
        write(*,'(A)') "Model            #param  fit_seconds  cumul_time  cumul_frac"
        write(*,'(A)') repeat("-", 62)
        total_seconds = sum(fit_seconds)
        cumul_seconds = 0.0_dp
        do i = 1, size(order)
            imod = order(i)
            cumul_seconds = cumul_seconds + fit_seconds(imod)
            if (total_seconds > 0.0_dp) then
                cumul_frac = cumul_seconds / total_seconds
            else
                cumul_frac = 0.0_dp
            end if
            write(*,'(A16,I8,3F13.4)') trim(model_names(imod)), param_count(imod), &
                fit_seconds(imod), cumul_seconds, cumul_frac
        end do
    end subroutine print_model_fit_times

    logical function comes_before(lhs, rhs, aic_wins, bic_wins)
        integer, intent(in) :: lhs, rhs
        integer, intent(in) :: aic_wins(:), bic_wins(:)

        if (aic_wins(lhs) /= aic_wins(rhs)) then
            comes_before = aic_wins(lhs) > aic_wins(rhs)
        else if (bic_wins(lhs) /= bic_wins(rhs)) then
            comes_before = bic_wins(lhs) > bic_wins(rhs)
        else
            comes_before = lhs < rhs
        end if
    end function comes_before

end module model_selection_mod

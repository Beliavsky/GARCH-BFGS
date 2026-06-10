! Small helpers for standalone main programs.

module program_utils_mod
    use kind_mod, only: dp
    implicit none
    private

    public :: read_integer_arg, elapsed_since

contains

    subroutine read_integer_arg(iarg, value)
        ! Read command argument iarg as a positive default integer.
        integer, intent(in) :: iarg
        integer, intent(out) :: value
        character(len=64) :: arg
        integer :: io

        call get_command_argument(iarg, arg)
        read(arg, *, iostat=io) value
        if (io /= 0 .or. value < 1) error stop "read_integer_arg: expected positive integer"
    end subroutine read_integer_arg

    real(dp) function elapsed_since(t_start)
        ! CPU seconds elapsed since t_start from cpu_time.
        real(dp), intent(in) :: t_start
        real(dp) :: t_now

        call cpu_time(t_now)
        elapsed_since = t_now - t_start
    end function elapsed_since

end module program_utils_mod

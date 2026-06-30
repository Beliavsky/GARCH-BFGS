module strings_mod
    implicit none
    private
    public :: uppercase, split_string, int_to_str

contains

    pure elemental function uppercase(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, code

        out = s
        do i = 1, len(s)
            code = iachar(out(i:i))
            if (code >= iachar('a') .and. code <= iachar('z')) out(i:i) = achar(code - 32)
        end do
    end function uppercase

    subroutine split_string(str, delim, tokens)
        character(len=*), intent(in) :: str
        character(len=*), intent(in) :: delim
        character(:), allocatable, intent(out) :: tokens(:)
        integer :: start, pos, i, count, n

        n = len_trim(str)
        if (n == 0) then
            allocate(character(len=0) :: tokens(1))
            tokens(1) = ""
            return
        end if

        count = 0
        start = 1
        do
            pos = index(str(start:), delim)
            if (pos == 0) then
                count = count + 1
                exit
            end if
            count = count + 1
            start = start + pos
        end do

        allocate(character(len=n) :: tokens(count))
        start = 1
        i = 1
        do
            pos = index(str(start:), delim)
            if (pos == 0) then
                tokens(i) = adjustl(str(start:))
                exit
            end if
            tokens(i) = adjustl(str(start:start+pos-2))
            start = start + pos
            i = i + 1
        end do
    end subroutine split_string


    function int_to_str(i) result(s)
        ! Convert integer to trimmed string with no leading blanks.
        integer, intent(in) :: i
        character(len=12) :: s
        write(s,'(I0)') i
    end function int_to_str

end module strings_mod

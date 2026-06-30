! Utility to collect input filenames from command-line arguments, a glob pattern, or a default list.

module input_files_mod
    use glob_mod, only: glob, MAX_PATH_LEN
    implicit none
    private
    public :: collect_input_filenames
    public :: MAX_PATH_LEN

contains

    ! Collect input file paths in priority order: command-line arguments, glob pattern, default list.
    !
    ! Command-line arguments containing '*' or '?' are expanded via glob; others are used as-is.
    ! With no command-line arguments: expands file_pattern (if present and non-empty) via glob;
    ! otherwise falls back to default_filenames (if present).  Errors if no source yields files.
    subroutine collect_input_filenames(filenames, file_pattern, default_filenames)
        character(len=MAX_PATH_LEN), allocatable, intent(out) :: filenames(:)
        character(len=*), optional, intent(in) :: file_pattern
        character(len=*), optional, intent(in) :: default_filenames(:)
        character(len=MAX_PATH_LEN) :: arg
        character(len=MAX_PATH_LEN), allocatable :: glob_matches(:), expanded(:)
        integer :: nargs, i

        nargs = command_argument_count()
        if (nargs > 0) then
            allocate(expanded(0))
            do i = 1, nargs
                call get_command_argument(i, arg)
                if (scan(arg, "*?") > 0) then
                    print '(A,A)', "Input file pattern: ", trim(arg)
                    call glob(trim(arg), glob_matches)
                    if (size(glob_matches) > 0) then
                        expanded = [character(len=MAX_PATH_LEN) :: expanded, glob_matches]
                    else
                        expanded = [character(len=MAX_PATH_LEN) :: expanded, arg]
                    end if
                else
                    expanded = [character(len=MAX_PATH_LEN) :: expanded, arg]
                end if
            end do
            filenames = expanded
        else if (present(file_pattern) .and. len_trim(file_pattern) > 0) then
            call glob(trim(file_pattern), filenames)
            print '(A,A)', "Input file pattern: ", trim(file_pattern)
        else if (present(default_filenames)) then
            if (size(default_filenames) < 1) error stop "collect_input_filenames: empty default_filenames"
            allocate(filenames(size(default_filenames)))
            do i = 1, size(default_filenames)
                filenames(i) = default_filenames(i)
            end do
        else
            error stop "collect_input_filenames: no command-line arguments, file pattern, or default filenames"
        end if
        print '(A,I0)', "# input files: ", size(filenames)
        print '(A,*(1X,A))', "Input files:", (trim(filenames(i)), i = 1, size(filenames))
    end subroutine collect_input_filenames

end module input_files_mod

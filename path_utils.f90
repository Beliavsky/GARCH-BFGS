! Small path-name helpers used by standalone programs.

module path_utils_mod
    use strings_mod, only: uppercase
    implicit none
    private
    public :: basename
    public :: dirname
    public :: basename_without_extension
    public :: basename_with_extension
    public :: has_extension
    public :: is_absolute_path
    public :: resolve_filename
    public :: files_with_extension_in_dir
    public :: csv_files_in_dir

contains

    ! Return the directory component of path including the trailing separator.
    ! Returns "." when path has no directory component.
    pure elemental function dirname(path) result(dir)
        character(len=*), intent(in) :: path
        character(len=512) :: dir
        integer :: slash1, slash2, slash_pos

        slash1 = scan(trim(path), "\", back=.true.)
        slash2 = scan(trim(path), "/", back=.true.)
        slash_pos = max(slash1, slash2)
        if (slash_pos > 0) then
            dir = path(1:slash_pos)
        else
            dir = "."
        end if
    end function dirname

    ! Return the final path component after either Windows or POSIX separators.
    pure elemental function basename(path) result(name)
        character(len=*), intent(in) :: path
        character(len=512) :: name
        integer :: slash1, slash2, slash_pos

        slash1 = scan(trim(path), "\", back=.true.)
        slash2 = scan(trim(path), "/", back=.true.)
        slash_pos = max(slash1, slash2)
        if (slash_pos > 0) then
            name = path(slash_pos + 1:)
        else
            name = path
        end if
    end function basename

    ! Return basename(path) with its final suffix removed.
    pure elemental function basename_without_extension(path) result(name)
        character(len=*), intent(in) :: path
        character(len=512) :: name
        character(len=512) :: base
        integer :: dot_pos

        base = basename(path)
        dot_pos = scan(trim(base), ".", back=.true.)
        if (dot_pos > 1) then
            name = base(:dot_pos - 1)
        else
            name = base
        end if
    end function basename_without_extension

    ! Return basename(path) with its final suffix replaced by extension.
    pure elemental function basename_with_extension(path, extension) result(name)
        character(len=*), intent(in) :: path, extension
        character(len=512) :: name
        character(len=512) :: base
        integer :: dot_pos

        base = basename(path)
        dot_pos = scan(trim(base), ".", back=.true.)
        if (dot_pos > 1) then
            name = trim(base(:dot_pos - 1)) // extension
        else
            name = trim(base) // extension
        end if
    end function basename_with_extension

    ! Return true when path ends in extension, comparing case-insensitively.
    pure elemental logical function has_extension(path, extension)
        character(len=*), intent(in) :: path, extension
        integer :: path_len, ext_len

        path_len = len_trim(path)
        ext_len = len_trim(extension)
        if (ext_len < 1 .or. path_len < ext_len) then
            has_extension = .false.
            return
        end if
        has_extension = uppercase(path(path_len - ext_len + 1:path_len)) == uppercase(extension(:ext_len))
    end function has_extension

    ! Return true for Windows drive/UNC paths or POSIX absolute paths.
    pure elemental logical function is_absolute_path(path)
        character(len=*), intent(in) :: path

        is_absolute_path = .false.
        if (len_trim(path) >= 1) then
            if (path(1:1) == "\" .or. path(1:1) == "/") is_absolute_path = .true.
        end if
        if (len_trim(path) >= 3) then
            if (path(2:2) == ":" .and. (path(3:3) == "\" .or. path(3:3) == "/")) is_absolute_path = .true.
        end if
    end function is_absolute_path

    ! Add data_dir to filename unless filename is already absolute.
    pure elemental function resolve_filename(filename, data_dir) result(path)
        character(len=*), intent(in) :: filename, data_dir
        character(len=512) :: path

        if (len_trim(data_dir) == 0 .or. is_absolute_path(filename)) then
            path = filename
        else if (data_dir(len_trim(data_dir):len_trim(data_dir)) == "\" .or. &
                 data_dir(len_trim(data_dir):len_trim(data_dir)) == "/") then
            path = trim(data_dir) // trim(filename)
        else
            path = trim(data_dir) // "\" // trim(filename)
        end if
    end function resolve_filename

    ! Return full paths of files with the requested extension in a directory using the host shell.
    subroutine files_with_extension_in_dir(dir, extension, files)
        character(len=*), intent(in) :: dir, extension
        character(len=512), allocatable, intent(out) :: files(:)
        character(len=512), allocatable :: tmp(:)
        character(len=512) :: list_file, line, pattern
        character(len=2048) :: command
        integer :: unit, io, n, i, exitstat

        if (len_trim(dir) == 0) error stop "files_with_extension_in_dir: empty directory"
        if (len_trim(extension) == 0) error stop "files_with_extension_in_dir: empty extension"
        list_file = "files_with_extension_in_dir.tmp"
        pattern = resolve_filename("*" // trim(extension), dir)
        command = 'cmd /c dir /b /a-d "' // trim(pattern) // '" > "' // trim(list_file) // '" 2> nul'
        call execute_command_line(command, wait=.true., exitstat=exitstat)
        if (exitstat /= 0) then
            allocate(files(0))
            open(newunit=unit, file=list_file, status='old', action='read', iostat=io)
            if (io == 0) close(unit, status='delete')
            return
        end if

        open(newunit=unit, file=list_file, status='old', action='read', iostat=io)
        if (io /= 0) error stop "files_with_extension_in_dir: could not read temporary listing"
        n = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0) exit
            if (len_trim(line) > 0) n = n + 1
        end do
        rewind(unit)
        allocate(tmp(max(n, 1)))
        i = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0) exit
            if (len_trim(line) == 0) cycle
            i = i + 1
            tmp(i) = resolve_filename(trim(line), dir)
        end do
        close(unit, status='delete')
        if (n < 1) then
            allocate(files(0))
        else
            allocate(files(n))
            files = tmp(1:n)
        end if
        deallocate(tmp)
    end subroutine files_with_extension_in_dir

    ! Return full paths of CSV files in a directory using the host shell.
    subroutine csv_files_in_dir(dir, files)
        character(len=*), intent(in) :: dir
        character(len=512), allocatable, intent(out) :: files(:)

        call files_with_extension_in_dir(dir, ".csv", files)
    end subroutine csv_files_in_dir

end module path_utils_mod

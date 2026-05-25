! Small path-name helpers used by standalone programs.

module path_utils_mod
    use strings_mod, only: uppercase
    implicit none
    private
    public :: basename
    public :: basename_with_extension
    public :: has_extension

contains

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

end module path_utils_mod

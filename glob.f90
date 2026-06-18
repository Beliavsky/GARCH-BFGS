module glob_mod
  use iso_c_binding, only: c_int32_t, c_char, c_intptr_t, c_null_char
  implicit none
  private
  public :: glob, MAX_PATH_LEN

  integer, parameter :: MAX_PATH_LEN = 260
  integer, parameter :: MAX_FILES    = 65536

  ! Mirror of WIN32_FIND_DATAA — each FILETIME is two DWORDs to avoid
  ! 8-byte alignment padding that would shift cFileName's offset.
  type, bind(C) :: WIN32_FIND_DATA
    integer(c_int32_t)     :: dwFileAttributes
    integer(c_int32_t)     :: ftCreationTimeLow
    integer(c_int32_t)     :: ftCreationTimeHigh
    integer(c_int32_t)     :: ftLastAccessTimeLow
    integer(c_int32_t)     :: ftLastAccessTimeHigh
    integer(c_int32_t)     :: ftLastWriteTimeLow
    integer(c_int32_t)     :: ftLastWriteTimeHigh
    integer(c_int32_t)     :: nFileSizeHigh
    integer(c_int32_t)     :: nFileSizeLow
    integer(c_int32_t)     :: dwReserved0
    integer(c_int32_t)     :: dwReserved1
    character(kind=c_char) :: cFileName(MAX_PATH_LEN)
    character(kind=c_char) :: cAlternateFileName(14)
  end type WIN32_FIND_DATA

  interface
    function FindFirstFileA(lpFileName, lpFindFileData) &
        bind(C, name='FindFirstFileA')
      import :: c_char, c_intptr_t, WIN32_FIND_DATA
      character(kind=c_char), intent(in)  :: lpFileName(*)
      type(WIN32_FIND_DATA),  intent(out) :: lpFindFileData
      integer(c_intptr_t)                 :: FindFirstFileA
    end function FindFirstFileA

    function FindNextFileA(hFindFile, lpFindFileData) &
        bind(C, name='FindNextFileA')
      import :: c_intptr_t, c_int32_t, WIN32_FIND_DATA
      integer(c_intptr_t), value, intent(in) :: hFindFile
      type(WIN32_FIND_DATA),      intent(out) :: lpFindFileData
      integer(c_int32_t)                      :: FindNextFileA
    end function FindNextFileA

    function FindClose(hFindFile) &
        bind(C, name='FindClose')
      import :: c_intptr_t, c_int32_t
      integer(c_intptr_t), value, intent(in) :: hFindFile
      integer(c_int32_t)                     :: FindClose
    end function FindClose
  end interface

contains

  ! Collect full paths for all entries matching a Windows glob pattern.
  ! Wildcards (* and ?) are supported only in the final path component,
  ! which is the same restriction the underlying API has.
  !
  ! On return:
  !   files(:)  allocated array of matching paths (caller must deallocate)
  subroutine glob(pattern, files)
    character(len=*),                         intent(in)  :: pattern
    character(len=MAX_PATH_LEN), allocatable, intent(out) :: files(:)

    type(WIN32_FIND_DATA)                     :: fd
    integer(c_intptr_t)                       :: handle
    integer(c_int32_t)                        :: rc
    character(kind=c_char)                    :: cpattern(len_trim(pattern)+1)
    character(len=MAX_PATH_LEN), allocatable  :: buf(:)
    character(len=MAX_PATH_LEN)               :: dir, name
    integer                                   :: i, n, sep, plen, nfiles

    plen = len_trim(pattern)

    do i = 1, plen
      cpattern(i) = pattern(i:i)
    end do
    cpattern(plen+1) = c_null_char

    ! Isolate the directory prefix (up to and including the last separator).
    sep = 0
    do i = plen, 1, -1
      if (pattern(i:i) == '\' .or. pattern(i:i) == '/') then
        sep = i
        exit
      end if
    end do
    dir = ''
    if (sep > 0) dir = pattern(1:sep)

    allocate(buf(MAX_FILES))
    n = 0

    handle = FindFirstFileA(cpattern, fd)
    if (handle == -1_c_intptr_t) then   ! INVALID_HANDLE_VALUE
      nfiles = 0
      allocate(files(0))
      return
    end if

    do
      name = c_to_f_str(fd%cFileName)
      if (trim(name) /= '.' .and. trim(name) /= '..') then
        n = n + 1
        if (n <= MAX_FILES) buf(n) = trim(dir) // trim(name)
      end if
      rc = FindNextFileA(handle, fd)
      if (rc == 0) exit
    end do

    rc = FindClose(handle)

    nfiles = min(n, MAX_FILES)
    allocate(files(nfiles))
    files(1:nfiles) = buf(1:nfiles)
    deallocate(buf)
  end subroutine glob

  pure function c_to_f_str(chars) result(s)
    character(kind=c_char), intent(in) :: chars(:)
    character(len=size(chars))         :: s
    integer :: i
    s = ' '
    do i = 1, size(chars)
      if (chars(i) == c_null_char) exit
      s(i:i) = chars(i)
    end do
  end function c_to_f_str

end module glob_mod

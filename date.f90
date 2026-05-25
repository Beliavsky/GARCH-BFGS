! Minimal date type and utilities.
!
! This module is based on the public-domain DataFrame date_mod implementation,
! with small conversion helpers for this repo's existing YYYYMMDD integer dates.

module date_mod
    implicit none
    private

    public :: date_t, date_time_t
    public :: seconds_per_minute, minutes_per_hour, seconds_per_hour
    public :: valid, date_from_iso, date_from_basic, date_from_yyyymmdd
    public :: yyyymmdd, date_label, day_of_year
    public :: datetime_from_iso
    public :: operator(+), operator(-), operator(==), operator(/=)
    public :: operator(<), operator(<=), operator(>), operator(>=)

    type :: date_t
        integer :: year = 0
        integer :: month = 0
        integer :: day = 0
    contains
        procedure :: to_str
        procedure :: to_yyyymmdd
    end type date_t

    integer, parameter :: seconds_per_minute = 60
    integer, parameter :: minutes_per_hour = 60
    integer, parameter :: seconds_per_hour = minutes_per_hour * seconds_per_minute

    type :: date_time_t
        type(date_t) :: date
        integer :: hour = 0
        integer :: minute = 0
        integer :: second = 0
        integer :: utc_offset_minutes = 0
    contains
        procedure :: seconds_since_midnight
        procedure :: datetime_to_str
        generic :: to_str => datetime_to_str
    end type date_time_t

    interface operator(+)
        module procedure add_days_right
        module procedure add_days_left
    end interface

    interface operator(-)
        module procedure subtract_days
        module procedure difference_days
    end interface

    interface operator(==)
        module procedure eq_date
    end interface

    interface operator(/=)
        module procedure ne_date
    end interface

    interface operator(<)
        module procedure lt_date
    end interface

    interface operator(<=)
        module procedure le_date
    end interface

    interface operator(>)
        module procedure gt_date
    end interface

    interface operator(>=)
        module procedure ge_date
    end interface

contains

    pure function to_str(this) result(s)
        class(date_t), intent(in) :: this
        character(len=10) :: s

        s = zero_pad_4(this%year) // '-' // zero_pad_2(this%month) // '-' // zero_pad_2(this%day)
    end function to_str

    pure integer function to_yyyymmdd(this)
        class(date_t), intent(in) :: this

        to_yyyymmdd = yyyymmdd(this)
    end function to_yyyymmdd

    pure elemental logical function valid(x)
        type(date_t), intent(in) :: x

        valid = .false.
        if (x%month < 1 .or. x%month > 12) return
        if (x%day < 1) return
        if (x%day > days_in_month(x%year, x%month)) return
        valid = .true.
    end function valid

    pure function date_from_iso(s) result(x)
        character(len=*), intent(in) :: s
        type(date_t) :: x
        character(len=len(s)) :: t
        integer :: y, m, d
        logical :: ok1, ok2, ok3

        x = date_t(0, 0, 0)
        t = adjustl(s)
        if (len_trim(t) /= 10) return
        if (t(5:5) /= '-' .or. t(8:8) /= '-') return
        call parse_uint(t(1:4), y, ok1)
        call parse_uint(t(6:7), m, ok2)
        call parse_uint(t(9:10), d, ok3)
        if (.not. (ok1 .and. ok2 .and. ok3)) return
        x = date_t(y, m, d)
        if (.not. valid(x)) x = date_t(0, 0, 0)
    end function date_from_iso

    pure function date_from_basic(s) result(x)
        character(len=*), intent(in) :: s
        type(date_t) :: x
        character(len=len(s)) :: t
        integer :: y, m, d
        logical :: ok1, ok2, ok3

        x = date_t(0, 0, 0)
        t = adjustl(s)
        if (len_trim(t) /= 8) return
        call parse_uint(t(1:4), y, ok1)
        call parse_uint(t(5:6), m, ok2)
        call parse_uint(t(7:8), d, ok3)
        if (.not. (ok1 .and. ok2 .and. ok3)) return
        x = date_t(y, m, d)
        if (.not. valid(x)) x = date_t(0, 0, 0)
    end function date_from_basic

    pure elemental function date_from_yyyymmdd(n) result(x)
        integer, intent(in) :: n
        type(date_t) :: x

        x = date_t(n / 10000, mod(n / 100, 100), mod(n, 100))
        if (.not. valid(x)) x = date_t(0, 0, 0)
    end function date_from_yyyymmdd

    pure elemental integer function yyyymmdd(x)
        type(date_t), intent(in) :: x

        yyyymmdd = 10000*x%year + 100*x%month + x%day
    end function yyyymmdd

    pure function date_label(n) result(label)
        integer, intent(in) :: n
        character(len=10) :: label
        type(date_t) :: x

        x = date_from_yyyymmdd(n)
        label = x%to_str()
    end function date_label

    pure elemental integer function day_of_year(x)
        type(date_t), intent(in) :: x
        integer :: m

        if (.not. valid(x)) then
            day_of_year = 0
            return
        end if
        day_of_year = x%day
        do m = 1, x%month - 1
            day_of_year = day_of_year + days_in_month(x%year, m)
        end do
    end function day_of_year

    function datetime_from_iso(s) result(dt)
        character(len=*), intent(in) :: s
        type(date_time_t) :: dt
        character(len=:), allocatable :: t
        integer :: pos_sign, hh, mm, ss, oh, om
        logical :: ok

        t = adjustl(s)
        if (len_trim(t) < 19) error stop "datetime_from_iso: timestamp too short"
        dt%date = date_from_iso(t(1:10))
        if (.not. valid(dt%date)) error stop "datetime_from_iso: invalid date"
        if (t(11:11) /= " " .and. t(11:11) /= "T") error stop "datetime_from_iso: invalid separator"
        call parse_uint(t(12:13), hh, ok)
        if (.not. ok) error stop "datetime_from_iso: invalid hour"
        call parse_uint(t(15:16), mm, ok)
        if (.not. ok) error stop "datetime_from_iso: invalid minute"
        call parse_uint(t(18:19), ss, ok)
        if (.not. ok) error stop "datetime_from_iso: invalid second"
        dt%hour = hh
        dt%minute = mm
        dt%second = ss
        dt%utc_offset_minutes = 0
        if (len_trim(t) >= 25) then
            pos_sign = 20
            if (t(pos_sign:pos_sign) == "+" .or. t(pos_sign:pos_sign) == "-") then
                call parse_uint(t(pos_sign+1:pos_sign+2), oh, ok)
                if (.not. ok) error stop "datetime_from_iso: invalid UTC offset hour"
                call parse_uint(t(pos_sign+4:pos_sign+5), om, ok)
                if (.not. ok) error stop "datetime_from_iso: invalid UTC offset minute"
                dt%utc_offset_minutes = 60*oh + om
                if (t(pos_sign:pos_sign) == "-") dt%utc_offset_minutes = -dt%utc_offset_minutes
            end if
        end if
        if (dt%hour < 0 .or. dt%hour > 23 .or. dt%minute < 0 .or. dt%minute > 59 .or. &
            dt%second < 0 .or. dt%second > 59) then
            error stop "datetime_from_iso: invalid time"
        end if
    end function datetime_from_iso

    pure integer function seconds_since_midnight(this)
        class(date_time_t), intent(in) :: this

        seconds_since_midnight = seconds_per_hour*this%hour + seconds_per_minute*this%minute + this%second
    end function seconds_since_midnight

    function datetime_to_str(this) result(s)
        class(date_time_t), intent(in) :: this
        character(len=25) :: s
        character(len=1) :: sign_char
        integer :: off_abs

        sign_char = "+"
        if (this%utc_offset_minutes < 0) sign_char = "-"
        off_abs = abs(this%utc_offset_minutes)
        write(s,'(A,1X,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') this%date%to_str(), &
            this%hour, ":", this%minute, ":", this%second, sign_char, off_abs/60, ":", mod(off_abs, 60)
    end function datetime_to_str

    pure elemental integer function days_in_month(year, month)
        integer, intent(in) :: year, month

        select case (month)
        case (1, 3, 5, 7, 8, 10, 12)
            days_in_month = 31
        case (4, 6, 9, 11)
            days_in_month = 30
        case (2)
            if ((mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. mod(year, 400) == 0) then
                days_in_month = 29
            else
                days_in_month = 28
            end if
        case default
            days_in_month = 0
        end select
    end function days_in_month

    pure elemental type(date_t) function add_days_right(x, n)
        type(date_t), intent(in) :: x
        integer, intent(in) :: n

        if (.not. valid(x)) then
            add_days_right = date_t(0, 0, 0)
        else
            add_days_right = from_day_number(day_number(x) + n)
        end if
    end function add_days_right

    pure elemental type(date_t) function add_days_left(n, x)
        integer, intent(in) :: n
        type(date_t), intent(in) :: x

        add_days_left = add_days_right(x, n)
    end function add_days_left

    pure elemental type(date_t) function subtract_days(x, n)
        type(date_t), intent(in) :: x
        integer, intent(in) :: n

        subtract_days = add_days_right(x, -n)
    end function subtract_days

    pure elemental integer function difference_days(x, y)
        type(date_t), intent(in) :: x, y

        if (.not. valid(x) .or. .not. valid(y)) then
            difference_days = 0
        else
            difference_days = day_number(x) - day_number(y)
        end if
    end function difference_days

    pure elemental logical function eq_date(x, y)
        type(date_t), intent(in) :: x, y

        eq_date = x%year == y%year .and. x%month == y%month .and. x%day == y%day
    end function eq_date

    pure elemental logical function ne_date(x, y)
        type(date_t), intent(in) :: x, y

        ne_date = .not. eq_date(x, y)
    end function ne_date

    pure elemental logical function lt_date(x, y)
        type(date_t), intent(in) :: x, y

        lt_date = x%year < y%year .or. &
                  (x%year == y%year .and. (x%month < y%month .or. &
                  (x%month == y%month .and. x%day < y%day)))
    end function lt_date

    pure elemental logical function le_date(x, y)
        type(date_t), intent(in) :: x, y

        le_date = lt_date(x, y) .or. eq_date(x, y)
    end function le_date

    pure elemental logical function gt_date(x, y)
        type(date_t), intent(in) :: x, y

        gt_date = .not. le_date(x, y)
    end function gt_date

    pure elemental logical function ge_date(x, y)
        type(date_t), intent(in) :: x, y

        ge_date = .not. lt_date(x, y)
    end function ge_date

    pure elemental integer function day_number(x)
        type(date_t), intent(in) :: x
        integer :: y, m, d, era, yoe, doy, doe, mp

        y = x%year
        m = x%month
        d = x%day
        if (m <= 2) y = y - 1
        era = floor_div(y, 400)
        yoe = y - era*400
        if (m > 2) then
            mp = m - 3
        else
            mp = m + 9
        end if
        doy = (153*mp + 2)/5 + d - 1
        doe = yoe*365 + yoe/4 - yoe/100 + doy
        day_number = era*146097 + doe - 719468
    end function day_number

    pure elemental type(date_t) function from_day_number(z)
        integer, intent(in) :: z
        integer :: zz, era, doe, yoe, y, doy, mp, m, d

        zz = z + 719468
        era = floor_div(zz, 146097)
        doe = zz - era*146097
        yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365
        y = yoe + era*400
        doy = doe - (365*yoe + yoe/4 - yoe/100)
        mp = (5*doy + 2)/153
        d = doy - (153*mp + 2)/5 + 1
        if (mp < 10) then
            m = mp + 3
        else
            m = mp - 9
        end if
        if (m <= 2) y = y + 1
        from_day_number = date_t(y, m, d)
    end function from_day_number

    pure elemental integer function floor_div(a, b)
        integer, intent(in) :: a, b

        floor_div = a / b
        if (mod(a, b) < 0) floor_div = floor_div - 1
    end function floor_div

    pure function zero_pad_2(n) result(s)
        integer, intent(in) :: n
        character(len=2) :: s

        if (n < 0 .or. n > 99) then
            s = '**'
            return
        end if
        s(1:1) = achar(iachar('0') + n/10)
        s(2:2) = achar(iachar('0') + mod(n, 10))
    end function zero_pad_2

    pure function zero_pad_4(n) result(s)
        integer, intent(in) :: n
        character(len=4) :: s
        integer :: m

        if (n < 0 .or. n > 9999) then
            s = '****'
            return
        end if
        m = n
        s(4:4) = achar(iachar('0') + mod(m, 10))
        m = m / 10
        s(3:3) = achar(iachar('0') + mod(m, 10))
        m = m / 10
        s(2:2) = achar(iachar('0') + mod(m, 10))
        m = m / 10
        s(1:1) = achar(iachar('0') + mod(m, 10))
    end function zero_pad_4

    pure subroutine parse_uint(s, n, ok)
        character(len=*), intent(in) :: s
        integer, intent(out) :: n
        logical, intent(out) :: ok
        integer :: i, m, digit
        character(len=len(s)) :: t

        n = 0
        ok = .false.
        t = adjustl(s)
        m = len_trim(t)
        if (m <= 0) return
        do i = 1, m
            if (t(i:i) < '0' .or. t(i:i) > '9') return
            digit = iachar(t(i:i)) - iachar('0')
            n = 10*n + digit
        end do
        ok = .true.
    end subroutine parse_uint

end module date_mod

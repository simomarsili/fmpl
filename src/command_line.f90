! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module command_line
  use constants, only: long_string
  implicit none
  private
  public :: read_args

contains

  subroutine read_args(data_file,wid,lambda,ignore_pivot,accuracy,scores_format,nerrs)
    use kinds
    character(long_string), intent(out) :: data_file
    real(kflt),             intent(out) :: wid
    real(kflt),             intent(out) :: lambda
    integer,                intent(out) :: ignore_pivot
    integer,                intent(out) :: accuracy
    character(1),           intent(out) :: scores_format
    integer,                intent(out) :: nerrs
    integer :: iarg,nargs
    integer :: err
    character(long_string) :: cmd,arg

    call get_command(cmd)
    nargs = command_argument_count()
    if (nargs == 0) then
       nerrs = -1
       return
    end if

    iarg = 1
    nerrs = 0
    ! set defaults
    data_file = ''
    lambda = 0.01_kflt
    wid = 0.0_kflt
    ignore_pivot = 0
    accuracy = 1
    scores_format = "r"
    ignore_pivot = 0
    do while(iarg <= nargs)
       call get_command_argument(iarg,arg)
       select case(trim(arg))
       case('-h')
          nerrs = -1
          return
       case('-i','--input')
          ! input file
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          data_file = arg
          if(len_trim(data_file) == 0) then
             nerrs = nerrs + 1
             write(0,*) 'error: check data file'
          end if
          if(data_file(1:1) == '-') then
             nerrs = nerrs + 1
             write(0,*) 'error: check data file'
          end if
       case('-w','--reweight')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) wid
          if ( err/=0 ) then
             write(0,*) "error ! check %id (reweight)"
             nerrs = nerrs + 1
          end if
       case('-l','--lambda')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) lambda
          if ( err/=0 ) then
             write(0,*) "error ! check lambda"
             nerrs = nerrs + 1
          end if
       case('-p','--ignore_pivot')
          ! remove contrib. from pivot state  to the final scores
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) ignore_pivot
          if ( err/=0 ) then
             write(0,*) "error ! check -p,--ignore_pivot"
             nerrs = nerrs + 1
          end if
       case('-a','--accuracy')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) accuracy
          if ( err/=0 ) then
             write(0,*) "error ! check accuracy"
             nerrs = nerrs + 1
          end if
       case('--scores_format')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          select case(trim(arg))
          case("rcv")
             scores_format = "r"
          case("array")
             scores_format = "a"
          case("coordinate")
             scores_format = "c"
          case default
             write(0,*) 'error ! possible formats for scores matrix format are {rcv, coordinate, array}'
             nerrs = nerrs + 1             
          end select
       case default
          write(0,'(a,1x,a)') 'error: unknown option',trim(arg)
          nerrs = nerrs + 1             
       end select
       iarg = iarg + 1
    end do

    ! check command line
    if (data_file == '') then
       write(0,*) 'error ! missing input data file'
       nerrs = nerrs + 1
    end if
    if (lambda < 1.e-5_kflt) then
       write(0,*) 'error ! regularization parameter too small (< 1.e-5)'
       nerrs = nerrs + 1
    end if
    if (.not. any([0,1,2] == accuracy)) then
       write(0,*) 'error ! possible accuracy values are 0,1,2'
       nerrs = nerrs + 1
    end if

  end subroutine read_args

end module command_line

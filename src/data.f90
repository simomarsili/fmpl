! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module data
  use kinds
  implicit none
  private

  public :: w

  public :: initialize_data, read_data

  integer :: nd   ! number of data samples
  integer :: nv   ! number of variables
  integer :: ns ! number of states per variable

  real(kflt), allocatable :: w(:)

contains

  subroutine initialize_data(udata,nd,nv,neff)
    use constants, only: long_string
    use parser, only: parser_nfields
    integer, intent(in) :: udata
    integer, intent(out) :: nv,nd
    real(kflt), intent(out) :: neff
    integer :: nfields
    integer :: err
    character(long_string) :: line,newline

    ! read the first string
    read(udata,'(a)',iostat=err) line
    call parser_nfields(line,newline,nfields)
    rewind(udata)

    ! set nv
    nv = nfields

    ! count lines for nd
    nd = 0
    do
       read(udata,'(a)',iostat=err) line
       if(err < 0) exit
       nd = nd + 1
    end do
    rewind(udata)

    ! allocate memory for the data

    allocate(w(nd),stat=err)

    neff = nd

  end subroutine initialize_data
  
  subroutine read_data(udata,w_id,ns,neff,data_samples)
    use constants, only: long_string
    use parser, only: parser_nfields
    integer,    intent(in) :: udata
    real(kflt), intent(in) :: w_id
    integer, intent(out) :: ns
    real(kflt), intent(out) :: neff
    integer, intent(out) :: data_samples(:,:)
    integer :: i
    integer :: nfields
    integer :: err
    character(long_string) :: line,newline

    nv = size(data_samples,dim=1)
    nd = size(data_samples,dim=2)
    do i = 1,nd
       if(mod(i,100) == 0)write(0,*) 'data: ', i
       read(udata,'(a)',iostat=err) line
       call parser_nfields(line,newline,nfields)
       read(newline,*,iostat=err) data_samples(:,i)
       if(err > 0) write(0,*) 'error: reading data'
    end do

    if (any(data_samples==0)) data_samples = data_samples + 1

    ns = maxval(data_samples)
    
    write(0,*) '(nd, nv, ns)', nd, nv, ns
    flush(0)
    
    if(w_id > 1.E-10_kflt) then
       write(0,*) 'computing weights...'
       call data_reweight(data_samples,w_id,neff)
    else
       w = 1.0_kflt / real(nd)
    end if

  end subroutine read_data

  subroutine data_reweight(data_samples,w_id,neff)
    integer, intent(in) :: data_samples(:,:)
    real(kflt), intent(in) :: w_id
    real(kflt), intent(out) :: neff
    integer :: id,jd
    integer :: err
    integer, allocatable :: x(:),y(:)

    nv = size(data_samples,dim=1)
    nd = size(data_samples,dim=2)
    allocate(x(nv),y(nv),stat=err)
    w = 1.0_kflt
    do id = 1,nd-1
       if(mod(id,1000) == 0) then
          write(0,'(f5.1,a)') 100.0 * real(id) / real(nd),' %'
       end if
       x = data_samples(:,id)
       do jd = id+1,nd
          y = data_samples(:,jd)
          if(count(x == y) >= nint(nv * w_id * 0.01_kflt)) then
             w(id) = w(id) + 1.0_kflt
             w(jd) = w(jd) + 1.0_kflt
          end if
       end do
    end do

    w = 1.0_kflt / w
    neff = sum(w)
    w = w / neff
   write(0,*) 'neff: ', neff

  end subroutine data_reweight

end module data

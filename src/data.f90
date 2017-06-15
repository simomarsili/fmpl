! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module data
  use kinds
  implicit none
  private

  public :: initialize_data, read_data, data_reweight

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

  end subroutine initialize_data
  
  subroutine read_data(udata,w_id,ns,neff,data_samples)
    use constants, only: long_string
    use parser, only: parser_nfields
    integer,    intent(in) :: udata
    real(kflt), intent(in) :: w_id
    integer, intent(out) :: ns
    real(kflt), intent(out) :: neff
    integer, intent(out) :: data_samples(:,:)
    integer :: nd,nv
    integer :: i
    integer :: nfields
    integer :: err
    character(long_string) :: line,newline
    real(kflt_single)      :: finish,start,t0,tm

    nv = size(data_samples,dim=1)
    nd = size(data_samples,dim=2)
    call cpu_time(t0)
    do i = 1,nd
       read(udata,'(a)',iostat=err) line
       call parser_nfields(line,newline,nfields)
       read(newline,*,iostat=err) data_samples(:,i)
       call cpu_time(finish)
       if (finish - t0 > 3.0) then
          write(0,'(i6,"/",i6)') i, nd
          t0 = finish
       else if (i == nd) then 
          write(0,'(i6,"/",i6)') i, nd         
       end if
       flush(0)
       if(err > 0) write(0,*) 'error ! reading data'
    end do

    if (any(data_samples==0)) data_samples = data_samples + 1

    ns = maxval(data_samples)
    
  end subroutine read_data

  subroutine data_reweight(data_samples,w_id,neff,w)
    integer, intent(in) :: data_samples(:,:)
    real(kflt), intent(in) :: w_id
    real(kflt), intent(out) :: neff
    real(kflt), intent(out) :: w(:)
    integer :: nd,nv
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

  end subroutine data_reweight

end module data

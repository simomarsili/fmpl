! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module optimize
  use kinds
  use model, only: etot
  ! wrapper to dvmlm subroutine 
  implicit none
  private 

  public :: fit

contains

  subroutine fit(nv,ns,nd,data_samples,w,prm,grd,accuracy,minimizer)
    integer, intent(in) :: nv,ns,nd
    integer, intent(in) :: data_samples(nv,nd)
    real(kflt), intent(in) :: w(nd)
    real(kflt), intent(inout) :: prm(ns+ns*ns*nv)
    real(kflt), intent(inout) :: grd(ns+ns*ns*nv)
    integer, intent(in) :: accuracy
    character(*) :: minimizer

    select case(trim(minimizer))
    case('dvmlm')
       call dvmlm_minimizer(nv,ns,nd,data_samples,w,prm,grd,accuracy)
    case default
       call dvmlm_minimizer(nv,ns,nd,data_samples,w,prm,grd,accuracy)
    end select
    
  end subroutine fit

  subroutine dvmlm_minimizer(nv,ns,nd,data_samples,w,prm,grd,accuracy)
    use model, only: update_gradient
    integer, intent(in) :: nv,ns,nd
    integer, intent(in) :: data_samples(nv,nd)
    real(kflt), intent(in) :: w(nd)
    real(kflt), intent(inout) :: prm(ns+ns*ns*nv)
    real(kflt), intent(inout) :: grd(ns+ns*ns*nv)
    integer, intent(in) :: accuracy
    integer :: niter,neval
    integer :: err
    integer :: ndim,mstep
    character(60) :: task
    integer :: lwa
    real(kflt), allocatable :: wa(:)
    integer       :: isave(5)
    real(kflt) :: dsave(24)
    real(kflt) :: f
    real(kflt) :: frtol,fatol,fmin
    external dvmlm

    ! set prms for minimization
    ! frtol: desired relative error
    ! fatol: desired absolute error
    select case(accuracy)
    case(0)
       ! low accuracy
       frtol = 1.0e-2_kflt
       fatol = 1.0e-3_kflt
    case(1)
       ! moderate accuracy
       frtol = 1.0e-4_kflt
       fatol = 1.0e-6_kflt
    case(2)
       ! high accuracy
       frtol = 1.0e-8_kflt
       fatol = 1.0e-12_kflt
    end select
    
    fmin = -1.9e30_kflt
    
    niter = 0
    neval = 0
    task = 'START'

    ndim = size(prm)
    ! this is the number of steps for hessian approximation
    mstep = 100
    lwa = 2*ndim*mstep + 2*ndim + mstep
    allocate(wa(lwa),stat=err)
    
    call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:))
    do 
       if(neval > 100) then 
          write(0,*) 'warning: neval > 100'
          flush(0)
       end if
       !call dvmlm_min(prm,grd,size(prm),accuracy,ndim,mstep,task,wa)
       f = - etot
       call dvmlm(ndim,prm,f,grd,frtol,fatol,fmin,task,mstep,&
            wa(1),wa(ndim*mstep+1),wa(2*ndim*mstep+1),&
            isave,dsave,wa(2*ndim*mstep+mstep+1),wa(2*ndim*mstep+mstep+ndim+1))
       if(task(1:2) == 'FG') then 
          ! update etot and gradient for line search
          neval = neval + 1
          call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:))
       elseif(task(1:4) == 'NEWX') then
          ! start new line search
          niter = niter + 1
       elseif(task(1:4) == 'WARN') then 
          write(0,*) 'warning ', niter
          flush(0)
       elseif(task(1:4) == 'CONV') then 
          ! compute final values for likelihood
          call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:))
          exit
       end if
    end do

    deallocate(wa)

  end subroutine dvmlm_minimizer

end module optimize
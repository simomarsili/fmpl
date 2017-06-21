! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module optimize
  use kinds
  use model, only: update_gradient
  ! wrapper to dvmlm subroutine 
  implicit none
  private 

  public :: fit

contains

  subroutine fit(nd,nv,ns,data_samples,w,prm,grd,lambda,accuracy,minimizer,ll,ereg,precision,fit_flag)
    integer, intent(in) :: nd,nv,ns
    integer, intent(in) :: data_samples(nv,nd)
    real(kflt), intent(in) :: w(nd)
    real(kflt), intent(inout) :: prm(ns+ns*ns*nv)
    real(kflt), intent(inout) :: grd(ns+ns*ns*nv)
    real(kflt), intent(in) :: lambda
    character(*), intent(in) :: minimizer
    character(12), intent(in) :: fit_flag
    real(kflt), intent(out) :: ll,ereg,precision
    integer, intent(in) :: accuracy

    select case(trim(fit_flag))
    case('optimization')
       select case(trim(minimizer))
       case('dvmlm')
          call dvmlm_minimizer(nv,ns,nd,data_samples,w,prm,grd,lambda,accuracy,ll,ereg,precision)
       case default
          call dvmlm_minimizer(nv,ns,nd,data_samples,w,prm,grd,lambda,accuracy,ll,ereg,precision)
       end select
    case('single-point')
       call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:),lambda,ll,ereg,precision)
    case default
       write(0,'("unknown flag for fit: ",a)') trim(fit_flag)
    end select
    
  end subroutine fit

  subroutine gd_minimizer(nv,ns,nd,data_samples,w,prm,grd,lambda,accuracy,ll,ereg,precision)
    integer, intent(in) :: nv,ns,nd
    integer, intent(in) :: data_samples(nv,nd)
    real(kflt), intent(in) :: w(nd)
    real(kflt), intent(inout) :: prm(ns+ns*ns*nv)
    real(kflt), intent(inout) :: grd(ns+ns*ns*nv)
    real(kflt), intent(in) :: lambda
    real(kflt), intent(out) :: ll,ereg,precision
    integer, intent(in) :: accuracy
    integer :: i

    do i = 1,1000
       call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:),lambda,ll,ereg,precision)
       prm = prm - 0.1*grd
    end do
    call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:),lambda,ll,ereg,precision)
  end subroutine gd_minimizer
  
  subroutine dvmlm_minimizer(nv,ns,nd,data_samples,w,prm,grd,lambda,accuracy,ll,ereg,precision)
    integer, intent(in) :: nv,ns,nd
    integer, intent(in) :: data_samples(nv,nd)
    real(kflt), intent(in) :: w(nd)
    real(kflt), intent(inout) :: prm(ns+ns*ns*nv)
    real(kflt), intent(inout) :: grd(ns+ns*ns*nv)
    real(kflt), intent(in) :: lambda
    integer, intent(in) :: accuracy
    real(kflt), intent(out) :: ll,ereg,precision
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
    logical :: verbose = .false.
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
    
    call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:),lambda,ll,ereg,precision)
    do 
       if(neval > 100) then 
          write(0,*) 'warning: neval > 100'
          flush(0)
       end if
       !call dvmlm_min(prm,grd,size(prm),accuracy,ndim,mstep,task,wa)
       f = - ll + ereg
       call dvmlm(ndim,prm,f,grd,frtol,fatol,fmin,task,mstep,&
            wa(1),wa(ndim*mstep+1),wa(2*ndim*mstep+1),&
            isave,dsave,wa(2*ndim*mstep+mstep+1),wa(2*ndim*mstep+mstep+ndim+1))
       if(task(1:2) == 'FG') then 
          ! update etot and gradient for line search
          neval = neval + 1
          call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:),lambda,ll,ereg,precision)
       elseif(task(1:4) == 'NEWX') then
          ! "A new iterate has been computed. Approximated solution, function value and gradient are available for examination."
          ! Start new line search
          niter = niter + 1
          if (verbose) write(*,*) niter, neval, real(ll), real(ereg), sqrt(sum(grd**2)/size(grd)),maxval(abs(grd))
       elseif(task(1:4) == 'WARN') then 
          write(0,*) 'warning ', niter
          flush(0)
       elseif(task(1:4) == 'CONV') then 
          ! compute final values for likelihood
          call update_gradient(nv,ns,nd,data_samples,w,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:),lambda,ll,ereg,precision)
          if (verbose) write(*,*) neval, real(ll), real(ereg)
          exit
       end if
    end do

    deallocate(wa)

  end subroutine dvmlm_minimizer

end module optimize

! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module dvmlm_wrapper
  use kinds
  use model, only: etot
  ! wrapper to dvmlm subroutine 
  implicit none
  private 

  public :: dvmlm_minimize

contains

  subroutine dvmlm_min(prm,grd,dim,accuracy,ndim,mstep,task,wa)
    integer :: dim
    real(kflt), intent(inout) :: prm(dim)
    real(kflt), intent(inout) :: grd(dim)
    integer,       intent(in) :: accuracy,ndim,mstep
    character(60), intent(inout) :: task
    real(kflt),    intent(in) :: wa(:)
    integer       :: isave(5)
    real(kflt) :: dsave(24)
    real(kflt) :: f
    real(kflt) :: frtol
    real(kflt) :: fatol
    real(kflt) :: fmin
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
    
    f = - etot
    call dvmlm(ndim,prm,f,grd,frtol,fatol,fmin,task,mstep,&
         wa(1),wa(ndim*mstep+1),wa(2*ndim*mstep+1),&
         isave,dsave,wa(2*ndim*mstep+mstep+1),wa(2*ndim*mstep+mstep+ndim+1))

  end subroutine dvmlm_min

  subroutine dvmlm_minimize(nv,ns,nd,dim,data_samples,prm,grd,accuracy,iter,totiter)
    use model, only: update_gradient
    integer, intent(in) :: nv,ns,nd,dim
    integer, intent(in) :: data_samples(nv,ns)
    real(kflt), intent(inout) :: prm(dim)
    real(kflt), intent(inout) :: grd(dim)
    integer, intent(in) :: accuracy
    integer, intent(out) :: iter,totiter
    integer :: err
    integer :: ndim,mstep
    character(60) :: task
    real(kflt), allocatable :: wa(:)
    integer :: lwa
    
    iter = 0
    totiter = 0
    task = 'START'

    ndim = size(prm)
    ! this is the number of steps for hessian approximation
    mstep = 100
    lwa = 2*ndim*mstep + 2*ndim + mstep
    allocate(wa(lwa),stat=err)
    
    call update_gradient(nv,ns,nd,data_samples,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:))
    do 
       if(totiter > 100) then 
          write(0,*) 'warning: totiter > 100'
          flush(0)
       end if
       call dvmlm_min(prm,grd,size(prm),accuracy,ndim,mstep,task,wa)
       if(task(1:2) == 'FG') then 
          ! update etot and gradient for line search
          totiter = totiter + 1
          call update_gradient(nv,ns,nd,data_samples,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:))
       elseif(task(1:4) == 'NEWX') then
          ! start new line search
          iter = iter + 1
       elseif(task(1:4) == 'WARN') then 
          write(0,*) 'warning ', iter
          flush(0)
       elseif(task(1:4) == 'CONV') then 
          ! compute final values for likelihood
          call update_gradient(nv,ns,nd,data_samples,prm(:ns),prm(ns+1:),grd(:ns),grd(ns+1:))
          exit
       end if
    end do

    deallocate(wa)

  end subroutine dvmlm_minimize

end module dvmlm_wrapper

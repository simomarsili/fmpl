! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module model
  use kinds
  implicit none
  private

  public :: initialize_model
  public :: model_reset
  public :: fix_gauge
  public :: update_gradient

  ! index of outcome variable (all the others nv - 1 variables are included in the set of explanatory variables)
  integer :: out_var
  
  ! arrays for fixed residue/sequence
  real(kflt), allocatable, save :: model_f1(:)       ! model single-variable frequencies (ns)
  real(kflt), allocatable, save :: model_f2(:,:,:)   ! model frequencies for pairs of variables (ns x nv x ns)
  real(kflt), allocatable, save :: data_f1(:)        ! data single-variable frequencies (ns)
  real(kflt), allocatable, save :: data_f2(:,:,:)    ! data frequencies for pairs of variables (ns x nv x ns)

  ! regularization parameters 
  real(kflt) :: regularization_strength=0.01_kflt ! regularization strength; the default is l2 with regularization_strength=0.01

contains

  subroutine initialize_model(nv,ns,lambda)
    integer, intent(in) :: nv,ns
    real(kflt), intent(in) :: lambda
    integer :: err

    regularization_strength = lambda

    allocate(data_f1(ns),stat=err)
    allocate(data_f2(ns,ns,nv),stat=err)
    allocate(model_f1(ns),stat=err)
    allocate(model_f2(ns,ns,nv),stat=err)

  end subroutine initialize_model

  subroutine update_data_averages(nd,nv,data_samples,w,err)
    integer, intent(in) :: nd,nv
    integer, intent(in) :: data_samples(:,:)
    real(kflt), intent(in) :: w(nd)
    integer :: id,jv,mys
    integer :: err
    integer, allocatable :: list(:)
    
    data_f1 = 0.0_kflt
    data_f2 = 0.0_kflt

    ! compute variable-specific arrays of frequencies
    allocate(list(nv),stat=err)
    data_f1 = 0.0_kflt
    data_f2 = 0.0_kflt
    do id = 1,nd
       list = data_samples(:,id)
       mys = list(out_var)
       data_f1(mys) = data_f1(mys) + w(id)
       do jv = 1,nv
          if(jv /= out_var) then 
             data_f2(mys,list(jv),jv) = data_f2(mys,list(jv),jv) + w(id)
          end if
       end do
    end do
    deallocate(list)
    
  end subroutine update_data_averages

  subroutine model_reset(nd,nv,iv,data_samples,w,vprm,err) ! couplings
    integer, intent(in) :: nd,nv
    integer, intent(in) :: iv
    integer, intent(in) :: data_samples(:,:)
    real(kflt), intent(in) :: w(nd)
    real(kflt), intent(out) :: vprm(:)
    ! make couplings given out_var 
    ! must be called before looping on data
    integer :: err

    out_var = iv
    vprm = 0.0_kflt

    ! compute variable-specific arrays of frequencies
    call update_data_averages(nd,nv,data_samples,w,err)

  end subroutine model_reset

  subroutine fix_gauge(nv,ns,fields,couplings)
    integer, intent(in) :: nv,ns
    real(kflt), intent(inout) :: fields(ns)
    real(kflt), intent(inout) :: couplings(ns,ns,nv)
    integer :: jv,is,js
    real(kflt) :: mat(ns,ns),arr(ns),marr
    real(kflt) :: rsum(ns),csum(ns),totsum

    arr = fields
    marr = sum(arr) / real(ns)
    arr = arr - marr
    do jv = 1,nv
       mat = couplings(:,:,jv)
       totsum = sum(mat)
       totsum = totsum / real(ns*ns)
       do is = 1,ns
          rsum(is) = sum(mat(is,:))
          csum(is) = sum(mat(:,is))
       end do
       rsum = rsum / real(ns)
       csum = csum / real(ns)
       if(jv /= out_var) arr = arr + csum - totsum
       do js = 1,ns
          do is = 1,ns
             mat(is,js) = mat(is,js) - rsum(is) - csum(js) + totsum
          end do
       end do
       couplings(:,:,jv) = mat
    end do
    fields = arr
    
  end subroutine fix_gauge

  subroutine update_model_averages(nv,ns,nd,data_samples,w,fields,couplings,cond_likelihood)
    integer, intent(in) :: nv,ns,nd
    integer, intent(in) :: data_samples(nv,nd)
    real(kflt), intent(in) :: w(nd)
    real(kflt), intent(in) :: fields(ns)
    real(kflt), intent(in) :: couplings(ns,ns,nv)
    real(kflt), intent(out) :: cond_likelihood
    integer :: list(nv)
    real(kflt) :: conp(ns)
    real(kflt) :: rsum
    integer :: mys
    integer :: id
    real(kflt) :: ww
    integer :: jv,js

    cond_likelihood = 0.0_kflt
    model_f1 = 0.0_kflt
    model_f2 = 0.0_kflt

    ! loop over data
    do id = 1,nd
       list = data_samples(:,id)
       ww = w(id)
       mys = list(out_var)
       
       ! loop over the states of out_var
       conp = fields
       do jv = 1,nv
          if(out_var /= jv) then
             js = list(jv)
             conp = conp + couplings(:,js,jv)
          end if
       end do
       conp = exp(conp)
       
       rsum = sum(conp)
       conp = conp / rsum
       
       cond_likelihood = cond_likelihood + ww * log(conp(mys))
       
       ! update histograms 
       ! loop over the states of out_var
       conp = conp * ww
       model_f1 = model_f1 + conp
       do jv = 1,nv
          if(out_var /= jv) then
             js = list(jv)
             model_f2(:,js,jv) = model_f2(:,js,jv) + conp
          end if
       end do
    end do
    
  end subroutine update_model_averages
  
  subroutine update_gradient(nv,ns,nd,data_samples,w,fields,couplings,grd1,grd2,cond_likelihood,ereg)
    ! update cost-related variables: cond_likelihood, ereg and gradient grd
    integer, intent(in) :: nv,ns,nd
    integer, intent(in) :: data_samples(nv,nd)
    real(kflt), intent(in) :: w(nd)
    real(kflt), intent(in) :: fields(ns)
    real(kflt), intent(in) :: couplings(ns,ns,nv)
    real(kflt), intent(out) :: grd1(ns)
    real(kflt), intent(out) :: grd2(ns,ns,nv)
    real(kflt), intent(out) :: cond_likelihood,ereg

    ! take averages over model distribution
    call update_model_averages(nv,ns,nd,data_samples,w,fields,couplings,cond_likelihood)
    ereg = regularization_strength * (sum(fields**2) + 0.5_kflt * sum(couplings**2))

    ! update gradient 
    grd1 = model_f1 - data_f1 + 2.0_kflt * regularization_strength * fields
    grd2 = model_f2 - data_f2 + 2.0_kflt * 0.5_kflt * regularization_strength * couplings
    
  end subroutine update_gradient

end module model
  

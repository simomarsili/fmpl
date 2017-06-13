! Copyright (C) 2015-2017, Simone Marsili
! All rights reserved.
! License: BSD 3 clause
! License: BSD 3 clause

module scrs
  use kinds
  implicit none
  private

  public :: compute_scores, go_scores
  public :: print_scores

  real(kflt), allocatable, save :: scores(:,:)

contains

  pure integer function kmap(v2,a1,a2) result(k)
    ! returns the index of element (b,v,a,u) in prm
    use data, only: nv,ns
    implicit none
    integer, intent(in) :: v2,a1,a2
    k = ns + (a1 - 1) * ns*nv + (v2 - 1) * ns + a2
  end function kmap

  real(kflt) function frobenius(a)
    real(kflt), intent(in) :: a(:,:)

    frobenius = sqrt(sum(a**2))

  end function frobenius

  subroutine apc_correction(a)
    real(kflt), intent(inout) :: a(:,:)
    real(kflt), allocatable :: sums(:)
    real(kflt) :: totsum
    integer :: err
    integer :: i,j
    integer :: n

    n = size(a,1)
    allocate(sums(n),stat=err)

    do i = 1,n
       sums(i) = sum(a(:,i))
    end do
    totsum = sum(sums)

    do j = 1,n
       do i = j,n
          a(i,j) = a(i,j) - sums(i)*sums(j)/totsum
          a(j,i) = a(i,j)
       end do
    end do

    deallocate(sums)

  end subroutine apc_correction

  subroutine symmetrize_prm(nv,ns,prm)
    ! compute scores (no temporary array)
    integer, intent(in) :: nv,ns
    real(kflt), intent(inout) :: prm(:,:)
    integer :: iv,jv,is,js,k,k1,k2
    integer :: err
    real(kflt) :: prmm

    do iv = 1,nv-1
       k = ns
       do is = 1,ns
          do jv = iv+1,nv
             do js = 1,ns
                k = k + 1
                k1 = kmap(jv,is,js)
                k2 = kmap(iv,js,is)
                prm(k1,iv) = prm(k1,iv) + prm(k2,jv)
                prm(k2,jv) = prm(k1,iv)
             end do
          end do
       end do
    end do
    prm = 0.5_kflt * prm

  end subroutine symmetrize_prm

  subroutine compute_scores(nv,ns,prm,skip_gaps,symmetrize)
    ! compute scores (no temporary array)
    integer, intent(in) :: nv,ns
    real(kflt), intent(inout) :: prm(:,:)
    logical, intent(in) :: skip_gaps,symmetrize
    integer :: iv,jv,is,js,k
    integer :: err
    real(kflt) :: prmm

    ! at the very end of the run
    allocate(scores(nv,nv),stat=err)
    scores = 0.0_kflt

    if (symmetrize) call symmetrize_prm(nv,ns,prm)

    do iv = 1,nv
       k = ns
       do is = 1,ns
          do jv = 1,nv
             do js = 1,ns
                k = k + 1
                if (skip_gaps .and. (is == 1 .or. js == 1)) cycle
                scores(iv,jv) = scores(iv,jv) + prm(k,iv)**2
             end do
          end do
       end do
    end do
    scores = sqrt(scores)

    call apc_correction(scores)

  end subroutine compute_scores

  subroutine print_scores(mat,uscrs,print_indices)
    real(kflt), intent(in) :: mat(:,:)
    integer, intent(in) :: uscrs
    logical, optional, intent(in) :: print_indices
    logical :: inds = .true.
    integer :: i,j,n1,n2

    n1 = size(mat,1)
    n2 = size(mat,2)
    if (present(print_indices)) then
       if (.not.print_indices) then
          inds = .false.
          write(uscrs,'(i5,1x,i5)') n1,n2
       end if
    end if
    
    do i = 1,n1-1
       do j = i+1,n2
          if (inds) then
             write(uscrs,'(i5,1x,i5,1x,f8.4)') i, j, mat(i,j)
          else
             write(uscrs,'(f8.4)') mat(i,j)
          end if
       end do
    end do

  end subroutine print_scores

  subroutine go_scores(nv,ns,couplings,scores)
    implicit none
    integer(kint), intent(in) :: nv,ns
    real(kflt), intent(in) :: couplings(ns,ns,nv,nv)
    real(kflt), intent(out) :: scores(nv,nv)
    integer :: iv,jv,err
    real(kflt) :: sij
    real(kflt), allocatable :: sums(:)
    real(kflt) :: totsum
    
    scores = 0.0_kflt
    do iv = 1,nv-1
       do jv = iv+1,nv
          scores(iv,jv) = sum((couplings(:,:,iv,jv) + transpose(couplings(:,:,jv,iv)))**2)
          scores(jv,iv) = scores(iv,jv)
       end do
    end do
    scores = 0.5_kflt * sqrt(scores)
    
    allocate(sums(nv),stat=err)
    do iv = 1,nv
       sums(iv) = sum(scores(:,iv))
    end do
    totsum = sum(sums)
    
    do jv = 1,nv
       do iv = jv,nv
          scores(iv,jv) = scores(iv,jv) - sums(iv)*sums(jv)/totsum
          scores(jv,iv) = scores(iv,jv)
       end do
    end do
    
    deallocate(sums)

!    call print_scores(uscrs)
!    do iv = 1,nv-1
!       do jv = iv+1,nv
!          sij = scores(iv,jv)
!          write(uscrs,'(i5,1x,i5,1x,f8.5)') iv,jv,sij
!       end do
!    end do
    
  end subroutine go_scores

end module scrs

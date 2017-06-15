! Copyright (C) 2015-2017, Simone Marsili
! All rights reserved.
! License: BSD 3 clause

module scrs
  use kinds
  implicit none
  private
  public :: compute_scores
  public :: print_mat

contains

  pure real(kflt) function frobenius(a)
    real(kflt), intent(in) :: a(:,:)
    
    frobenius = sqrt(sum(a**2))
  end function frobenius

  pure subroutine apc_correction(a)
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

  subroutine print_mat(mat,uscrs,frmt,nerr)
    ! https://people.sc.fsu.edu/~jburkardt/data/mm/mm.html
    real(kflt), intent(in) :: mat(:,:)
    integer, intent(in) :: uscrs
    character(1), intent(in) :: frmt
    integer, intent(out) :: nerr
    logical :: inds = .true.
    integer :: i,j,n1,n2

    nerr = 0
    ! check format
    if (.not. any(["r","c","a"] == frmt)) then
       write(0,*) "error ! unknown format for scores matrix file"
       nerr = nerr + 1
       return
    end if
    n1 = size(mat,1)
    n2 = size(mat,2)
    ! write header line
    select case(frmt)
    case("r")
       do j = 1,n1
          do i = j,n2
             write(uscrs,'(i5,1x,i5,1x,f8.4)') i, j, mat(i,j)
          end do
       end do
    case("a")
       write(uscrs,'(a)') "%%MatrixMarket matrix array real symmetric"
       write(uscrs,'(a)') "%Lower triangular, column-major order"
       write(uscrs,'(i5,1x,i5)') n1,n2
       do j = 1,n1
          do i = j,n2
             write(uscrs,'(f8.4)') mat(i,j)
          end do
       end do
    case("c")
       write(uscrs,'(a)') "%%MatrixMarket matrix coordinate real symmetric"
       write(uscrs,'(a)') "%Lower triangular, column-major order"
       write(uscrs,'(i5,1x,i5,1x,i5)') n1,n2,n1*n2
       do j = 1,n1
          do i = j,n2
             write(uscrs,'(i5,1x,i5,1x,f8.4)') i, j, mat(i,j)
          end do
       end do
    end select

  end subroutine print_mat

  subroutine compute_scores(nv,ns,couplings,ignore_pivot,scores)
    implicit none
    integer(kint), intent(in) :: nv,ns
    real(kflt), intent(in) :: couplings(ns,ns,nv,nv)
    integer, intent(in) :: ignore_pivot
    real(kflt), intent(out) :: scores(nv,nv)
    integer :: iv,jv,err,is,js
    real(kflt) :: sums(nv),rsum(ns),csum(ns),matij(ns,ns),matji(ns,ns)
    real(kflt) :: totsum
    integer :: inds(ns-1)
    
    scores = 0.0_kflt
    do iv = 1,nv-1
       do jv = iv+1,nv
          if (ignore_pivot >= 1 .and. ignore_pivot <= ns) then
             is = 0
             do js = 1,ns
                if (js /= ignore_pivot) then
                   is = is + 1
                   inds(is) = js
                end if
             end do
             matij = couplings(inds,inds,iv,jv)
             matji = couplings(inds,inds,jv,iv)
          else
             matij = couplings(:,:,iv,jv)
             matji = couplings(:,:,jv,iv)
          end if
          call double_center(matij)
          call double_center(matji)
          scores(iv,jv) = sum((matij + transpose(matji))**2)
          scores(jv,iv) = scores(iv,jv)
       end do
    end do
    scores = 0.5_kflt * sqrt(scores)
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
  end subroutine compute_scores

  pure subroutine double_center(a)
    real(kflt), intent(inout) :: a(:,:)
    integer :: i,j,n,m
    real(kflt) :: totsum
    real(kflt), allocatable :: rsum(:),csum(:)
    integer :: err
    n = size(a,dim=1)
    m = size(a,dim=2)
    allocate(rsum(n),csum(m),stat=err)
    totsum = sum(a) / real(n*m)
    csum = sum(a,dim=2) / real(n)
    rsum = sum(a,dim=1) / real(m)
    do j = 1,m
       do i = 1,n
          a(i,j) =  a(i,j) - csum(i) - rsum(j) + totsum
       end do
    end do
    deallocate(rsum,csum)
  end subroutine double_center

end module scrs

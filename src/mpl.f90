! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program mpl
  use kinds
  use constants,     only: long_string
  use command_line,  only: read_args
  use data,          only: initialize_data,read_data,data_reweight
  use model,         only: initialize_model, model_set_myv,fix_gauge
  use scrs,          only: print_mat, compute_scores
  use dvmlm_wrapper, only: dvmlm_minimize
  implicit none
  ! input variables
  character(long_string) :: data_file, scores_file
  real(kflt)             :: w_id, lambda
  integer                :: ignore_pivot, accuracy
  character(1)           :: scores_format
  ! data dimensions
  integer :: nd ! n. of samples
  integer :: nv ! n. of variables
  integer :: ns ! n. of classes
  real(kflt) :: neff
  ! main arrays for the run
  integer,    allocatable :: data_samples(:,:)  ! nv x nd
  real(kflt), allocatable :: w(:)
  real(kflt), allocatable :: prm(:,:)           ! (ns + ns x ns x nv) x nv
  real(kflt), allocatable :: grd(:)             ! (ns + ns x ns x nv) x nv
  real(kflt), allocatable :: fields(:,:)        ! ns x nv
  real(kflt), allocatable :: couplings(:,:,:,:) ! ns x ns x nv x nv
  real(kflt), allocatable :: scores(:,:)        ! nv x nv
  !
  integer                :: udata,uscrs
  integer                :: err,iv
  integer                :: niter,neval
  real(kflt_single)      :: finish,start,start_min,tpv,elapsed_time,expected_time
  character(4) :: time_unit
  

  ! get command line 
  call read_args(data_file,w_id,lambda,ignore_pivot,accuracy,scores_format,err)
  if(err /= 0) then 
     write(0,100) 
     stop
  end if
  ! print an header for the run
  write(0,101) trim(data_file), lambda, accuracy
  if (w_id > 0.0_kflt) write(0,102) w_id

  ! open data file
  open(newunit=udata,file=data_file,status='old',action='read',iostat=err)
  if(err /= 0) then 
     write(0,'("error ! cannot access ",a,": file not found")') trim(data_file)
     stop
  end if
  
  ! open scores file
  scores_file = trim(data_file)//'.scores'
  open(newunit=uscrs,file=scores_file,status='replace',action='write',iostat=err)

  call initialize_data(udata,nd,nv,neff)
  allocate(data_samples(nv,nd),stat=err)
  allocate(w(nd),stat=err)

  write(0,'(/a)') "Reading data.."
  call read_data(udata,w_id,ns,neff,data_samples)
  write(0,103) nd, nv, ns
  flush(0)
  if(w_id > 1.E-10_kflt) then
     write(0,'(a)') "Computing weights.."
     call data_reweight(data_samples,w_id,neff,w)
     write(0,'(a,f8.1)') "Neff: ",neff
  else
     neff = nd
     w = 1.0_kflt / neff
  end if

  ! allocate parameters and gradient
  allocate(prm(ns + ns*ns*nv,nv),stat=err)
  allocate(grd(ns + ns*ns*nv),stat=err)
  allocate(scores(nv,nv),stat=err)
  call initialize_model(nv,ns,lambda)

  write(0,'(/,a)') 'Running..'
  call cpu_time(start_min)
  tpv = 0.0_kflt
  ! loop over features
  do iv = 1,nv
     call model_set_myv(nd,nv,iv,data_samples,w,prm(:,iv),grd,err)
     niter = 0
     call cpu_time(start)
     call dvmlm_minimize(nv,ns,nd,size(prm(:,iv)),data_samples,w,prm(:,iv),grd,accuracy,niter,neval)
     call cpu_time(finish)
     elapsed_time = finish - start_min
     tpv = elapsed_time / real(iv)
     expected_time = tpv * (nv - iv)
     time_unit = "secs"
     if (elapsed_time + expected_time > 2*60.) then
        elapsed_time = elapsed_time / 60.0_kflt
        expected_time = expected_time / 60.0_kflt
        time_unit = "mins"
     end if
     if (mod(iv,int(3.0/tpv)+1)==0) &
          write(0,'(i4,"/",i4," completed in ",f5.1,1x,a,"; ",f5.1," to end")')&
          iv, nv, elapsed_time, time_unit, expected_time
     if (iv == nv) &
          write(0,'("*** ",i4,"/",i4," completed in ",f9.3,1x,a," ***")')&
          iv, nv, elapsed_time, "secs"
     !write(0,'(a,i5,a,2i5,a,f8.3,a)') ' variable ', iv, &
     !'  converged (niter,neval) ', niter, neval, ' in ', finish-start, ' secs'
     !call fix_gauge(nv,ns,prm(:ns,iv),prm(ns+1:,iv))
  end do
  flush(0)
  
  ! reorder prm array into fields and couplings 
  allocate(fields(ns,nv),couplings(ns,ns,nv,nv),stat=err)
  call reshape_prm(nv,ns,prm,fields,couplings)
  deallocate(prm,grd)

  ! compute scores and print
  call compute_scores(nv,ns,couplings,ignore_pivot,scores)
  call print_mat(scores,uscrs,scores_format,err)
  if (err /= 0) stop

101 format(&
         '# mpl                                                              '/&
         '#                                                                  '/&
         '# input data file: ',a,'                                           '/&
         '# regularization: ',f6.4,'                                         '/&
         '# accuracy level: ',i1,'                                           ')

102 format(&
         '# reweight samples; % id threshold: ',f5.1,'                       ')
 
103 format(&
         'Sample size         : ',i6,'                                       '/&
         'Dimensionality      : ',i6,'                                       '/&
         'Classes per variable: ',i6,'                                       ')

100 format(/&
         'mpl                                                            '/&
         '                                                               '/&
         'Usage:                                                         '/&
         '    mpl [options] -i <file>                                    '/&
         '                                                               '/&
         'Options:                                                       '/&
         '-h, --help                                                     '/&
         '    Display this help and exit.                                '/&
         '                                                               '/&
         '-l, --lambda <regularization_parameter> float, optional        '/&
         '    Controls L_2 regularization strength.                      '/&
         '    [default: 0.01]                                            '/&
         '                                                               '/&
         '-w, --reweight <percentage_identity> float, optional           '/&
         '    Reweight data using a percentage identity threshold.       '/&
         '    [default: no reweight.]                                    '/&
         '                                                               '/&
         '-a, --accuracy <accuracy_level> integer, optional              '/&
         '    Larger values correspond to stricter convergence criteria  '/&
         '    and longer covergence time. Possible values are {0, 1, 2}. '/&
         '    [default: 1.]                                              '/&
         '                                                               '/&
         '--ignore_pivot <pivot_state> integer, optional                 '/&
         '    Ignore the contribution of <pivot_state> to final scores.  '/&
         '    [default: include all states.]                             '/&
         '                                                               '/&
         '--scores_format <scores_matrix_format> - string, optional      '/&
         '    Possible values are {"rcv", "array", "coordinate"}.        '/&
         '    "rcv": (row_index, column_index, value) format             '/&
         '    "coordinate", "array": Matrix Market (MM) formats. See:    '/&
         '    https://people.sc.fsu.edu/~jburkardt/data/mm/mm.html       '/&
         '    [default: "rcv"].                                          '/)
  
  
contains

  subroutine reshape_prm(nv,ns,prm,fields,couplings)
    ! reorder prm array into fields and couplings 
    implicit none
    integer, intent(in) :: nv,ns
    real(kflt), intent(in) :: prm(:,:)
    real(kflt), intent(out) :: fields(ns,nv)
    real(kflt), intent(out) :: couplings(ns,ns,nv,nv)
    integer :: err
    integer :: iv,jv,is,js,k1,k2,k
    
    k = 0
    do iv = 1,nv
       fields(:,iv) = prm(:ns,iv)
       k = ns
       do jv = 1,nv
          do js = 1,ns
             do is = 1,ns
                k = k + 1
                couplings(js,is,jv,iv) = prm(k,iv)
             end do
          end do
       end do
    end do
  end subroutine reshape_prm

end program mpl

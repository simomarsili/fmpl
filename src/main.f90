! Copyright (C) 2015-2017, Simone Marsili
! All rights reserved.
! License: BSD 3 clause

program fmpl
  use kinds
  use constants,     only: long_string
  use command_line,  only: read_args
  use data,          only: initialize_data,read_data,data_reweight
  use model,         only: initialize_model, model_reset,fix_gauge
  use scrs,          only: print_mat, compute_scores
  use optimize,      only: fit
  implicit none
  ! input variables
  character(long_string) :: data_file, prm_file, scores_file
  real(kflt)             :: w_id, lambda
  integer                :: ignore_pivot, accuracy
  character(1)           :: scores_format
  logical                :: dump_prm
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
  integer                :: udata,uscrs,uprm
  integer                :: err,iv
  real(kflt_single)      :: finish,start,start_min,tpv,elapsed_time,expected_time
  character(4) :: time_unit
  logical :: exist
  character(long_string) :: string
  character(long_string) :: minimizer='dvmlm'

  ! get command line
  call read_args(data_file,prm_file,w_id,lambda,ignore_pivot,accuracy,scores_format,dump_prm,err)
  if(err /= 0) then
     write(0,100)
     stop
  end if
  
  if (len_trim(prm_file) > 0) then
     ! no fitting; read prm from prm file
     ! open unit
     open(newunit=uprm,file=prm_file,status='old',action='read',&
          access='stream',form='unformatted',iostat=err)
     if (err /= 0) then
        write(0,'("error ! cannot access ",a,": file not found")') trim(prm_file)
        stop
     end if

     ! read prm from binary file
     read(uprm) nv
     read(uprm) ns
     allocate(prm(ns + ns*ns*nv,nv),stat=err)
     allocate(scores(nv,nv),stat=err)
     read(uprm) prm
     close(uprm)
     
     ! print an header for the run
     write(0,104) trim(prm_file), nv, ns
  else
     ! print an header for the run
     write(0,101) trim(data_file), lambda, accuracy
     if (w_id > 0.0_kflt) write(0,102) w_id

     ! open data file
     open(newunit=udata,file=data_file,status='old',action='read',iostat=err)
     if(err /= 0) then
        write(0,'("error ! cannot access ",a,": file not found")') trim(data_file)
        stop
     end if

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
        call model_reset(nd,nv,iv,data_samples,w,prm(:,iv),err)
        call cpu_time(start)
        call fit(nv,ns,nd,data_samples,w,prm(:,iv),grd,accuracy,minimizer)
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
        !call fix_gauge(nv,ns,prm(:ns,iv),prm(ns+1:,iv))
     end do
     flush(0)
  end if

  ! dump prm file
  if (dump_prm) then
     open(newunit=uprm,file=trim(data_file)//'.prm',status='unknown',action='write',&
          access='stream',form='unformatted',iostat=err)
     write(uprm) nv
     write(uprm) ns
     write(uprm) prm
     close(uprm)
  end if

  ! reorder prm array into fields and couplings
  allocate(fields(ns,nv),couplings(ns,ns,nv,nv),stat=err)
  call reshape_prm(nv,ns,prm,fields,couplings)
  deallocate(prm)
  if (allocated(grd)) deallocate(grd)

  ! compute scores and print
  call compute_scores(nv,ns,couplings,ignore_pivot,scores)

  ! open scores file
  if (len_trim(prm_file) > 0) then
     scores_file = trim(prm_file)//'.scores'
  else
     scores_file = trim(data_file)//'.scores'
  end if
  open(newunit=uscrs,file=scores_file,status='replace',action='write',iostat=err)

  call print_mat(scores,uscrs,scores_format,err)
  if (err /= 0) stop

101 format(&
       '# fmpl                                                             '/&
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

104 format(&
       '# fmpl                                                             '/&
       '#                                                                  '/&
       '# prm file: ',a,'                                                  '/&
       '                                                                   '/&
       'Dimensionality      : ',i6,'                                       '/&
       'Classes per variable: ',i6,'                                       ')

100 format(&
       'fmpl                                                           '/&
       '                                                               '/&
       'Usage:                                                         '/&
       '    fmpl [options] -i <data_file>                              '/&
       '    fmpl [options] -p <prm_file>                               '/&
       '                                                               '/&
       'Description:                                                   '/&
       '    Either read a matrix of data (-i <file>) and fit a set of  '/&
       '    parameters, or read them directly (-p <file>).             '/&
       '    Output a file of scores (<file>.scores).                   '/&
       '                                                               '/&
       'Options:                                                       '/&
       '-h, --help                                                     '/&
       '    Display this help and exit.                                '/&
       '                                                               '/&
       '-l, --lambda <regularization_parameter>, float                 '/&
       '    Controls L_2 regularization strength.                      '/&
       '    [default: 0.01]                                            '/&
       '                                                               '/&
       '-w, --reweight <percentage_identity>, float                    '/&
       '    Reweight data using a percentage identity threshold.       '/&
       '    [default: no reweight.]                                    '/&
       '                                                               '/&
       '-a, --accuracy <accuracy_level>, integer                       '/&
       '    Larger values correspond to stricter convergence criteria  '/&
       '    and longer covergence times. Possible values are {0, 1, 2}.'/&
       '    [default: 1.]                                              '/&
       '                                                               '/&
       '--dump_prm, bool                                               '/&
       '    Dump parameters to binary file <data_file>.                '/&
       '    [default: false]                                           '/&
       '                                                               '/&
       '--ignore_pivot <pivot_state>, integer                          '/&
       '    Ignore the contribution of <pivot_state> to final scores.  '/&
       '    [default: include all states]                              '/&
       '                                                               '/&
       '--scores_format <scores_matrix_format>, string                 '/&
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

end program fmpl

! Femsim
!
! Simulates fluctuation electron microscopy V(k)
! from an input model.
!
! Voyles research group: Jinwoo Hwang, Feng Yi, and
! Paul Voyles, begun 12/16/08

program femsim

    use omp_lib
    use rmc_global
    use readinputs
    use model_mod
    use fem_mod
    implicit none
    include 'mpif.h'
    type(model) :: m
    character (len=256) :: model_filename
    character (len=256) :: param_filename  
    character (len=512) :: comment
    real :: Q, res
    real, pointer, dimension(:) :: vk, vk_exp, k, vk_exp_err, v_background
    real, pointer, dimension(:,:) :: scatfact_e
    integer :: i
    integer :: nk
    integer :: ntheta, nphi, npsi
    integer :: istat, status2
    integer :: iseed2
    integer :: ipvd
    doubleprecision :: t0, t1 !timers

    call mpi_init_thread(MPI_THREAD_MULTIPLE, ipvd, mpierr) !http://www.open-mpi.org/doc/v1.5/man3/MPI_Init_thread.3.php
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call mpi_comm_size(mpi_comm_world, numprocs, mpierr)
    
    t0 = mpi_wtime()

    call mpi_comm_size(mpi_comm_world, numprocs, mpierr)

    !model_filename = 'model.xyz'
    param_filename = '/export/home/group/femsim/param_file.in'

    call get_command_argument(1, model_filename)
    if(len_trim(model_filename) == 0) then
        write(*,*) "You must input a model file as the first parameter"
        stop
    endif

    ! Read input model
    call read_model(model_filename, comment, m, istat)
    call check_model(m, istat)
    call recenter_model(0.0, 0.0, 0.0, m)

    ! Read input parameters
    ! Many of these parameters you can ignore for femsim.
    call read_inputs(param_filename, vk_exp, k, vk_exp_err, v_background, ntheta, nphi, npsi, Q, status2)

    res = 0.61/Q
    nk = size(k)

    iseed2 = 104756

    call fem_initialize(m, res, k, nk, ntheta, nphi, npsi, scatfact_e, istat)
    allocate(vk(size(vk_exp)))
    write(*,*) "Calculating V(k)... will be output to vk.txt"
    call fem(m, res, k, vk, v_background, scatfact_e, mpi_comm_world, istat)

    t1 = mpi_wtime()

    if(myid.eq.0)then
        ! Write initial vk 
        open(unit=52,file="vk.txt",form='formatted',status='unknown')
            do i=1, nk
                write(52,*)k(i),vk(i)
            enddo
        close(52)
    endif

    call mpi_finalize(mpierr)

end program femsim


! This module contains functions for femsim.

module fem_mod
    use  model_mod
    use  RMC_Global
    use  scattering_factors 

    implicit none
    private
    public :: fem_initialize, fem
    integer, save :: nk, nrot  ! number of k points, pixels, and rotations
    type pix_array
        real, dimension(:,:), pointer :: pix ! npix x 2 list of pixel positions
        integer :: npix, npix_1D ! number of pixels and number of pixels in 1 dimension
        real :: phys_diam
        real :: dr ! Distance between pixels. Note, therefore, that there is half this distance between the pixels and the world edge. This is NOT the distance between the pixel centers. This is the distance between the edges of two different pixels. dr + phys_diam is the distance between the pixel centers!
    end type pix_array
    real, save, dimension(:,:), pointer :: rot ! nrot x 3 list of (phi, psi, theta) rotation angles
    real, save, dimension(:,:,:), pointer :: int_i, int_sq  ! nk x npix x nrot.  int_sq == int_i**2
    real, save, dimension(:), pointer :: int_sum, int_sq_sum  ! nk long sums of int and int_sq arrays for calculating V(k)
    real, save, allocatable, dimension(:) :: j0, A1                                               
    type(model), save, dimension(:), pointer :: mrot  ! array of rotated models
    type(pix_array), save :: pa

contains

    subroutine fem_initialize(m, res, k, nki, ntheta, nphi, npsi, scatfact_e, istat)
        type(model), intent(in) :: m 
        real, intent(in) :: res
        real, dimension(:), intent(in) :: k 
        integer, intent(in) :: nki, ntheta, nphi, npsi 
        real, dimension(:,:), pointer :: scatfact_e
        integer, intent(out) :: istat
        !real :: dr ! Distance between pixels
        real r_max, const1, const2, const3
        integer bin_max
        integer i, j 
        integer const4 
        double precision b_x, b_j0 , b_j1 

        r_max = 2*res     !assuming resolution=radius

        bin_max = int(r_max/fem_bin_width)+1

        const1 = twopi*(0.61/res)/fem_bin_width  !(0.61/res = Q) 
        const2 = 1/fem_bin_width
        const3 = (const1/(0.61/res))/const2
        const4 = int(bin_max*const3*CEILING(k(SIZE(k))))+1

        allocate(j0(0:const4),a1(0:const4), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Failed to allocate memory for Bessel and Airy functions.'
            return
        endif

        !calculating bessel function
        j0=0.0
        do i=1, const4
            b_x = i*fem_bin_width
            call bessel_func(b_x, b_j0,b_j1)
            j0(i) = sngl(b_j0)
            a1(i) = 2*sngl(b_j1)/b_x
        enddo
        a1(0)=1.0
        j0(0)=1.0

        nk = nki

        call init_rot(ntheta, nphi, npsi, nrot, istat)
        call init_pix(m, res, istat, .true.)

        allocate(int_i(nk, pa%npix, nrot), int_sq(nk, pa%npix, nrot), int_sum(nk), int_sq_sum(nk), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Cannot allocate memory in fem_initialize.'
            return
        endif

        if( mod(m%lx,pa%phys_diam) >= 0.001 ) then
            ! The following warning is for reverse monte carlo - not femsim.
            !write(*,*) "WARNING! Your world size should be an integer multiple of the resolution. Pixel diameter = ", pa%phys_diam, ". World size = ", m%lx
        endif

        call read_f_e
        allocate(scatfact_e(m%nelements,nk), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Allocation of electron scattering factors table failed.'
            return
        endif

        do j=1,m%nelements
            do i=1, nk
                scatfact_e(j,i)=f_e(m%atom_type(j),k(i))
            enddo
        enddo

    end subroutine fem_initialize


    subroutine init_rot(ntheta, nphi, npsi, num_rot, istat)
    ! Calculates the rotation angles and initializes them into the global
    ! rotation array rot. The rot_temp variable is probably unnecessary.
        integer, intent(in) :: ntheta, nphi, npsi
        integer, intent(out) :: istat
        integer :: i,j, k, jj
        real,dimension(3) :: step_size
        integer :: ntheta_w(nphi*npsi)
        integer, intent(out) :: num_rot
        real, dimension(:,:), allocatable :: rot_temp
        real :: psi_temp
        integer :: pp

        allocate(rot_temp(ntheta*nphi*npsi,3), stat=istat)
        if (istat /= 0) then
           write (*,*) 'Cannot allocate temporary rotations array.'
           return
        endif

        !phi runs from 0 to 2 PI
        !psi runs from 0 to 2 PI
        !theta runs from 0 to PI   !not sure any more after weighting by psi angle - JWH 09/03/09
        !step_size(1) for phi step    
        !step_size(2) for psi step
        !step_size(3) for theta step
        step_size(1) = TWOPI / nphi
        step_size(2) = TWOPI / npsi
        step_size(3) = PI / ntheta  !not used any more after weighting by psi angle - JWH 09/03/09

        jj = 1
        do i=1, nphi
            do j=1, npsi/2
                psi_temp = (j-1)*step_size(2)
                ntheta_w(j) = int(sin(psi_temp)*ntheta)
                if(ntheta_w(j).ge.0)then
                    if(ntheta_w(j).gt.0)then
                        pp = 2*(ntheta_w(j)-1)
                    endif
                    if(ntheta_w(j).eq.0)then
                        pp = 1
                    endif
                    do k=1, pp
                        if(k*(pi/(ntheta_w(j)-1)).lt.pi)then
                            rot_temp(jj,1) = (i-1)*step_size(1)
                            rot_temp(jj,2) = (j-1)*step_size(2)
                            rot_temp(jj,3) = k*(pi/(ntheta_w(j)-1))
                            jj = jj + 1
                        endif
                    enddo
                endif
            enddo
        enddo

        num_rot = jj - 1

        allocate(rot(num_rot, 3), stat=istat)
        if (istat /= 0) then
           write (*,*) 'Cannot allocate rotations array.'
           return
        endif

        do i=1, num_rot
            rot(i,1) = rot_temp(i,1)
            rot(i,2) = rot_temp(i,2)
            rot(i,3) = rot_temp(i,3)
        enddo

        deallocate(rot_temp)
    end subroutine init_rot


    subroutine init_pix(m, res, istat, square_pixel)
        type(model), intent(in) :: m
        real, intent(in) :: res ! Pixel width.
        integer, intent(out) :: istat
        logical, optional, intent(in) :: square_pixel
        integer :: i, j, k
        logical :: pixel_square

        if(present(square_pixel)) then
            pixel_square = square_pixel
        else
            pixel_square = .FALSE.
        endif

        if(pixel_square) then
            pa%phys_diam = res * sqrt(2.0)
        else
            pa%phys_diam = res
        endif
        pa%npix_1D = floor( m%lx / pa%phys_diam )
        pa%npix = pa%npix_1D**2

        pa%dr = m%lx/pa%npix_1D - pa%phys_diam

        allocate(pa%pix(pa%npix, 2), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Cannot allocate pixel position array.'
            return
        endif

        k=1
        do i=1, pa%npix_1D
            do j=1, pa%npix_1D
                pa%pix(k,1) = -m%lx/2.0 + (pa%phys_diam+pa%dr)/2.0 + (pa%phys_diam+pa%dr)*(i-1)
                pa%pix(k,2) = -m%ly/2.0 + (pa%phys_diam+pa%dr)/2.0 + (pa%phys_diam+pa%dr)*(j-1)
                k = k + 1
            enddo
        enddo

        ! Uncomment below lines if you want to see pixel information.
        !if(myid.eq.0)then
        !    write(*,*)"pixels=", pa%npix_1D, "by", pa%npix_1D
        !    write(*,*) "They are centered at:"
        !    k=1
        !    do i=1, pa%npix_1D
        !        do j=1, pa%npix_1D
        !            write(*,*)"(", pa%pix(k,1), ",", pa%pix(k,2), ")"
        !            k=k+1
        !        enddo
        !    enddo
        !    write(*,*) "with a distance between pixels of", pa%dr
        !endif
    end subroutine init_pix


    subroutine fem(m, res, k, vk, v_background, scatfact_e, comm, istat, rot_begin, rot_end)
        use mpi
        implicit none
        type(model), intent(in) :: m
        real, intent(in) :: res
        real, dimension(:), intent(in) :: k
        real, dimension(:), INTENT(OUT) :: Vk
        real, dimension(:), intent(in) :: v_background
        real, dimension(:,:), pointer :: scatfact_e
        integer, intent(out) :: istat
        integer, optional, intent(in) :: rot_begin, rot_end
        real, dimension(:), allocatable :: psum_int, psum_int_sq, sum_int, sum_int_sq  !mpi
        integer :: comm
        integer :: i, j
        integer begin_rot, end_rot

        if(present(rot_begin)) then
            begin_rot = rot_begin
        else
            begin_rot  = 1
        endif

        if(present(rot_end)) then
            end_rot = rot_end
        else
            end_rot = nrot ! This is set in fem_initialize->init_rot
        endif

        allocate (psum_int(size(k)), psum_int_sq(size(k)), sum_int(size(k)), sum_int_sq(size(k)), stat=istat)
        sum_int = 0.0
        sum_int_sq = 0.0
        psum_int = 0.0
        psum_int_sq = 0.0

        ! Initialize the rotated models
        allocate(mrot(nrot), stat=istat)
        call check_allocation(istat, 'Cannot allocate rotated model array.')

        ! Calculate all the rotated models and save them in mrot.
        do i=myid+1, nrot, numprocs
            call rotate_model(rot(i, 1), rot(i, 2), rot(i, 3), m, mrot(i), istat)
        call check_allocation(istat, 'Failed to rotate model')
        enddo

        ! Calculate intensities for every single pixel in every single model. This is very expensive.
        !write(*,*); write(*,*) "Calculating intensities over the models: nrot = ", nrot; write(*,*)
        do i=myid+1, nrot, numprocs
            do j=1, pa%npix
                !write(*,*) "Calling intensity on pixel (", pa%pix(j,1), ",",pa%pix(j,2), ") in rotated model ", i
                call intensity(mrot(i), res, pa%pix(j, 1), pa%pix(j, 2), k, int_i(1:nk, j, i), scatfact_e, istat)
                int_sq(1:nk, j, i) = int_i(1:nk, j, i)**2
                psum_int(1:nk) = psum_int(1:nk) + int_i(1:nk, j, i)
                psum_int_sq(1:nk) = psum_int_sq(1:nk) + int_sq(1:nk, j, i)
            enddo
        enddo

        call mpi_reduce (psum_int, sum_int, size(k), mpi_real, mpi_sum, 0, comm, mpierr)
        call mpi_reduce (psum_int_sq, sum_int_sq, size(k), mpi_real, mpi_sum, 0, comm, mpierr)

        if(myid.eq.0)then
            do i=1, nk
                Vk(i) = (sum_int_sq(i)/(pa%npix*nrot))/((sum_int(i)/(pa%npix*nrot))**2)-1.0
                Vk(i) = Vk(i) - v_background(i)  ! background subtraction   052210 JWH
            end do
        endif

        deallocate(psum_int, psum_int_sq, sum_int, sum_int_sq)
    end subroutine fem


    subroutine intensity(m_int, res, px, py, k, int_i, scatfact_e, istat)
    ! Calculates int_i.
        use  omp_lib
        type(model), intent(in) :: m_int
        real, intent(in) :: res, px, py
        real, dimension(nk), intent(in) :: k
        real, dimension(nk), intent(out) :: int_i
        real, dimension(:,:), pointer :: scatfact_e
        integer, intent(out) :: istat
        real, dimension(:,:,:), allocatable :: gr_i   ! unneeded 'save' keyword removed pmv 03/18/09  !tr re-ok -jwh
        real, dimension(:), allocatable ::x1, y1, rr_a
        real, dimension(:,:), allocatable :: sum1
        real :: x2, y2, rr, t1, t2, const1, const2, const3, pp, r_max
        integer, pointer, dimension(:) :: pix_atoms, znum_r
        integer :: i,j,ii,jj,kk
        integer :: bin_max, size_pix_atoms
        real, allocatable, dimension(:) :: rr_x, rr_y
        real :: sqrt1_2_res
        real :: k_1
        real :: timer1, timer2
        integer :: nthr, thrnum

        call cpu_time(timer1)

        sqrt1_2_res = SQRT(0.5) * res
        r_max = 2*res !small pixel inscribed in airy circle
        call hutch_list_pixel_sq(m_int, px, py, pa%phys_diam, pix_atoms, istat)
        allocate( rr_x(size(pix_atoms)),rr_y(size(pix_atoms)), stat=istat)

        size_pix_atoms = size(pix_atoms)
        bin_max = int(r_max/fem_bin_width)+1

        allocate(gr_i(m_int%nelements,m_int%nelements, 0:bin_max), stat=istat)
        allocate(x1(size_pix_atoms),y1(size_pix_atoms),rr_a(size_pix_atoms), stat=istat)
        allocate(sum1(m_int%nelements,size_pix_atoms), stat=istat)
        allocate(znum_r(size_pix_atoms), stat=istat)

        do i=1, size_pix_atoms
            znum_r(i) = m_int%znum_r%ind(pix_atoms(i))
        enddo

        gr_i = 0.0; int_i = 0.0; x1 = 0.0; y1 = 0.0; rr_a = 0.0

        x2 = 0.0; y2 = 0.0
        const1 = twopi*(0.61/res)/fem_bin_width  !(0.61/res = Q)
        const2 = 1/fem_bin_width
        const3 = TWOPI

        ! Uncomment the following lines if you want to see how many threads
        ! openmp is using.
        !!$omp parallel do
        !do i=1,1
        !    nthr = omp_get_num_threads() !omp_get_max_threads()
        !    thrnum = omp_get_thread_num()
        !    write(*,*) "We are using", nthr, " thread(s) in Intensity."
        !enddo
        !!$omp end parallel do

        ! Calculate sum1 for gr_i calculation in next loop.
        !$omp parallel do private(i, j, ii, jj, kk, rr, t1, t2, pp, r_max, x2, y2) shared(pix_atoms, A1, rr_a, const1, const2, const3, x1, y1, gr_i, int_i, znum_r, sum1, rr_x, rr_y)
        do i=1,size_pix_atoms
            x2=m_int%xx%ind(pix_atoms(i))-px
            y2=m_int%yy%ind(pix_atoms(i))-py
            x2=x2-m_int%lx*anint(x2/m_int%lx)
            y2=y2-m_int%ly*anint(y2/m_int%ly)
            rr_x(i) = ABS(x2)
            rr_y(i) = ABS(y2)
            rr_a(i)=sqrt(x2*x2 + y2*y2)
            !if((rr_x(i).le.res) .AND. (rr_y(i) .le. res))then
            if((rr_x(i) .le. sqrt1_2_res) .AND. (rr_y(i) .le.  sqrt1_2_res))then !small pixel inscribed in Airy circle
                k_1=0.82333
                x1(i)=x2
                y1(i)=y2
                j=int(const1*rr_a(i))
                sum1(znum_r(i),i)=A1(j)
            endif
        enddo
        !$omp end parallel do

        ! Calculate gr_i for int_i in next loop.
        !$omp parallel do private(i, j, ii, jj, kk, rr, t1, t2, pp, r_max, x2, y2) shared(pix_atoms, A1, rr_a, const1, const2, const3, x1, y1, gr_i, int_i, znum_r, sum1, rr_x, rr_y)
        do i=1,size_pix_atoms
            if((rr_x(i).le.sqrt1_2_res) .and. (rr_y(i) .le.  sqrt1_2_res))then
                do j=i,size_pix_atoms
                    if((rr_x(j).le.sqrt1_2_res) .and. (rr_y(j) .le. sqrt1_2_res))then
                        x2=x1(i)-x1(j)
                        y2=y1(i)-y1(j)
                        rr=sqrt(x2*x2 + y2*y2)
                        kk=int(const2*rr)
                        if(i == j)then
                            t1=sum1(znum_r(i),i)
                            gr_i(znum_r(i),znum_r(j),kk)=gr_i(znum_r(i),znum_r(j),kk)+t1*t1
                        else
                            t1=sum1(znum_r(i),i)
                            t2=sum1(znum_r(j),j)
                            gr_i(znum_r(i),znum_r(j),kk)=gr_i(znum_r(i),znum_r(j),kk)+2.0*t1*t2 !changed by FY on 05/04/2009
                        endif
                    endif
                enddo
            endif
        enddo
        !$omp end parallel do

        !$omp parallel do private(i, j, ii, jj, kk, rr, t1, t2, pp, r_max, x2, y2) shared(pix_atoms, A1, rr_a, const1, const2, const3, x1, y1, gr_i, int_i, znum_r, sum1, rr_x, rr_y, k)
        do i=1,nk
            do j=0,bin_max
                do ii=1,m_int%nelements
                    do jj=1,m_int%nelements
                        pp=const3*j*k(i)
                        int_i(i)=int_i(i)+scatfact_e(ii,i)*scatfact_e(jj,i)*J0(INT(pp))*gr_i(ii,jj,j)
                    enddo
                enddo
            end do
        end do
        !$omp end parallel do

        if(allocated(gr_i)) deallocate(gr_i)
        if(allocated(x1)) deallocate(x1,y1, rr_a, znum_r)
        if(size(pix_atoms) .gt. 0) deallocate(pix_atoms)
        if(allocated(sum1)) deallocate(sum1)
        if(allocated(rr_x)) deallocate(rr_x, rr_y)

        call cpu_time(timer2)
        !write ( *, * ) 'Elapsed CPU time in Intensity call = ', timer2 - timer1
    end subroutine intensity


    subroutine bessel_func(x,bj0,bj1)
        IMPLICIT none 
        doubleprecision A,B,A1,B1,BJ0,BJ1,BY0,BY1,DY0,DY1,X,X2,RP2,DJ0,DJ1,R,DABS
        doubleprecision EC,CS0,W0,R0,CS1,W1,R1,T1,T2,P0,P1,Q0,Q1,CU,DCOS,DSIN
        !integer i,j,k,l,m,n,k0
        integer k,k0 
        DIMENSION A(12),B(12),A1(12),B1(12)

        RP2=0.63661977236758D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN 
            BJ0=1.0D0
            BJ1=0.0D0
            DJ0=0.0D0
            DJ1=0.5D0
            BY0=-1.0D+300
            BY1=-1.0D+300
            DY0=1.0D+300
            DY1=1.0D+300
            RETURN
        ENDIF
        IF (X.LE.12.0D0) THEN 
            BJ0=1.0D0
            R=1.0D0
            DO 5 K=1,30
                R=-0.25D0*R*X2/(K*K)
                BJ0=BJ0+R
                IF (DABS(R).LT.DABS(BJ0)*1.0D-15) GO TO 10
5           CONTINUE
10          BJ1=1.0D0
            R=1.0D0
            DO 15 K=1,30
                R=-0.25D0*R*X2/(K*(K+1.0D0))
                BJ1=BJ1+R
                IF (DABS(R).LT.DABS(BJ1)*1.0D-15) GO TO 20
15          CONTINUE
20          BJ1=0.5D0*X*BJ1
            EC=DLOG(X/2.0D0)+0.5772156649015329D0
            CS0=0.0D0
            W0=0.0D0
            R0=1.0D0
            DO 25 K=1,30
                W0=W0+1.0D0/K
                R0=-0.25D0*R0/(K*K)*X2
                R=R0*W0
                CS0=CS0+R
                IF (DABS(R).LT.DABS(CS0)*1.0D-15) GO TO 30
25          CONTINUE
30          BY0=RP2*(EC*BJ0-CS0)
            CS1=1.0D0
            W1=0.0D0
            R1=1.0D0
            DO 35 K=1,30
                W1=W1+1.0D0/K
                R1=-0.25D0*R1/(K*(K+1))*X2
                R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
                CS1=CS1+R
                IF (DABS(R).LT.DABS(CS1)*1.0D-15) GO TO 40
35          CONTINUE
40          BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
        ELSE

            DATA A/-.7031250000000000D-01,.1121520996093750D+00, &
            -.5725014209747314D+00,.6074042001273483D+01, &
            -.1100171402692467D+03,.3038090510922384D+04, &
            -.1188384262567832D+06,.6252951493434797D+07, &
            -.4259392165047669D+09,.3646840080706556D+11, &
            -.3833534661393944D+13,.4854014686852901D+15/
            DATA B/ .7324218750000000D-01,-.2271080017089844D+00, &
            .1727727502584457D+01,-.2438052969955606D+02, &
            .5513358961220206D+03,-.1825775547429318D+05, &
            .8328593040162893D+06,-.5006958953198893D+08, &
            .3836255180230433D+10,-.3649010818849833D+12, &
            .4218971570284096D+14,-.5827244631566907D+16/
            DATA A1/.1171875000000000D+00,-.1441955566406250D+00, &
            .6765925884246826D+00,-.6883914268109947D+01, &
            .1215978918765359D+03,-.3302272294480852D+04, &
            .1276412726461746D+06,-.6656367718817688D+07, &
            .4502786003050393D+09,-.3833857520742790D+11, &
            .4011838599133198D+13,-.5060568503314727D+15/
            DATA B1/-.1025390625000000D+00,.2775764465332031D+00, &
            -.1993531733751297D+01,.2724882731126854D+02, &
            -.6038440767050702D+03,.1971837591223663D+05, &
            -.8902978767070678D+06,.5310411010968522D+08, &
            -.4043620325107754D+10,.3827011346598605D+12, &
            -.4406481417852278D+14,.6065091351222699D+16/

            K0=12
            IF (X.GE.35.0) K0=10
            IF (X.GE.50.0) K0=8
            T1=X-0.25D0*PI
            P0=1.0D0
            Q0=-0.125D0/X
            DO 45 K=1,K0
                P0=P0+A(K)*X**(-2*K)
45              Q0=Q0+B(K)*X**(-2*K-1)
            CU=DSQRT(RP2/X)
            BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
            BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
            T2=X-0.75D0*PI
            P1=1.0D0
            Q1=0.375D0/X
            DO 50 K=1,K0
                P1=P1+A1(K)*X**(-2*K)
50              Q1=Q1+B1(K)*X**(-2*K-1)
            CU=DSQRT(RP2/X)
            BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
            BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
        ENDIF
        DJ0=-BJ1
        DJ1=BJ0-BJ1/X
        DY0=-BY1
        DY1=BY0-BY1/X
        RETURN
    end subroutine bessel_func

end module fem_mod


!
! This module loads the experiment data.
!
module ReadInputs

contains

    subroutine read_inputs(param_filename, V, k, V_err, V_background, ntheta, nphi, npsi, Q, status2)
    ! param_filename=input file name containing initilizing parameters,such as temperature
    ! V=FEM intensity variance
    ! k=scattering vector
    ! V_err=FEM measurement variance error
    ! V_background=
    ! ntheta, nphi, npsi=rotation number for calculating simulated V
    ! Q=resolution for FEM
    ! status2=whether param_filename is opened with success or not

        implicit none
        character (len=*), intent(in) :: param_filename  !Assume size array, be careful
        real, pointer, dimension(:) :: v, k, v_err, v_background
        integer, intent(out) :: ntheta, nphi, npsi
        real, intent(out) :: q
        integer, intent(out) :: status2
        character (len=80) comment1  
        character (len=80) femfile ! The fem data file name, at most 80 characters
        integer filenamelength !The file name length in scatteringfile or femfile
        real indicator_end !Indicator_end=-1 means the end of file reaches
        integer status1 !Indicate the status of opening file in this subroutine
        logical file_end !Indicate fileend has reached
        integer num_line !Number of lines in each data file except comment line
        real, pointer, dimension(:) :: tempdata
        integer stat_allocate1, stat_allocate2, stat_allocate3,stat_allocate4 !Allocate status, 0 means success

        open(20, file=param_filename,iostat=status2, status='old')
        if(status2 .ne. 0) then ! Open fails
            print *, 'Cannot open file with name:  ', param_filename
            return
        endif
        read(20, '(a80)') comment1 ! Comment line

        num_line=0

        read(20, '(a)') femfile
        femfile = adjustl(femfile)
        read(20, *) nphi, npsi, ntheta
        read(20, *) q

        filenamelength=len_trim(femfile)
        file_end=.false.
        open(30,file=femfile(1:filenamelength),iostat=status1,status='old') 
        if(status1 .eq. 0) then ! Open succeeds
            read(30, '(a80)') comment1 ! First line is comment
            ! Count how many data pairs are in the file. -1 denotes EOF
            do while( .not. file_end)
                read(30, *) indicator_end
                if(abs(indicator_end+1.0) .le. 1e-6) then ! EOF is reached
                    exit ! Exit this loop
                else
                    num_line = num_line + 1
                endif
            enddo
            rewind(30) ! Go to the beginning of the file
            read(30, '(a80)') comment1
            allocate(tempdata(4*num_line),stat=stat_allocate1)
            ! Read k, V, and V_err data.
            ! Read k first, then V, and last V_err.
            ! First line is comment
            ! Count how many data pairs are in the file first.
            ! -1 denotes the last point
            if(stat_allocate1 .eq. 0) then
                read(30, *) tempdata
                allocate(k(num_line), stat=stat_allocate2)
                allocate(v(num_line), stat=stat_allocate3)
                allocate(v_err(num_line), stat=stat_allocate4)
                allocate(v_background(num_line))

                if ((stat_allocate2 .eq. 0) .and. (stat_allocate3 .eq. 0) .and. (stat_allocate4 .eq. 0)) then
                    k=tempdata(1:4*num_line:4)
                    v=tempdata(2:4*num_line:4)
                    v_err=tempdata(3:4*num_line:4)
                    v_background=tempdata(4:4*num_line:4)
                else
                    print *, 'An error occured reading the fem file.'
                    return
                endif ! allocate2, 3 and 4
            else
                print *, 'Fem allocation fails.'
                return
            endif
        deallocate(tempdata)
        else
            print *, 'Open fem file fails ', femfile(1:filenamelength)
        endif
        close(30)
        close(20)
    end subroutine read_inputs

end module readinputs


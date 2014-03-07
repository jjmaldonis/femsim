!
!This module contains the main functions for the model.
!

module model_mod
    use RMC_Global  ! Global variables
    implicit none
    ! derived data type for the hutches.  at is a pointer to a list
    ! of the indices
    ! of the atoms in this hutch.  nat is the number of atoms in
    ! this hutch
    type hutch
        integer, dimension(:), allocatable :: at
        integer :: nat
    end type hutch

    ! derived data type for the hutch array: contains an array of hutches,
    ! plus supporting information
    type hutch_array
        ! array of hutch objects
        type(hutch), dimension(:,:,:), pointer :: h
        ! number of hutches in x, y, and z
        integer :: nhutch_x, nhutch_y, nhutch_z
        ! physical size of a hutch in Angstroms
        real :: hutch_size
        ! list of the hutch indices for every atom. we don't use it so im
        ! getting rid of it. Jason 20130729
        integer, pointer, dimension(:,:) :: atom_hutch
    end type hutch_array

    ! adjustable-size list of atom indices
    type index_list
        integer :: nat
        integer, allocatable, dimension(:) :: ind
    end type index_list

    type real_index_list
        integer :: nat
        real, allocatable, dimension(:) :: ind
    end type real_index_list


    ! list for a rotated model that is the original models natoms long.  each atom in the list
    ! list is itself a list of the indices in the new, rotated model of the atoms corresponding
    ! to the the original atom in the unrotated model.
    ! for a rotated model, a list of the indices in the rotated model corresponding to each atom
    ! index in the original, unrotated model.  the original, unrotated model's index long.
    ! Defined type for a structural model with atoms positions and a bunch of metadata
    type model
        integer :: natoms                              ! number of atoms in the model
        type(real_index_list) :: xx, yy, zz      ! atom positions in Angstroms
        type(index_list) :: znum, znum_r         ! atom atomic numbers, and reduced z numbners
        real :: lx, ly, lz                             ! box size, in Angstroms
        integer :: nelements                           ! # of elements in the model
        integer, allocatable, dimension(:) :: atom_type    ! array listing atomic numbers present
        real, allocatable, dimension(:) :: composition     ! fractional composition in the order of atom_type
        type(hutch_array) :: ha                        ! hutch data structure
        logical :: rotated                             ! TRUE if model has been rotated, FALSE otherwise
        integer :: unrot_natoms
        type(index_list), dimension(:), allocatable :: rot_i ! list of which atoms in the rotated model correspond
        ! to the index i in the unrotated model
    end type model

contains

    subroutine read_model(model_filename, comment, m, istat)
    ! Reads a model in the Kirkland .xyz file format from the file model_filename.
    ! Puts the first-line comment in "comment", the model in m, and returns 0 in
    ! istat if the file cant be opened or the memory allocation fails.
        implicit none
        character (len=*),intent(in) :: model_filename
        character (LEN=*),intent(out) :: comment
        type(model), intent(out) :: m
        integer, intent(out) :: istat      !0 for successful open, others for failure.
        integer :: i, j, atom_count=0, nat=0, atom_temp
        integer, dimension(103) :: elements=0
        real :: comp_temp

        ! Open file that contains the model information.
        open(1,file=model_filename,iostat=istat,status='old')
        call check_allocation(istat, "Error in opening flie, "//model_filename)

        read(1,*) ! Comment line.
        read(1,*) ! This line contains the box size (lx, ly, lz).
        ! Count how many atoms there are in the model.
        do while( atom_count .ne. -1)
            read(1,*) atom_count ! Read line.
            nat=nat+1.0
        enddo
        nat=nat-1.0

        rewind(1)

        ! Set the number of atoms in the model m and allocate space for each
        ! coordinate.
        m%natoms = nat
        allocate(m%xx%ind(nat), m%yy%ind(nat), m%zz%ind(nat), m%znum%ind(nat), stat=istat)
        m%xx%nat = nat
        m%yy%nat = nat
        m%zz%nat = nat
        m%znum%nat = nat
        m%znum_r%nat = nat
        call check_allocation(istat, 'Unable to allocate memory for the model being read.')

        ! Read in the first 80 characters of the comment
        read(1,'(a80)') comment
        ! Read in the box size.
        read(1,*) m%lx,m%ly,m%lz
        ! If the model is not a perfect cube then the rest of the calculations
        ! wont work, so we really should check that.
        if((m%lx /= m%ly) .or. (m%lx /= m%lz)) then
            write(*,*) "The model is not a cube and will work correctly. Exiting."
            return
        endif
        ! Read the atomic numbers and atom positions directly into the model.
        do i=1,nat
            read(1,*) m%znum%ind(i),m%xx%ind(i),m%yy%ind(i),m%zz%ind(i)
            ! If this atom has atomic number z, then increment the z position in
            ! the array elements. This counts the number of each atom type we have.
            elements(m%znum%ind(i)) = elements(m%znum%ind(i)) + 1
        enddo
        close(1)

        ! Count the number of elements we have in our model.
        m%nelements=0
        do i=1, 103
            if(elements(i) /= 0) then
                m%nelements = m%nelements + 1
            end if
        end do

        ! Note: nelements is usually between 1 and 5 so these loops are tiny.
        ! Set m%atom_type to contain the atom types; set m%composition to
        ! contain the percent composition (as a number between 0 and 1).
        allocate(m%atom_type(m%nelements), m%composition(m%nelements), stat=istat)
        call check_allocation(istat, 'Unable to allocate memory for m%atom_type and m%composition.')
        ! Initialize the composition array to 0.0
        m%composition = 0.0
        ! i corresponds to the atomic number.
        ! j is the next open position in composition and atom_type.
        j = 1
        do i=1, 103
            if(elements(i) /= 0) then
                ! If we reach a non-zero element in elements then there are
                ! atoms with atomic number i in the model. Append this atomc
                ! number to atom_types and calculate the fractional composition
                ! of this element in the model, storing it in m%composition.
                ! Increment j to move to the next open spot.
                m%atom_type(j) = i
                m%composition(j) = real(elements(i)) / real(m%natoms)
                j = j + 1
            end if
        end do
        ! Note that m%atom_type and m%composition are now linked. You can look
        ! at m%composition, but you still won't know which element has this
        ! fractional composition. You need to look at the same index in
        ! m%atom_type in order to figure this out.

        ! Sort atom_type by increasing atomic order. Re-order composition in the
        ! same way so the indices stay in sync. (Insertion sort)
        do i=1, m%nelements
            do j=1, i
                if( m%atom_type(i) < m%atom_type(j) ) then
                    atom_temp = m%atom_type(i)
                    comp_temp = m%composition(i)
                    m%atom_type(i) = m%atom_type(j)
                    m%composition(i) = m%composition(j)
                    m%atom_type(j) = atom_temp
                    m%composition(j) = comp_temp
                end if
            end do
        end do

        ! For each atom i, add a parameter znum_r(i) that corresponds to
        ! m%atom_type and m%composition for fast lookup.
        allocate(m%znum_r%ind(m%natoms), stat=istat)
        call check_allocation(istat, 'Unable to allocate memory for m%znum_r.')
        m%znum_r%ind = 0.0
        do i=1, m%natoms
            do j=1, m%nelements
                if(m%znum%ind(i) .eq. m%atom_type(j)) then
                    m%znum_r%ind(i) = j
                end if
            end do
        end do

        m%rotated = .FALSE.

        call recenter_model(0.0, 0.0, 0.0, m)

        call model_init_hutches(m, istat)
    end subroutine read_model

    subroutine recenter_model(xc, yc, zc, m)
    ! Shifts the atom positions in model so that the mid-point between the maximum and
    ! and minimum atom positions in each dimensions sits at the position (xc, yc, zc),
    ! measured in units of the model supercell.
        real, intent(in) :: xc, yc, zc
        type(model), intent(inout) :: m
        real :: xshift, yshift, zshift
        ! maxval calculates the maximum value in the array. There are some
        ! nice parameters for it described online by the way.
        xshift = xc*m%lx - (maxval(m%xx%ind) + minval(m%xx%ind))/2.0
        yshift = yc*m%ly - (maxval(m%yy%ind) + minval(m%yy%ind))/2.0
        zshift = zc*m%lz - (maxval(m%zz%ind) + minval(m%zz%ind))/2.0

        m%xx%ind = m%xx%ind+xshift
        m%yy%ind = m%yy%ind+yshift
        m%zz%ind = m%zz%ind+zshift
    end subroutine recenter_model

    subroutine model_init_hutches(m, status)
    ! Initializes the hutch_array ha within the model m. It calcualtes the
    ! hutch_size (based on the model box size and the parameter
    ! ATOMS_PER_HUTCH) and the number of hutches in the array (nhutch_x, nhutch_y,
    ! and hhutch_z). It then assigns all the atoms in the current model atom
    ! position arrays xa, ya, and za to the appropriate hutches.  It does NOT
    ! check whether ha has already been initialized, so this routine should
    ! NEVER be called more than once for the same hutch_array.
        type(model), intent(inout) :: m
        integer, intent(out) :: status
        integer :: istat, numhutches, hx, hy, hz, i
        ! Note: numhutches is not the total number of hutches, it is the number
        ! of hutches in each dimension. So numhutches^3 is the total.

        status = 0

        numhutches = anint( (m%natoms/ATOMS_PER_HUTCH)**(1./3.) )
        m%ha%hutch_size = m%lx / numhutches 
        !write (*,*) 'Hutch size is ',m%ha%hutch_size,' Angstroms.'
        !write (*,*) 'Number of hutch in each dimension is: ', numhutches

        m%ha%nhutch_x = numhutches 
        m%ha%nhutch_y = numhutches 
        m%ha%nhutch_z = numhutches 

        allocate(m%ha%h(m%ha%nhutch_x, m%ha%nhutch_y, m%ha%nhutch_z), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Cannot allocate memory for hutch algorithm.  Exiting.'
            status = 1
            return
        end if

        allocate(m%ha%atom_hutch(m%natoms, 3), stat=istat)
        if (istat /= 0) then
            write(*,*) 'Cannot allocate memory for hutch algorithm.  Exiting.'
            status = 1
            return
        end if
        m%ha%atom_hutch = 0

        ! These hutch atom arrays are allocated and initialized in
        ! hutch_add_atom. We just need to initialize them to empty 
        ! and nat to 0 so that we can add atoms to them correctly.
        do hx = 1, m%ha%nhutch_x
            do hy = 1, m%ha%nhutch_y
                do hz = 1, m%ha%nhutch_z
                    if(allocated(m%ha%h(hx,hy,hz)%at)) deallocate(m%ha%h(hx,hy,hz)%at)
                    m%ha%h(hx, hy, hz)%nat = 0
                end do
            end do
        end do

        ! Calculate which hutch each atom should be in and add it to that hutch.
        do i=1, m%natoms
            call hutch_position(m, m%xx%ind(i), m%yy%ind(i), m%zz%ind(i), hx, hy, hz)
            call hutch_add_atom(m, i, hx, hy, hz)
        end do
    end subroutine model_init_hutches


    subroutine hutch_position(m, xx, yy, zz, hx, hy, hz)
    ! returns the indices of the hutch that encompasses position (xx, yy, zz) in
    ! the hutch_array in the integers (hx, hy, hz).  It assumes that the model 
    ! extends from -lx/2 to lx/2, -ly/2 to ly/2 and -lz/2 to lz/2 and does no
    ! error checking.
        type(model), intent(in) :: m 
        real, intent(in) :: xx, yy, zz
        integer, intent(out) :: hx, hy, hz

        ! This makes the range of hx, hy, and hz from 0 to nhutch_i, however
        ! the only time one of them will be 0 is if the position is exactly on
        ! the left edge. Thats what the next set of if statements is for: If 
        ! they are on an edge just move them a little bit over. Technically you 
        ! can get statistically more atoms in hutches on the 3 "left" edges but 
        ! it will happen extremely rarely so it wont matter. By the time we 
        ! are done hx, hy, and hz are restrained from 1 to nhutch_i so we can 
        ! convienently put them in an array.
        hx = ceiling( (xx + 0.5*m%lx) / m%ha%hutch_size )
        hy = ceiling( (yy + 0.5*m%ly) / m%ha%hutch_size )
        hz = ceiling( (zz + 0.5*m%lz) / m%ha%hutch_size )

        if (hx == 0) hx = 1
        if (hy == 0) hy = 1
        if (hz == 0) hz = 1
    end subroutine hutch_position


    subroutine hutch_add_atom(m, atom, hx, hy, hz)
    ! Adds the atom with index atom to the hutch_array in hutch hx, hy, hz.
        type(model), target, intent(inout) :: m
        integer, intent(in) :: atom, hx, hy, hz
        integer :: nat
        integer, dimension(m%ha%h(hx, hy, hz)%nat+1) :: scratch_atoms
        integer, dimension(:,:), allocatable :: temp_atom_hutch
        type(hutch_array), pointer :: ha
        ha => m%ha
        ! ha%h(hx,hy,hz)%nat is set to 0 in a do loop in model_init_hutches,
        ! slightly before this function is called for each atom.
        nat = ha%h(hx,hy,hz)%nat
        if(nat > 0) then
            scratch_atoms(1:nat) = ha%h(hx, hy, hz)%at
            scratch_atoms(nat+1) = atom
            ! Reallocate with new size
            deallocate(ha%h(hx,hy,hz)%at)
            allocate(ha%h(hx,hy,hz)%at(1:nat+1)) ! +1 for extra atom
            ha%h(hx,hy,hz)%at = scratch_atoms
        else
            allocate(ha%h(hx,hy,hz)%at(1:1))
            ha%h(hx,hy,hz)%at(1) = atom
        end if

        ha%h(hx,hy,hz)%nat = nat+1
        ! Create space if there isnt already.
        ! I am lucky that this array allocation works.
        if( size(ha%atom_hutch) / 3 < atom ) then
            allocate(temp_atom_hutch( size(ha%atom_hutch) / 3, 3))
            temp_atom_hutch = ha%atom_hutch
            deallocate(ha%atom_hutch)
            allocate(ha%atom_hutch(m%natoms, 3))
            ha%atom_hutch = temp_atom_hutch
            deallocate(temp_atom_hutch)
        endif
        ha%atom_hutch(atom, 1) = hx
        ha%atom_hutch(atom, 2) = hy
        ha%atom_hutch(atom, 3) = hz
    end subroutine hutch_add_atom


    subroutine check_model(m, istat)
    ! simple error checking on a model.  Currently checks: are all the  
    ! atoms in the box?  Are all the atomic numbers between 1 and 103
    ! (the range for which Kirkland calculated electron scattering factors)
    ! More should be added as we think of it.
        type(model), intent(in) :: m 
        integer, intent(out) :: istat
        real xlen, ylen, zlen 
        istat = 0

        xlen = maxval(m%xx%ind) - minval(m%xx%ind)
        ylen = maxval(m%yy%ind) - minval(m%yy%ind)
        zlen = maxval(m%zz%ind) - minval(m%zz%ind)

        if ( xlen > m%lx ) then 
            write (*,*) 'Maximum x distance of ',xlen,' Ang exceeds box size ',m%lx,' Ang.'
            istat = 1
        end if

        if ( ylen > m%ly ) then 
            write (*,*) 'Maximum y distance of ',ylen,' Ang exceeds box size ',m%ly,' Ang.'
            istat = 1
        end if

        if ( zlen > m%lz ) then 
            write (*,*) 'Maximum z distance of ',zlen,' Ang exceeds box size ',m%lz,' Ang.'
            istat = 1
        end if

        if (minval(m%znum%ind) < 1) then 
            write (*,*) 'Minimum atomic number of ', minval(m%znum%ind, 1), 'is less than zero.'
            istat = 1
        end if

        if (maxval(m%znum%ind) > 103) then 
            write (*,*) 'Maximum atomic number of ', maxval(m%znum%ind, 1), 'is greater than 103.'
            istat = 1
        end if
    end subroutine check_model

    subroutine rotate_model(phi, psi, theta, min, mrot, istat)
        ! rotates model min by angles phi, psi, theta and puts the results in mrot. min is unchanged.
        real, intent(in) :: theta, phi, psi
        type(model), intent(in) :: min
        type(model), intent(out) :: mrot 
        integer, intent(out) :: istat
        real, dimension(3,3) :: r                         ! rotation matrix
        real :: cpsi, cphi, ctheta, sphi, spsi, stheta    ! sines and cosines of the angles
        integer :: i, j                                   ! loop counters
        real :: x, y, z                                   ! temporary positions
        real :: lx2, ly2, lz2                             ! half box sizes
        type(model) :: mt                                 ! temporary oversize model
        integer, dimension(:), allocatable :: orig_indices

        ! periodic continue mt to 3x3x3 of the original model
        istat = 0
        call periodic_continue_model(3, 3, 3, min, mt, .FALSE., istat)
        if (istat /= 0) return

        allocate(orig_indices(mt%natoms), stat=istat)
        if (istat /= 0) then 
           write (*,*) 'Memory allocation failure in rotate_model.'
           return
        endif

        do i=1,mt%natoms   !Array loop was temporarily changed due to error in visual fortran - Jinwoo Hwang
            if(mod(i,min%natoms) .eq. 0)then
                orig_indices(i) = min%natoms
            else
                orig_indices(i) = mod(i,min%natoms)
            endif
            if((orig_indices(i) .gt. min%natoms) .or. (orig_indices(i) .lt. 1)) then
                write(*,*) 'wrong here', i, orig_indices(i)
            endif
        enddo
        orig_indices = (/ (mod(i,min%natoms)+1, i=1,mt%natoms) /)

        ! generate the members of a 3x3 rotation matrix.  Use the Goldstein "x-convention"
        ! and Euler angles phi theta, psi.
        cpsi = cos(psi)
        cphi = cos(phi)
        ctheta = cos(theta)
        sphi = sin(phi)
        spsi = sin(psi)
        stheta = sin(theta)

        !phi ignored - JWH 09/02/09
        r(1,1) = cpsi
        r(1,2) = spsi
        r(1,3) = 0.0
        r(2,1) = -ctheta*spsi
        r(2,2) = ctheta*cpsi
        r(2,3) = stheta
        r(3,1) = stheta*spsi
        r(3,2) = -stheta*cpsi
        r(3,3) = ctheta

        ! Rotate the position vectors in mt (the temporary 3x3x3 model).
        do i=1,mt%natoms
            if(abs(mt%xx%ind(i)).le.1.2*sqrt(2.0)*min%lx/2)then
                if(abs(mt%yy%ind(i)).le.1.2*sqrt(2.0)*min%ly/2)then
                    if(abs(mt%zz%ind(i)).le.1.2*sqrt(2.0)*min%lz/2)then
                        x = mt%xx%ind(i)*r(1,1) + mt%yy%ind(i)*r(1,2) + mt%zz%ind(i)*r(1,3)
                        y = mt%xx%ind(i)*r(2,1) + mt%yy%ind(i)*r(2,2) + mt%zz%ind(i)*r(2,3)
                        z = mt%xx%ind(i)*r(3,1) + mt%yy%ind(i)*r(3,2) + mt%zz%ind(i)*r(3,3)
                        mt%xx%ind(i) = x
                        mt%yy%ind(i) = y
                        mt%zz%ind(i) = z
                        !write(1008,*)i, mt%znum_r%ind(i), mt%xx%ind(i), mt%yy%ind(i), mt%zz%ind(i)
                    endif
                endif
            endif
        end do

        ! Cut the temporary model back to the original box size.
        ! First count the atoms in the box.
        mrot%natoms = 0
        lx2 = min%lx / 2.0
        ly2 = min%ly / 2.0
        lz2 = min%lz / 2.0
        do i=1, mt%natoms
            if((mt%xx%ind(i) <= lx2 .AND. mt%xx%ind(i) >= -1.0*lx2) .and. &
               (mt%yy%ind(i) <= ly2 .AND. mt%yy%ind(i) >= -1.0*ly2) .and. &
               (mt%zz%ind(i) <= lz2 .AND. mt%zz%ind(i) >= -1.0*lz2)) then
                mrot%natoms = mrot%natoms + 1
            endif
        enddo
        ! Allocate memory for the new atoms.
        mrot%unrot_natoms = min%natoms
        allocate(mrot%xx%ind(mrot%natoms), mrot%yy%ind(mrot%natoms), mrot%zz%ind(mrot%natoms), &
            mrot%znum%ind(mrot%natoms),  mrot%rot_i(mrot%unrot_natoms), mrot%znum_r%ind(mrot%natoms), stat=istat) !add mrot%znum_r here by Feng Yi on 03/19/2009
        mrot%xx%nat = mrot%natoms
        mrot%yy%nat = mrot%natoms
        mrot%zz%nat = mrot%natoms
        mrot%znum%nat = mrot%natoms
        mrot%znum_r%nat = mrot%natoms
        call check_allocation(istat, 'Problem allocating memory in rotate_model.')

        do i=1,mrot%unrot_natoms
           mrot%rot_i(i)%nat = 0
           if(allocated(mrot%rot_i(i)%ind)) deallocate(mrot%rot_i(i)%ind)
        enddo

        ! now copy just the atoms inside the original box size 
        ! from the temp model to the rotated one.
        j=1
        do i=1, mt%natoms
            if (mt%xx%ind(i) <= lx2 .AND. mt%xx%ind(i) >= -1.0*lx2) then
                if (mt%yy%ind(i) <= ly2 .AND. mt%yy%ind(i) >= -1.0*ly2) then
                    if (mt%zz%ind(i) <= lz2 .AND. mt%zz%ind(i) >= -1.0*lz2) then
                        mrot%xx%ind(j) = mt%xx%ind(i)
                        mrot%yy%ind(j) = mt%yy%ind(i)
                        mrot%zz%ind(j) = mt%zz%ind(i)
                        mrot%znum%ind(j) = mt%znum%ind(i)
                        mrot%znum_r%ind(j) = mt%znum_r%ind(i) !Added by Feng Yi on 03/19/2009   !Bug fixed : j to i -JWH 09/03/09
                        ! add_index is basically just the general 
                        ! "append(list, element)" function except 
                        ! it takes a type object containing the list 
                        ! and an int equal to its size.
                        call add_index(mrot%rot_i(orig_indices(i)), j)
                        j = j+1
                    endif
                endif
            endif
        enddo

        !release the memory allocated to mt
        deallocate(mt%atom_type, mt%composition)
        deallocate(mt%znum%ind,mt%znum_r%ind, mt%xx%ind, mt%yy%ind, mt%zz%ind)

        ! set the rest of of the rotated model paramters
        mrot%lx = min%lx
        mrot%ly = min%ly
        mrot%lz = min%lz
        mrot%rotated = .TRUE.
        mrot%unrot_natoms = min%natoms

        if(mrot%natoms .ne. 0) then !added by jwh 03/26/2009
            call composition_model(mrot) ! have to recalculate this because the # of atoms may have changed a little
        endif

        call model_init_hutches(mrot, istat)

        if(allocated(orig_indices)) then !added by feng yi on 3/14/2009
            deallocate(orig_indices)
        endif
    end subroutine rotate_model


    subroutine destroy_model(m)
    ! Deallocates all the various allocatable arrays and sub-arrays in a model.
        type(model), intent(inout) :: m 
        deallocate(m%xx%ind, m%yy%ind, m%zz%ind, m%znum%ind, m%znum_r%ind)
        if(allocated(m%atom_type)) deallocate(m%atom_type)
        if(allocated(m%composition)) deallocate(m%composition)
        call destroy_hutch(m%ha)
        call destroy_rot_indices(m%unrot_natoms, m%rot_i)
    end subroutine destroy_model

    subroutine destroy_hutch(ha)
    ! Deallocates the hutch_array ha and the atom lists inside it.
    ! Used by destroy_model.
        type(hutch_array), intent(inout) :: ha
        integer i, j, k
        if(associated(ha%h)) then
            do i=1,ha%nhutch_x
                do j=1,ha%nhutch_y
                    do k=1,ha%nhutch_z
                        if(ha%h(i,j,k)%nat .gt. 0) then
                            deallocate(ha%h(i,j,k)%at)
                        endif
                    enddo
                enddo
            enddo
            deallocate(ha%h, ha%atom_hutch)
        endif !if associated(ha%h)
        ! I wonder if there is a memory leak because we dont actually delete the
        ! rest of the ha variables and ha itself. TODO
        ! This would occur for all the rot_atom models in fem_update.
    end subroutine destroy_hutch

    subroutine destroy_rot_indices(unrot_natoms, ri)
    ! Deallocates all of the allocatable arrays and sub-arrays in an index list.
    ! Used by destroy_model.
        integer, intent(in) :: unrot_natoms
        type(index_list), allocatable, dimension(:) :: ri
        integer :: i
        do i=1, unrot_natoms
            if(ri(i)%nat .gt. 0) then  !added by feng yi
                deallocate(ri(i)%ind)
            endif
        enddo
        if(allocated(ri)) deallocate(ri)
    end subroutine destroy_rot_indices


    subroutine add_index(il, i)
        type(index_list), intent(inout) :: il
        integer, intent(in) :: i
        integer, dimension(:), allocatable :: scratch
        if( il%nat >= 1 ) then
            ! If there is space no need to reallocate. If not, reallocate.
            if(size(il%ind) .ge. il%nat+1) then
                if(il%nat == -1) il%nat = 0 ! We set old_index(i) to -1 sometimes
                il%nat = il%nat + 1
                il%ind(il%nat) = i
            else
                allocate(scratch(il%nat))
                scratch = il%ind
                il%nat = il%nat + 1
                deallocate(il%ind)
                allocate(il%ind( il%nat + 1 ))
                il%ind(1:il%nat-1) = scratch
                il%ind(il%nat) = i
            endif
        else
            il%nat = 1
            allocate(il%ind(1))
            il%ind(1) = i
        endif
        if(allocated(scratch)) then
            deallocate(scratch)
        endif
    end subroutine add_index


    subroutine composition_model(m)
    ! Calculates the composition of the model and fills in nelements, atom_type,
    ! and composition.
        type(model), intent(inout) :: m
        integer, dimension(103) :: znum_list
        integer :: i, j, isnew
        integer temp1
        real temp2 ! for temporary storage

        m%nelements=1
        znum_list(1) = m%znum%ind(1)
        do i=1,m%natoms
            isnew = 1
            do j=1,m%nelements
                ! If atom i's atomic number is already in the list, don't add
                ! its atomic number to the list again.
                if(m%znum%ind(i) == znum_list(j)) isnew = 0
            enddo
            if (isnew /= 0) then
                m%nelements=m%nelements+1
                znum_list(m%nelements) = m%znum%ind(i)
            endif
        enddo

        if( .not. allocated(m%atom_type) ) then
            allocate(m%atom_type(m%nelements))
        endif
        if( .not. allocated(m%composition) ) then
            allocate(m%composition(m%nelements))
        endif
        m%atom_type = znum_list(1:m%nelements)
        m%composition = 0.0

        do i = 1, m%natoms
            do j=1,m%nelements
                if(m%atom_type(j) == m%znum%ind(i)) then
                    m%composition(j) = m%composition(j) + 1.0
                    cycle
                endif
            enddo
        enddo

        m%composition = m%composition / real(m%natoms)

        ! Floating bubble method
        do i=1, (size(m%atom_type)-1)
            do j=1, (size(m%atom_type)-i)
                if(m%atom_type(j) .gt. m%atom_type(j+1)) then
                    temp1 = m%atom_type(j)
                    m%atom_type(j) = m%atom_type(j+1)
                    m%atom_type(j+1) = temp1

                    temp2 = m%composition(j)
                    m%composition(j) = m%composition(j+1)
                    m%composition(j+1) = temp2
                endif
            enddo
        enddo
    end subroutine composition_model


    subroutine periodic_continue_model(xp, yp, zp, min, mout, init_hutch, istat)
    ! Makes (xp, yp, zp) copies of the input model min and puts them in the output model
    ! mout.  Returns non-zero in istat is the memory can't be allocated.
        integer, intent(in):: xp, yp, zp
        type(model), intent(in) :: min
        type(model), intent(out) :: mout
        logical, intent(in) :: init_hutch
        integer, intent(out) :: istat
        integer :: i, j, k, c
        real :: shift_x, shift_y, shift_z

        mout%natoms = min%natoms*xp*yp*zp
        allocate(mout%xx%ind(mout%natoms), mout%yy%ind(mout%natoms), mout%zz%ind(mout%natoms), &
             mout%znum%ind(mout%natoms), mout%znum_r%ind(mout%natoms),stat=istat) !modified by Feng Yi on 03/19/2009
        mout%xx%nat = mout%natoms
        mout%yy%nat = mout%natoms
        mout%zz%nat = mout%natoms
        mout%znum%nat = mout%natoms
        mout%znum_r%nat = mout%natoms
        call check_allocation(istat, 'Error allocating memory for the periodic continued model.')

        mout%lx = min%lx*real(xp)
        mout%ly = min%ly*real(yp)
        mout%lz = min%lz*real(zp)

        c=0
        do i = -(xp-1)/2, (xp-1)/2     !jwh fyi 040809
            shift_x = real(i)*min%lx
            do j = -(yp-1)/2, (yp-1)/2
                shift_y = real(j)*min%ly
                do k = -(zp-1)/2, (zp-1)/2
                    shift_z = real(k)*min%lz
                    mout%xx%ind(c*min%natoms+1:(c+1)*min%natoms) = min%xx%ind + shift_x
                    mout%yy%ind(c*min%natoms+1:(c+1)*min%natoms) = min%yy%ind + shift_y
                    mout%zz%ind(c*min%natoms+1:(c+1)*min%natoms) = min%zz%ind + shift_z
                    mout%znum%ind(c*min%natoms+1:(c+1)*min%natoms) = min%znum%ind
                    mout%znum_r%ind(c*min%natoms+1:(c+1)*min%natoms) = min%znum_r%ind  !added by Feng Yi on 03/19/2009
                    c = c+1
                end do
            end do
        end do

        mout%nelements = min%nelements
        allocate(mout%atom_type(mout%nelements), mout%composition(mout%nelements), stat=istat)
        call check_allocation(istat, 'Problem allocating memory in periodic_continue_model.')
        mout%atom_type = min%atom_type
        mout%composition = mout%composition

        if(init_hutch) then
            call model_init_hutches(mout, istat)
            call check_allocation(istat, 'Cannot allocate memeory for the new hutch_array.')
        endif
    end subroutine periodic_continue_model


    subroutine hutch_list_pixel_sq(m, px, py, diameter, atoms, istat)
        type(model), target, intent(in) :: m
        real, intent(in) :: px, py, diameter
        integer, pointer, dimension(:) :: atoms !output of atom indices
        integer, intent(out) :: istat
        integer :: nh           ! number of hutches corresponding to diameter
        integer :: nlist        ! number of atoms in list
        integer :: i, j, k      ! counting variables
        integer, dimension(:), allocatable, target :: temp_atoms
        integer :: i_start, i_end, j_start, j_end, trash

        !write(*,*) "Number of hutches in the x, y, and z directions:", m%ha%nhutch_x, m%ha%nhutch_y, m%ha%nhutch_z
        allocate(temp_atoms(m%natoms), stat=istat)
        if (istat /= 0) then
            write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel'
            return
        end if

        ! Calculate the start and end hutch positions.
        call hutch_position(m, px-diameter/2.0, py-diameter/2.0, 0.0, i_start, j_start, trash)
        call hutch_position(m, px+diameter/2.0, py+diameter/2.0, 0.0, i_end, j_end, trash)
        nh = (i_end-i_start+1)*(j_end-j_start+1)*(m%ha%nhutch_z)
        
        ! Fill in the list.
        nlist = 1
        do i = i_start, i_end
            do j = j_start, j_end
                do k = 1, m%ha%nhutch_z
                    if(m%ha%h(i, j, k)%nat /= 0) then
                        temp_atoms(nlist:nlist+m%ha%h(i, j, k)%nat-1) = m%ha%h(i, j, k)%at(1:m%ha%h(i, j, k)%nat)
                        nlist = nlist + m%ha%h(i, j, k)%nat
                    endif
                enddo
            enddo
        enddo

        ! Assign atoms to the subset of temp_atoms that was filled in.
        if( nlist > 1 ) then
            allocate(atoms(nlist-1), stat=istat)
            if (istat /= 0) then
                write (*,*) 'Unable to allocate memory for atom indices in hutch_list_pixel.'
                return
            endif
            atoms = temp_atoms
        else
            nullify(atoms)
            istat = -1
        endif

        !write(*,*) "pixel (", px,py, ") has diameter", diameter, "and contains", nlist, "atoms and ", nh, &
            !"hutches !<= ", ( (ceiling(diameter/m%ha%hutch_size)+1) * (ceiling(diameter/m%ha%hutch_size)+1) * 11 ) ! debug

        if(allocated(temp_atoms)) deallocate(temp_atoms)
    end subroutine hutch_list_pixel_sq


    subroutine check_allocation(istat, message)
        integer, intent(in) :: istat
        character(len=*), intent(in) :: message
        if (istat /= 0) then
            write (*,*) message
            return
        endif
    end subroutine check_allocation

end module model_mod

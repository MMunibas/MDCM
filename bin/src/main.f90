!///////////////////////////////////////////////////////////////////////////////
!
!      main.f90
!      Created: 9 May 2016 at 09:54
!      Author: Oliver Unke
!      Updated and extended by Mike Devereux
!
!///////////////////////////////////////////////////////////////////////////////
program cubefit
use differential_evolution
use symmetry
implicit none

real(rp), parameter :: bohr2angstrom      = 0.52917721067_rp
real(rp), parameter :: angstrom2bohr      = 1._rp/bohr2angstrom
real(rp), parameter :: hartree2kcal       = 627.509_rp
real(rp), parameter :: vbig               = 1.e99_rp ! huge() not working well
!https://de.wikipedia.org/wiki/Van-der-Waals-Radius
real(rp), parameter :: vdW_scaling = 1.0_rp/3.0_rp ! r must be smaller than scaling*vdW_radius to be feasible
real(rp), dimension(55), parameter :: vdW_radius = (/ &
                                     1.20_rp*angstrom2bohr, & ! H
                                     1.40_rp*angstrom2bohr, & ! He
                                     1.82_rp*angstrom2bohr, & ! Li
                                     1.53_rp*angstrom2bohr, & ! Be
                                     1.92_rp*angstrom2bohr, & ! B
                                     1.70_rp*angstrom2bohr, & ! C
                                     1.55_rp*angstrom2bohr, & ! N
                                     1.52_rp*angstrom2bohr, & ! O
                                     1.47_rp*angstrom2bohr, & ! F
                                     1.54_rp*angstrom2bohr, & ! Ne
                                     2.27_rp*angstrom2bohr, & ! Na
                                     1.73_rp*angstrom2bohr, & ! Mg
                                     1.84_rp*angstrom2bohr, & ! Al
                                     2.10_rp*angstrom2bohr, & ! Si
                                     1.80_rp*angstrom2bohr, & ! P
                                     1.80_rp*angstrom2bohr, & ! S
                                     1.75_rp*angstrom2bohr, & ! Cl
                                     1.88_rp*angstrom2bohr, & ! Ar
                                     2.75_rp*angstrom2bohr, & ! K
                                     2.31_rp*angstrom2bohr, & ! Ca
                                     0.00_rp*angstrom2bohr, & ! Sc
                                     0.00_rp*angstrom2bohr, & ! Ti
                                     0.00_rp*angstrom2bohr, & ! V
                                     0.00_rp*angstrom2bohr, & ! Cr
                                     0.00_rp*angstrom2bohr, & ! Mn
                                     0.00_rp*angstrom2bohr, & ! Fe
                                     0.00_rp*angstrom2bohr, & ! Co
                                     1.63_rp*angstrom2bohr, & ! Ni
                                     1.40_rp*angstrom2bohr, & ! Cu
                                     1.39_rp*angstrom2bohr, & ! Zn
                                     1.87_rp*angstrom2bohr, & ! Ga
                                     2.11_rp*angstrom2bohr, & ! Ge
                                     1.85_rp*angstrom2bohr, & ! As
                                     1.90_rp*angstrom2bohr, & ! Se
                                     1.85_rp*angstrom2bohr, & ! Br
                                     2.02_rp*angstrom2bohr, & ! Kr
                                     3.03_rp*angstrom2bohr, & ! Rb
                                     2.49_rp*angstrom2bohr, & ! Sr
                                     0.00_rp*angstrom2bohr, & ! Y
                                     0.00_rp*angstrom2bohr, & ! Zr
                                     0.00_rp*angstrom2bohr, & ! Nb
                                     0.00_rp*angstrom2bohr, & ! Mo
                                     0.00_rp*angstrom2bohr, & ! Tc
                                     0.00_rp*angstrom2bohr, & ! Ru
                                     0.00_rp*angstrom2bohr, & ! Rh
                                     1.63_rp*angstrom2bohr, & ! Pd
                                     1.72_rp*angstrom2bohr, & ! Ag
                                     1.58_rp*angstrom2bohr, & ! Cd
                                     1.93_rp*angstrom2bohr, & ! In
                                     2.17_rp*angstrom2bohr, & ! Sn
                                     2.06_rp*angstrom2bohr, & ! Sb
                                     2.06_rp*angstrom2bohr, & ! Te
                                     1.98_rp*angstrom2bohr, & ! I
                                     2.16_rp*angstrom2bohr, & ! Xe
                                     0.00_rp*angstrom2bohr  &
                                     /)

! program parameters
logical  :: verbose           = .false. ! toggles printing mode
logical  :: use_greedy_fit    = .false. ! greedily fit individual multipoles first
logical  :: greedy_only_multi = .false. ! stops greedy fit after fitting atomic multipoles
logical  :: simplex_only      = .false. ! perform only simplex refinement, no D.E.
logical  :: use_symmetry      = .false. ! toggles symmetry constraint mode
logical  :: fit_multipoles    = .false. ! fit multipoles instead of charges
logical  :: generate_atomic   = .false. ! for visualizing charge fits
logical  :: refine_solution   = .false. ! when used to refine a solution
integer, parameter :: lmax    = 5       ! where is the multipole expansion truncated? maximum = 5!
integer  :: lstart            = 0       ! where do we start the fitting of the multipole expansion?
integer  :: lstop             = lmax    ! where do we start the fitting of the multipole expansion?
integer  :: lcur              = 0       ! where is the multipole expansion truncated? maximum = 5!
integer  :: num_charges       = 1       ! how many charges to fit (current)
integer  :: num_charges_min   = 2       ! maximum number of charges to fit
integer  :: num_charges_max   = 2       ! maximum number of charges to fit
integer  :: num_charges_min_multipole = 1 
integer  :: num_charges_max_multipole = 5
integer  :: num_trials        = 1       ! maximum number of trials per number of charges
real(rp) :: total_charge      = 0._rp   ! total charge
real(rp) :: total_charge2     = 0._rp   ! total charge
real(rp) :: max_charge        = 1._rp  ! maximum allowed absolute value of a charge (in search range)
real(rp) :: max_extend        = 5._rp   ! gets calculated automatically from vdW radii

character(len=1024) :: input_esp_cubefile = '', &
                       compare_esp_cubefile = '', &
                       input_multipolefile = '', &
                       input_xyzfile = '', &
                       input_density_cubefile = '', &
                       prefix = '', &
                       dummystring = '', dummystring2 = '', dummystring3 = '' ! input files
integer :: ppos !for use with scan()

! other program "modes", only used for analysis or generating cube files at better resolution
logical :: analysis_mode = .false. ! analysis mode: compare two cube files and compute RMSDs
logical :: generate_mode = .false. ! for generating esp cube file from charges (at higher resolution)


!how to "prune" grid points
logical  :: use_vdW_grid_cutoff = .true.
real(rp) :: vdw_grid_min_cutoff = 1.20_rp      !radius smaller than this*vdw_radius gets ignored
real(rp) :: vdw_grid_max_cutoff = 2.20_rp      !radius larger than this*vdw_radius gets ignored
logical  :: use_density_grid_cutoff = .false.
real(rp) :: density_grid_min_cutoff = 3.162277660e-4_rp !density smaller than this is ignored
real(rp) :: density_grid_max_cutoff = 2.0e-3 ! 1.0e-3_rp !density larger than this is ignored


! information about the ESP
integer :: Natom ! number of atoms
real(rp), dimension(3)                :: origin, axisX, axisY, axisZ   ! coordinate system of grid
integer                               :: NgridX, NgridY, NgridZ, Ngrid ! number of grid points in each direction
real(rp)                              :: Ngridr                        ! total number of grid points as real (for mean etc.)
real(rp), dimension(:),   allocatable :: esp_grid, esp_grid2           ! values of the electrostatic potential
real(rp), dimension(:,:), allocatable :: gridval                       ! stores at which gridpoints the ESP is interesting (outside atoms)
integer,  dimension(:),   allocatable :: atom_num                      ! stores the atomic numbers of the atoms
real(rp), dimension(:,:), allocatable :: atom_pos                      ! stores the atomic positions 
real(rp), dimension(3)                :: atom_com                      ! stores the atomic center of mass (for checking charge symmetry)
real(rp), dimension(:), allocatable :: multipole                       ! stores multipole parameters

real(rp), dimension(:,:,:), allocatable :: multipole_solutions         ! stores charge positions for the multipole fits (num_charges_max_multipole*4,num_charges_max_multipole,Natom)
real(rp), dimension(:), allocatable :: mapsol                          ! temporary charge array
real(rp), dimension(:,:,:), allocatable :: sym_multipole_solutions     ! stores charge positions for the multipole fits with symmetry constraints (num_charges_max_multipole*4,num_charges_max_multipole,Natom)
real(rp), dimension(:,:), allocatable   :: multipole_solutions_rmse    ! stores RMSE of the multipole solutions (needed for greedy mode)
integer, dimension(:,:), allocatable    :: num_sym_solution_ops                        ! stores number of symmetry ops for best fit to each number of atomic charges
integer, dimension(:,:,:), allocatable  :: sym_solution_ops                          ! stores symmetry ops for best fit to each number of atomic charges
real(rp), dimension(:), allocatable :: multipole_best, multipole_save  ! current best charge fit in the multipole fit
!real(rp), dimension(:,:,:), allocatable, save :: symmetry_ops          ! stores the symmetry operations in matrix form

! these are only used for visualizing the grid with R in the end
real(rp), dimension(:,:), allocatable :: sliceXY,sliceXZ,sliceYZ       ! cuts along planes
real(rp), dimension(:,:), allocatable :: sliceXY2,sliceXZ2,sliceYZ2    ! cuts along planes
logical,  dimension(:,:), allocatable :: usedXY, usedXZ, usedYZ        ! is the grid point at that cut used in the fit?

! fitting variables
integer ::  natmfit = 0                               ! defines number of atoms included in fit
integer,  dimension(:),   allocatable :: fitatoms     ! contains list of atom indices to be fitted
integer :: fitAtm   ! contains atom currently being fitted to multipolar ESP
integer ::  num_symFitAtms = 0                        ! defines number of symmetry-unique atoms included in symmetry-constratined fit
integer,  dimension(:),   allocatable :: symFitAtms   ! contains list of symmetry-unique atom indices to be fitted in symmetry-constrained fit
real(rp), dimension(:),   allocatable :: charges      ! contains position (x,y,z) and magnitude (q) of charges
                                                      ! x1: charge(1), y1: charge(2), z1: charge(3), q1: charge(4)
                                                      ! x2: charge(5), ...
real(rp), dimension(:),   allocatable :: sym_charges  ! contains fitted parameters (distances along axes or in symmetry planes) for symmetry-constrained fits
real(rp), dimension(:),   allocatable :: bestcharges  ! best solution found so far
integer, dimension(:), allocatable   :: best_symatm_combo ! num chgs per atom for best symmetric charge guess
real(rp), dimension(:,:), allocatable :: search_range ! has a minimum and a maximum value for every entry of charges


! for error handling
integer :: ios
logical :: file_exists

integer :: i,j,k,a,a1,a2,b,try,qdim,sqdim,flag,lcheck !current dimensionality of charge

real(rp) :: tmp, RMSE_best, RMSE_tmp, MAE_tmp, maxAE_tmp, lrmse_best


call system("mkdir -p slices")

! read command line arguments
call read_command_line_arguments()

!make a folder to store output for cleanliness
if(trim(prefix) /= '') then
    call execute_command_line("mkdir -p "//trim(prefix))
end if 

! initialize RNG (for drawing the population)
call init_random_seed(0)
! read in the cube file
call read_cube_file(trim(input_esp_cubefile),trim(input_density_cubefile))

! analysis mode: evaluate quality of fit only
if(analysis_mode) then
    !save old values
    if(.not.allocated(esp_grid2)) allocate(esp_grid2(NgridX*NgridY*NgridZ))
    if(.not.allocated(sliceXY2)) allocate(sliceXY2(NgridX,NgridY))
    if(.not.allocated(sliceXZ2)) allocate(sliceXZ2(NgridX,NgridZ))
    if(.not.allocated(sliceYZ2)) allocate(sliceYZ2(NgridY,NgridZ))
    sliceXY2 = sliceXY
    sliceXZ2 = sliceXZ
    sliceYZ2 = sliceYZ
    esp_grid2 = esp_grid
    
    !read in the second esp file
    call read_cube_file(trim(compare_esp_cubefile),trim(input_density_cubefile))
    call do_analysis(vdw_grid_min_cutoff,vdw_grid_max_cutoff) !calculate all the errors in different regions
    call write_error_cube_file() !writes the difference between true and fittes esp file to a cube file 
    
    ! plot with R
    call write_image_slice_data_analysis() !for visualization with R
!    call execute_command_line("./visualize.r",wait=.true.)
    if(trim(prefix) /= '') then
        dummystring = trim(prefix)//"/analysis_comparison.png"
    else
        dummystring = "analysis_comparison.png"
    end if
    call execute_command_line("mv comparison.png "//trim(dummystring),wait=.true.)
    if(verbose) write(*,'(A)') "Did analysis and wrote slice data"
    stop
end if
! generate mode: generate new cube file only
if(generate_mode) then
    if(fit_multipoles) then
        !allocate memory
        if(.not.allocated(multipole))      allocate(multipole((lmax+1)**2*Natom))
        if(.not.allocated(multipole_best)) allocate(multipole_best((lmax+1)**2*Natom))
        call read_multipole_file(input_xyzfile)
        
        if(generate_atomic) then
            do a = 1,Natom
                do num_charges = 1,num_charges_max_multipole
                    !calculate esp grid and slice data for the multipole of atom a
                    call calc_multipole_grid_and_slice_data(multipole,a) 
                    !read the charges
                    write(dummystring, '(I0)') a
                    write(dummystring2,'(I0)') num_charges
                    input_xyzfile = "multipole"//trim(dummystring)//"_"//trim(dummystring2)//"charges.xyz"
                    if(trim(prefix) /= '') input_xyzfile  = trim(prefix)//"/"//trim(input_xyzfile)
                    call read_xyz_file()
                    ! plot with R
                    call write_image_slice_data(charges) !for visualization with R
!                    call execute_command_line("./visualize.r", wait=.true.)
                    dummystring3 = "multipole"//trim(dummystring)//"_"//trim(dummystring2)//"charges_comparison.png"
                    if(trim(prefix) /= '') dummystring3 = trim(prefix)//"/"//trim(dummystring3)
                    call execute_command_line("mv comparison.png "//trim(dummystring3),wait=.true.) 
                end do
            end do
        else                                     
            !decide filename
            select case (lcur)
                case(0) 
                    write(dummystring,'(A)') "monopole"
                case(1)
                    write(dummystring,'(A)') "dipole"
                case(2) 
                    write(dummystring,'(A)') "quadrupole"
                case(3)
                    write(dummystring,'(A)') "octopole"
                case(4) 
                    write(dummystring,'(A)') "hexadecapole"
                case(5)
                    write(dummystring,'(A)') "ditriantapole"
            end select       
            call write_cube_file_multipole(multipole)
            ! plot with R
            call write_image_slice_data_multipole(multipole) !for visualization with R
!            call execute_command_line("./visualize.r",wait=.true.)
            if(trim(prefix) /= '') then
                dummystring = trim(prefix)//"/"//trim(dummystring)//"_expansion_comparison.png"
            else
                dummystring = trim(dummystring)//"_expansion_comparison.png"
            end if
            call execute_command_line("mv comparison.png "//trim(dummystring),wait=.true.) 
            if(verbose) write(*,'(A)') "Written cubefile and slice data"  
        end if
    else
        call read_xyz_file()
        call write_cube_file(charges)
        ! plot with R
        call write_image_slice_data(charges) !for visualization with R
!        call execute_command_line("./visualize.r",wait=.true.)
        write(dummystring,'(I0)') num_charges 
        if(trim(prefix) /= '') then
            dummystring = trim(prefix)//"/"//trim(dummystring)//"charges_comparison.png"
        else
            dummystring = trim(dummystring)//"charges_comparison.png"
        end if
        call execute_command_line("mv comparison.png "//trim(dummystring),wait=.true.)
        if(verbose) write(*,'(A)') "Written cubefile and slice data"
    end if
    stop
end if

! initialize symmetry module
if(use_symmetry) then
        if(verbose) write(*,'(A)')
        if(verbose) write(*,'(A)') "Initializing symmetry module..."
        if(verbose) write(*,'(A)')
        call symmetry_init(Natom, transpose(atom_pos), real(atom_num(:),rp))
        atom_com = centerOfMass(transpose(atom_pos))
!        if(.not.allocated(symmetry_ops)) allocate(symmetry_ops(3,3,size(unique_ops,dim=1)))
!        do i = 1,size(unique_ops,dim=1)
!            symmetry_ops(:,:,i) = unique_ops(i)%M
!        end do
        if(verbose) write(*,'(A)')
        if(verbose) write(*,'(A)') "...done!"    
end if

if(verbose) write(*,*) 
if(verbose) write(*,*) 
if(verbose) write(*,*) 

!when fitting the atomic multipoles instead of charges
if(fit_multipoles) then
    if(use_symmetry) then ! should have been caught earlier during argument parsing!
      call throw_error('sorry, symmetry module does not yet work with multipole fits')
    end if
    write(*,'(A)') "FITTING ATOMIC MULTIPOLES"
    write(*,'(A)') "WARNING: THIS ROUTINE IS NOT EFFICIENT!"
    write(*,'(A)') "USE FOR TESTING PURPOSES ONLY!!"
    !allocate memory
    if(.not.allocated(multipole))      allocate(multipole((lmax+1)**2*Natom))
    if(.not.allocated(multipole_best)) allocate(multipole_best((lmax+1)**2*Natom))
    if(.not.allocated(multipole_save)) allocate(multipole_save((lmax+1)**2*Natom))
    if(.not.allocated(search_range))   allocate(search_range(2,size(multipole,dim=1)))
    multipole = 0._rp
    multipole_best = 0._rp
    lcheck = 0
    
    !initialize search_range
    do i = 1,size(multipole,dim=1)
        search_range(:,i) = (/-0.1_rp, 0.1_rp/)
    end do
    !initialize differential evolution
    !LOOK AT THESE VALUES AGAIN! POP SIZE AND MAXGENS MIGHT NOT REALLY BE OPTIMAL
    call DE_init(set_range            = search_range,    &
                 set_popSize          = 100,             &
                 set_maxGens          = 10,              &
                 set_maxChilds        = 1,               &
                 set_forceRange       = .false.,         &
                 set_mutationStrategy = DEtargettobest1, &
                 set_crossProb        = 1.d0,            &
                 set_verbose          = verbose,         &
                 set_Nprint     = 10)  
    do lcur = 0,lmax
        if(lcur > lstop) exit
        !decide filename
        select case (lcur)
            case(0) 
                write(dummystring,'(A)') "monopole"
            case(1)
                write(dummystring,'(A)') "dipole"
            case(2) 
                write(dummystring,'(A)') "quadrupole"
            case(3)
                write(dummystring,'(A)') "octopole"
            case(4) 
                write(dummystring,'(A)') "hexadecapole"
            case(5)
                write(dummystring,'(A)') "ditriantapole"
        end select
        
        !load RMSE_best from files
        dummystring2 = trim(dummystring)//"_expansion.txt"
        if(trim(prefix) /= '') dummystring2 = trim(prefix)//"/"//trim(dummystring2)
        !print*, trim(dummystring2)
        open(30, file=trim(dummystring2), action='read', status='old', iostat=ios)
        if(ios == 0) then
            read(30,*) dummystring2, dummystring3, RMSE_best
            read(30,*)
            i = 1
            do while(.true.)
                read(30,*,iostat=ios) dummystring2, dummystring3, multipole_best(i)
                if(ios /= 0) then
                    exit
                else
                    i = i + 1
                end if
            end do
            close(30)   
        else
            RMSE_best = vbig
        end if   
        
        if(lcur < lstart) cycle !skip multipoles we dont want to fit
    
        do try = 1,num_trials
            if(verbose) write(*,'(A,I0)') "Fitting "//trim(dummystring)//" expansion, try ", try
            call DE_optimize(rmse_multipole,feasible_multipole,sum_constr_multipole,multipole,&
                             guess=multipole_best) 
            RMSE_tmp = rmse_multipole(multipole)
            if(verbose) write(*,'(A,ES23.9,A)') "RMSE ", RMSE_tmp, " Hartree"
            if(RMSE_tmp < RMSE_best) then
                if(verbose) write(*,'(A)') "NEW BEST!"   
                RMSE_best = RMSE_tmp
                multipole_best = multipole
                        call write_multipole_file(multipole_best) !write output file
                ! plot with R
                call write_image_slice_data_multipole(multipole_best) !for visualization with R
!                call execute_command_line("./visualize.r",wait=.true.)
                dummystring2 = trim(dummystring)//"_expansion_comparison.png"
                if(trim(prefix) /= '') dummystring2 = trim(prefix)//"/"//trim(dummystring2)
                call execute_command_line("mv comparison.png "//trim(dummystring2),wait=.true.)   
            end if
        end do
    end do 
    call DE_exit()
    stop
end if



!when fitting atomic charges to the atomic multipoles
if(use_greedy_fit) then

    if(use_symmetry) then 
      if(.not.allocated(symFitAtms)) allocate(symFitAtms(Natom))
      ! first remove atoms that are redundant due to symmetry
      b=1
      do while (b.le.natmfit)
        a = fitatoms(b)
        flag=0
        do i=1,Natom 
          if(count(atom_sea(i,:) /= 0) == 0) exit
          if(a.eq.atom_sea(i,1)) then
            flag=1
            num_symFitAtms=num_symFitAtms+1
            symFitAtms(num_symFitAtms)=i
            b=b+1
          endif
        end do
        if(flag.eq.0)then
          if(verbose) then
            write(*,'(A,I0,A)') 'Atom ',a,' is redundant by symmetry, ignoring.'
            do i=b,natmfit-1
              fitatoms(i)=fitatoms(i+1)
            enddo
            natmfit=natmfit-1
            write(*,'(I0,A)') natmfit,' atoms remain'
          endif
        endif
      end do
      ! of the remaining symmetry-unique atoms, check which symmetry ops leave
      ! nuclear coordinates untouched as these need to be applied to charges
      ! during the fit
!      if(.not.allocated(num_sym_ops)) allocate(num_sym_ops(Natom))
      do b = 1,num_symFitAtms
        a = symFitAtms(b)
        i=getChargeOps(a)
      enddo ! b
      ! now initialize how many charges are spawned by each of these sym ops
      ! for charges placed along rotational axes, in mirror planes etc.
      do b = 1,num_symFitAtms
        a = symFitAtms(b)
        call get_chgs_spawned(a)
      enddo ! b

    endif ! use_symmetry

    write(*,'(A)') "USING GREEDY MODE"
    !save the values for the full ESP (individual multipoles will be saved into there)
    if(.not.allocated(esp_grid2)) allocate(esp_grid2(NgridX*NgridY*NgridZ))
    if(.not.allocated(sliceXY2))  allocate(sliceXY2(NgridX,NgridY))
    if(.not.allocated(sliceXZ2))  allocate(sliceXZ2(NgridX,NgridZ))
    if(.not.allocated(sliceYZ2))  allocate(sliceYZ2(NgridY,NgridZ))
    sliceXY2 = sliceXY
    sliceXZ2 = sliceXZ
    sliceYZ2 = sliceYZ
    esp_grid2 = esp_grid
    total_charge2 = total_charge
    
    !allocate memory for the multipole information and read it
    if(.not.allocated(multipole))                 allocate(multipole((lmax+1)**2*Natom))
    if(.not.allocated(multipole_solutions))       allocate( &
       multipole_solutions(num_charges_max_multipole*4,num_charges_max_multipole,Natom))
    if(.not.allocated(sym_multipole_solutions))       allocate( &
       sym_multipole_solutions(num_charges_max_multipole*4,num_charges_max_multipole,Natom))
    if(.not.allocated(multipole_solutions_rmse))  allocate( &
       multipole_solutions_rmse(0:num_charges_max_multipole,Natom))
    if(.not.allocated(multipole_best))            allocate( &
       multipole_best(num_charges_max_multipole*4))

    !read in multipole data
    call read_multipole_file(input_multipolefile)

    ! allocate memory
    if(.not.allocated(search_range)) allocate(search_range(2,4*num_charges_max_multipole), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    
    !loop over the individual atoms selected for fitting (all atoms by default)
    do b = 1,natmfit
        a = fitatoms(b)
        fitAtm = a
        !calculate the total charge of this multipole
        total_charge = multipole(a)
        
        !calculate esp grid and slice data for the multipole of atom a
        call calc_multipole_grid_and_slice_data(multipole,a) 

!        if(num_charges_max < 5) then
!            num_charges_max_multipole = num_charges_max
!        end if
        
        
        !calculate the RMSE of using no charge at all (needed later in sampling the population)        
        multipole_solutions_rmse(0,a) = 0._rp
        do k = 1,Ngrid
            multipole_solutions_rmse(0,a) = multipole_solutions_rmse(0,a) + esp_grid(k)**2
        end do
        multipole_solutions_rmse(0,a) = sqrt(multipole_solutions_rmse(0,a)/Ngridr)

        outer: do num_charges = num_charges_min_multipole,num_charges_max_multipole
            !current dimensionality
            qdim = 4*num_charges-1
        
            !check whether a fit exists already, if yes, we do nothing
            if(use_symmetry)then !read existing symmetry-constrained parameters
              if(.not.allocated(sym_solution_ops)) then
                allocate(num_sym_solution_ops(num_charges_max_multipole,Natom))
                allocate(sym_solution_ops(num_charges_max_multipole,Natom,&
                         num_charges_max_multipole))
              endif
              write(dummystring,'(I0)') a
              write(dummystring2,'(I0)') num_charges
              if(trim(prefix) /= '') then
                  dummystring = trim(prefix)//"/symatm"//trim(dummystring)//"_"//trim(dummystring2)//"charges.fit"
              else
                  dummystring = "symatm"//trim(dummystring)//"_"//trim(dummystring2)//"charges.fit"
              end if
              
              file_exists=read_sym_file(dummystring, &
                  multipole_solutions(:,num_charges,a),sym_solution_ops,&
                  num_sym_solution_ops,sqdim,num_charges,a,num_charges_max_multipole)
              if(file_exists)then
                !calculate the RMSE of the loaded solution
                multipole_solutions_rmse(num_charges,a) = sym_atm_rmse_qtot( &
                                  multipole_solutions(1:sqdim,num_charges,a))
                !print*, multipole_solutions_rmse(num_charges,a)
                if(verbose) then
                    write(*,'(A)') 'File "'//trim(dummystring)//'" already exists.'//&
                                                   " Fitting procedure is skipped."
  
                end if
                cycle
              endif
            endif
            dummystring=''
            dummystring2=''
            write(dummystring,'(I0)') a
            write(dummystring2,'(I0)') num_charges
            if(trim(prefix) /= '') then
                dummystring = trim(prefix)//"/multipole"//trim(dummystring)//"_"//trim(dummystring2)//"charges.xyz"
            else
                dummystring = "multipole"//trim(dummystring)//"_"//trim(dummystring2)//"charges.xyz"
            end if
            open(30,file=trim(dummystring), status='old', action='read',iostat = ios)
            if(ios == 0) then 
                !actually read in the solution
                read(30,*,iostat=ios) i
                if(ios /= 0 .or. i /= num_charges) call throw_error('Could not read "'//trim(dummystring)// &
                '" or the file contains the wrong number of charges.')
                read(30,*,iostat=ios) !skip comment line
                !read charge positions and charges
                do i = 1,4*num_charges,4
                    read(30,*,iostat=ios) dummystring2, multipole_solutions(i:i+3,num_charges,a)
                    if(ios /= 0) call throw_error('Could not read "'//trim(dummystring)// &
                    '" or the file contains the wrong number of charges.')
                    multipole_solutions(i:i+2,num_charges,a) = multipole_solutions(i:i+2,num_charges,a)*angstrom2bohr
                end do                          
                close(30)
                !calculate the RMSE of the loaded solution
                multipole_solutions_rmse(num_charges,a) = rmse_qtot(multipole_solutions(1:qdim,num_charges,a))
                !print*, multipole_solutions_rmse(num_charges,a)
                if(verbose) then
                    write(*,'(A)') 'File "'//trim(dummystring)//'" already exists.'//&
                                                   " Fitting procedure is skipped."
                                                   
                end if  
                cycle
            end if        
                
            ! initialize RMSE    
            multipole_solutions_rmse(num_charges,a) = vbig

            if(.not.use_symmetry)then
              ! general coding strategy: create a new parameter array based on axes and planes and fit that, then convert back to overwrite qdim with the best solution
  
              ! initialize search_range
              call init_search_range()
  
              call DE_init(set_range            = search_range(:,1:qdim), &
                           set_popSize          = 10*qdim,                &
                           set_maxGens          = 2000*num_charges,       &
                           set_crossProb        = 1.00_rp,                &
                           set_maxChilds        = 1,                      &
                           set_forceRange       = .false.,                &
                           set_mutationStrategy = DEtargettobest1,        &
                           set_verbose          = verbose,                &
                           set_Nprint           = 100)  
            endif
            do try = 1,num_trials
                if(verbose) write(*,'(3(A,I0))') "Starting fitting procedure for multipole expansion of atom ",a,&
                                                   " with ",num_charges," charges, trial ",try
                if(use_symmetry)then
                  ! (re)initialize symmetry-constrained search
                  call DE_exit()
                  if(.not.init_atm_sym_search(num_charges,a)) then
                    cycle outer !skip if no fit exists with this combination of number of charges and symmetry constraints
                  endif
                  call init_sym_search_range(search_range,symFitAtms,num_symFitAtms,&
                       max_extend,max_charge,sqdim)
    
                  call DE_init(set_range            = search_range(:,1:sqdim), &
                               set_popSize          = 10*sqdim,                &
                               set_maxGens          = 2000*num_charges,        &
                               set_crossProb        = 1.00_rp,                 &
                               set_maxChilds        = 1,                       &
                               set_forceRange       = .false.,                 &
                               set_mutationStrategy = DEtargettobest1,         &
                               set_verbose          = verbose,                 &
                               set_Nprint           = 100)
                  call DE_optimize(sym_atm_rmse_qtot,sym_atm_feasible,sym_sum_atm_constr,&
                        sym_multipole_solutions(1:sqdim,num_charges,a),&
                        init_pop=sym_init_pop_multipole)
                  call spawn_sym_chgs(sym_multipole_solutions(1:sqdim,num_charges,a),&
                        multipole_solutions(1:qdim,num_charges,a),&
                        num_charges,symFitAtms,num_symFitAtms,total_charge,.false.)
                else
                  call DE_optimize(rmse_qtot,feasible,sum_constr,&
                        multipole_solutions(1:qdim,num_charges,a),init_pop=init_pop_multipole)
                endif

                ! measure the quality of the fit
                RMSE_tmp = rmse_qtot(multipole_solutions(1:qdim,num_charges,a)) 
                if(verbose) write(*,'(A,ES23.9,A)') "RMSE ", RMSE_tmp, "Hartree"
                if(RMSE_tmp < multipole_solutions_rmse(num_charges,a)) then
                    if(verbose) write(*,'(A)') "NEW BEST!"                       
                    multipole_solutions_rmse(num_charges,a) = RMSE_tmp
                    multipole_best = multipole_solutions(:,num_charges,a)
                    MAE_tmp   = mae_qtot(multipole_solutions(1:qdim,num_charges,a))
                    maxAE_tmp = max_ae_qtot(multipole_solutions(1:qdim,num_charges,a))
                    ! write results to file        
                    call write_xyz_file(multipole_solutions(1:qdim,num_charges,a),a)
                    if(use_symmetry)then
                      if(allocated(mapsol)) deallocate(mapsol)
                      allocate(mapsol(num_charges*4))
                      !create xyz files for symmetry equivalent atoms
                      do i=1,Natom
                        if(count(atom_sea(i,:) /= 0) == 0) exit
                        a1=atom_sea(i,1)
                        if(a1.ne.a) cycle ! cycle if this is the wrong sea group for atom a
                        do j=2,Natom 
                          if(atom_sea(i,j) == 0) exit
                          a2=atom_sea(i,j)
                          call sym_map_sea_q_coords(&
                             multipole_solutions(1:qdim,num_charges,a),mapsol,a1,a2,&
                             num_charges)
                          call write_xyz_file(mapsol,a2)
                        enddo
                      enddo
                      call write_sym_file(sym_multipole_solutions(1:sqdim,num_charges,a),&
                                          a,num_charges)
                    endif
                    !call write_cube_file(multipole_solutions(1:qdim,num_charges,a),a)
                    ! plot with R
                    call write_image_slice_data(multipole_solutions(1:qdim,num_charges,a)) !for visualization with R
!                    call execute_command_line("./visualize.r",wait=.true.)
                    write(dummystring,'(I0)') a
                    write(dummystring2,'(I0)') num_charges 
                    if(trim(prefix) /= '') then
                        dummystring = trim(prefix)//"/multipole"//trim(dummystring)//&
                                      "_"//trim(dummystring2)//"charges_comparison.png"
                    else
                        dummystring = "multipole"//trim(dummystring)//"_"//trim(dummystring2)//"charges_comparison.png"
                    end if
                    call execute_command_line("mv comparison.png "&
                    //trim(dummystring),wait=.true.)
                    if(verbose) write(*,*)
                end if
            end do
            
            !load the best solution
            multipole_solutions(:,num_charges,a) = multipole_best
            
            !calculate the magnitude of the last charge
            multipole_solutions(qdim+1,num_charges,a) = 0._rp
            do i = 1,(num_charges-1)*4,4
                multipole_solutions(qdim+1,num_charges,a) = multipole_solutions(qdim+1,num_charges,a) + &
                                                            multipole_solutions(i+3,num_charges,a)            
            end do
            multipole_solutions(qdim+1,num_charges,a) = total_charge - multipole_solutions(qdim+1,num_charges,a)         

            RMSE_tmp  = rmse_qtot(multipole_solutions(1:qdim,num_charges,a))
            MAE_tmp   = mae_qtot(multipole_solutions(1:qdim,num_charges,a)) 
            maxAE_tmp = max_ae_qtot(multipole_solutions(1:qdim,num_charges,a))
            ! write results to console
            if(verbose) then
                write(*,'(A,I0,A)') "Best found solution for ",num_charges &
                  ," charges:"
                write(*,'(A,ES23.9,A)') "        RMSE ", &
                  RMSE_tmp*hartree2kcal," kcal/mol"
                write(*,'(A,ES23.9,A)') "         MAE ", &
                  MAE_tmp*hartree2kcal," kcal/mol"
                write(*,'(A,ES23.9,A)') "     max. AE ", &
                  maxAE_tmp*hartree2kcal," kcal/mol"
                if(.not.feasible(multipole_solutions(1:qdim,num_charges,a))) write(*,'(A)') "WARNING, SOLUTION IS NOT FEASIBLE!"
                write(*,*)
                write(*,'(4A14)') "x[bohr]","y[bohr]","z[bohr]","q[e]"
                tmp = 0._rp
                do i = 1,qdim,4
                    if(i+3 <= qdim) then
                        write(*,'(4F14.7)') multipole_solutions(i,num_charges,a), multipole_solutions(i+1,num_charges,a),&
                                            multipole_solutions(i+2,num_charges,a), multipole_solutions(i+3,num_charges,a)
                        tmp = tmp + multipole_solutions(i+3,num_charges,a)
                    else
                        write(*,'(4F14.7)') multipole_solutions(i,num_charges,a), multipole_solutions(i+1,num_charges,a),&
                                            multipole_solutions(i+2,num_charges,a), total_charge-tmp
                                            
                    end if
                end do
            end if
            if(verbose) write(*,*)
            
            ! clean up
            call DE_exit()                        
        end do outer
    end do   

    !load back the full grid
    sliceXY = sliceXY2
    sliceXZ = sliceXZ2
    sliceYZ = sliceYZ2
    esp_grid = esp_grid2
    total_charge = total_charge2
    
    if(allocated(esp_grid2))    deallocate(esp_grid2)
    if(allocated(sliceXY2))     deallocate(sliceXY2)
    if(allocated(sliceXZ2))     deallocate(sliceXZ2)
    if(allocated(sliceYZ2))     deallocate(sliceYZ2)
    if(allocated(search_range)) deallocate(search_range)    
    if(greedy_only_multi) stop 
end if

!Change reference MEP to the multipolar MEP if we're fitting a fragment rather
!than the whole molecule

!allocate memory for the multipole information and read it
if(.not.allocated(multipole))                 allocate(multipole((lmax+1)**2*Natom))
if(.not.allocated(multipole_solutions))       allocate( &
   multipole_solutions(num_charges_max_multipole*4,num_charges_max_multipole,Natom))
if(.not.allocated(multipole_solutions_rmse))  allocate( &
   multipole_solutions_rmse(0:num_charges_max_multipole,Natom))
if(.not.allocated(multipole_best))            allocate( &
   multipole_best(num_charges_max_multipole*4))

!read in multipole data
call read_multipole_file(input_multipolefile)

!!subtract ESP at each grid point due to multipoles of atoms excluded from fit
!outer: do a=1,Natom
!  do b=1, natmfit  !skip atoms in fragment to be fitted
!    i=fitatoms(b)
!    if(i == a) cycle outer
!  enddo
!  call subtract_atom_multipole_ESP_from_grid(multipole,a)
!enddo outer

!if atoms are excluded from the fit then fit to multipolar ESP of remaining
!atoms

if(natmfit < Natom.and..not.use_symmetry) then ! some atoms are excluded
  write(*,'(A)') 'Fragment fit detected'
  write(*,'(A)') 'Fitting to MTP grid, discarding reference MEP cube data'
  write(*,'(A)') 'Subsequent stats etc. are compared to MTP data, not the'
  write(*,'(A)') 'MEP from the reference cube file'
  ! first zero reference grid
  do i = 1,Ngrid
     esp_grid(i) = 0.d0
  end do
  ! now repopulate with multipolar ESP
  do b=1, natmfit  
    a=fitatoms(b)
    call add_atom_multipole_ESP_to_grid(multipole,a)
  enddo
  !subtract charge from atomic multipoles of atoms that are not being fitted from total
  !charge of molecule (to allow fitting of fragments)
  do b=1, natmfit
    a=fitatoms(b)
    total_charge=total_charge + multipole(a)
  enddo
  if(verbose) write(*,'(A,ES23.9)') "Total fragment charge (a.u.): ",total_charge
endif

!refining the solution
if(refine_solution) then
    ! symmetry not yet implemented for refinement
    if(use_symmetry)then
      call throw_error('Sorry, symmetry is not yet implemented for model refinement')
    endif
    open(30, file=trim(input_xyzfile), status="old", action="read", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "'//trim(input_xyzfile)//'" for reading')
    read(30,*,iostat=ios) num_charges
    allocate(charges(4*num_charges-1), bestcharges(4*num_charges-1), search_range(2,4*num_charges-1), stat=ios)
    qdim = 4*num_charges-1
    if((ios /= 0).or.(num_charges < 1)) call throw_error('"'//trim(input_xyzfile)//'" has the wrong format.')
    read(30,*,iostat=ios) !skip comment line
    total_charge  = 0._rp
    !read charge coordinates and magnitude (also automatically determines total charge)
    do i = 1,qdim,4
        if(i+3 < qdim) then
            read(30,*) dummystring, charges(i:i+3)
            total_charge  = total_charge  + charges(i+3)
        else
            read(30,*) dummystring, charges(i:i+2), tmp
            total_charge = total_charge + tmp
        end if
        charges(i:i+2) = charges(i:i+2)*angstrom2bohr
        !print*, charges(i:i+2)
    end do
    read(30,*,iostat=ios) !skip blank line
!    read(30,*) dummystring, RMSE_best !read RMSE (commented as we may be using
!    fragments, in which case RMSE in reconstructed file is incorrect)
    close(30)

    !initialize best RMSE
    RMSE_best = rmse_qtot(charges(1:qdim))

    if(verbose) then
      write(*,'(A,A)') " Refining solution from file",trim(input_xyzfile)
      write(*,'(A,I6,A)') " Contains ",num_charges," charges"
      write(*,'(A,ES23.9,A)') "      Initial RMSE = ",RMSE_best*hartree2kcal," kcal/mol"
    endif

    bestcharges = charges 
    
    !initialize DE
    call DE_init(set_range            = search_range(:,1:qdim), &
                 set_popSize          = 10*qdim,                &
                 set_maxGens          = 2000*num_charges,       &
                 set_crossProb        = 1.00_rp,                &
                 set_maxChilds        = 1,                      &
                 set_forceRange       = .false.,                &
                 set_mutationStrategy = DEtargettobest1,        &
                 set_verbose          = verbose,                &
                 set_Nprint           = 100)  
    if(simplex_only .and. num_trials > 1) then
      num_trials=1
      write(*,'(A)') "Warning: number of trials should be 1 for simplex optimization"
      write(*,'(A)') "Setting num_trials=1 and continuing"
    endif
    do try = 1,num_trials
        if(simplex_only) then
          if(verbose) then
             write(*,'(A)')
             write(*,'(A)') "User requested refinement using simplex algorithm (no DE)"
             write(*,'(A)')
          endif
          call DE_simplex(rmse_qtot,feasible,charges(1:qdim))
        else
          call DE_optimize(rmse_qtot,feasible,sum_constr,charges(1:qdim),guess=charges(1:qdim))
          if(verbose) write(*,'(A,I0,A,I0)') "Starting refinement for ",num_charges," charges, trial ",try
        endif
        ! measure the quality of the fit
        RMSE_tmp = rmse_qtot(charges(1:qdim)) 
        if(verbose) write(*,'(A,ES23.9,A)') "RMSE ", RMSE_tmp, " Hartree"
        if(RMSE_tmp < RMSE_best) then
            if(verbose) write(*,'(A)') "NEW BEST!"   
            RMSE_best = RMSE_tmp
            MAE_tmp   = mae_qtot(charges(1:qdim))
            maxAE_tmp = max_ae_qtot(charges(1:qdim))
            bestcharges(1:qdim) = charges(1:qdim) 
            ! write results to file
            ppos=scan(trim(input_xyzfile),".",BACK=.true.)
            write(dummystring,'(I0,A20)') num_charges,"_charges_refined.xyz"
            call write_xyz_file(charges(1:qdim),filename=dummystring)
            !call write_cube_file(charges(1:qdim))
            ! plot with R
            call write_image_slice_data(charges(1:qdim)) !for visualization with R
!            call execute_command_line("./visualize.r",wait=.true.)
            dummystring = trim(input_xyzfile)//"_comparison.png"
            call execute_command_line("mv comparison.png "//trim(dummystring),wait=.true.)
        end if
    end do                
    !write the results for the best solution to console
    !charges(1:qdim) = bestcharges(1:qdim)
    !RMSE_tmp        = rmse_qtot(charges(1:qdim))
    !MAE_tmp         = mae_qtot(charges(1:qdim)) 
    !maxAE_tmp       = max_ae_qtot(charges(1:qdim))
    if(verbose) then
        write(*,'(A,I0,A)') "Best found solution for ",num_charges," charges:"
        write(*,'(A,ES23.9,A)') "        RMSE ", &
                  RMSE_best*hartree2kcal," kcal/mol"
        write(*,'(A,ES23.9,A)') "         MAE ", &
                  MAE_tmp*hartree2kcal," kcal/mol"
        write(*,'(A,ES23.9,A)') "     max. AE ", &
                  maxAE_tmp*hartree2kcal," kcal/mol"
        if(.not.feasible(charges(1:qdim))) write(*,'(A)') "WARNING, SOLUTION IS NOT FEASIBLE!"
        write(*,*)
        write(*,'(4A14)') "x[bohr]","y[bohr]","z[bohr]","q[e]"
        tmp = 0._rp
        do i = 1,qdim,4
            if(i+3 <= qdim) then
                write(*,'(4F14.7)') bestcharges(i), bestcharges(i+1), bestcharges(i+2), bestcharges(i+3)
                tmp = tmp + bestcharges(i+3)
            else
                write(*,'(4F14.7)') bestcharges(i), bestcharges(i+1), bestcharges(i+2), total_charge-tmp
            end if
        end do
    end if
    if(verbose) write(*,*)
    call DE_exit()
    stop
end if


! allocate memory
allocate(charges(4*num_charges_max-1), bestcharges(4*num_charges_max-1), search_range(2,4*num_charges_max-1), stat=ios)
if(ios /= 0) call throw_error('Could not allocate memory.')
if(use_symmetry) allocate(sym_charges(4*num_charges_max-1)) !max. possible size

!other fits are not meaningful for now
!check that the requested number of charges is possible with the atom fits available:
if(num_charges_min < 2)       num_charges_min = 2
if(.not.use_symmetry.and.num_charges_max > num_charges_max_multipole*natmfit) then
  write(*,*)
  write(*,*) 'WARNING:'
  write(*,'(A,I4,A,I4,A,I4,A)') "Can't create model with ",num_charges_max, &
    " charges using",num_charges_max_multipole," charges per atom for ",natmfit,&
    " atoms"
  if(num_charges_max_multipole*natmfit < num_charges_min) then
    write(*,'(A)') 'Exiting!'
    write(*,*)
    stop
  else
    write(*,'(A,I4,A)') "Reducing to maximum ",num_charges_max_multipole*natmfit,&
      " charges"
  end if
  num_charges_max = num_charges_max_multipole*natmfit
  write(*,*)
endif

! start fit for full molecule (or fragment defined by -atom flag)
do num_charges = num_charges_min,num_charges_max
    !calculate dimensionality
    qdim = 4*num_charges-1 

    ! load RMSE_best from file
    write(dummystring,'(I0)') num_charges
    dummystring = trim(dummystring)//"charges.xyz"
    if(trim(prefix) /= '') dummystring = trim(prefix)//"/"//trim(dummystring)
    open(30, file=dummystring, status='old', action='read', iostat=ios)
    if(ios == 0) then
        !skip unnecessary lines
        read(30,*) !Natom line
        read(30,*) !comment line
        !read charge coordinates and magnitude
        do i = 1,qdim,4
            if(i+3 < qdim) then
                read(30,*) dummystring, bestcharges(i:i+3)
                !print*, trim(dummystring), bestcharges(i:i+3)
            else
                read(30,*) dummystring, bestcharges(i:i+2), tmp
                !print*, trim(dummystring), bestcharges(i:i+2), tmp 
            end if
            charges(i:i+2) = charges(i:i+2)*angstrom2bohr
        end do
        read(30,*) !empty line
        read(30,*) dummystring, RMSE_best
        read(30,*) dummystring, MAE_tmp
        read(30,*) dummystring, dummystring2, maxAE_tmp
        close(30)
    endif
    RMSE_best = vbig ! reset RMSE

    ! initialize search_range
    if(.not.use_symmetry)then
      if(allocated(search_range)) deallocate(search_range)
      allocate(search_range(2,4*num_charges))
      call init_search_range()
      call DE_init(set_range            = search_range(:,1:qdim), &
                     set_popSize          = 10*qdim,                &
                     set_maxGens          = 2000*num_charges,       &
                     set_crossProb        = 1.00_rp,                &
                     set_maxChilds        = 1,                      &
                     set_forceRange       = .false.,                &
                     set_mutationStrategy = DEtargettobest1,        &
                     set_verbose          = verbose,                &
                     set_Nprint           = 100)  
    endif
    do try = 1,num_trials
        if(verbose) write(*,'(A,I0,A,I0)') "Starting fitting procedure for ",num_charges," charges, trial ",try
        if(use_symmetry)then
          ! symmetry-constrained search
          call DE_exit()
          if(allocated(search_range)) deallocate(search_range)
          allocate(search_range(2,4*num_charges))
          if(.not.init_sym_search(num_charges)) then
             cycle !skip if no fit exists with this combination of number of charges and symmetry constraints
          endif
          call init_sym_search_range(search_range,symFitAtms,num_symFitAtms,&
                           max_extend,max_charge,sqdim)
          call DE_init(set_range            = search_range(:,1:sqdim), &
                       set_popSize          = 10*sqdim,                &
                       set_maxGens          = 2000*num_charges,        &
                       set_crossProb        = 1.00_rp,                 &
                       set_maxChilds        = 1,                       &
                       set_forceRange       = .false.,                 &
                       set_mutationStrategy = DEtargettobest1,         &
                       set_verbose          = verbose,                 &
                       set_Nprint           = 100)
          call DE_optimize(sym_rmse_qtot,sym_feasible,sym_sum_constr,&
                sym_charges(1:sqdim),init_pop=sym_init_pop)
          call spawn_sym_chgs(sym_charges(1:sqdim),&
                charges,&
                num_charges,symFitAtms,num_symFitAtms,total_charge,.true.)
        else
          if(use_greedy_fit) then
              call DE_optimize(rmse_qtot,feasible,sum_constr,charges(1:qdim),init_pop=init_pop_greedy)
          else
              call DE_optimize(rmse_qtot,feasible,sum_constr,charges(1:qdim),init_pop=init_pop)
          end if
        endif
        ! measure the quality of the fit
        RMSE_tmp = rmse_qtot(charges(1:qdim)) 
        if(verbose) write(*,'(A,ES23.9,A)') "RMSE ", RMSE_tmp, " Hartree"
        if(RMSE_tmp < RMSE_best) then
            if(verbose) write(*,'(A)') "NEW BEST!"   
            RMSE_best = RMSE_tmp
            MAE_tmp   = mae_qtot(charges(1:qdim)) 
            maxAE_tmp = max_ae_qtot(charges(1:qdim))
            bestcharges(1:qdim) = charges(1:qdim) 
            ! write results to file
            call write_xyz_file(charges(1:qdim))
            !call write_cube_file(charges(1:qdim))
            ! plot with R
            call write_image_slice_data(charges(1:qdim)) !for visualization with R
!            call execute_command_line("./visualize.r",wait=.true.)
            write(dummystring,'(I0)') num_charges 
            dummystring = trim(dummystring)//"charges_comparison.png"
            if(trim(prefix) /= '') dummystring = trim(prefix)//"/"//trim(dummystring)
            call execute_command_line("mv comparison.png "//trim(dummystring),wait=.true.)
        end if
        if(verbose) write(*,*)
    end do
    !write the results for the best solution to console
    !charges(1:qdim) = bestcharges(1:qdim)
    !RMSE_tmp        = rmse_qtot(charges(1:qdim))
    !MAE_tmp         = mae_qtot(charges(1:qdim)) 
    !maxAE_tmp       = max_ae_qtot(charges(1:qdim))
    if(verbose) then
        write(*,'(A,I0,A)') "Best found solution for ",num_charges," charges:"
        write(*,'(A,ES23.9,A)') "        RMSE ", &
                  RMSE_tmp*hartree2kcal," kcal/mol"
        write(*,'(A,ES23.9,A)') "         MAE ", &
                  MAE_tmp*hartree2kcal," kcal/mol"
        write(*,'(A,ES23.9,A)') "     max. AE ", &
                  maxAE_tmp*hartree2kcal," kcal/mol"
        if(.not.feasible(charges(1:qdim))) write(*,'(A)') "WARNING, SOLUTION IS NOT FEASIBLE!"
        write(*,*)
        write(*,'(4A14)') "x[bohr]","y[bohr]","z[bohr]","q[e]"
        tmp = 0._rp
        do i = 1,qdim,4
            if(i+3 <= qdim) then
                write(*,'(4F14.7)') bestcharges(i), bestcharges(i+1), bestcharges(i+2), bestcharges(i+3)
                tmp = tmp + bestcharges(i+3)
            else
                write(*,'(4F14.7)') bestcharges(i), bestcharges(i+1), bestcharges(i+2), total_charge-tmp
            end if
        end do
    end if
    if(verbose) write(*,*)
    flush(6)
    ! clean up
    call DE_exit()
end do
stop
call dealloc()

contains
!-------------------------------------------------------------------------------
! checks whether a value is feasible (in constraints)
logical function feasible_multipole(m)
    real(rp), dimension(:) :: m
    feasible_multipole = .true.;
end function feasible_multipole
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes the amount of constraint violation
real(rp) function sum_constr_multipole(m)
    real(rp), dimension(:) :: m
    sum_constr_multipole = 0._rp
end function sum_constr_multipole
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! checks feasibility of fitting solution with symmetry constraints for atom fits
logical function sym_atm_feasible(sq)
    implicit none
    real(rp), dimension(:) :: sq ! input charges
    real(rp), dimension(num_charges*4) :: q   ! complete charges
    real(rp) :: qsum
    qsum = 0._rp

    !if one atom meets constraints, so must all its sea's, so we don't need to check
    call spawn_sym_chgs(sq,q,num_charges,symFitAtms,num_symFitAtms,total_charge,.false.) ! apply sym ops to sqin to populate q

    sym_atm_feasible = feasible(q)
end function sym_atm_feasible
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! checks feasibility of fitting solution with symmetry constraints
logical function sym_feasible(sq)
    implicit none
    real(rp), dimension(:) :: sq ! input charges
    real(rp), dimension(num_charges*4) :: q   ! complete charges
    real(rp) :: qsum
    qsum = 0._rp

    !if one atom meets constraints, so must all its sea's, so we don't need to check
    call spawn_sym_chgs(sq,q,num_charges,symFitAtms,num_symFitAtms,total_charge,.true.) ! apply sym ops to sqin to populate q

    sym_feasible = feasible(q)
end function sym_feasible
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! checks whether a value is feasible (in constraints)
logical function feasible(q)
    implicit none
    real(rp), dimension(:) :: q
    real(rp), parameter :: max_charge_magnitude = 1._rp
    real(rp) :: qtot
    real(rp) :: r, rmin ! shortest distance to atom
    integer  :: i, a, b
    qtot = 0._rp
    ! loop over the charge coordinates 
    do i = 1,size(q,dim=1),4
        !check the charge magnitude
        if(i+3 <= size(q,dim=1)) then
            if(abs(q(i+3)) >  max_charge_magnitude) then
                feasible = .false.
                !print*,'feasible: charge ',i/4,i,' too large: ',abs(q(i+3))
                return
            end if
            qtot = qtot + q(i+3)
        else
            qtot = total_charge - qtot
            if(abs(qtot) >  max_charge_magnitude) then
                feasible = .false.
                !print*,'feasible: final charge too large: ',abs(qtot)
                return
            end if
        end if
    
    
        !find atom with minimal relative distance to atom included in fit
        rmin = vbig
        if(use_symmetry) then
          do a = 1,Natom
              r = sqrt(sum((atom_pos(:,a)-q(i:i+2))**2))/(vdW_scaling*vdW_radius(atom_num(a)))
              !print*, atom_num(a), vdW_radius(atom_num(a))*bohr2angstrom
              if(r < rmin) rmin = r
          end do
        else
          do b = 1,natmfit
              a = fitatoms(b)
              r = sqrt(sum((atom_pos(:,a)-q(i:i+2))**2))/(vdW_scaling*vdW_radius(atom_num(a)))
              !print*, atom_num(a), vdW_radius(atom_num(a))*bohr2angstrom
              if(r < rmin) rmin = r
          end do
        endif
        if(rmin > 1._rp) then !this means that rmin is larger than the vdW radius
                feasible = .false.
                !print*,'feasible: charge ',(i+3)/4,i,' outside spatial bounds ',rmin
!do b=1,size(q,dim=1),4
!     print*, 'H ',q(b:b+3)
!enddo
!print*,''
                return
        end if
    end do
    
!    if(use_symmetry) then
!        if(.not.is_symmetric(q)) then
!            feasible = .false.
!            return
!        end if
!    end if
    
    feasible = .true.
!    print*,'feasible!'
    return
end function feasible
!-------------------------------------------------------------------------------

!!-------------------------------------------------------------------------------
!! checks whether the input charges are symmetric
!logical function is_symmetric(q)
!    implicit none
!    real(rp), dimension(:) :: q
!    real(rp), dimension(num_charges,3) :: pos       !stores all the position
!    real(rp), dimension(num_charges)   :: magnitude !these are used as identifier
!    real(rp) :: qsum
!    integer :: counter    
!    
!    ! loop over the charges, read out their positions and magnitude
!    qsum = 0d0
!    counter = 1
!    do i = 1,size(q,dim=1),4
!        pos(counter,1:3) = q(i:i+2) - atom_com !referenced to center of mass
!        if(i+3 <= size(q,dim=1)) then
!            qsum = qsum + q(i+3)
!            magnitude(counter) = q(i+3)
!        else
!            magnitude(counter) = total_charge - qsum
!        end if
!        counter = counter + 1
!    end do
!    
!    !loop over the unique symmetry operations and check whether it is symmetric or not
!    do i = 1,size(symmetry_ops,dim=3)
!        if(isSymmetryOperation(pos,magnitude,symmetry_ops(:,:,i))) then !one positive hit is sufficient
!            is_symmetric = .true.
!            return
!        end if
!    end do
!    is_symmetric = .false.
!    return
!end function is_symmetric
!!-------------------------------------------------------------------------------

!!-------------------------------------------------------------------------------
!! returns the smallest symmetry violation, needed for sum_constr
!real(rp) function symmetry_violation(q)
!    implicit none
!    real(rp), dimension(:) :: q
!    real(rp), dimension(num_charges,3) :: pos       !stores all the position
!    real(rp), dimension(num_charges)   :: magnitude !these are used as identifier
!    real(rp) :: qsum, tmperr, minerr
!    integer :: counter
!    ! loop over the charges, read out their positions and magnitude
!    qsum = 0._rp
!    counter = 1
!    do i = 1,size(q,dim=1),4
!        pos(counter,1:3) = q(i:i+2) - atom_com !referenced to center of mass
!        if(i+3 <= size(q,dim=1)) then
!            qsum = qsum + q(i+3)
!            magnitude(counter) = q(i+3)
!        else
!            magnitude(counter) = total_charge - qsum
!        end if
!        counter = counter + 1
!    end do
!    
!    !loop over the unique symmetry operations and record the minimal error
!    minerr = huge(0._rp)
!    do i = 1,size(symmetry_ops,dim=3)
!        tmperr = asymmetryMagnitude(pos,magnitude,symmetry_ops(:,:,i))
!        if(tmperr < minerr) minerr = tmperr
!    end do
!    symmetry_violation = minerr
!end function symmetry_violation
!!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
real(rp) function sym_sum_atm_constr(sq)
    implicit none
    real(rp), dimension(:) :: sq ! input charges
    real(rp), dimension(num_charges*4) :: q   ! complete charges
    real(rp) :: qsum
    qsum = 0._rp

    call spawn_sym_chgs(sq,q,num_charges,symFitAtms,num_symFitAtms,total_charge,.false.) ! apply sym ops to sqin to populate q

    sym_sum_atm_constr = sum_constr(q)
end function sym_sum_atm_constr
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
real(rp) function sym_sum_constr(sq)
    implicit none
    real(rp), dimension(:) :: sq ! input charges
    real(rp), dimension(num_charges*4) :: q   ! complete charges
    real(rp) :: qsum
    qsum = 0._rp

    call spawn_sym_chgs(sq,q,num_charges,symFitAtms,num_symFitAtms,total_charge,.true.) ! apply sym ops to sqin to populate q

    sym_sum_constr = sum_constr(q)
end function sym_sum_constr
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! returns the sum of constraint violations
real(rp) function sum_constr(q)
    implicit none
    real(rp), dimension(:) :: q
    real(rp) :: r, rmin ! distance to atom
    integer  :: i, a, b
    sum_constr = 0._rp
    ! loop over the charge coordinates 
    do i = 1,size(q,dim=1),4
        !find atom with minimal relative distance
        rmin = vbig
        do b = 1,natmfit
            a = fitatoms(b)
            r = sqrt(sum((atom_pos(:,a)-q(i:i+2))**2))/(vdW_scaling*vdW_radius(atom_num(a)))
            if(r < rmin) rmin = r
        end do
        if(rmin > 1._rp) then !this means that rmin is larger than the vdW radius
            sum_constr = sum_constr + (rmin-1._rp)
        end if
    end do
    
!    !check symmetry violation
!    if(use_symmetry) then
!        sum_constr = sum_constr + symmetry_violation(q) 
!    end if
end function sum_constr
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! computes root mean squared error of the current fit to the true esp, using constraint charges
! (this means that the total charge must add up to a specific value)
real(rp) function rmse_qtot(qin)
    implicit none
    real(rp), dimension(:) :: qin ! input charges
    real(rp), dimension(size(qin,dim=1)+1) :: q   ! complete charges
    real(rp) :: qsum
    integer :: i
    qsum = 0._rp
    do i = 1,size(qin,dim=1)-3,4
        qsum = qsum + qin(i+3)
    end do
    q(1:size(qin,dim=1)) = qin(:)
    q(size(q,dim=1)) = total_charge - qsum
    rmse_qtot = rmse(q)
end function rmse_qtot
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! computes root mean squared error of the current fit to the true esp, using constraint
! charges (this means that the total charge must add up to a specific value) and 
! applyng symmetry operations to spawn additional charges until symmetry constraints
! are satisfied (i.e. charge arrangement possesses same symmetry as parent molecule)
real(rp) function sym_rmse_qtot(sqin)
    implicit none
    real(rp), dimension(:) :: sqin ! input charges
    real(rp), dimension(num_charges*4) :: q   ! complete charges
    real(rp) :: qsum
    qsum = 0._rp

    call spawn_sym_chgs(sqin,q,num_charges,symFitAtms,num_symFitAtms,total_charge,.true.) ! apply sym ops to sqin to populate q

    sym_rmse_qtot = rmse(q)
end function sym_rmse_qtot
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! computes root mean squared error of a given atom fit to the true esp, using constraint
! charges (this means that the total charge must add up to a specific value) and 
! applyng symmetry operations to spawn additional charges until symmetry constraints
! are satisfied (i.e. charge arrangement possesses same symmetry as parent molecule)
real(rp) function sym_atm_rmse_qtot(sqin)
    implicit none
    real(rp), dimension(:) :: sqin ! input charges
    real(rp), dimension(num_charges*4) :: q   ! complete charges
    real(rp) :: qsum
    integer,dimension(1) :: atms
    qsum = 0._rp

    atms(1)=fitAtm
    call spawn_sym_chgs(sqin,q,num_charges,atms,1,total_charge,.false.) ! apply sym ops to sqin to populate q

    sym_atm_rmse_qtot = rmse(q)
end function sym_atm_rmse_qtot
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! computes root mean squared error of the current fit to the true esp
real(rp) function rmse(q)
    implicit none
    real(rp), dimension(:) :: q ! input charges
    integer :: idx
    rmse = 0._rp
    do idx = 1,Ngrid
        rmse = rmse + (coulomb_potential(gridval(:,idx),q) - esp_grid(idx))**2
    end do
    rmse = sqrt(rmse/Ngridr)
end function rmse
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes root mean squared error of the current fit to the true esp
real(rp) function rmse_multipole(m)
    implicit none
    real(rp), dimension(:) :: m ! input multipoles
    integer :: idx
    rmse_multipole = 0._rp
    do idx = 1,Ngrid
        rmse_multipole = rmse_multipole + (coulomb_potential_multipole(gridval(:,idx),m) - esp_grid(idx))**2
    end do
    rmse_multipole = sqrt(rmse_multipole/Ngridr)
end function rmse_multipole
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! computes the RMSE 
subroutine do_analysis(vdw_grid_min_cutoff,vdw_grid_max_cutoff)
    implicit none
    real(rp), parameter :: hartree2kcalmol = 627.509469_rp
    integer :: idx, Nclose, Nmedium, Nfar, npts
    real(rp) :: rmse_tot, rmse_close, rmse_medium, rmse_far, maxerror, currerror
    real(rp) :: vdw_grid_min_cutoff,  vdw_grid_max_cutoff
    real(rp) :: tmp_min_cutoff, tmp_max_cutoff

     
    write(*,'(A30,A23,A23,A10)')  "Measure","RMSE[kcal/mol/e]", "MaxError[kcal/mol/e]","NPoints"
    !compute total RMSE
    rmse_tot = 0._rp
    maxerror = 0._rp
    npts = 0
    tmp_min_cutoff = vdw_grid_min_cutoff
    vdw_grid_min_cutoff = 1.2_rp
    do idx = 1,Ngrid
      if(in_interaction_belt(gridval(:,idx),vdw_grid_min_cutoff, &
          vdw_grid_max_cutoff)) then
        npts = npts + 1
        currerror = (esp_grid2(idx) - esp_grid(idx))**2 
        rmse_tot = rmse_tot + currerror  
        if(currerror > maxerror) maxerror = currerror
      end if
    end do    
    write(*,'(A30,2ES23.9,I10)') "Total", sqrt(rmse_tot/npts)*hartree2kcalmol, &
        sqrt(maxerror)*hartree2kcalmol,npts
    vdw_grid_min_cutoff = tmp_min_cutoff
    
    !compute close RMSE (r < 1.66)
    Nclose = 0
    rmse_close = 0._rp
    maxerror = 0._rp
    tmp_min_cutoff = vdw_grid_min_cutoff
    tmp_max_cutoff = vdw_grid_max_cutoff
    vdw_grid_min_cutoff = 1.2_rp
    vdw_grid_max_cutoff = 1.66_rp
    do idx = 1,Ngrid
        if(in_interaction_belt(gridval(:,idx),vdw_grid_min_cutoff, &
           vdw_grid_max_cutoff)) then
            Nclose = Nclose + 1
            currerror = (esp_grid2(idx) - esp_grid(idx))**2 
            rmse_close = rmse_close + currerror 
            if(currerror > maxerror) maxerror = currerror
        end if
    end do
    if(Nclose > 0)  write(*,'(A7,F6.2,A7,F6.2,A4,2ES23.9,I10)') "Close (", &
       vdw_grid_min_cutoff," < r < ",vdw_grid_max_cutoff,")  ", &
         sqrt(rmse_close/real(Nclose,rp))*hartree2kcalmol, &
         sqrt(maxerror)*hartree2kcalmol, Nclose
    vdw_grid_min_cutoff = tmp_min_cutoff
    vdw_grid_max_cutoff = tmp_max_cutoff
    
    !compute medium RMSE (1.66 < r < 2.20)
    Nmedium = 0
    rmse_medium = 0._rp
    maxerror = 0._rp
    tmp_min_cutoff = vdw_grid_min_cutoff
    tmp_max_cutoff = vdw_grid_max_cutoff
    vdw_grid_min_cutoff = 1.66_rp
    vdw_grid_max_cutoff = 2.20_rp
    do idx = 1,Ngrid
        if(in_interaction_belt(gridval(:,idx),vdw_grid_min_cutoff, &
           vdw_grid_max_cutoff)) then
            Nmedium = Nmedium + 1
            currerror = (esp_grid2(idx) - esp_grid(idx))**2 
            rmse_medium = rmse_medium + currerror 
            if(currerror > maxerror) maxerror = currerror  
        end if
    end do
    if(Nmedium > 0)  write(*,'(A11,F5.2,A7,F5.2,A2,2ES23.9,I10)') "Mid-range (",&
         vdw_grid_min_cutoff," < r < ", &
         vdw_grid_max_cutoff,") ", &
         sqrt(rmse_medium/real(Nmedium,rp))*hartree2kcalmol,&
         sqrt(maxerror)*hartree2kcalmol, Nmedium
    vdw_grid_min_cutoff = tmp_min_cutoff
    vdw_grid_max_cutoff = tmp_max_cutoff
    
    !compute far RMSE (2.20 < r < infinity)
    Nfar = 0
    rmse_far = 0._rp
    maxerror = 0._rp
    tmp_min_cutoff = vdw_grid_min_cutoff
    vdw_grid_min_cutoff = 2.20_rp
    do idx = 1,Ngrid
        if(in_interaction_belt(gridval(:,idx),vdw_grid_min_cutoff, &
           vdw_grid_max_cutoff)) then
            Nfar = Nfar + 1
            currerror = (esp_grid2(idx) - esp_grid(idx))**2 
            rmse_far = rmse_far + currerror  
            if(currerror > maxerror) maxerror = currerror
        end if
    end do
    if(Nfar > 0)  write(*,'(A11,F6.2,A13,2ES23.9,I10)') "Far-range (", &
           vdw_grid_min_cutoff," < r)", &
           sqrt(rmse_far/real(Nfar,rp))*hartree2kcalmol, &
           sqrt(maxerror)*hartree2kcalmol, Nfar
    vdw_grid_min_cutoff = tmp_min_cutoff
    

end subroutine do_analysis
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! computes root mean squared error of the current fit to the true esp, using constraint charges
! (this means that the total charge must add up to a specific value)
real(rp) function mae_qtot(qin)
    implicit none
    real(rp), dimension(:) :: qin ! input charges
    real(rp), dimension(size(qin,dim=1)+1) :: q   ! complete charges
    real(rp) :: qsum
    integer :: i
    qsum = 0._rp
    do i = 1,size(qin,dim=1)-3,4
        qsum = qsum + qin(i+3)
    end do
    q(1:size(qin,dim=1)) = qin(:)
    q(size(q,dim=1)) = total_charge - qsum
    mae_qtot = mae(q)
end function mae_qtot
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes mean absolute error of the current fit to the true esp
real(rp) function mae(q)
    implicit none
    real(rp), dimension(:) :: q ! input charges
    integer :: idx
    mae = 0._rp
    do idx = 1,Ngrid
        mae = mae + abs(coulomb_potential(gridval(:,idx),q) - esp_grid(idx))
    end do
    mae = mae/Ngridr
end function mae
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes the maximum squared error of the current fit to the true esp, using constraint charges
! (this means that the total charge must add up to a specific value)
real(rp) function max_ae_qtot(qin)
    implicit none
    real(rp), dimension(:) :: qin ! input charges
    real(rp), dimension(size(qin,dim=1)+1) :: q   ! complete charges
    real(rp) :: qsum
    integer :: i
    qsum = 0._rp
    do i = 1,size(qin,dim=1)-3,4
        qsum = qsum + qin(i+3)
    end do
    q(1:size(qin,dim=1)) = qin(:)
    q(size(q,dim=1)) = total_charge - qsum
    max_ae_qtot = max_ae(q)
end function max_ae_qtot
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes the maximum squared error of the current fit to the true esp
real(rp) function max_ae(q)
    implicit none
    real(rp), dimension(:) :: q ! input charges
    real(rp) :: ae
    integer :: idx
    max_ae = -vbig
    do idx = 1,Ngrid
        ae = abs(coulomb_potential(gridval(:,idx),q) - esp_grid(idx))
        if(ae > max_ae) max_ae = ae
    end do
end function max_ae
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes Coulomb potential for the given charges q at position x, using constraint charges
! (this means that the total charge must add up to a specific value)
real(rp) function coulomb_potential_qtot(x,qin)
    implicit none
    real(rp), dimension(3), intent(in)     :: x ! position
    real(rp), dimension(:) :: qin ! input charges
    real(rp), dimension(size(qin,dim=1)+1) :: q   ! complete charges
    real(rp) :: qsum
    integer :: i
    qsum = 0._rp
    !do i = 1,size(qin,dim=1),4
    do i = 1,size(qin,dim=1)-3,4
        qsum = qsum + qin(i+3)
    end do
    q(1:size(qin,dim=1)) = qin(:)
    q(size(q,dim=1)) = total_charge - qsum
    coulomb_potential_qtot = coulomb_potential(x,q)
end function coulomb_potential_qtot
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes Coulomb potential for the given charges q at position x
real(rp) function coulomb_potential(x,q)
    implicit none
    real(rp), dimension(3), intent(in) :: x ! position
    real(rp), dimension(:) :: q ! charges
    real(rp) :: r
    integer :: i
    coulomb_potential = 0._rp
    do i=1,size(q,dim=1),4
        r = sqrt(sum((q(i:i+2)-x)**2))  ! compute distance
        if(r < 1.e-9_rp) r = 1.e-9_rp ! prevent division by 0
        coulomb_potential = coulomb_potential + q(i+3)/r
    end do
end function coulomb_potential
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes Coulomb potential for the given mutipole expansion at position x, 
real(rp) function coulomb_potential_multipole(x,mp)
    implicit none
    real(rp), dimension(3), intent(in)    :: x ! position
    real(rp), dimension(:) :: mp ! multipoles
    real(rp) :: rmag(Natom), r(3,Natom), qsum
    integer :: idx, l, m, a, load
    
    !precalculate distances and vectors for all atoms, as well as the total charge
    qsum = 0._rp
    do a = 1,Natom
        r(:,a) = x - atom_pos(:,a)
        rmag(a) = sqrt(sum(r(:,a)**2))
        r(:,a) = r(:,a)/rmag(a)
        if((lcur == 0).and.(a < Natom)) qsum = qsum + mp(a) ! add charge to total charge
    end do
    
    if(lcur == 0 .and. .not. generate_mode) then
        mp(Natom) = total_charge - qsum ! constrains the total charge to 0
    else !load the previously found solution
        load = lcur**2*Natom
        mp(1:load) = multipole_best(1:load)
    end if
    

    idx = 0
    coulomb_potential_multipole = 0._rp
    
    do l = 0,lcur    
        do a = 1,Natom
            do m = -l,l
               idx = idx + 1
                   coulomb_potential_multipole = coulomb_potential_multipole + mp(idx)*mESP(l,m,rmag(a),r(:,a))
            end do
        end do
    end do    
         
end function coulomb_potential_multipole
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! computes Coulomb potential for a single given atomic mutipole expansion at position x
real(rp) function coulomb_potential_single_multipole(x,mp,a)
    implicit none
    integer, intent(in) :: a
    real(rp), dimension(3), intent(in)    :: x ! position    
    real(rp), dimension(:) :: mp ! multipoles
    
    real(rp)  :: rmag, r(3)
    integer :: idx, l, m, i
    
    !precalculate distance
    r = x - atom_pos(:,a)
    rmag = sqrt(sum(r**2))
    r = r/rmag

    
    idx = 0
    coulomb_potential_single_multipole = 0._rp
    do l = 0,lcur    
        do i = 1,Natom
            do m = -l,l
                idx = idx + 1
                if(i /= a) cycle
                coulomb_potential_single_multipole = coulomb_potential_single_multipole + mp(idx)*mESP(l,m,rmag,r)
            end do
        end do
    end do    
         
end function coulomb_potential_single_multipole
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! checks whether a point is inside the interaction belt as defined by the vdW radii
! Atomic Multipoles: Electrostatic Potential Fit, Local Reference Axis Systems,
! and Conformational Dependence
logical function in_interaction_belt(q,mincut,maxcut)
    implicit none
    real(rp), dimension(:) :: q
    real(rp) :: r, rmin ! shortest distance to atom
    real(rp) :: mincut, maxcut ! minimum and maximum grid cut-offs passed
    integer  :: i, a, b
    ! loop over the coordinates of the point q
    do i = 1,size(q,dim=1),4
        !first check we're not inside the molecule
        !find atom with minimal relative distance
        rmin = vbig
        do a = 1,Natom
            r = sqrt(sum((atom_pos(:,a)-q(i:i+2))**2))/(vdW_radius(atom_num(a)))
            if(r .lt. rmin) then
              rmin = r
            endif
        end do
        ! this means the radius is not in the defined interaction cutoff
        if(.not.(rmin >= mincut) ) then
                in_interaction_belt = .false.
                return
        end if
        !now check we're not too far from any atoms in the fragment (or molecule)
        rmin = vbig
        do b = 1,natmfit
            a = fitatoms(b)
            r = sqrt(sum((atom_pos(:,a)-q(i:i+2))**2))/(vdW_radius(atom_num(a)))
            if(r .lt. rmin) then
              rmin = r
            endif
        end do
        ! this means the radius is not in the defined interaction cutoff
        if(.not.(rmin <= maxcut)) then
                in_interaction_belt = .false.
                return
        end if
    end do
    in_interaction_belt = .true.
    return
end function in_interaction_belt
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine init_search_range()
    implicit none
    real(rp) :: x_min, x_max, y_min, y_max, z_min, z_max, q_min, q_max
    real(rp), dimension(3,natmfit) :: fitatom_pos 
    integer :: i,a,b

    do b = 1,natmfit
      a = fitatoms(b)
      fitatom_pos(:,b) = atom_pos(:,a)
    end do
    
    x_min = minval(fitatom_pos(1,:)) - max_extend
    x_max = maxval(fitatom_pos(1,:)) + max_extend
    y_min = minval(fitatom_pos(2,:)) - max_extend
    y_max = maxval(fitatom_pos(2,:)) + max_extend
    z_min = minval(fitatom_pos(3,:)) - max_extend
    z_max = maxval(fitatom_pos(3,:)) + max_extend
    q_min = -max_charge
    q_max =  max_charge
    
    do i = 1,qdim,4
        search_range(1,i  ) = x_min
        search_range(2,i  ) = x_max
        search_range(1,i+1) = y_min
        search_range(2,i+1) = y_max
        search_range(1,i+2) = z_min
        search_range(2,i+2) = z_max
        if(i+3 < size(charges,dim=1)) then ! we use constraint charges
            search_range(1,i+3) = q_min
            search_range(2,i+3) = q_max
        end if
    end do
    
end subroutine init_search_range
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! initializes the population to only feasible solutions
subroutine init_pop(pop)
    implicit none
    real(rp), dimension(:,:), intent(out) :: pop 
    real(rp), dimension(3) :: ranvec
    real(rp) :: ran
    integer  :: popSize, a, b, d, p
    
    popSize = size(pop,dim=2)
    
    ! loop over population
    do p = 1,popSize
        do d = 1,qdim,4
            ! determine atom number
            if(d/4 < (num_charges/natmfit)*natmfit) then 
                ! choose atom through equal sampling
                b = mod(d/4,natmfit) + 1   
                a = fitatoms(b)
            else
                ! determine a random atom
                call random_number(ran)
                b = ceiling(ran*natmfit)
                a = fitatoms(b)
            end if
            !write(*,'(A,I0,A,I0,A,L)') "pop ", p, " atom# ", a ," random? ", .not.(d/4 < (num_charges/Natom)*Natom)
            
            ! draw a random vector and normalize (random direction)
            call random_number(ranvec)
            ranvec = ranvec - 0.5_rp
            ranvec = ranvec/sqrt(sum(ranvec**2))
            
            ! draw a random scaling factor and scale vector
            call random_number(ran)
            ranvec = ran*vdW_scaling*vdW_radius(atom_num(a))*ranvec
            
            ! set charge position
            pop(d:d+2,p) = atom_pos(:,a) + ranvec
            
            ! set charge magnitude
            if(d+3 < qdim) then ! in case we use constraint charges
                call random_number(ran)
                pop(d+3,p) = -max_charge + 2*ran*max_charge
            end if
        end do
    end do     
end subroutine init_pop
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! initializes the population to only feasible solutions generated using the greedy fit
subroutine init_pop_greedy(pop)
    implicit none
    real(rp), dimension(:,:), intent(out) :: pop 
    real(rp), dimension(size(pop,dim=1))  :: prototype
    real(rp) :: ran, sumP, charge_correction
    real(rp), dimension(Natom) :: improvementProbability
    integer, dimension(Natom)  :: chargesPerAtom
    integer  :: popSize, a, p, dimStart, dimEnd
    
    popSize = size(pop,dim=2)    
    
    ! loop over population        
    do p = 1,popSize
        !decide how to distribute the charges
        chargesPerAtom = 0 !initialize to 0
        improvementProbability = 1._rp/real(natmfit,rp) !initialize
        
        do while(sum(chargesPerAtom) < num_charges .and. &
                 sum(chargesPerAtom) < num_charges_max_multipole*natmfit) !we randomly decide where to put charges
            
            !determine how likely it is that placing a charge will improve the solution    
            do b = 1,natmfit
                a = fitatoms(b)
                if(chargesPerAtom(a) < num_charges_max_multipole) then
                    !compare RMSE with a chgs to RMSE with a+1 for atomic solution
                    improvementProbability(a) = multipole_solutions_rmse(chargesPerAtom(a),a) &
                                              - multipole_solutions_rmse(chargesPerAtom(a)+1,a)
                    if(improvementProbability(a) < 0._rp) improvementProbability(a) = 0._rp
                else
                    improvementProbability(a) = 0._rp
                end if
            end do   
            !normalize
            if(sum(improvementProbability) > epsilon(0._rp)) then
                improvementProbability = improvementProbability/sum(improvementProbability)
            else
                improvementProbability = 1._rp/real(natmfit,rp)
            end if
            call random_number(ran) !draw a random number
            sumP = 0._rp
            do b = 1,natmfit
                a = fitatoms(b)
                sumP = sumP + improvementProbability(a)
                if(sumP > ran) then
                    chargesPerAtom(a) = chargesPerAtom(a) + 1
                    exit 
                end if
            end do
        end do
        
        !in case we have 0 charges on at least 1 atom, the total charge will be messed up.
        !So, all charges need to be scaled accordingly. Here we determine by what they are scaled.
        charge_correction = 0._rp
        do b = 1,natmfit
            a = fitatoms(b)
            if(chargesPerAtom(a) == 0) then
                charge_correction = charge_correction + multipole(a)
            end if
        end do
        charge_correction = charge_correction/real(num_charges,rp)
        
        !write(*,*) "charges per atom:", chargesPerAtom(:)
        !write(*,*) "probabilites:", improvementProbability
        !write(*,*)
        
        !build the prototype vector
        dimStart = 1
        do b = 1,natmfit
            a = fitatoms(b)
            if(chargesPerAtom(a) == 0) cycle
            !determine exactly what to load
            dimEnd = chargesPerAtom(a)*4-1
            if(dimStart+dimEnd > qdim) dimEnd = dimEnd-1
            prototype(dimStart:dimStart+dimEnd) = multipole_solutions(1:1+dimEnd,chargesPerAtom(a),a)
            dimStart = dimStart+dimEnd+1
            if(dimStart > qdim) exit
            !correct charge, if necessary
            prototype(dimStart-1) = prototype(dimStart-1) + charge_correction
        end do
        
        !calculate a factor with which to multiply the prototype
        !also, this factor gets further scaled down for the initial vectors in the population
        call random_number(pop(:,p))
        ran = 2*real(p,rp)/real(popSize,rp)
        if(ran > 1._rp) ran = 1._rp !clamp to 1
        pop(:,p) = 1._rp + (2*pop(:,p)-1) * ran

        !set vector
        pop(:,p) = pop(:,p) * prototype
    end do     
end subroutine init_pop_greedy
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! initializes the population to only feasible solutions
subroutine init_pop_multipole(pop)
    implicit none
    real(rp), dimension(:,:), intent(out) :: pop 
    real(rp), dimension(3) :: ranvec
    real(rp) :: ran
    integer  :: popSize, d, p
    
    popSize = size(pop,dim=2)
    
    ! loop over population
    do p = 1,popSize
        do d = 1,qdim,4            
            ! draw a random vector and normalize (random direction)
            call random_number(ranvec)
            ranvec = ranvec - 0.5_rp
            ranvec = ranvec/sqrt(sum(ranvec**2))
            
            ! draw a random scaling factor and scale vector
            call random_number(ran)
            ranvec = ran*vdW_scaling*vdW_radius(atom_num(a))*ranvec
            
            ! set charge position
            pop(d:d+2,p) = atom_pos(:,a) + ranvec
            
            ! set charge magnitude
            if(d+3 < qdim) then ! in case we use constraint charges
                call random_number(ran)
                pop(d+3,p) = -abs(multipole(a)) + 2*ran*abs(multipole(a))
            end if
        end do
    end do     
end subroutine init_pop_multipole
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! initializes the population to only feasible symmetry solutions
subroutine sym_init_pop_multipole(pop)
    implicit none
    real(rp), dimension(:,:), intent(out) :: pop
    integer  :: popSize, p

    popSize = size(pop,dim=2)

    ! loop over population
    do p = 1,popSize
      call sym_init_pars(pop(:,p),symFitAtms,num_symFitAtms,multipole,vdW_scaling,&
           vdW_radius)
    end do
end subroutine sym_init_pop_multipole
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! initializes the population to only feasible symmetry solutions
subroutine sym_init_pop(pop)
    implicit none
    real(rp), dimension(:,:), intent(out) :: pop
    integer  :: popSize, p

    popSize = size(pop,dim=2)

    ! loop over population
    do p = 1,popSize
      call sym_init_pars(pop(:,p),symFitAtms,num_symFitAtms,multipole,vdW_scaling,&
           vdW_radius)
    end do
end subroutine sym_init_pop
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! initializes the search parameters to only feasible symmetry-constrianed solutions
function init_sym_search(num_charges)
    implicit none
    integer :: num_charges !number of charges to fit for molecule / fragment
    logical :: init_sym_search
    integer :: i,nchg,nchg_tot
    integer, dimension(Natom) :: num_seas
    integer,dimension(num_symFitAtms) :: combo
    real(rp) :: rmse_tot

    init_sym_search=.false.
    nchg=0
    if(allocated(best_symatm_combo)) deallocate(best_symatm_combo)
    allocate(best_symatm_combo(num_symFitAtms))
    
    !find how many sea's each symmetry-unique atom has
    do i=1,Natom ! loop over sea sets
      if(count(atom_sea(i,:) /= 0) == 0) exit
      num_seas(atom_sea(i,1))=count(atom_sea(i,:) /= 0)
    enddo
    !find permutation of charges for symmetry-unique atoms with fewer than num_charges
    !charges in total and lowest estimated RMSE
    lrmse_best=vbig
    best_symatm_combo=combo
    rmse_tot=0.d0
    nchg_tot=0
    combo(1:num_symFitAtms)=1
    call sym_permutate_rmse(1,rmse_tot,nchg_tot,num_charges,num_seas,combo)
    if(lrmse_best.lt.vbig) then
      init_sym_search=.true.
      write(*,'(/,A,I0,A)') 'Best initial guess with ',num_charges,' charges found for:'
      do i=1,num_symFitAtms
        write(*,'(2(A,I0),A)') 'Atom ',symFitAtms(i),' with ',best_symatm_combo(i),&
          ' charges'
      enddo
    else
      write(*,'(/,A,I0,A,/)') 'No suitable solution found with ',num_charges,' charges'
      stop
    endif
    !now we chose a set of fitting operations, initialize arrays in symmetry routines
    call sym_init_fit_ops(symFitAtms,num_symFitAtms,best_symatm_combo,num_charges,&
             sym_solution_ops,num_sym_solution_ops)
end function init_sym_search

!-------------------------------------------------------------------------------
!find permutation of atomic solutions with up to maxchg charges that provides
!lowest RMSE. Solutions with under maxchg charges are allowed as many charge
!numbers aren't possible within symmetry constraints
recursive subroutine sym_permutate_rmse(N,rmse_tot,nchg_tot,nChgFit,num_seas,combo)
implicit none
integer  :: N,i,nchg_tot,lnchg_tot,nChgFit,atm
integer, dimension(:) :: num_seas
integer,dimension(num_symFitAtms) :: combo
real(8)  :: rmse_tot,lrmse_tot

! for each atom in symFitAtms
atm=symFitAtms(N)
lrmse_tot=rmse_tot
lnchg_tot=nchg_tot

! start with case where we have zero charges for this atom
if(lnchg_tot.le.nChgFit)then
  combo(N)=0
  lrmse_tot=rmse_tot+multipole_solutions_rmse(0,atm)*num_seas(atm)
  if(N<num_symFitAtms) then
    call sym_permutate_rmse(N+1,lrmse_tot,lnchg_tot,nChgFit,num_seas,combo)
  else !reached last atom in fit
    if(lrmse_tot.lt.lrmse_best.and.lnchg_tot.eq.nChgFit) then
      lrmse_best=lrmse_tot
      best_symatm_combo(:)=combo(:)
    endif
  endif
endif
! loop over best fitting results for each number of charges of this atom
do i=num_charges_min_multipole,num_charges_max_multipole
  ! add total number of charges spawned by this atom across all sea's
  lnchg_tot=nchg_tot+i*num_seas(atm)
  ! track which solutions we currently test
  combo(N)=i
  ! if total number of chgs is less than molecular target:
  if(lnchg_tot.le.nChgFit)then
    !check there was actually a solution for this many charges:
    if(multipole_solutions_rmse(i,atm).ne.vbig)then
      lrmse_tot=rmse_tot+multipole_solutions_rmse(i,atm)*num_seas(atm)
      if(N<num_symFitAtms) then
        call sym_permutate_rmse(N+1,lrmse_tot,lnchg_tot,nChgFit,num_seas,combo)
      else !reached last atom in fit
        if(lrmse_tot.lt.lrmse_best.and.lnchg_tot.eq.nChgFit) then
          lrmse_best=lrmse_tot
          best_symatm_combo(:)=combo(:)
        endif
      endif
    endif
  endif
end do

end subroutine sym_permutate_rmse

!-------------------------------------------------------------------------------
subroutine read_multipole_file(inpfile)
implicit none
character(len=*), intent(in) :: inpfile
integer :: ios,idx,l,a,m
character(len=1024) :: dummy1,dummy2

open(30, file=trim(inpfile), status="old", action="read", iostat = ios)
if(ios /= 0) call throw_error('Could not open "'//trim(inpfile)//'" for reading')
read(30,*,iostat=ios) !skip total charge
read(30,*,iostat=ios) !skip comment
read(30,*) lcur !read l
idx = 0
do l = 0,lcur
    do a = 1,Natom
        do m = -l,l
            idx = idx + 1
            read(30,*,iostat=ios) dummy1, dummy2, multipole(idx)
        end do
    end do
end do
multipole_best = multipole
close(30)



end subroutine read_multipole_file

!-------------------------------------------------------------------------------
!recursive subroutine permutate(E, P)
!implicit none
!integer, intent(in)  :: E(:)       ! array of objects
!integer, intent(out) :: P(:,:)     ! permutations of E
!integer  :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1)
!N = size(E); Nfac = size(P,1);
!do i=1,N                           ! cases with E(i) in front
!  if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S)
!  forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/)
!end do
!end subroutine permutate


!-------------------------------------------------------------------------------
subroutine read_xyz_file()
implicit none
integer :: i,ios
character(len=1024) :: dummy
real(rp) :: tmp

open(30, file=trim(input_xyzfile), status="old", action="read", iostat = ios)
if(ios /= 0) call throw_error('Could not open "'//trim(input_xyzfile)//'" for reading')
read(30,*,iostat=ios) num_charges
qdim = 4*num_charges-1
if((ios /= 0).or.(num_charges < 1)) call throw_error('"'//trim(input_xyzfile)//'" hast the wrong format.')
if(.not.allocated(charges)) allocate(charges(4*num_charges-1))
read(30,*,iostat=ios) !skip comment line
total_charge  = 0._rp
!read charge coordinates and magnitude
do i = 1,qdim,4
    if(i+3 < qdim) then
        read(30,*) dummy, charges(i:i+3)
        total_charge  = total_charge  + charges(i+3)
    else
        read(30,*) dummy, charges(i:i+2), tmp
        total_charge = total_charge + tmp
    end if
    charges(i:i+2) = charges(i:i+2)*angstrom2bohr
    !print*, charges(i:i+2)
end do
end subroutine read_xyz_file
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! writes a xyz file containing the results
subroutine write_xyz_file(charges,a,filename)
implicit none 
integer :: ios,i
integer, optional :: a
real(rp), dimension(qdim), intent(in) :: charges
character(len=*), intent(in), optional :: filename
character(len=1024) :: outfile, dummy

write(outfile,'(I0)') num_charges
if(present(a)) then
    write(dummy,'(I0)') a
    outfile = "multipole"//trim(dummy)//"_"//trim(outfile)//"charges.xyz"
else
    outfile = trim(outfile)//"charges.xyz"
end if
if(trim(prefix) /= '') outfile = trim(prefix)//"/"//trim(outfile)

if(present(filename)) outfile = filename


open(30, file=trim(outfile), status="replace", action="write", iostat = ios)
if(ios /= 0) call throw_error('Could not open "'//trim(outfile)//'" for writing')
write(30,'(I0)') num_charges
write(30,'(A,4A26)')  "s","x[A]","y[A]","z[A]","q[e]"
tmp = 0._rp
do i = 1,qdim,4
    if(i+3 <= qdim) then
        if(charges(i+3) > 0._rp) then
            write(30,'(A,1X,4(F25.16,1X))') "N",charges(i:i+2)*bohr2angstrom,charges(i+3)
        else
            write(30,'(A,1X,4(F25.16,1X))') "O",charges(i:i+2)*bohr2angstrom,charges(i+3)
        end if
        tmp = tmp + charges(i+3)
    else
        tmp = total_charge-tmp
        if(tmp > 0._rp) then
            write(30,'(A,1X,4(F25.16,1X))') "N",charges(i:i+2)*bohr2angstrom,tmp
        else
            write(30,'(A,1X,4(F25.16,1X))') "O",charges(i:i+2)*bohr2angstrom,tmp 
        end if
    end if
end do
write(30,*)
if(.not.feasible(charges(1:qdim))) write(30,'(A)') "WARNING, SOLUTION IS NOT FEASIBLE!"
write(30,'(A,ES23.9,A)') "        RMSE ", &
                  RMSE_tmp*hartree2kcal," kcal/mol"
write(30,'(A,ES23.9,A)') "         MAE ", &
                  MAE_tmp*hartree2kcal," kcal/mol"
write(30,'(A,ES23.9,A)') "     max. AE ", &
                  maxAE_tmp*hartree2kcal," kcal/mol"
write(30,*)
write(30,'(A)') "Coordinates in bohr"
write(30,'(A,4A26)')  "s","x[bohr]","y[bohr]","z[bohr]","q[e]"
tmp = 0._rp
do i = 1,qdim,4
    if(i+3 <= qdim) then
        if(charges(i+3) > 0._rp) then
            write(30,'(A,1X,4(F25.16,1X))') "+",charges(i:i+2),charges(i+3)
        else
            write(30,'(A,1X,4(F25.16,1X))') "-",charges(i:i+2),charges(i+3)
        end if
        tmp = tmp + charges(i+3)
    else
        tmp = total_charge-tmp
        if(tmp > 0._rp) then
            write(30,'(A,1X,4(F25.16,1X))') "+",charges(i:i+2),tmp
        else
            write(30,'(A,1X,4(F25.16,1X))') "-",charges(i:i+2),tmp 
        end if
    end if
end do
close(30)

end subroutine write_xyz_file
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! write multipole solution to a file
subroutine write_multipole_file(mp)
    implicit none 
    integer :: ios,l,m,a,idx
    real(rp), dimension(:) :: mp !multipoles
    character(len=1024) :: outfile
    
    !decide filename
    select case (lcur)
        case(0) 
            write(outfile,'(A)') "monopole"
        case(1)
            write(outfile,'(A)') "dipole"
        case(2) 
            write(outfile,'(A)') "quadrupole"
        case(3)
            write(outfile,'(A)') "octopole"
        case(4) 
            write(outfile,'(A)') "hexadecapole"
        case(5)
            write(outfile,'(A)') "ditriantapole"
    end select    
    outfile = trim(outfile)//"_expansion.txt"
    if(trim(prefix) /= '')  outfile = trim(prefix)//"/"//trim(outfile)
    open(30, file=trim(outfile), status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "'//trim(outfile)//'" for writing')
    
    write(30,'(ES23.9)') total_charge
    write(30,'(A,ES23.9,A)') "# RMSE: ", rmse_multipole(mp), "      Qa(l,m) [a = atom index]"
    write(30,'(I0)') lcur
    idx = 0
    do l = 0,lcur
        do a = 1,Natom
            do m = -l,l
                idx = idx + 1
                write(30,'(A,I0,AI0,A,I0,AF25.16)') "Q",a,"(",l,",",m,")",mp(idx)
            end do
        end do
    end do
    close(30)
end subroutine write_multipole_file
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! writes files for plotting in R, used in analysis mode
subroutine write_image_slice_data_analysis()
    implicit none
    integer :: i,j,ios
    
    !write XY slices
    open(30, file="slices/sliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceXY.csv" for writing')
    open(31, file="slices/truesliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceXY.csv" for writing')
    open(32, file="slices/fullsliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceXY.csv" for writing')
    open(33, file="slices/truefullsliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceXY.csv" for writing')
    do i = 1,NgridX
        do j = 1,NgridY
            !determine whether to write value or not
            if(usedXY(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') sliceXY(i,j),','
                write(31,'(ES13.5,A)', advance='no') sliceXY2(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') sliceXY(i,j),','
            write(33,'(ES13.5,A)', advance='no') sliceXY2(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)
    
    !write XZ slices
    open(30, file="slices/sliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceXZ.csv" for writing')
    open(31, file="slices/truesliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceXZ.csv" for writing')
    open(32, file="slices/fullsliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceXZ.csv" for writing')
    open(33, file="slices/truefullsliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceXZ.csv" for writing')
    do i = 1,NgridX
        do j = 1,NgridZ
            !determine whether to write value or not
            if(usedXZ(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') sliceXZ(i,j),','
                write(31,'(ES13.5,A)', advance='no') sliceXZ2(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') sliceXZ(i,j),','
            write(33,'(ES13.5,A)', advance='no') sliceXZ2(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)  
    
    !write YZ slices
    open(30, file="slices/sliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceYZ.csv" for writing')
    open(31, file="slices/truesliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceYZ.csv" for writing')
    open(32, file="slices/fullsliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceYZ.csv" for writing')
    open(33, file="slices/truefullsliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceYZ.csv" for writing')
    do i = 1,NgridY
        do j = 1,NgridZ
            !determine whether to write value or not
            if(usedYZ(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') sliceYZ(i,j),','
                write(31,'(ES13.5,A)', advance='no') sliceYZ2(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') sliceYZ(i,j),','
            write(33,'(ES13.5,A)', advance='no') sliceYZ2(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)  

end subroutine write_image_slice_data_analysis
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! writes files for plotting in R
subroutine write_image_slice_data(charges)
    implicit none
    real(rp), dimension(qdim) :: charges
    real(rp), dimension(3) :: x ! position
    integer :: i,j,ios
    
    !write XY slices
    open(30, file="slices/sliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceXY.csv" for writing')
    open(31, file="slices/truesliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceXY.csv" for writing')
    open(32, file="slices/fullsliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceXY.csv" for writing')
    open(33, file="slices/truefullsliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceXY.csv" for writing')
    do i = 1,NgridX
        do j = 1,NgridY
            x = origin + (i-1)*axisX + (j-1)*axisY + (NgridZ/2-1)*axisZ
            !determine whether to write value or not
            if(usedXY(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') coulomb_potential_qtot(x,charges),','
                write(31,'(ES13.5,A)', advance='no') sliceXY(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') coulomb_potential_qtot(x,charges),','
            write(33,'(ES13.5,A)', advance='no') sliceXY(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)
    
    !write XZ slices
    open(30, file="slices/sliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceXZ.csv" for writing')
    open(31, file="slices/truesliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceXZ.csv" for writing')
    open(32, file="slices/fullsliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceXZ.csv" for writing')
    open(33, file="slices/truefullsliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceXZ.csv" for writing')
    do i = 1,NgridX
        do j = 1,NgridZ
            x = origin + (i-1)*axisX + (NgridY/2-1)*axisY + (j-1)*axisZ
            !determine whether to write value or not
            if(usedXZ(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') coulomb_potential_qtot(x,charges),','
                write(31,'(ES13.5,A)', advance='no') sliceXZ(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') coulomb_potential_qtot(x,charges),','
            write(33,'(ES13.5,A)', advance='no') sliceXZ(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)  
    
    !write YZ slices
    open(30, file="slices/sliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceYZ.csv" for writing')
    open(31, file="slices/truesliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceYZ.csv" for writing')
    open(32, file="slices/fullsliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceYZ.csv" for writing')
    open(33, file="slices/truefullsliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceYZ.csv" for writing')
    do i = 1,NgridY
        do j = 1,NgridZ
            x = origin + (NgridX/2-1)*axisX + (i-1)*axisY + (j-1)*axisZ
            !determine whether to write value or not
            if(usedYZ(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') coulomb_potential_qtot(x,charges),','
                write(31,'(ES13.5,A)', advance='no') sliceYZ(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') coulomb_potential_qtot(x,charges),','
            write(33,'(ES13.5,A)', advance='no') sliceYZ(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)  

end subroutine write_image_slice_data
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! writes files for plotting in R
subroutine write_image_slice_data_multipole(multipole)
    implicit none
    real(rp), dimension(:) :: multipole
    real(rp), dimension(3) :: x ! position
    integer :: i,j,ios
    
    !write XY slices
    open(30, file="slices/sliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceXY.csv" for writing')
    open(31, file="slices/truesliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceXY.csv" for writing')
    open(32, file="slices/fullsliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceXY.csv" for writing')
    open(33, file="slices/truefullsliceXY.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceXY.csv" for writing')
    do i = 1,NgridX
        do j = 1,NgridY
            x = origin + (i-1)*axisX + (j-1)*axisY + (NgridZ/2-1)*axisZ
            !determine whether to write value or not
            if(usedXY(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') coulomb_potential_multipole(x,multipole),','
                write(31,'(ES13.5,A)', advance='no') sliceXY(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') coulomb_potential_multipole(x,multipole),','
            write(33,'(ES13.5,A)', advance='no') sliceXY(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)
    
    !write XZ slices
    open(30, file="slices/sliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceXZ.csv" for writing')
    open(31, file="slices/truesliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceXZ.csv" for writing')
    open(32, file="slices/fullsliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceXZ.csv" for writing')
    open(33, file="slices/truefullsliceXZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceXZ.csv" for writing')
    do i = 1,NgridX
        do j = 1,NgridZ
            x = origin + (i-1)*axisX + (NgridY/2-1)*axisY + (j-1)*axisZ
            !determine whether to write value or not
            if(usedXZ(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') coulomb_potential_multipole(x,multipole),','
                write(31,'(ES13.5,A)', advance='no') sliceXZ(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') coulomb_potential_multipole(x,multipole),','
            write(33,'(ES13.5,A)', advance='no') sliceXZ(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)  
    
    !write YZ slices
    open(30, file="slices/sliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/sliceYZ.csv" for writing')
    open(31, file="slices/truesliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truesliceYZ.csv" for writing')
    open(32, file="slices/fullsliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/fullsliceYZ.csv" for writing')
    open(33, file="slices/truefullsliceYZ.csv", status="replace", action="write", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "slices/truefullsliceYZ.csv" for writing')
    do i = 1,NgridY
        do j = 1,NgridZ
            x = origin + (NgridX/2-1)*axisX + (i-1)*axisY + (j-1)*axisZ
            !determine whether to write value or not
            if(usedYZ(i,j)) then 
                write(30,'(ES13.5,A)', advance='no') coulomb_potential_multipole(x,multipole),','
                write(31,'(ES13.5,A)', advance='no') sliceYZ(i,j),','
            else
                write(30,'(A)', advance='no') 'NA,'
                write(31,'(A)', advance='no') 'NA,'
            end if
            !full slices are always written
            write(32,'(ES13.5,A)', advance='no') coulomb_potential_multipole(x,multipole),','
            write(33,'(ES13.5,A)', advance='no') sliceYZ(i,j),','
        end do
        write(30,*) !line break
        write(31,*) !line break
        write(32,*) !line break
        write(33,*) !line break
    end do
    close(30)
    close(31)
    close(32)
    close(33)  

end subroutine write_image_slice_data_multipole
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! writes a cube file containing the "error surface" (used in analysis mode)
subroutine write_error_cube_file()
implicit none 
character(len=1024) outfile
integer :: i,j,k,lcount,ios,idx

outfile = trim(prefix)//"error.cube"

open(30, file=trim(outfile), status="replace", action="write", iostat = ios)
if(ios /= 0) call throw_error('Could not open "'//trim(outfile)//'" for writing')
!write header
write(30,'(1X,A)')       "Error surface"
write(30,'(1X,A,I0,A)')  "Difference of true and fitted electrostatic potential"
write(30,'(I5,3(F12.6),I5)') Natom, origin, 1
write(30,'(I5,3(F12.6))')    NgridX, axisX
write(30,'(I5,3(F12.6))')    NgridY, axisY
write(30,'(I5,3(F12.6))')    NgridZ, axisZ
do i = 1,Natom
    write(30,'(I5,4(F12.6))') atom_num(i),real(atom_num(i),rp),atom_pos(:,i)
end do
!write grid data
lcount = 0
idx = 0
do i = 1,NgridX
    do j = 1,NgridY
        do k = 1,NgridZ
            idx = idx + 1
            lcount = lcount + 1
            write(30,'(ES13.5)', advance='no') esp_grid2(idx)-esp_grid(idx)
            !to accomodate weird format, we sometimes have to skip lines
            if(mod(lcount,6) == 0 .or. mod(lcount,NgridZ) == 0) then
                write(30,*) ! line break
                if(mod(lcount,NgridZ) == 0) lcount = 0
            end if
        end do
    end do
end do
close(30)

end subroutine write_error_cube_file
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! writes a cube file containing the esp for direct comparison
subroutine write_cube_file(charges,multipole)
implicit none 
character(len=1024) outfile, dummy
real(rp), dimension(qdim) :: charges
real(rp) :: x(3) ! dummy vector
integer, optional :: multipole
integer :: i,j,k,lcount,ios

if(present(multipole)) then
    write(dummy,'(I0)') multipole
    write(outfile,'(I0)') num_charges
    outfile = trim(prefix)//"multipole"//trim(dummy)//"_"//trim(outfile)//"charges.cube"
else
    write(outfile,'(I0)') num_charges
    if(trim(prefix) /= '') then
        outfile = trim(prefix)//"/"//trim(outfile)//"charges.cube"
    else
        outfile = trim(outfile)//"charges.cube"
    end if
end if

print*,'writing cube file ',outfile
open(30, file=trim(outfile), status="replace", action="write", iostat = ios)
if(ios /= 0) call throw_error('Could not open "'//trim(outfile)//'" for writing')
!write header
write(30,'(1X,A)')       "Fitted ESP"
write(30,'(1X,A,I0,A)')  "Electrostatic potential fitted from ",num_charges," point charges."
write(30,'(I5,3(F12.6),I5)') Natom, origin, 1
write(30,'(I5,3(F12.6))')    NgridX, axisX
write(30,'(I5,3(F12.6))')    NgridY, axisY
write(30,'(I5,3(F12.6))')    NgridZ, axisZ
do i = 1,Natom
    write(30,'(I5,4(F12.6))') atom_num(i),real(atom_num(i),rp),atom_pos(:,i)
end do
!write grid data
lcount = 0
do i = 1,NgridX
    do j = 1,NgridY
        do k = 1,NgridZ
            lcount = lcount + 1
            x = origin + (i-1)*axisX + (j-1)*axisY + (k-1)*axisZ
            write(30,'(ES13.5)', advance='no') coulomb_potential_qtot(x,charges(1:qdim))
            !to accomodate weird format, we sometimes have to skip lines
            if(mod(lcount,6) == 0 .or. mod(lcount,NgridZ) == 0) then
                write(30,*) ! line break
                if(mod(lcount,NgridZ) == 0) lcount = 0
            end if
        end do
    end do
end do
close(30)

end subroutine write_cube_file
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! writes a cube file containing the esp for direct comparison
subroutine write_cube_file_multipole(multipole)
implicit none 
character(len=1024) outfile, dummy
real(rp), dimension(:) :: multipole
real(rp) :: x(3) ! dummy vector
real(rp) vtot
integer :: i,j,k,lcount,ios,a,b

!decide filename
select case (lcur)
    case(0) 
        write(dummy,'(A)') "monopole"
    case(1)
        write(dummy,'(A)') "dipole"
    case(2) 
        write(dummy,'(A)') "quadrupole"
    case(3)
        write(dummy,'(A)') "octopole"
    case(4) 
        write(dummy,'(A)') "hexadecapole"
    case(5)
        write(dummy,'(A)') "ditriantapole"
end select       

if(trim(prefix) /= '') then
    outfile = trim(prefix)//"/"//trim(dummy)//"_expansion.cube"
else
    outfile = trim(dummy)//"_expansion.cube"
end if

print*,'writing cube file ',outfile
open(30, file=trim(outfile), status="replace", action="write", iostat = ios)
if(ios /= 0) call throw_error('Could not open "'//trim(outfile)//'" for writing')
!write header
write(30,'(1X,A)')       "Fitted ESP"
write(30,'(1X,A,I0,A)')  "Electrostatic potential fitted from "//trim(dummy)//" expansion."
write(30,'(I5,3(F12.6),I5)') Natom, origin, 1
write(30,'(I5,3(F12.6))')    NgridX, axisX
write(30,'(I5,3(F12.6))')    NgridY, axisY
write(30,'(I5,3(F12.6))')    NgridZ, axisZ
do i = 1,Natom
    write(30,'(I5,4(F12.6))') atom_num(i),real(atom_num(i),rp),atom_pos(:,i)
end do
!write grid data
lcount = 0
do i = 1,NgridX
    do j = 1,NgridY
        do k = 1,NgridZ
            lcount = lcount + 1
            x = origin + (i-1)*axisX + (j-1)*axisY + (k-1)*axisZ
            if(natmfit < Natom) then
              vtot=0._rp
              do b=1,natmfit
                a=fitatoms(b)
                vtot=vtot + coulomb_potential_single_multipole(x,multipole,a)
              end do
            else 
              vtot=coulomb_potential_multipole(x,multipole)
            end if
            write(30,'(ES13.5)', advance='no') vtot
            !to accomodate weird format, we sometimes have to skip lines
            if(mod(lcount,6) == 0 .or. mod(lcount,NgridZ) == 0) then
                write(30,*) ! line break
                if(mod(lcount,NgridZ) == 0) lcount = 0
            end if
        end do
    end do
end do
close(30)

end subroutine write_cube_file_multipole
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine calc_multipole_grid_and_slice_data(mp,a)
    implicit none
    integer :: a !input atom
    real(rp), dimension(:) :: mp !input multipole
    integer :: i,j
    
    !calculate the esp grid generated from just this multipole
    do i = 1,Ngrid
        esp_grid(i) = coulomb_potential_single_multipole(gridval(:,i),mp,a) 
    end do
    
    !calculate the slice data from just this multipole
    do i = 1,NgridX
        do j = 1,NgridY
                sliceXY(i,j) = coulomb_potential_single_multipole(origin + (i-1)*axisX + (j-1)*axisY + (NgridZ/2-1)*axisZ,mp,a)             
        end do
    end do
    
    do i = 1,NgridX
        do j = 1,NgridZ
                sliceXZ(i,j) = coulomb_potential_single_multipole(origin + (i-1)*axisX + (NgridY/2-1)*axisY + (j-1)*axisZ,mp,a)             
        end do
    end do
    
    do i = 1,NgridY
        do j = 1,NgridZ
                sliceYZ(i,j) = coulomb_potential_single_multipole(origin + (NgridX/2-1)*axisX + (i-1)*axisY + (j-1)*axisZ,mp,a)             
        end do
    end do
    
end subroutine calc_multipole_grid_and_slice_data
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine subtract_atom_multipole_ESP_from_grid(mp,a)
    implicit none
    integer :: a !input atom
    real(rp), dimension(:) :: mp !input multipole
    integer :: i,j

    !calculate the esp grid generated from just this multipole and subtract from
    !total
    do i = 1,Ngrid
        esp_grid(i) = esp_grid(i) - coulomb_potential_single_multipole(gridval(:,i),mp,a)
    end do

    !calculate the slice data from just this multipole
    do i = 1,NgridX
        do j = 1,NgridY
                sliceXY(i,j) = sliceXY(i,j) - coulomb_potential_single_multipole( &
                     origin + (i-1)*axisX + (j-1)*axisY + (NgridZ/2-1)*axisZ,mp,a)
        end do
    end do

    do i = 1,NgridX
        do j = 1,NgridZ
                sliceXZ(i,j) = sliceXZ(i,j) - coulomb_potential_single_multipole( &
                     origin + (i-1)*axisX + (NgridY/2-1)*axisY + (j-1)*axisZ,mp,a)
        end do
    end do

    do i = 1,NgridY
        do j = 1,NgridZ
                sliceYZ(i,j) = sliceYZ(i,j) - coulomb_potential_single_multipole( &
                     origin + (NgridX/2-1)*axisX + (i-1)*axisY + (j-1)*axisZ,mp,a)
        end do
    end do

end subroutine subtract_atom_multipole_ESP_from_grid
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine add_atom_multipole_ESP_to_grid(mp,a)
    implicit none
    integer :: a !input atom
    real(rp), dimension(:) :: mp !input multipole
    integer :: i,j

    !calculate the esp grid generated from just this multipole and subtract from
    !total
    do i = 1,Ngrid
        esp_grid(i) = esp_grid(i) + coulomb_potential_single_multipole(gridval(:,i),mp,a)
    end do

    !calculate the slice data from just this multipole
    do i = 1,NgridX
        do j = 1,NgridY
                sliceXY(i,j) = sliceXY(i,j) + coulomb_potential_single_multipole( &
                     origin + (i-1)*axisX + (j-1)*axisY + (NgridZ/2-1)*axisZ,mp,a)
        end do
    end do

    do i = 1,NgridX
        do j = 1,NgridZ
                sliceXZ(i,j) = sliceXZ(i,j) + coulomb_potential_single_multipole( &
                     origin + (i-1)*axisX + (NgridY/2-1)*axisY + (j-1)*axisZ,mp,a)
        end do
    end do

    do i = 1,NgridY
        do j = 1,NgridZ
                sliceYZ(i,j) = sliceYZ(i,j) + coulomb_potential_single_multipole( &
                     origin + (NgridX/2-1)*axisX + (i-1)*axisY +(j-1)*axisZ,mp,a)
        end do
    end do

end subroutine add_atom_multipole_ESP_to_grid
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! read Gaussian cube files
subroutine read_cube_file(filepath,density_filepath)
    implicit none
    character(len=*), intent(in)           :: filepath
    character(len=*), intent(in), optional :: density_filepath !for density cutoff
    character(len=128) :: ctmp ! dummy character variable
    real(rp) :: x(3) ! dummy vector
    real(rp) :: rtmp,rtmp2 ! dummy real variable
    integer :: ios ! keeps track of io status
    integer :: i,j,k,idx,lcount
    
    if(use_vdW_grid_cutoff .and. use_density_grid_cutoff) &
        call throw_error('Only one grid cutoff scheme can be selected at once.')
        
    if(use_density_grid_cutoff.and..not.present(density_filepath)) &
        call throw_error('Density grid cutoff scheme was selected, but no density cube file was given.')
    
    
    
    if(verbose) write(*,'(A)') 'Reading "'//trim(filepath)//'"...'
    if(verbose) write(*,*)
    
    !open cube files
    open(30, file=trim(filepath), status= "old", action= "read", iostat = ios)
    if(ios /= 0) call throw_error('Could not open "'//trim(filepath)//'".')
    if(use_density_grid_cutoff) then
        open(31, file=trim(density_filepath), status= "old", action= "read", iostat = ios)
        if(ios /= 0) call throw_error('Could not open "'//trim(density_filepath)//'".')
    end if
    
    ! skip the title in the header
    read(30,'(A128)',iostat = ios) ctmp
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Bad Format.')
    if(verbose) write(*,*) trim(ctmp)
    read(30,'(A128)',iostat = ios) ctmp
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Bad Format.')
    if(verbose) write(*,*) trim(ctmp)
    if(use_density_grid_cutoff) then
        read(31,'(A128)',iostat = ios) ctmp
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Bad Format.')
        read(31,'(A128)',iostat = ios) ctmp
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Bad Format.')
    end if
    
    
    
    ! read information about coordinate system and the number of atoms
    read(30,*,iostat = ios) Natom, origin(:)
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Missing Header?')
    read(30,*,iostat = ios) NgridX, axisX(:)
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Missing Header?')
    read(30,*,iostat = ios) NgridY, axisY(:)
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Missing Header?')
    read(30,*,iostat = ios) NgridZ, axisZ(:)
    if(ios /= 0) call throw_error('Could not read "'//trim(filepath)//'". Missing Header?')
    if(use_density_grid_cutoff) then ! skip this information in density cube file
        read(31,*,iostat = ios) 
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Missing Header?')
        read(31,*,iostat = ios) 
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Missing Header?')
        read(31,*,iostat = ios) 
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Missing Header?')
        read(31,*,iostat = ios) 
        if(ios /= 0) call throw_error('Could not read "'//trim(density_filepath)//'". Missing Header?')
    end if
    
    ! initialize list of atoms to be fitted if nothing was defined with -atom flag
    if(natmfit == 0) then
      natmfit = Natom
      allocate(fitatoms(natmfit))
      do i=1,Natom
        fitatoms(i)=i
      enddo
    endif

    if(verbose) write(*,'(I5,3F12.6)') Natom, origin(:)
    if(verbose) write(*,'(I5,3F12.6)') NgridX, axisX(:)
    if(verbose) write(*,'(I5,3F12.6)') NgridY, axisY(:)
    if(verbose) write(*,'(I5,3F12.6)') NgridZ, axisZ(:) 
    
    ! allocate memory to store slices through the potential (for visualization with R)
    if(.not.allocated(sliceXY)) allocate(sliceXY(NgridX,NgridY), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for sliceXY.')
    if(.not.allocated(sliceXZ)) allocate(sliceXZ(NgridX,NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for sliceXZ.')
    if(.not.allocated(sliceYZ)) allocate(sliceYZ(NgridY,NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for sliceYZ.')
    if(.not.allocated(usedXY))  allocate(usedXY (NgridX,NgridY), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for usedXY.')
    if(.not.allocated(usedXZ))  allocate(usedXZ (NgridX,NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for usedXZ.')
    if(.not.allocated(usedYZ))  allocate(usedYZ (NgridY,NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory for usedYZ.')
    usedXY = .true.; usedXZ = .true.; usedYZ = .true.; !initialize to true
    
    
    ! allocate memory to store atom information
    if(.not.allocated(atom_num)) allocate(atom_num(Natom),   stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    if(.not.allocated(atom_pos)) allocate(atom_pos(3,Natom), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    
    ! read atom information
    do i = 1,Natom
        read(30,*,iostat=ios) atom_num(i), rtmp, atom_pos(:,i)
        if(ios /= 0) call throw_error('Could not read atom information.')
        if(verbose) write(*,'(I5,4F12.6)') atom_num(i), real(atom_num(i),rp), atom_pos(:,i)
    end do
    if(use_density_grid_cutoff) then ! skip atom information in density cube file
        do i = 1,Natom
            read(31,*,iostat=ios) 
            if(ios /= 0) call throw_error('Density cube file has unmatching format.')
        end do
    end if
        
    ! allocate memory to store grid information
    if(.not.allocated(esp_grid)) allocate(esp_grid(NgridX*NgridY*NgridZ), stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    if(.not.allocated(gridval))  allocate(gridval(3,NgridX*NgridY*NgridZ),stat=ios)
    if(ios /= 0) call throw_error('Could not allocate memory.')
    
    
    ! find the "interesting" grid points
    if(use_vdW_grid_cutoff) then !based on vdW radii                
        ! read grid information (the grid is stored in a flat array)
        ! for cache friendliness (looping over these values is
        ! performance critical for cost function evaluation)
        idx = 0
        lcount = 0
        do i = 1,NgridX
            do j = 1,NgridY
                do k = 1,NgridZ
                    lcount = lcount + 1
                    x = origin + (i-1)*axisX + (j-1)*axisY + (k-1)*axisZ
                    if(in_interaction_belt(x,vdw_grid_min_cutoff, &
                       vdw_grid_max_cutoff)) then !value gets added
                        idx = idx + 1
                        read(30,'(ES13.5)',advance='no',iostat = ios) esp_grid(idx)
                        if(ios /= 0) call throw_error('Could not read ESP grid.')
                        gridval(:,idx) = x
                    else
                        read(30,'(ES13.5)',advance='no',iostat = ios) rtmp
                        if(ios /= 0) call throw_error('Could not read ESP grid.')
                    end if
                    
                    !to accomodate weird format, we sometimes have to skip lines
                    if(mod(lcount,6) == 0 .or. mod(lcount,NgridZ) == 0) then
                        read(30,*) ! line break
                        if(mod(lcount,NgridZ) == 0) lcount = 0
                    end if
                    
                    ! everything below is for later visualization with R (very useful to assess quality)
                    if( i == NgridX/2 ) then  
                        if(.not.in_interaction_belt(x,vdw_grid_min_cutoff, &
                           vdw_grid_max_cutoff)) then !value is NOT used for fitting
                            usedYZ(j,k)  = .false.
                            sliceYZ(j,k) = rtmp
                        else
                            sliceYZ(j,k) = esp_grid(idx)
                        end if
                    end if
                    if( j == NgridY/2 ) then
                        if(.not.in_interaction_belt(x,vdw_grid_min_cutoff, &
                           vdw_grid_max_cutoff)) then !value is NOT used for fitting
                            usedXZ(i,k)  = .false.
                            sliceXZ(i,k) = rtmp
                        else
                            sliceXZ(i,k) = esp_grid(idx)
                        end if
                    end if
                    if( k == NgridZ/2 ) then
                        if(.not.in_interaction_belt(x,vdw_grid_min_cutoff, &
                           vdw_grid_max_cutoff)) then !value is NOT used for fitting
                            usedXY(i,j)  = .false.
                            sliceXY(i,j) = rtmp
                        else
                            sliceXY(i,j) = esp_grid(idx)
                        end if
                    end if
                end do
            end do
        end do  
    else if(use_density_grid_cutoff) then ! based on density cube file        
        idx = 0
        lcount = 0
        do i = 1,NgridX
            do j = 1,NgridY
                do k = 1,NgridZ
                    lcount = lcount + 1
                    read(31,'(ES13.5)',advance='no',iostat = ios) rtmp
                    if(ios /= 0) call throw_error('Density cube file has wrong format.')
                    if((rtmp < density_grid_max_cutoff) .and. (rtmp > density_grid_min_cutoff)) then !add point
                        x = origin + (i-1)*axisX + (j-1)*axisY + (k-1)*axisZ
                        idx = idx + 1
                        read(30,'(ES13.5)',advance='no',iostat = ios) esp_grid(idx)
                        if(ios /= 0) call throw_error('Could not read ESP grid.')
                        gridval(:,idx) = x
                    else !skip point
                        read(30,'(ES13.5)',advance='no',iostat = ios) rtmp2
                        if(ios /= 0) call throw_error('Could not read ESP grid.')
                    end if
                    !to accomodate weird format, we sometimes have to skip lines
                    if(mod(lcount,6) == 0 .or. mod(lcount,NgridZ) == 0) then
                        read(30,*) ! line break
                        read(31,*)
                        if(mod(lcount,NgridZ) == 0) lcount = 0
                    end if
                    
                    if( i == NgridX/2 ) then  
                        if(.not.((rtmp < density_grid_max_cutoff) .and. (rtmp > density_grid_min_cutoff))) then !value is NOT used for fitting
                            usedYZ(j,k)  = .false.
                            sliceYZ(j,k) = rtmp2
                        else
                            sliceYZ(j,k) = esp_grid(idx)
                        end if
                    end if
                    if( j == NgridY/2 ) then
                        if(.not.((rtmp < density_grid_max_cutoff) .and. (rtmp > density_grid_min_cutoff))) then !value is NOT used for fitting
                            usedXZ(i,k)  = .false.
                            sliceXZ(i,k) = rtmp2
                        else
                            sliceXZ(i,k) = esp_grid(idx)
                        end if
                    end if
                    if( k == NgridZ/2 ) then
                        if(.not.((rtmp < density_grid_max_cutoff) .and. (rtmp > density_grid_min_cutoff))) then !value is NOT used for fitting
                            usedXY(i,j)  = .false.
                            sliceXY(i,j) = rtmp2
                        else
                            sliceXY(i,j) = esp_grid(idx)
                        end if
                    end if
                    
                end do
            end do
        end do  
        close(31)
    end if
    !store number of interesting grid points
    Ngrid = idx  
    Ngridr = real(Ngrid,rp) 
    
    if(verbose) write(*,'(6ES13.5)') esp_grid(1:6)
    if(verbose) write(*,'(6(A8,5X))') ".",".",".",".",".","."
    if(verbose) write(*,'(6(A8,5X))') ".",".",".",".",".","."
    if(verbose) write(*,'(6(A8,5X))') ".",".",".",".",".","."
    if(verbose) write(*,'(6ES13.5)') esp_grid(Ngrid-5:Ngrid)
    if(verbose) write(*,*)
    if(verbose) write(*,'(I0,A,I0,A,F5.1,A)') Ngrid, ' out of ', NgridX*NgridY*NgridZ, &
            ' gridpoints are considered (',100*Ngrid/real(NgridX*NgridY*NgridZ,rp),'%).'
    
    close(30)
    
    if(verbose) write(*,*)
    if(verbose) write(*,'(A)') '...done!'
    if(verbose) write(*,*)

    return
end subroutine read_cube_file
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! multipole-charge interaction (aka ESP) given in the appendix of
! The Theory of Intermolecular Forces by A. J. Stone
real(rp) function mESP(l,m,rmag,r)
    implicit none
    integer, intent(in)  :: l,m !l and m index of spherical harmonic
    real(rp), intent(in) :: rmag, r(3)   !distance and unit vector in distance direction
    !just for staying closer to the maths formula
    integer, parameter :: x = 1, y = 2, z = 3 
    !for saving square root calculations
    real(rp), parameter :: sqrt3   = sqrt(3._rp)   ,&
                           sqrt5   = sqrt(5._rp)   ,&
                           sqrt6   = sqrt(6._rp)   ,&
                           sqrt10  = sqrt(10._rp)  ,&
                           sqrt14  = sqrt(14._rp)  ,&
                           sqrt15  = sqrt(15._rp)  ,&
                           sqrt35  = sqrt(35._rp)  ,&
                           sqrt70  = sqrt(70._rp)  ,&
                           sqrt105 = sqrt(105._rp)  
    real(rp) :: Rpow    
    Rpow = 1._rp/rmag**(l+1) ! calculate the appropriate power law

    select case(l)
        ! monopole     
        case(0) !00
            mESP = Rpow
        ! dipole  
        case(1)
            select case(m)
                case( 0) !10
                    mESP = Rpow * r(z)
                case( 1) !11c
                    mESP = Rpow * r(x)
                case(-1) !11s   
                    mESP = Rpow * r(y)
                case default
                    call throw_error("multipole_interaction must be called with m between -l and l")
            end select
        ! quadrupole
        case(2)
            select case(m)
                case( 0) !20
                    mESP = Rpow * 0.5_rp * (3*r(z)**2 - 1)
                case( 1) !21c
                    mESP = Rpow * sqrt3 * r(x) * r(z)
                case(-1) !21s
                    mESP = Rpow * sqrt3 * r(y) * r(z) 
                case( 2) !22c
                    mESP = Rpow * 0.5_rp * sqrt3 * (r(x)**2 - r(y)**2)
                case(-2) !22s
                    mESP = Rpow * sqrt3 * r(x) * r(y)
                case default
                    call throw_error("multipole_interaction must be called with m between -l and l")
            end select 
        ! octopole
        case(3)
            select case(m)
                case( 0) !30
                    mESP = Rpow * 0.5_rp * (5*r(z)**3 - 3*r(z))
                case( 1) !31c
                    mESP = Rpow * 0.25_rp * sqrt6 * r(x) * (5*r(z)**2 - 1)
                case(-1) !31s
                    mESP = Rpow * 0.25_rp * sqrt6 * r(y) * (5*r(z)**2 - 1) 
                case( 2) !32c
                    mESP = Rpow * 0.5_rp * sqrt15 * r(z) * (r(x)**2 - r(y)**2)
                case(-2) !32s
                    mESP = Rpow * sqrt15 * r(x) * r(y) * r(z)
                case( 3) !33c
                    mESP = Rpow * 0.25_rp * sqrt10 * r(x) * (r(x)**2 - 3*r(y)**2)
                case(-3) !33s
                    mESP = Rpow * 0.25_rp * sqrt10 * r(y) * (3*r(x)**2 - r(y)**2)
                case default
                    call throw_error("multipole_interaction must be called with m between -l and l")
            end select 
        ! hexadecapole        
        case(4)
            select case(m)
                case( 0) !40
                    mESP = Rpow * 0.125_rp * (35*r(z)**4 - 30*r(z)**2 + 3)
                case( 1) !41c
                    mESP = Rpow * 0.25_rp * sqrt10 * (7*r(x)*r(z)**3 - 3*r(x)*r(z))
                case(-1) !41s
                    mESP = Rpow * 0.25_rp * sqrt10 * (7*r(y)*r(z)**3 - 3*r(y)*r(z)) 
                case( 2) !42c
                    mESP = Rpow * 0.25_rp * sqrt5 * (7*r(z)**2 -1) * (r(x)**2 - r(y)**2)
                case(-2) !42s
                    mESP = Rpow * 0.5_rp * sqrt5 * (7*r(z)**2 -1) * r(x) * r(y)
                case( 3) !43c
                    mESP = Rpow * 0.25_rp * sqrt70 * r(x) * r(z) * (r(x)**2 - 3*r(y)**2)
                case(-3) !43s
                    mESP = Rpow * 0.25_rp * sqrt70 * r(y) * r(z) * (3*r(x)**2 - r(y)**2)
                case( 4) !44c
                    mESP = Rpow * 0.125_rp * sqrt35 * (r(x)**4 - 6*r(x)**2*r(y)**2 + r(y)**4)
                case(-4) !44s
                    mESP = Rpow * 0.5_rp * sqrt35 * r(x) * r(y) * (r(x)**2 - r(y)**2)
                case default
                    call throw_error("multipole_interaction must be called with m between -l and l")
            end select
        ! diatriantapole
        case(5)
            select case(m)
                case( 0) !50
                    mESP = Rpow * 0.125_rp * (63*r(z)**5 - 70*r(z)**3 + 15*r(z))
                case( 1) !51c
                    mESP = Rpow * 0.125_rp * sqrt15 * (21*r(x)*r(z)**4 - 14*r(x)*r(z)**2 + r(x))
                case(-1) !51s
                    mESP = Rpow * 0.125_rp * sqrt15 * (21*r(y)*r(z)**4 - 14*r(y)*r(z)**2 + r(y))
                case( 2) !52c
                    mESP = Rpow * 0.25_rp * sqrt105 * (3*r(x)**2*r(z)**3 - 3*r(y)**2*r(z)**3 - r(x)**2*r(z) + r(y)**2*r(z))
                case(-2) !52s
                    mESP = Rpow * 0.5_rp * sqrt105 * (3*r(x)*r(y)*r(z)**3 - r(x)*r(y)*r(z))
                case( 3) !53c
                    mESP = Rpow * 0.0625_rp * sqrt70 * (9*r(x)**3*r(z)**2 - 27*r(x)*r(y)**2*r(z)**2 - r(x)**3 + 3*r(x)*r(y)**2)
                case(-3) !53s
                    mESP = Rpow * 0.0625_rp * sqrt70 * (27*r(x)**2*r(y)*r(z)**2 - 9*r(y)**3*r(z)**2 - 3*r(x)**2*r(y) + r(y)**3)
                case( 4) !54c
                    mESP = Rpow * 0.375_rp * sqrt35 * (r(x)**4*r(z) - 6*r(x)**2*r(y)**2*r(z) + r(y)**4*r(z))
                case(-4) !54s
                    mESP = Rpow * 1.5_rp * sqrt35 * (r(x)**3*r(y)*r(z) - r(x)*r(y)**3*r(z))
                case( 5) !55c
                    mESP = Rpow * 0.1875_rp * sqrt14 * (r(x)**5 - 10*r(x)**3*r(y)**2 + 5*r(x)*r(y)**4)
                case(-5) !55s
                    mESP = Rpow * 0.1875_rp * sqrt14 * (5*r(x)**4*r(y) - 10*r(x)**2*r(y)**3 + r(y)**5)
                case default
                    call throw_error("multipole_interaction must be called with m between -l and l")
            end select 
        !something must be wrong
        case default
            call throw_error("multipole_interaction must be called with l between 0 and 5")
    end select


end function mESP
!-------------------------------------------------------------------------------





!-------------------------------------------------------------------------------
! throws an error message and terminates the code
subroutine throw_error(message)
    implicit none
    character(len=*), intent(in) :: message
    write(*,'(A)') "ERROR: "//message
    call dealloc()
    stop   
end subroutine throw_error
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read the command line arguments
subroutine read_command_line_arguments()
    implicit none
    integer :: i, l, cmd_count, ios
    character(len=1024) :: arg, bla
    
    !read argument count
    cmd_count = command_argument_count()

    natmfit = 0 ! initialize number of atoms to fit
    
    if(cmd_count == 0) then
        write(*,'(A)') "You have not provided any command line arguments"
        write(*,'(A)') "Expected usage (flags in [] are optional, arguments in <> need to be replaced by the user): "
        write(*,*)
        write(*,'(A)') "When running in analysis mode (-analysis)"
        write(*,'(A)',advance='no') "./"
        !$ write(*,'(A)',advance='no') "p"
        write(*,'(A)') "cubefit.x -analysis -esp <filepath> -esp2 <filepath> -dens <filepath> [-prefix <string>] [-v]"
        write(*,*)
        write(*,'(A)') "-esp     <string>    filepath of the cube file containing the reference ESP data" 
        write(*,'(A)') "-dens    <string>    filepath of the cube file containing the density data" 
        write(*,'(A)') "-esp2    <string>    filepath of the cube file containing the fitted ESP data" 
        write(*,'(A)') "-prefix  <string>    identifier prefix for output files"      
        write(*,'(A)') "-v                   verbose output" 
        write(*,*)
        write(*,'(A)') "When running in generate mode (-generate)"
        write(*,'(A)',advance='no') "./"
        !$ write(*,'(A)',advance='no') "p"
        write(*,'(A)') "cubefit.x -generate [-multipole] [-ofatoms] -esp <filepath> -dens <filepath>"//&
                       " -xyz <filepath> [-prefix <string>] [-v]"
        write(*,*)
        write(*,'(A)') "-multipole           if present, expects multipole input, else expects charge input" 
        write(*,'(A)') "-ofatoms             if present, visualizes charge fits to individual atomic multipoles instead"
        write(*,'(A)') "-esp     <string>    filepath of the cube file containing the desired output format" 
        write(*,'(A)') "-dens    <string>    filepath of the cube file containing the density data" 
        write(*,'(A)') "-xyz     <string>    filepath of the file containing the charges/multipoles from which to generate cubefile"  
        write(*,'(A)') "-prefix  <string>    identifier prefix for output files"      
        write(*,'(A)') "-v                   verbose output" 
        write(*,*)
        write(*,'(A)') "When running in multipole fitting mode"
        write(*,'(A)',advance='no') "./"
        !$ write(*,'(A)',advance='no') "p"
        write(*,'(A)') "cubefit.x -multipole -esp <filepath> -dens <filepath> [-qtot <real>] [-ntry <int>]"// &
                       "[-lstart <int>] [-prefix <string>] [-v]"
        write(*,*)
        write(*,'(A)') "-esp     <string>    filepath of the cube file containing the ESP data" 
        write(*,'(A)') "-dens    <string>    filepath of the cube file containing the density data" 
        write(*,'(A)') "-qtot    <real>      total charge (default = 0)"
        write(*,'(A)') "-lstart  <int>       where to start the multipole fitting (lower orders must exist already!)"
        write(*,'(A)') "-lstop   <int>       where to stop  the multipole fitting"
        write(*,'(A)') "-ntry    <int>       maximum number of trials per atomic multipole (default = 1)"   
        write(*,'(A)') "-prefix  <string>    identifier prefix for output files (needed in parallel mode to prevent overwrite)" 
        write(*,'(A)') "-v                   verbose output" 
        write(*,*)
        write(*,'(A)') "When running in charge fitting mode (default mode)"
        write(*,'(A)',advance='no') "./"
        !$ write(*,'(A)',advance='no') "p"
        write(*,'(A)') "cubefit.x -esp <filepath> -dens <filepath> [-qtot <real>] [-ncmin <int>]"//&
                                " [-ncmax <int>] [-ntry <int>] [-prefix <string>] [-greedy <filepath>]"//&
                                " [-onlymultipoles] [-xyz <filepath>] [-sym] [-v]"
        write(*,*)
        write(*,'(A)') "-atom    <int<,int>> index of atom(s) to fit"
        write(*,'(A)') "-esp     <string>    filepath of the cube file containing the ESP data" 
        write(*,'(A)') "-dens    <string>    filepath of the cube file containing the density data" 
        write(*,'(A)') "-qtot    <real>      total charge (default = 0)"
        write(*,'(A)') "-ncmin   <int>       minimum number of charges to fit (default = 2)" 
        write(*,'(A)') "-ncmax   <int>       maximum number of charges to fit (default = 2)" 
        write(*,'(A)') "-nacmin   <int>      minimum number of charges to fit to multipoles per atom (default = 1)"
        write(*,'(A)') "-nacmax   <int>      maximum number of charges to fit to multipoles per atom (default = 5)"
        write(*,'(A)') "-ntry    <int>       maximum number of trials per number of charges (default = 1)"   
        write(*,'(A)') "-prefix  <string>    identifier prefix for output files (needed in parallel mode to prevent overwrite)" 
        write(*,'(A)') "-greedy  <string>    filepath of the file containing the multipoles for the greedy fit" 
        write(*,'(A)') "-onlymultipoles      stops the greedy fit after fitting atomic multipoles" 
        write(*,'(A)') "-xyz <str1> <str2>   filepaths containing the guess (incompatible with greedy mode)", &
          " and previously-fitted multipoles" 
        write(*,'(A)') "-simplex             for refining charge models only: uses",&
          " simplex algorithm rather than D.E. (useful for final refinement with too",&
          " many charges for D.E.)"
        write(*,'(A)') "-sym                 turn on symmetry constrained mode (default is off)"      
        write(*,'(A)') "-v                   verbose output" 
        write(*,*)
        call throw_error("Could not read command line arguments")
    end if
        
    ! loop through command line arguments
    do i = 1,cmd_count
        call get_command_argument(i, arg, l)
        ! for debugging
        !write(*,'(A)'), arg(1:l)    
        
        ! read verbose flag
        if(arg(1:l) == '-v') verbose = .true. 
        
        ! read flag used in greedyfit (stops greedy fit after fitting atomic multipoles)
        if(arg(1:l) == '-onlymultipoles') greedy_only_multi = .true. 

        
        ! read symmetry flag
        if(arg(1:l) == '-sym') use_symmetry = .true. 

        ! read list of atoms to fit (fortran is terrible at this...)
        if(arg(1:l) == '-atom') then
            call get_command_argument(i+1, arg, l)
            natmfit=1
            k=1
            bla=''
            do j = 1,len_trim(arg)
              if(arg(j:j) .eq. ',') then
                natmfit=natmfit+1
                cycle
              endif
            enddo
            allocate(fitatoms(natmfit))
            natmfit=1
            do j = 1,len_trim(arg)
              if(arg(j:j) == ',') then
                read(bla,*) fitatoms(natmfit)
                bla=''
                natmfit=natmfit+1
                k=1
                cycle
              endif
              bla(k:k)=arg(j:j)
              k=k+1
            enddo
            read(bla,*) fitatoms(natmfit)
            if(ios /= 0) call throw_error('Could not read command line argument "-atom"')
        end if

        ! read analysis flag
        if(arg(1:l) == '-analysis') analysis_mode = .true. 
        
        ! read generate flag
        if(arg(1:l) == '-generate') generate_mode = .true. 
        
        ! read additional generate flag "ofatoms"
        if(arg(1:l) == '-ofatoms') generate_atomic = .true. 
      
        ! read multipole flag (if this is set, multipoles are fitted instead of charges)
        if(arg(1:l) == '-multipole') fit_multipoles = .true.
                
        ! lstart
        if(arg(1:l) == '-lstart') then
            call get_command_argument(i+1, arg, l)
            read(arg,*,iostat = ios) lstart
            if(ios /= 0) call throw_error('Could not read command line argument "-lstart"')
        end if
        
        ! lstop
        if(arg(1:l) == '-lstop') then
            call get_command_argument(i+1, arg, l)
            read(arg,*,iostat = ios) lstop
            if(ios /= 0) call throw_error('Could not read command line argument "-lstop"')
        end if
        
        ! input esp cube file
        if(arg(1:l) == '-esp') then
            call get_command_argument(i+1, arg, l)
            read(arg,'(A)',iostat = ios) input_esp_cubefile
            if(ios /= 0) call throw_error('Could not read command line argument "-esp"')
        end if
        
        ! read second input esp cube file
        if(arg(1:l) == '-esp2') then
            call get_command_argument(i+1, arg, l)
            read(arg,'(A)',iostat = ios) compare_esp_cubefile
            if(ios /= 0) call throw_error('Could not read command line argument "-esp2"')
        end if
        
        ! input xyz cube file
        if(arg(1:l) == '-xyz') then
            call get_command_argument(i+1, arg, l)
            read(arg,'(A)',iostat = ios) input_xyzfile
            if(ios /= 0) call throw_error('Could not read command line argument "-xyz"')
            call get_command_argument(i+2, arg, l)
            read(arg,'(A)',iostat = ios) input_multipolefile
        end if
       
        ! simplex refinement only
        if(arg(1:l) == '-simplex') simplex_only=.true.
 
        ! input dens cube file
        if(arg(1:l) == '-dens') then
            call get_command_argument(i+1, arg, l)
            read(arg,'(A)',iostat = ios) input_density_cubefile
            if(ios /= 0) call throw_error('Could not read command line argument "-dens"')
        end if
        
        ! total charge
        if(arg(1:l) == '-qtot') then
            call get_command_argument(i+1, arg, l)
            read(arg,*,iostat = ios) total_charge
            if(ios /= 0) call throw_error('Could not read command line argument "-qtot"')
        end if
        
        ! minimum number of charges
        if(arg(1:l) == '-ncmin') then
            call get_command_argument(i+1, arg, l)
            read(arg,*,iostat = ios) num_charges_min
            if(ios /= 0) call throw_error('Could not read command line argument "-ncmin"')
        end if
        
        ! minimum number of charges
        if(arg(1:l) == '-ncmax') then
            call get_command_argument(i+1, arg, l)
            read(arg,*,iostat = ios) num_charges_max
            if(ios /= 0) call throw_error('Could not read command line argument "-ncmax"')
        end if

        ! minimum number of charges per atom for fitting to atomic multipoles
        if(arg(1:l) == '-nacmin') then
            call get_command_argument(i+1, arg, l)
            read(arg,*,iostat = ios) num_charges_min_multipole
            if(ios /= 0) call throw_error('Could not read command line argument "-nacmin"')
        end if
        
        ! maximum number of charges per atom for fitting to atomic multipoles
        if(arg(1:l) == '-nacmax') then
            call get_command_argument(i+1, arg, l)
            read(arg,*,iostat = ios) num_charges_max_multipole
            if(ios /= 0) call throw_error('Could not read command line argument "-nacmax"')
        end if
        
        ! maximum number of trials
        if(arg(1:l) == '-ntry') then
            call get_command_argument(i+1, arg, l)
            read(arg,*,iostat = ios) num_trials
            if(ios /= 0) call throw_error('Could not read command line argument "-ntry"')
        end if
            
        ! greedy mode
        if(arg(1:l) == '-greedy') then
            use_greedy_fit = .true.
            call get_command_argument(i+1, arg, l)
            read(arg,'(A)',iostat = ios) input_multipolefile
            if(ios /= 0) call throw_error('Could not read command line argument "-greedy"')
        end if
        
        ! prefix file identifier
        if(arg(1:l) == '-prefix') then
            call get_command_argument(i+1, arg, l)
            read(arg,'(A)',iostat = ios) prefix
            if(ios /= 0) call throw_error('Could not read command line argument "-prefix"')
        end if
        
    end do
    
    if(analysis_mode.and.generate_mode) then
        call throw_error('Only one flag allowed: either "-analysis" or "-generate"')
    end if

    if(use_symmetry.and.fit_multipoles) then
        call throw_error('Only one flag allowed: either "-multipole" or "-sym"')
    end if
    
    ! check for the presence of required arguments
    if(analysis_mode) then !analysis mode
        vdw_grid_min_cutoff = 0.0_rp
        vdw_grid_max_cutoff = vbig
        if((trim(input_esp_cubefile) == '').or.(trim(compare_esp_cubefile) == '') &
           .or.(trim(input_density_cubefile) == '')) then
            call throw_error('Missing required arguments "-esp" and/or "-esp2" and/or "-dens"')
        end if
        if(verbose) write(*,'(A)') "Running in analysis mode"
    else if (generate_mode) then !generate mode
        if((trim(input_esp_cubefile) == '').or.(trim(input_xyzfile) == '') &
           .or.(trim(input_density_cubefile) == '')) then
            call throw_error('Missing required arguments "-esp" and/or "-xyz" and/or "-dens"')
        end if
        if(verbose) write(*,'(A)') "Running in generate mode"
    else !normal mode
        if((trim(input_esp_cubefile) == '').or.(trim(input_density_cubefile) == '')) then
            call throw_error('Missing required arguments "-esp" and/or "-dens"')
        end if
        if(trim(input_xyzfile) /= '') refine_solution = .true.
        if(use_greedy_fit .and. refine_solution) then
            call throw_error('Refining a solution is incompatible with greedy mode!')
        end if
        if(simplex_only .and. .not. refine_solution) then
            call throw_error('simplex algorithm should be used for refinement only!')
        end if
        if(verbose) write(*,'(A)') "Running in fitting mode"
    end if
    
    if(verbose) write(*,*)
   
end subroutine read_command_line_arguments
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! This subroutine is used to initialize the seed of the RNG to a given seed
subroutine init_random_seed(seed_number)
    integer, intent(in) :: seed_number
    integer :: i, n, iseed
    integer, dimension(:), allocatable :: seed
    
    if(seed_number == 0) then
        call system_clock(count=iseed)   
    else
        iseed = seed_number
    end if

    call random_seed(size = n)
    allocate(seed(n))
    
    seed = iseed + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
    
    return
end subroutine init_random_seed
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! deallocate memory
subroutine dealloc()
    implicit none
    if(allocated(atom_num))     deallocate(atom_num)
    if(allocated(atom_pos))     deallocate(atom_pos)
    if(allocated(esp_grid))     deallocate(esp_grid)
    if(allocated(charges))      deallocate(charges)
    if(allocated(bestcharges))  deallocate(bestcharges)
    if(allocated(search_range)) deallocate(search_range)
    if(allocated(sliceXY))      deallocate(sliceXY)
    if(allocated(usedXY))       deallocate(usedXY)
    if(allocated(sliceXZ))      deallocate(sliceXZ)
    if(allocated(usedXZ))       deallocate(usedXZ)
    if(allocated(sliceYZ))      deallocate(sliceYZ)
    if(allocated(usedYZ))       deallocate(usedYZ)   
end subroutine dealloc
!-------------------------------------------------------------------------------
end program cubefit


!!-------------------------------------------------------------------------------
!! computes relative root mean squared error of the current fit to the true esp, using constraint charges
!! (this means that the total charge must add up to a specific value)
!real(rp) function rel_rmse_qtot(qin)
!    implicit none
!    real(rp), dimension(:), intent(in)     :: qin ! input charges
!    real(rp), dimension(size(qin,dim=1)+1) :: q   ! complete charges
!    real(rp) :: qsum
!    integer :: i
!    qsum = 0._rp
!    do i = 1,size(qin,dim=1)-3,4
!        qsum = qsum + qin(i+3)
!    end do
!    q = qin
!    q(size(q,dim=1)) = total_charge - qsum
!    rel_rmse_qtot = rel_rmse(q)
!end function rel_rmse_qtot
!!-------------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------------
!! computes relative root mean squared error of the current fit to the true esp
!real(rp) function rel_rmse(q)
!    implicit none
!    real(rp), dimension(:), intent(in) :: q ! input charges
!    real(rp), dimension(3) :: x ! position
!    integer :: i,j,k,idx
!    rel_rmse = 0._rp
!    do idx = 1,Ngrid
!        rel_rmse = rel_rmse + ((coulomb_potential(gridval(:,idx),q) - esp_grid(idx))/esp_grid(idx))**2
!    end do
!    rel_rmse = sqrt(rel_rmse/Ngridr)
!end function rel_rmse
!!-------------------------------------------------------------------------------


!!-------------------------------------------------------------------------------
!! computes the maximum squared error of the current fit to the true esp, using constraint charges
!! (this means that the total charge must add up to a specific value)
!real(rp) function max_se_qtot(qin)
!    implicit none
!    real(rp), dimension(:), intent(in)     :: qin ! input charges
!    real(rp), dimension(size(qin,dim=1)+1) :: q   ! complete charges
!    real(rp) :: qsum
!    integer :: i
!    qsum = 0._rp
!    do i = 1,size(qin,dim=1)-3,4
!        qsum = qsum + qin(i+3)
!    end do
!    q = qin
!    q(size(q,dim=1)) = total_charge - qsum
!    max_se_qtot = max_se(q)
!end function max_se_qtot
!!-------------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------------
!! computes the maximum squared error of the current fit to the true esp
!real(rp) function max_se(q)
!    implicit none
!    real(rp), dimension(:), intent(in) :: q ! input charges
!    real(rp), dimension(3) :: x ! position
!    real(rp) :: se
!    integer :: idx
!    max_se = -huge(0._rp)
!    do idx = 1,Ngrid
!        se = (coulomb_potential(gridval(:,idx),q) - esp_grid(idx))**2
!        if(se > max_se) max_se = se
!    end do
!end function max_se
!!-------------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------------
!! computes the maximum relative squared error of the current fit to the true esp, using constraint charges
!! (this means that the total charge must add up to a specific value)
!real(rp) function max_rel_se_qtot(qin)
!    implicit none
!    real(rp), dimension(:), intent(in)     :: qin ! input charges
!    real(rp), dimension(size(qin,dim=1)+1) :: q   ! complete charges
!    real(rp) :: qsum
!    integer :: i
!    qsum = 0._rp
!    do i = 1,size(qin,dim=1)-3,4
!        qsum = qsum + qin(i+3)
!    end do
!    q = qin
!    q(size(q,dim=1)) = total_charge - qsum
!    max_rel_se_qtot = max_rel_se(q)
!end function max_rel_se_qtot
!!-------------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------------
!! computes the maximum relative squared error of the current fit to the true esp
!real(rp) function max_rel_se(q)
!    implicit none
!    real(rp), dimension(:), intent(in) :: q ! input charges
!    real(rp), dimension(3) :: x ! position
!    real(rp) :: rel_se
!    integer :: idx
!    max_rel_se = -huge(0._rp)
!    do idx = 1,Ngrid
!        rel_se = ((coulomb_potential(gridval(:,idx),q) - esp_grid(idx))/esp_grid(idx))**2
!        if(rel_se > max_rel_se) max_rel_se = rel_se
!    end do
!end function max_rel_se
!!-------------------------------------------------------------------------------

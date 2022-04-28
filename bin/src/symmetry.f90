!///////////////////////////////////////////////////////////////////////////////
!
!      symmetry.f90
!      Created: 31 May 2016 at 10:50
!      Authors: Oliver Unke, Mike Devereux
!
!///////////////////////////////////////////////////////////////////////////////
module symmetry

private
public :: symmetry_init, equal, isSymmetryOperation, asymmetryMagnitude, centerOfMass
public :: symmetryOperation, getShiftedPos, getChargeOps, chgsSpawned
public :: rotax_to_cartesian, refplane_to_cartesian, get_atm_sym_op, point_to_cartesian
public :: init_sym_search_range, get_chgs_spawned, spawn_sym_chgs
public :: getNumSymOps, sym_init_pars, read_sym_file, write_sym_file, sym_map_sea_q_coords
public :: init_atm_sym_search, sym_init_fit_ops
          !calc_eigenvalues, calc_eigenvectors, rotation !only for test reasons public
          
integer, parameter :: rp = kind(0d0)
!type definition of symmetry operations
type symmetryOperation
    character(len=4)        :: label !identifies the symmetryOperation, e.g. Cn, sig, i, etc.
    character(len=3)        :: typ   !type of operation (rot(ation),ref(lection),imp(roper),inv(ersion))
    real(rp), dimension(3,3) :: M     !matrix that describes the symmetry operation
    real(rp), dimension(3)   :: axis  !vector that describe rotation axis or plane normal
    integer                 :: n     !symmetry order
end type symmetryOperation          
integer, parameter :: maxOperations = 20000
          

real(rp), save :: tolerance = 0.03d0 !tolerance to consider 2 reals "equal" 0.1
real(rp), save :: tolerance_b = 0.002d0 !percentage tolerance to consider 2 reals "equal" 0.01
logical, save :: verbose = .true.

integer, save :: Natom ! number of atoms
integer, dimension(:), allocatable, save :: num_atm_sym_ops
integer, dimension(:), allocatable, save :: num_atm_fit_ops
integer, dimension(:,:), allocatable,save :: sea_ops ! array of sym ops linking sea's
integer, dimension(:,:), allocatable,save :: num_spawned ! array of charges spawned by symmetry operations, plus charges spawned by randomly placed charge plus charges spawned by a charge at the nuclear position


integer, dimension(:,:), allocatable, save :: atm_sym_ops ! array of symmetry operations to apply to charges
integer, dimension(:,:), allocatable, save :: atm_fit_ops ! array of symmetry operations for each atom to use during charge fitting

real(rp), dimension(:,:), allocatable, save :: atom_pos ! stores the atomic positions [Natom,3]
real(rp), dimension(3),                save :: atom_com, atom_axisA, atom_axisB, atom_axisC ! centre of mass of the atoms and axes
real(rp), dimension(3,3),              save :: atom_inertia !inertia tensor of atoms
real(rp),                              save :: atom_Ia, atom_Ib, atom_Ic !principal moments of inertia
integer,                              save :: atom_rotor_type ! rotor type based on moments of inertia
integer, dimension(:,:), allocatable, save :: atom_sea        ! symmetry equivalent atoms
real(rp), dimension(:,:), allocatable, save :: atom_dist       ! interdistance matrix
real(rp), dimension(:),   allocatable, save :: atom_identifier ! unique identifier for atom type, e. g. order number
                                                              ! (for charges, this could be the charge magnitude!)

public :: atom_sea
                                                     

! storage of symmetry operations and their number   
!rotations                                                         
type(symmetryOperation), dimension(maxOperations), save :: rotations 
integer , save :: num_rotations = 0   
!reflections        
type(symmetryOperation), dimension(maxOperations), save :: reflections
integer , save :: num_reflections = 0        
!improper rotations
type(symmetryOperation), dimension(maxOperations), save :: impropers
integer , save :: num_impropers = 0      
!rotations                                                         
type(symmetryOperation), dimension(maxOperations), save :: allSymOps
integer , save :: num_SymOps = 0


! this array will contain only all the unique symmetry operations
!type(symmetryOperation), dimension(:), allocatable, save :: unique_ops

public :: num_rotations,num_reflections,num_impropers       
public :: rotations,reflections,impropers
                                                              
                                                              
!enum to distinguish between rotor types
enum, bind(C)
    enumerator :: LINEAR_TOP, SPHERICAL_TOP, SYMMETRIC_TOP, ASYMMETRIC_TOP
end enum

!enum to distinguish arrangements of symmetry equivalent atoms
enum, bind(C)
    enumerator :: SINGLE_ATOM, LINEAR_ARRANGEMENT, POLYGONAL_ARRANGEMENT_REGULAR,&
                  POLYGONAL_ARRANGEMENT_IRREGULAR, POLYHEDRAL_ARRANGEMENT_SYMMETRIC, &
                  POLYHEDRAL_ARRANGEMENT_ASYMMETRIC, POLYHEDRAL_ARRANGEMENT_SPHERIC
end enum


contains 
!-------------------------------------------------------------------------------
subroutine symmetry_init(set_Natom, set_atom_pos, set_atom_identifier)
    implicit none 
    integer, intent(in) :: set_Natom
    real(rp), dimension(set_Natom,3) :: set_atom_pos
    real(rp), dimension(set_Natom)   :: set_atom_identifier
    real(rp), dimension(3) :: v1, v2, axis, zero, tpos, tpos2
    real(rp) :: tmp
    integer :: arrangement
    integer :: i,j,k,l,m
    logical :: op_found
    type(symmetryOperation) :: tmp_sym_op

    Natom = set_Natom 
    zero(1:3) = 0.D0
    
    !allocate space for atom positions
    if(     allocated(atom_pos)) deallocate(atom_pos)
    if(.not.allocated(atom_pos)) allocate  (atom_pos(Natom,3))
    
    !allocate symmetry equivalent atoms (maximum Natom in each dimension)
    if(     allocated(atom_sea)) deallocate(atom_sea)
    if(.not.allocated(atom_sea)) allocate  (atom_sea(Natom,Natom+1))
    
    !allocate interatomic distance matrix (Natom in each dimension)
    if(     allocated(atom_dist)) deallocate(atom_dist)
    if(.not.allocated(atom_dist)) allocate  (atom_dist (Natom,Natom))
    
    !allocate atomic identifiers (Natom)
    if(     allocated(atom_identifier))  deallocate(atom_identifier )
    if(.not.allocated(atom_identifier )) allocate  (atom_identifier(Natom))
    
    atom_pos        = set_atom_pos
    atom_identifier = set_atom_identifier
    
    !write atom coordinates
    if(verbose) then
        write(*,'(A)') "Coordinates"
        do i = 1,Natom
            if(verbose) write(*,'(3F12.6)') atom_pos(i,:)
        end do
        write(*,'(A)')
    end if
    
    !calculate center of mass
    atom_com = centerOfMass(atom_pos)
    
    !shift atom positions such that they are in the centre of mass
    if(verbose) write(*,'(A)')  "Shifted coordinates (centre of mass in origin)"
    do i = 1,Natom
        atom_pos(i,:) = atom_pos(i,:) - atom_com
        if(verbose) write(*,'(3F12.6)') atom_pos(i,:)
    end do
    if(verbose) write(*,'(A)') 
    
    !calculate the inertia tensor
    atom_inertia = inertiaTensor(atom_pos)
    if(verbose) then
        write(*,'(A)') "Inertia tensor"
        do i = 1,3
             write(*,'(3F12.6)') atom_inertia(i,:)
        end do
        write(*,'(A)')
    end if
    
    !calculate the principal moments of inertia
    call calc_eigenvalues(atom_inertia,atom_Ia,atom_Ib,atom_Ic)
    if(verbose) then
        write(*,'(A)') "Principal moments of inertia"
        write(*,'(A7,F12.6)') "Ia ",atom_Ia
        write(*,'(A7,F12.6)') "Ib ",atom_Ib
        write(*,'(A7,F12.6)') "Ic ",atom_Ic
        write(*,'(A)')
    end if
    
    !determine the rotortype based on the principal moments of inertia
    atom_rotor_type = rotorType(atom_Ia,atom_Ib,atom_Ic)
    if(verbose) then
        if      (atom_rotor_type == LINEAR_TOP)        then
            write(*,'(A)') "Rotortype: Linear"
            write(*,'(A)') "Possible point groups: Cinfv, Dinfh"
        else if (atom_rotor_type == SPHERICAL_TOP) then
            write(*,'(A)') "Rotortype: Spherical top" 
            write(*,'(A)') "Possible point groups: Td, Th, Oh, Ih, K" 
        else if (atom_rotor_type == SYMMETRIC_TOP) then
            write(*,'(A)') "Rotortype: Symmetric top"  
            write(*,'(A)') "Possible point groups: Cn, Cnh, Cnv, Dn, Dnh, Dnd, Sn"            
        else
            write(*,'(A)') "Rotortype: Asymmetric top"
            write(*,'(A)') "Possible point groups: C1, Cs, Ci, Cn, Cnh, Cnv, Dn, Dnh, Dnd, Sn"              
        end if 
        write(*,'(A)') 
    end if
    
    !calculate the interdistance matrix
    atom_dist = interdistanceMatrix(atom_pos)
    if(verbose) then
        write(*,'(A)') "Interdistance matrix"
        do i = 1,Natom
            do j = 1,Natom
                write(*,'(F12.6)', advance='no') atom_dist(i,j)
            end do
            write(*,'(A)')
        end do
        write(*,'(A)')
    end if
    
    !find the symmetrically equivalent atoms
    call find_sea(atom_sea,atom_dist,atom_identifier)
    if(verbose) then
        write(*,'(A)') "Symmetrically equivalent atoms"
        do i = 1,Natom 
            if(count(atom_sea(i,:) /= 0) == 0) exit 
            if(verbose) write(*,'(A,I0,A)', advance='no') "set ",i,": "
            do j = 1,Natom
                if(atom_sea(i,j) == 0) exit
                write(*,'(2X,I0)',advance='no') atom_sea(i,j)   
            end do
            write(*,'(A)')
        end do
        write(*,'(A)')
    end if
    
    !check all sets of symmetry equivalent atoms for their symmetry elements
    if(verbose) write(*,'(A)') "Found symmetry operations"
    do i = 1,Natom
        if(count(atom_sea(i,:) /= 0) == 0) exit ! no sets left
        call find_symmetry_elements_of_sea_set(atom_sea(i,1:minloc(atom_sea(i,:),dim=1)-1), &
                                               arrangement)
        if(verbose) then
            if      (arrangement == SINGLE_ATOM)                       then
                write(*,'(A,I0,A)') "in set ",i,": single atom"
            else if (arrangement == LINEAR_ARRANGEMENT)                then
                write(*,'(A,I0,A)') "in set ",i,": linear arrangement"
            else if (arrangement == POLYGONAL_ARRANGEMENT_REGULAR)     then
                write(*,'(A,I0,A)') "in set ",i,": regular polygonal arrangement"
            else if (arrangement == POLYGONAL_ARRANGEMENT_IRREGULAR)   then
                write(*,'(A,I0,A)') "in set ",i,": irregular polygonal arrangement"
            else if (arrangement == POLYHEDRAL_ARRANGEMENT_SPHERIC)  then
                write(*,'(A,I0,A)') "in set ",i,": spherical arrangement"
            else if (arrangement == POLYHEDRAL_ARRANGEMENT_SYMMETRIC)  then
                write(*,'(A,I0,A)') "in set ",i,": symmetric polyhedral arrangement"
            else if (arrangement == POLYHEDRAL_ARRANGEMENT_ASYMMETRIC) then
                write(*,'(A,I0,A)') "in set ",i,": asymmetric polyhedral arrangement"
            end if
        end if                                     
    end do

    !detect special case of nonlinear triatomic with only Cs symmetry (e.g. SCN-)
    if(Natom.eq.3 .and..not. (equal(atom_dist(1,2),atom_dist(1,3)).or.equal(  &
      atom_dist(1,2),atom_dist(2,3)).or.equal(atom_dist(1,3),atom_dist(2,3))) &
      .and.atom_rotor_type /= LINEAR_TOP)then
      if(verbose) write(*,'(A)') "Cs nonlinear triatomic detected:"
      v1 = atom_pos(2,:)-atom_pos(1,:)
      v2 = atom_pos(3,:)-atom_pos(1,:)
      !check whether the vectors are colinear
      axis(1) = v1(2)*v2(3) - v1(3)*v2(2)
      axis(2) = v1(3)*v2(1) - v1(1)*v2(3)
      axis(3) = v1(1)*v2(2) - v1(2)*v2(1)
      if(.not.equal(sum(axis**2),0d0)) then
        tmp = sqrt(sum(axis**2))
        axis = axis/tmp
        tmp_sym_op%M = reflection(axis)
        tmp_sym_op%typ = "ref"
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = 1
        write(tmp_sym_op%label,'(A)') "sigh"
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
            call add_symmetry_operation_to_list(tmp_sym_op,reflections,num_reflections,.true.)
      end if
    end if
    
    !this is a rare special case, but it is needed to identify symmetry in some molecules like trans-1,2-dichloroethene
    if(verbose) write(*,'(A)') "Searching other symmetry elements:"
    do i = 1,Natom
        if(count(atom_sea(i,:) /= 0) <  2) cycle
        v1 = atom_pos(atom_sea(i,1),:)-atom_pos(atom_sea(i,2),:)
        exit
    end do
    do i = 1,Natom
        if(count(atom_sea(i,:) /= 0) == 0) exit ! no sets left
        if(count(atom_sea(i,:) /= 0) /= 2) cycle ! only if 2 elements in set
        v2 = atom_pos(atom_sea(i,1),:)-atom_pos(atom_sea(i,2),:)
        !check whether the vectors are colinear
        axis(1) = v1(2)*v2(3) - v1(3)*v2(2)
        axis(2) = v1(3)*v2(1) - v1(1)*v2(3)
        axis(3) = v1(1)*v2(2) - v1(2)*v2(1)  
        if(equal(sum(axis**2),0d0)) cycle
        tmp = sqrt(sum(axis**2))
        axis = axis/tmp
        tmp_sym_op%label = "C2"
        tmp_sym_op%typ = "rot"
        tmp_sym_op%M     = rotation(axis,n=2)
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = 1
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
            call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
            !also check the reflection
            tmp_sym_op%M = reflection(axis)
            tmp_sym_op%typ = "ref"
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            write(tmp_sym_op%label,'(A)') "sigh"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,reflections,num_reflections,.true.)
            !also check the impropers
            tmp_sym_op%M = improper(axis, n = 2)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S2"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            tmp_sym_op%M = improper(axis, n = 4)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S4"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
        end if
    end do    
    do i = 1,Natom
        if(count(atom_sea(i,:) /= 0) == 0) exit ! no sets left

    end do
    !build a proper list of symmetry operations, also check inversion
    tmp_sym_op%M = inversion()
    tmp_sym_op%axis  = zero
    tmp_sym_op%n     = 1
    tmp_sym_op%typ = "inv"
    write(tmp_sym_op%label,'(A)') "i"
    if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
!        call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
        call add_symmetry_operation_to_list(tmp_sym_op,allSymOps,num_SymOps,.true.)
    end if
    if(verbose) then
        write(*,'(A)') 
        write(*,'(A)') "Adding unique impropers and reflections to allSymOps list"
    end if
    do i = 1,num_rotations
        call add_symmetry_operation_to_list(rotations(i),allSymOps,num_SymOps,.false.)
    end do
    do i = 1,num_impropers
        call add_symmetry_operation_to_list(impropers(i),allSymOps,num_SymOps,.false.)  
    end do
    do i = 1,num_reflections
        call add_symmetry_operation_to_list(reflections(i),allSymOps,num_SymOps,.false.)  
    end do
    if(verbose) write(*,'(2A)') "Appending nuclear charge and free charge ",&
      "operations for fitting purposes"
    !free charge
    tmp_sym_op%M = identity() !not really, will handle separately...
    tmp_sym_op%axis  = zero
    tmp_sym_op%n     = 1
    tmp_sym_op%typ = "fre"
    write(tmp_sym_op%label,'(A)') "free"
    call add_symmetry_operation_to_list(tmp_sym_op,allSymOps,num_SymOps,.true.)
    !nuclear charge
    tmp_sym_op%M = identity() !for our purposes (atom fitting) the identity matrix
    tmp_sym_op%axis  = zero
    tmp_sym_op%n     = 1
    tmp_sym_op%typ = "nuc"
    write(tmp_sym_op%label,'(A)') "nuc"
    call add_symmetry_operation_to_list(tmp_sym_op,allSymOps,num_SymOps,.true.)
    if(verbose) write(*,'(A)') 
    if(verbose) write(*,'(A,I0,A)') "A total of ",num_SymOps," unique symmetry operations was found (no redundancy)."
!    if(.not.allocated(unique_ops)) allocate(unique_ops(num_SymOps))
    do i = 1,num_SymOps
        if(verbose) write(*,'(I0,A)') i, " "//allSymOps(i)%label
    end do
    
    !finally, determine symmetry ops required to transform the first sea in each group to each of the others
    allocate(sea_ops(Natom,4))
    do i = 1,Natom ! loop over sea sets
      if(count(atom_sea(i,:) /= 0) == 0) exit
      do j=2,Natom ! loop over sea's in set
        if(atom_sea(i,j) == 0) exit
        op_found=.false.
        do k=1,num_SymOps
          tpos(:)=atom_pos(atom_sea(i,1),:)
          if(allSymOps(k)%typ.eq.'rot'.or.allSymOps(k)%typ.eq.'imp')then
            do l=1,allSymOps(k)%n
              tpos=matmul(tpos(:),allSymOps(k)%M)
              if(equal(sqrt(sum((atom_pos(atom_sea(i,j),:)-tpos(:))**2)),0d0))then
                !store sym op to get from first sea in set to jth sea
                sea_ops(atom_sea(i,j),1)=i !set atom is in
                sea_ops(atom_sea(i,j),2)=k !index of sym op to use
                sea_ops(atom_sea(i,j),3)=0 !index of 2nd sym op to use
                sea_ops(atom_sea(i,j),4)=l !how many times to apply sym op
                op_found=.true.
              endif
            enddo
          else
            tpos=matmul(tpos(:),allSymOps(k)%M)
            if(equal(sqrt(sum((atom_pos(atom_sea(i,j),:)-tpos(:))**2)),0d0))then
              !store sym op to get from first sea in set to jth sea
              sea_ops(atom_sea(i,j),1)=i !set atom is in
              sea_ops(atom_sea(i,j),2)=k !index of sym op to use
              sea_ops(atom_sea(i,j),3)=0 !index of 2nd sym op to use
              sea_ops(atom_sea(i,j),4)=1 !how many times to apply sym op
              op_found=.true.
            endif
          endif
        enddo
        ! if we didn't find any single operation, try a combination of two operations
        if(.not.op_found)then
          do k=1,num_SymOps
            tpos(:)=atom_pos(atom_sea(i,1),:)
            if(allSymOps(k)%typ.eq.'rot'.or.allSymOps(k)%typ.eq.'imp')then
              do l=1,allSymOps(k)%n
                tpos=matmul(tpos(:),allSymOps(k)%M)
                do m=1,num_SymOps
                  tpos2=matmul(tpos(:),allSymOps(m)%M)
                  if(equal(sqrt(sum((atom_pos(atom_sea(i,j),:)-tpos2(:))**2)),0d0))then
                    !store sym op to get from first sea in set to jth sea
                    sea_ops(atom_sea(i,j),1)=i !set atom is in
                    sea_ops(atom_sea(i,j),2)=k !index of sym op to use
                    sea_ops(atom_sea(i,j),3)=m !index of sym op to use
                    sea_ops(atom_sea(i,j),4)=l !how many times to apply sym op
                    op_found=.true.
                  endif
                enddo
              enddo
            else
              tpos=matmul(tpos(:),allSymOps(k)%M)
              do m=1,num_SymOps
                tpos2=matmul(tpos(:),allSymOps(m)%M)
                if(equal(sqrt(sum((atom_pos(atom_sea(i,j),:)-tpos2(:))**2)),0d0))then
                  !store sym op to get from first sea in set to jth sea
                  sea_ops(atom_sea(i,j),1)=i !set atom is in
                  sea_ops(atom_sea(i,j),2)=k !index of sym op to use
                  sea_ops(atom_sea(i,j),3)=m !index of sym op to use
                  sea_ops(atom_sea(i,j),4)=1 !how many times to apply sym op
                  op_found=.true.
                endif
              enddo
            endif
          enddo
        endif
        ! thrown error if no sym op transformation was found for this atom
        if(.not.op_found)then
          write(*,'(2(A,I0))') 'No symmetry transormation found to link atom ',&
            atom_sea(i,1),' with atom ',atom_sea(i,j)
          call throw_error('No transformation found in symmetry_init')
        endif
      enddo
    enddo
flush(6)
end subroutine symmetry_init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! eigenvalues of a real symmetric 3x3 matrix
! Smith, Oliver K. (April 1961), "Eigenvalues of a symmetric 3 × 3 matrix.", 
! Communications of the ACM 4 (4): 168, doi:10.1145/355578.366316
subroutine calc_eigenvalues(A,eig1,eig2,eig3)
    implicit none
    real(rp), parameter                  :: pi = acos(-1d0)
    real(rp), dimension(3,3), parameter  :: I = reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3,3/))
    real(rp), dimension(3,3), intent(in) :: A
    real(rp), dimension(3,3)             :: B
    real(rp), intent(out) :: eig1, eig2, eig3
    real(rp) :: p1, p2, p, q, r, phi
    
    p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2
    if(equal(p1,0d0)) then !A is diagonal, the task is easy
        if      (A(1,1) <= A(2,2) .and. A(2,2) <= A(3,3)) then 
            eig1 = A(1,1)
            eig2 = A(2,2)
            eig3 = A(3,3)
        else if (A(1,1) <= A(3,3) .and. A(3,3) <= A(2,2)) then 
            eig1 = A(1,1)
            eig3 = A(2,2)
            eig2 = A(3,3)
        else if (A(2,2) <= A(1,1) .and. A(1,1) <= A(3,3)) then
            eig2 = A(1,1)
            eig1 = A(2,2)
            eig3 = A(3,3)
        else if (A(2,2) <= A(3,3) .and. A(3,3) <= A(1,1)) then
            eig3 = A(1,1)
            eig1 = A(2,2)
            eig2 = A(3,3)
        else if (A(3,3) <= A(1,1) .and. A(1,1) <= A(2,2)) then
            eig2 = A(1,1)
            eig3 = A(2,2)
            eig1 = A(3,3)
        else if (A(3,3) <= A(2,2) .and. A(2,2) <= A(1,1)) then
            eig3 = A(1,1)
            eig2 = A(2,2)
            eig1 = A(3,3)
        end if
    else
        q  = (A(1,1) + A(2,2) + A(3,3))/3d0 !trace
        p2 = (A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + 2*p1
        p  = sqrt(p2/6d0)
        B  = (1d0/p) * (A - q*I)
        r  = ( B(1,1)*B(2,2)*B(3,3) + B(1,2)*B(2,3)*B(3,1) + B(1,3)*B(2,1)*B(3,2) &
              -B(1,3)*B(2,2)*B(3,1) - B(1,2)*B(2,1)*B(3,3) - B(1,1)*B(2,3)*B(3,2))/2d0 ! determinant/2
        if( r <= -1d0) then
            phi = pi/3d0
        else if (r >= 1d0) then
            phi = 0d0
        else
            phi = acos(r)/3d0
        end if
        eig3 = q + 2*p*cos(phi)
        eig1 = q + 2*p*cos(phi + (2d0*pi/3d0))
        eig2 = 3*q - eig1 - eig3
    end if
    return
end subroutine calc_eigenvalues
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Eigenvectors of a 3x3 matrix, given the eigenvalues (by Cayley-Hamilton theorem)
subroutine calc_eigenvectors(A,eig1,eig2,eig3,vec1,vec2,vec3)
    implicit none
    real(rp), dimension(3,3), parameter   :: I = reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3,3/))
    real(rp), dimension(3,3), intent(in)  :: A
    real(rp), intent(in)                  :: eig1, eig2, eig3
    real(rp), dimension(3),   intent(out) :: vec1, vec2, vec3
    real(rp), dimension(3,3)              :: Am1, Am2, Am3, P
    real(rp) :: tmp
    integer :: j
    
    Am1 = A - eig1*I
    Am2 = A - eig2*I
    Am3 = A - eig3*I
    vec1 = 0d0
    vec2 = 0d0
    vec3 = 0d0
    
    if(.not.(equal(eig1,eig2).or.equal(eig1,eig3).or.equal(eig2,eig3))) then !distinct eigenvalues        
        !first eigenvector
        P = matmul(Am2,Am3)
        do j = 1,3
            tmp = sum(P(:,j)**2)
            if(tmp < tolerance) cycle !vector is close to zero
            vec1 = P(:,j)/sqrt(tmp) !normalized vector
            exit
        end do
        
        !second eigenvector
        P = matmul(Am1,Am3)
        do j = 1,3
            tmp = sum(P(:,j)**2)
            if(tmp < tolerance) cycle !vector is close to zero
            vec2 = P(:,j)/sqrt(tmp) !normalized vector
            exit
        end do 
        
        !third eigenvector
        P = matmul(Am1,Am2)
        do j = 1,3
            tmp = sum(P(:,j)**2)
            if(tmp < tolerance) cycle !vector is close to zero
            vec3 = P(:,j)/sqrt(tmp) !normalized vector
            exit
        end do        
    else !degenerate eigenvalues
        if(equal(eig1,eig2)) then
            !third eigenvector
            P = matmul(Am1,Am2)
            do j = 1,3
                tmp = sum(P(:,j)**2)
                if(tmp < tolerance) cycle !vector is close to zero
                vec3 = P(:,j)/sqrt(tmp) !normalized vector
                exit
            end do
            
            P = Am3
            do j = 1,3
                tmp = sum(P(:,j)**2)
                if(tmp < tolerance) cycle !vector is close to zero
                vec1 = P(:,j)/sqrt(tmp) !normalized vector
                exit
            end do
            
            vec2 = crossProduct(vec1, vec3)
        else if (equal(eig2,eig3)) then    
            !first eigenvector
            P = matmul(Am2,Am3)
            do j = 1,3
                tmp = sum(P(:,j)**2)
                if(tmp < tolerance) cycle !vector is close to zero
                vec1 = P(:,j)/sqrt(tmp) !normalized vector
                exit
            end do
            
            P = Am1
            do j = 1,3
                tmp = sum(P(:,j)**2)
                if(tmp < tolerance) cycle !vector is close to zero
                vec2 = P(:,j)/sqrt(tmp) !normalized vector
                exit
            end do
            
            vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
            vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
            vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)  
        else
            !second eigenvector
            P = matmul(Am1,Am3)
            do j = 1,3
                tmp = sum(P(:,j)**2)
                if(tmp < tolerance) cycle !vector is close to zero
                vec2 = P(:,j)/sqrt(tmp) !normalized vector
                exit
            end do 
            
            P = Am2
            do j = 1,3
                tmp = sum(P(:,j)**2)
                if(tmp < tolerance) cycle !vector is close to zero
                vec1 = P(:,j)/sqrt(tmp) !normalized vector
                exit
            end do
            
            vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
            vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
            vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)  
        end if      
    end if
end subroutine calc_eigenvectors
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! finds all sets of symmetrically equivalent atoms and stores their indices
subroutine find_sea(sea,dist,identifier)
    implicit none
    integer, dimension(:,:), intent(out) :: sea  !output set of SEAs
    real(rp), dimension(:,:), intent(in)  :: dist !interdistance matrix
    real(rp), dimension(:),   intent(in)  :: identifier !"labels" of the atoms
    integer, dimension(size(sea,dim=1))  :: numInSet
    integer :: i,j,numSets
    logical :: foundMatch
    !initialize to 0
    sea = 0 
    numInSet = 0
   
    !the first atom HAS to be in some set (there is at least 1 set)
    sea(1,1) = 1
    numInSet(1) = 1
    numSets = 1
    
    !loop over atoms
    do i = 2,size(dist,dim=1)
        foundMatch = .false.
        !loop over sets
        do j = 1,numSets
            !it is always sufficient to check the first element of a set (symmetry-equivalent atoms must be same atom-type!)
            if(identifier(i) == identifier(sea(j,1))) then !if identifiers are not equal, we don't need to bother
                if(sameElementsInSet(dist(:,i),dist(:,sea(j,1)))) then
                    numInSet(j) = numInSet(j) + 1 !increment element count
                    sea(j,numInSet(j)) = i !store index
                    foundMatch = .true.
                    exit
                end if
            end if
        end do
        !no match was found, so we add the atom to a new set
        if(.not.foundMatch) then
            numSets = numSets + 1
            numInSet(numSets) = 1
            sea(numSets,numInSet(numSets)) = i
        end if
    end do
end subroutine find_sea

!-------------------------------------------------------------------------------
! adds a symmetry operation to a specified list of symmetry operations, but only
! if it is not redundant!
subroutine add_symmetry_operation_to_list(sym_op,list,num,override)
    implicit none
    type(symmetryOperation), intent(in) :: sym_op
    type(symmetryOperation), dimension(maxOperations), intent(out) :: list
    logical :: override !force addition of an operation
    integer, intent(out) :: num !number of operations in list
    integer :: i
    
    if(num == 0) then
        num = num + 1
        list(num) = sym_op
        if(verbose) write(*,'(A)') sym_op%label
        return
    else
        !loop through symmetry operations and check whether they are equal
        do i = 1,num
            if(equal(sum((list(i)%M-sym_op%M)**2)/9d0,0d0).and..not.override) return ! operations are equal, so we return         
        end do
        !symmetry operation is not found in list, so we add it
        num = num + 1
        list(num) = sym_op
        if(verbose) write(*,'(A)') sym_op%label
        return    
    end if
end subroutine add_symmetry_operation_to_list

!-------------------------------------------------------------------------------
!finds symmetry elements for each set of sea separately and checks whether the molecule
!also has these symmetry elements (if not, they are useless)
subroutine find_symmetry_elements_of_sea_set(set,arrangement)
    implicit none
    real(rp), parameter :: pi = acos(-1d0)
    real(rp), parameter :: goldenangle = pi*(3d0-sqrt(5d0)) !for checking C infinity
    integer, dimension(:), intent(in) :: set
    integer, intent(out) :: arrangement
    real(rp), dimension(3,3) :: inertia !atom inertia tensor 
    real(rp), dimension(3) :: com !center of mass
    real(rp), dimension(size(set, dim=1),3) :: pos !shifted coordinates
    real(rp) :: Ia, Ib, Ic, tmp !principal moments of inertia
    
    !helper variables to detect symmetry operation
    real(rp), dimension(3)   :: axis, v1, v2, v3   !an axis for rotation
    
    type(symmetryOperation) :: tmp_sym_op
    
       
    integer :: k !number of elements in set
    integer :: i,j,l
    
    k = size(set, dim = 1)
    !calculate centre of mass
    com = centerOfMass(atom_pos(set,:))
    !positions referenced to the set's center of mass
    do i = 1,k
        pos(i,:) = atom_pos(set(i),:)-com
    end do
        
    !check what arrangement the set has
    if(k > 2) then
        inertia = inertiaTensor(pos)
        call calc_eigenvalues(inertia, Ia, Ib, Ic)
        arrangement = arrangementType(Ia, Ib, Ic)
    else
        if(k == 2) then
            arrangement = LINEAR_ARRANGEMENT
        else
            arrangement = SINGLE_ATOM
        end if
    end if

    !depending on the arrangement, different operations are performed
    if      (arrangement == SINGLE_ATOM)                       then
        if(Natom == 3) then !so we don't miss a Cinfinity axis
            axis = (atom_pos(1,:)-atom_pos(2,:))/sqrt(sum((atom_pos(1,:)-atom_pos(2,:))**2))
            tmp_sym_op%M = rotation(axis, alpha = goldenangle)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 99 ! infinity would break the code...
            tmp_sym_op%typ = "rot"
            write(tmp_sym_op%label,'(A)') "C%"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
        end if
    else if (arrangement == LINEAR_ARRANGEMENT)                then
        !calculate the axis (normalized)
        axis = (pos(1,:)-pos(2,:))/sqrt(sum((pos(1,:)-pos(2,:))**2))
        !check the Cinfinity rotation
        tmp_sym_op%M = rotation(axis, alpha = goldenangle)
        tmp_sym_op%axis  = axis
        tmp_sym_op%typ = "rot"
        tmp_sym_op%n     = 99 ! infinity would break the code...
        write(tmp_sym_op%label,'(A)') "C%"
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M )) &
            call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
        !check the existence of a C2 perpendicular to the Cinfinity axis
        v1 = axis
        do i = 1,3 !construct some non-colinear axis
            if(equal(v1(i),0d0)) v1(i) = 1d0
            v1(i) = -v1(i)
        end do
        axis = crossProduct(axis,v1)
        tmp_sym_op%M = rotation(axis, n = 2)
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = 1
        tmp_sym_op%typ = "rot"
        write(tmp_sym_op%label,'(A)') "C2"
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
            call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
            !also check the impropers
            tmp_sym_op%M = improper(axis, n = 2)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S2"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            tmp_sym_op%M = improper(axis, n = 4)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S4"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
        end if
    
        !check other C2 rotations that are perpendicular to the Cinfinity axis
        do i = 1,Natom !loop over all possible atoms
            if(any(set == i)) cycle !atoms from this set are not interesting
            axis = (atom_pos(i,:)-com)/sqrt(sum((atom_pos(i,:)-com)**2))
            if(.not.equal(sum((atom_pos(set(1),:)-atom_pos(set(2),:))*axis),0d0)) cycle !axis is not perpendicular
            tmp_sym_op%M = rotation(axis, n = 2)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "rot"
            write(tmp_sym_op%label,'(A)') "C2"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
                call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                !also check the impropers
                tmp_sym_op%M = improper(axis, n = 2)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A)') "S2"
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                tmp_sym_op%M = improper(axis, n = 4)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A)') "S4"
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            end if
        end do
    else if (arrangement == POLYGONAL_ARRANGEMENT_REGULAR)     then
        !calculate an axis perpendicular to the plane of atoms in set
        axis = crossProduct(pos(1,:)-pos(2,:),pos(1,:)-pos(3,:))
        !check the Ck rotation
        tmp_sym_op%M  = rotation(axis, n = k)
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = k-1
        tmp_sym_op%typ = "rot"
        write(tmp_sym_op%label,'(A,I0)') "C",k
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M )) then
            call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
            !also check the impropers
            tmp_sym_op%M = improper(axis, n = k)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A,I0)') "S",k
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            tmp_sym_op%M = improper(axis, n = 2*k)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A,I0)') "S",2*k
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
        end if
        !check divisors of k (Ci rotations)
        do i = 2,k/2
            tmp_sym_op%M = rotation(axis, n = i)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = i-1
            tmp_sym_op%typ = "rot"
            write(tmp_sym_op%label,'(A,I0)') "C",i
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
                call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                !also check the impropers
                tmp_sym_op%M = improper(axis, n = i)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A,I0)') "S",i
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                tmp_sym_op%M = improper(axis, n = 2*i)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A,I0)') "S",2*i
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            end if
        end do
        !check the reflection along the principal axis
        tmp_sym_op%M  = reflection(axis)
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = 1
        tmp_sym_op%typ = "ref"
        write(tmp_sym_op%label,'(A)') "sigh"
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
            call add_symmetry_operation_to_list(tmp_sym_op,reflections,num_reflections,.true.)
    else if (arrangement == POLYGONAL_ARRANGEMENT_IRREGULAR)   then
        !calculate an axis perpendicular to the plane of atoms in set
        axis = crossProduct(pos(1,:)-pos(2,:),pos(1,:)-pos(3,:))
        !check divisors of k (Ci rotations) (Ck does not exist in irregular)
        do i = 2,k/2
            tmp_sym_op%M = rotation(axis, n = i)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = i-1
            tmp_sym_op%typ = "rot"
            write(tmp_sym_op%label,'(A,I0)') "C",i
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
                call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                !also check the impropers
                tmp_sym_op%M = improper(axis, n = i)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A,I0)') "S",i
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                tmp_sym_op%M = improper(axis, n = 2*i)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A,I0)') "S",2*i
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            end if
        end do
        !check the reflection along the principal axis
        tmp_sym_op%M  = reflection(axis)
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = 1
        tmp_sym_op%typ = "ref"
        write(tmp_sym_op%label,'(A)') "sigh"
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
            call add_symmetry_operation_to_list(tmp_sym_op,reflections,num_reflections,.true.)
    else if (arrangement == POLYHEDRAL_ARRANGEMENT_SPHERIC) then
        !the possible rotation axis will always pass through one of the atoms and the center of mass,
        !so we try all atoms
        do i = 1,k
            !catch divide by zero
            if(sum((pos(i,:)-com)**2) == 0.d0) then
              axis(:)=0.d0
            else
              axis = (pos(i,:)-com)/sqrt(sum((pos(i,:)-com)**2))
            endif
            do j = 2,k-1
                tmp_sym_op%M = rotation(axis, n = j)  
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = j-1
                tmp_sym_op%typ = "rot"
                write(tmp_sym_op%label,'(A,I0)') "C",j
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M )) then
                    call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                    !also check the impropers
                    tmp_sym_op%M = improper(axis, n = j)
                    tmp_sym_op%axis  = axis
                    tmp_sym_op%n     = 1
                    tmp_sym_op%typ = "imp"
                    write(tmp_sym_op%label,'(A,I0)') "S",j
                    if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                        call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                    tmp_sym_op%M = improper(axis, n = 2*j)
                    tmp_sym_op%axis  = axis
                    tmp_sym_op%n     = 1
                    tmp_sym_op%typ = "imp"
                    write(tmp_sym_op%label,'(A,I0)') "S",2*j
                    if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                        call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                end if
            end do
            !or it passes through the midpoint of atoms
            do j = i+1,k
                v1 = (pos(i,:)+pos(j,:))/2d0
                if(sum(v1**2) == 0.d0) then
                  axis(:) = 0.d0
                else
                  axis = v1/sqrt(sum(v1**2))  
                endif
                do l = 2,k-1
                    tmp_sym_op%M = rotation(axis, n = l)  
                    tmp_sym_op%axis  = axis
                    tmp_sym_op%n     = l-1
                    tmp_sym_op%typ = "rot"
                    write(tmp_sym_op%label,'(A,I0)') "C",l
                    if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M )) then
                        call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                        !also check the impropers
                        tmp_sym_op%M = improper(axis, n = l)
                        tmp_sym_op%axis  = axis
                        tmp_sym_op%n     = 1
                        tmp_sym_op%typ = "imp"
                        write(tmp_sym_op%label,'(A,I0)') "S",l
                        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                            call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                        tmp_sym_op%M = improper(axis, n = 2*l)
                        tmp_sym_op%axis  = axis
                        tmp_sym_op%n     = 1
                        tmp_sym_op%typ = "imp"
                        write(tmp_sym_op%label,'(A,I0)') "S",2*l
                        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                            call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                    end if
                end do  
            end do
        end do
    else if (arrangement == POLYHEDRAL_ARRANGEMENT_SYMMETRIC)  then
        !first, it is necessary to locate 3 atoms lying in the same plane
        !in order to define the rotation axis
        outerloop: do i = 1,k
            do j = 1,k
                do l = 1,k
                    if(equal(sum(crossProduct(pos(i,:)-pos(j,:),pos(i,:)-pos(l,:))*(pos(j,:)-pos(l,:))), 0d0)) &
                        exit outerloop
                end do
            end do
        end do outerloop
        !calculate an axis perpendicular to the plane of atoms in set
        axis = crossProduct(pos(i,:)-pos(j,:),pos(i,:)-pos(l,:))
        !check divisors of k (Ci rotations) 
        do i = 2,k/2
            tmp_sym_op%M = rotation(axis, n = i)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = i-1
            tmp_sym_op%typ = "rot"
            write(tmp_sym_op%label,'(A,I0)') "C",i
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
                call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                !also check the impropers
                tmp_sym_op%M = improper(axis, n = i)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A,I0)') "S",i
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                tmp_sym_op%M = improper(axis, n = 2*i)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A,I0)') "S",2*i
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            end if
        end do
        !check the existence of C2 axes perpendicular to the principal axis
        do i = 1,k
            do j = 1,k
                if(equal(dot_product(axis,pos(i,:)-pos(j,:)),0d0)) then
                    v1 = (pos(i,:)-pos(j,:))/sqrt(sum((pos(i,:)-pos(j,:))**2))
                    tmp_sym_op%M = rotation(v1, n = 2)
                    tmp_sym_op%axis  = v1
                    tmp_sym_op%n     = 1
                    tmp_sym_op%typ = "rot"
                    write(tmp_sym_op%label,'(A)') "C2"
                    if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
                        call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                        !also check the impropers
                        tmp_sym_op%M = improper(v1, n = 2)
                        tmp_sym_op%axis  = v1
                        tmp_sym_op%n     = 1
                        tmp_sym_op%typ = "imp"
                        write(tmp_sym_op%label,'(A)') "S2"
                        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                            call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                        tmp_sym_op%M = improper(v1, n = 4)
                        tmp_sym_op%axis  = v1
                        tmp_sym_op%n     = 1
                        tmp_sym_op%typ = "imp"
                        write(tmp_sym_op%label,'(A)') "S4"
                        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                            call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                    end if
                end if
            end do
        end do       
    else if (arrangement == POLYHEDRAL_ARRANGEMENT_ASYMMETRIC) then
        call calc_eigenvectors(inertia,Ia,Ib,Ic,v1,v2,v3)
        write(tmp_sym_op%label,'(A)') "C2"
        axis = v1
        tmp_sym_op%M  = rotation(axis, n = 2)
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = 1
        tmp_sym_op%typ = "rot"
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
            call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
            !also check the impropers
            tmp_sym_op%M = improper(axis, n = 2)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S2"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            tmp_sym_op%M = improper(axis, n = 4)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S4"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
        end if
        axis = v2
        tmp_sym_op%M = rotation(axis, n = 2)
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = 1
        tmp_sym_op%typ = "rot"
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
            call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
            !also check the impropers
            tmp_sym_op%M = improper(axis, n = 2)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S2"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            tmp_sym_op%M = improper(axis, n = 4)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S4"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
        end if
        axis = v3
        tmp_sym_op%M = rotation(axis, n = 2)
        tmp_sym_op%axis  = axis
        tmp_sym_op%n     = 1
        tmp_sym_op%typ = "rot"
        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
            call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
            !also check the impropers
            tmp_sym_op%M = improper(axis, n = 2)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S2"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            tmp_sym_op%M = improper(axis, n = 4)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "imp"
            write(tmp_sym_op%label,'(A)') "S4"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
        end if
    end if
    
    !a non-linear molecule with not more than 2 SEAs can have more C2 axes, which need
    !to be located differently
    if(atom_rotor_type /= LINEAR_TOP) then
        write(tmp_sym_op%label,'(A)') "C2"
        !A)
        do i = 1,k
            do j = i+1,k
                v1 = (pos(i,:)-pos(j,:))/2d0
                tmp = sqrt(sum(v1**2))
                if(.not.equal(tmp,0d0)) then !v1 is not in the center of mass
                    axis = v1/tmp
                    tmp_sym_op%M = rotation(axis, n = 2)
                    tmp_sym_op%axis  = axis
                    tmp_sym_op%n     = 1
                    tmp_sym_op%typ = "rot"
                    if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
                        call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                        !also check the impropers
                        tmp_sym_op%M = improper(axis, n = 2)
                        tmp_sym_op%axis  = axis
                        tmp_sym_op%n     = 1
                        tmp_sym_op%typ = "imp"
                        write(tmp_sym_op%label,'(A)') "S2"
                        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                            call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                        tmp_sym_op%M = improper(axis, n = 4)
                        tmp_sym_op%axis  = axis
                        tmp_sym_op%n     = 1
                        tmp_sym_op%typ = "imp"
                        write(tmp_sym_op%label,'(A)') "S4"
                        if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                            call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                    end if
                end if
            end do
        end do
        !B)
        do i = 1,k
            !catch divide by zero
            if(sum(pos(i,:)**2) == 0) then
              axis(:)=0.d0
            else
              axis = pos(i,:)/sqrt(sum(pos(i,:)**2))
            endif
            tmp_sym_op%M = rotation(axis, n = 2)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "rot"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
                call add_symmetry_operation_to_list(tmp_sym_op,rotations,num_rotations,.true.)
                !also check the impropers
                tmp_sym_op%M = improper(axis, n = 2)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A)') "S2"
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
                tmp_sym_op%M = improper(axis, n = 4)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "imp"
                write(tmp_sym_op%label,'(A)') "S4"
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,impropers,num_impropers,.true.)
            end if
        end do
    end if
    
    
    
    !look for possible reflections
    write(tmp_sym_op%label,'(A)') "sigv"
    do i = 1,k
        do j = i+1,k   
            axis = (pos(i,:)-pos(j,:))/sqrt(sum((pos(i,:)-pos(j,:))**2))
            tmp_sym_op%M = reflection(axis)
            tmp_sym_op%axis  = axis
            tmp_sym_op%n     = 1
            tmp_sym_op%typ = "ref"
            if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) then
                call add_symmetry_operation_to_list(tmp_sym_op,reflections,num_reflections,.true.)
            endif
        end do
    end do
    
    !look for special reflections
    if(arrangement == LINEAR_ARRANGEMENT) then
        write(tmp_sym_op%label,'(A)') "sigv"
        do i = 1,Natom
            if (count(atom_sea(i,:) /= 0) == 0) exit !none left
            if (count(atom_sea(i,:) /= 0) == 1) then !another atom_sea is a single atom
                v1 = pos(1,:)-pos(2,:)
                v2 = pos(1,:)-atom_pos(atom_sea(i,1),:)
                axis = crossProduct(v1,v2)
                tmp_sym_op%M = reflection(axis)
                tmp_sym_op%axis  = axis
                tmp_sym_op%n     = 1
                tmp_sym_op%typ = "ref"
                if(isSymmetryOperation(atom_pos,atom_identifier,tmp_sym_op%M)) &
                    call add_symmetry_operation_to_list(tmp_sym_op,reflections,num_reflections,.true.)
            end if
        end do    
    end if
    
    

end subroutine find_symmetry_elements_of_sea_set
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!checks whether two sets contain the same elements (this is needed to get symmetry equivalent atoms from distance matrix)
logical function sameElementsInSet(s1,s2)
    implicit none
    real(rp), dimension(:), intent(in)  :: s1,s2
    logical, dimension(size(s2,dim=1)) :: alreadyMatched !keeps track of whether an element of s2 is already matched
    logical :: foundMatch
    integer :: i,j
    
    !if the sets are not of the same size, they cannot contain the same elements
    if(size(s1,dim=1) /= size(s2,dim=1)) then
        sameElementsInSet = .false.
        return
    end if
    
    alreadyMatched = .false.
    do i = 1,size(s1,dim=1)
        foundMatch = .false.
        do j = 1,size(s2,dim=1)
            if(alreadyMatched(j)) cycle !only consider values not matched previously
            if(equal(s1(i),s2(j))) then
               alreadyMatched(j) = .true.
               foundMatch = .true.
               exit 
            end if
        end do
        !no matching element was found
        if(.not.foundMatch) then
            sameElementsInSet = .false.
            return
        end if
    end do
    !all elements could be matched
    sameElementsInSet = .true.
end function sameElementsInSet
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! determines whether 2 reals are equal (up to a certain tolerance threshold)
logical function equal(r1,r2)
    implicit none
    real(rp), intent(in) :: r1,r2
    if(abs(r1-r2) <= tolerance) then
        equal = .true.
    elseif(abs(r1).gt.10.d0.and.abs((r1-r2)/r1) <= tolerance_b) then
        equal = .true.
    else
        equal = .false.
    end if
end function equal
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! determines whether a given operation op is a symmetry operation for a set of
! positions and identifiers (atom types)
logical function isSymmetryOperation(pos,ident,op)
    implicit none
    real(rp), dimension(:,:),               intent(in)   :: pos
    real(rp), dimension(size(pos,dim=1)),   intent(in)   :: ident
    real(rp), dimension(3,3),               intent(in)   :: op
    real(rp), dimension(size(pos,dim=1),size(pos,dim=2)) :: posXop
    logical :: foundMatch
    integer :: i,j
    posXop = matmul(pos,op)
    do i = 1,size(pos,dim=1) !loop over positions
        foundMatch = .false.
        do j = 1,size(posXop,dim=1) !loop over transformed positions
            !overlapping positions were found
            if(equal(sqrt(sum((pos(i,:)-posXop(j,:))**2)),0d0)) then
                if(equal(ident(i),ident(j))) then !the identifiers also match
                    foundMatch = .true.
                    exit
                else !identifiers do not match
                    exit ! (foundMatch is false already)
                end if
            end if
        end do
        if(.not.foundMatch) then !no match was found => no symmetry operation
            isSymmetryOperation = .false.
            return
        end if
    end do
    !all checks were passed
    isSymmetryOperation = .true.
end function isSymmetryOperation
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! returns the magnitude of assymetry for a given symmetry operation op
real(rp) function asymmetryMagnitude(pos,ident,op)
    implicit none
    real(rp), dimension(:,:),               intent(in)   :: pos
    real(rp), dimension(size(pos,dim=1)),   intent(in)   :: ident
    real(rp), dimension(3,3),               intent(in)   :: op
    real(rp), dimension(size(pos,dim=1),size(pos,dim=2)) :: posXop
    real(rp) :: diff
    integer :: i,j   
    asymmetryMagnitude = 0d0
    posXop = matmul(pos,op)
    do i = 1,size(pos,dim=1) !loop over positions
        do j = 1,size(posXop,dim=1) !loop over transformed positions
            diff = sqrt(sum((pos(i,:)-posXop(j,:))**2))
            if(.not.equal(diff,0d0)) then
                asymmetryMagnitude = asymmetryMagnitude + (diff-tolerance)
            end if
        end do
    end do
end function asymmetryMagnitude
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! determines the rotor type of the molecule based on principal moments of inertia
integer function rotorType(Ia,Ib,Ic)
    implicit none
    real(rp), intent(in) :: Ia,Ib,Ic
    !determine the rotortype based on the principal moments of inertia
    if(      (equal(Ia,0d0).and.equal(Ib,Ic)).or.&
             (equal(Ib,0d0).and.equal(Ia,Ic)).or.&
             (equal(Ic,0d0).and.equal(Ia,Ib))) then
        rotorType = LINEAR_TOP
    else if ( equal(Ia,Ib).and.equal(Ib,Ic)) then
        rotorType = SPHERICAL_TOP
    else if ((equal(Ia,Ib).and..not.equal(Ib,Ic)).or.&
             (equal(Ia,Ic).and..not.equal(Ic,Ib)).or.&
             (equal(Ib,Ic).and..not.equal(Ic,Ia))) then
        rotorType = SYMMETRIC_TOP         
    else
        rotorType = ASYMMETRIC_TOP      
    end if    
end function rotorType
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! determines the type of arangement given the principal moments of inertia of a set
! of symmetry equivalent atoms
integer function arrangementType(Ia,Ib,Ic)
    implicit none
    real(rp), intent(in) :: Ia,Ib,Ic
    if      ( equal(Ia,0d0).and.equal(Ib,0d0).and.equal(Ic,0d0) ) then
        arrangementType = SINGLE_ATOM
    else if ((equal(Ia,0d0).and.equal(Ib,Ic)).or.&
             (equal(Ib,0d0).and.equal(Ia,Ic)).or.&
             (equal(Ic,0d0).and.equal(Ia,Ib))) then
        arrangementType = LINEAR_ARRANGEMENT
    else if ( equal(Ia+Ib,Ic).or.equal(Ia+Ic,Ib).or.equal(Ib+Ic,Ia)) then
        if  ( equal(Ia,Ib).or.equal(Ia,Ic).or.equal(Ib,Ic)) then
            arrangementType = POLYGONAL_ARRANGEMENT_REGULAR
        else
            arrangementType = POLYGONAL_ARRANGEMENT_IRREGULAR
        end if
    else
        if(equal(Ia,Ib).or.equal(Ia,Ic).or.equal(Ib,Ic)) then
            if(equal(Ia,Ib).and.equal(Ib,Ic)) then
                arrangementType = POLYHEDRAL_ARRANGEMENT_SPHERIC 
            else
                arrangementType = POLYHEDRAL_ARRANGEMENT_SYMMETRIC
            end if
        else
            arrangementType = POLYHEDRAL_ARRANGEMENT_ASYMMETRIC
        end if
    end if    
end function arrangementType
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!returns the center of mass of a set of coordinates
function centerOfMass(pos)
    implicit none
    real(rp), dimension(:,:), intent(in) :: pos !dimension (?,3)
    real(rp), dimension(3) :: centerOfMass
    integer :: i
    centerOfMass = 0d0
    do i = 1,size(pos,dim=1)
        centerOfMass = centerOfMass + pos(i,:)
    end do
    centerOfMass = centerOfMass/real(size(pos,dim=1),8)
end function centerOfMass
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!returns the inertia tensor of a set of coordinates
function inertiaTensor(pos)
    implicit none
    real(rp), dimension(:,:), intent(in)  :: pos
    real(rp), dimension(3,3) :: inertiaTensor
    !Ixx
    inertiaTensor(1,1) = sum(pos(:,2)**2 + pos(:,3)**2)
    !Ixy
    inertiaTensor(1,2) = sum(-pos(:,1)*pos(:,2))
    !Ixz
    inertiaTensor(1,3) = sum(-pos(:,1)*pos(:,3))
    !Iyx
    inertiaTensor(2,1) = inertiaTensor(1,2)
    !Iyy
    inertiaTensor(2,2) = sum(pos(:,1)**2 + pos(:,3)**2)
    !Iyz
    inertiaTensor(2,3) = sum(-pos(:,2)*pos(:,3))
    !Izx
    inertiaTensor(3,1) = inertiaTensor(1,3)
    !Izy
    inertiaTensor(3,2) = inertiaTensor(2,3)
    !Izz
    inertiaTensor(3,3) = sum(pos(:,1)**2 + pos(:,2)**2)
end function inertiaTensor
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!returns the normalized crossProduct of two vectors v1 and v2
function crossProduct(v1,v2)
    implicit none
    real(rp), dimension(3), intent(in) :: v1, v2
    real(rp), dimension(3) :: crossProduct
    crossProduct(1) = v1(2)*v2(3) - v1(3)*v2(2)
    crossProduct(2) = v1(3)*v2(1) - v1(1)*v2(3)
    crossProduct(3) = v1(1)*v2(2) - v1(2)*v2(1)  
    crossProduct = crossProduct/sqrt(sum(crossProduct**2)) !normalize
end function crossProduct
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! calculates the interdistance matrix of a set of coordinates
function interdistanceMatrix(pos)
    implicit none
    real(rp), dimension(:,:), intent(in)  :: pos
    real(rp), dimension(size(pos,dim=1),size(pos,dim=1)) :: interdistanceMatrix
    integer :: i,j,dims
    interdistanceMatrix = 0d0
    dims = size(pos,dim=1)
    do i = 1,dims
        do j = i+1,dims
            interdistanceMatrix(i,j) = sqrt(sum((pos(i,:)-pos(j,:))**2))
            interdistanceMatrix(j,i) = interdistanceMatrix(i,j)
        end do
    end do
end function interdistanceMatrix
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! calculates the interdistance matrix of a set of coordinates
function getShiftedPos(atom)
    implicit none
    real(rp), dimension(3) :: getShiftedPos
    integer :: atom

    getShiftedPos = atom_pos(atom,:)
end function getShiftedPos

!-------------------------------------------------------------------------------
! find and store the symmetry ops that affect off-center charges for a given
! atom (those that don't affect nuclear positions). These will potenitally be
! needed for atomic charge fitting.
function getChargeOps(atom)
    implicit none
    integer i,atom,getChargeOps
    logical toadd

    real(rp), dimension(3) :: posXop, shift_pos    ! atomic positions shifted to CoM and transformed by symmetry ops

    if(.not.allocated(atm_sym_ops)) allocate(atm_sym_ops(Natom,num_SymOps))
    if(.not.allocated(num_atm_sym_ops)) allocate(num_atm_sym_ops(Natom))
    num_atm_sym_ops(atom)=0
    do i=1,num_SymOps
      toadd=.false.
      shift_pos = getShiftedPos(atom)
      posXop = matmul(shift_pos(:),allSymOps(i)%M) ! multiply atom position by sym op

      if(equal(sqrt(sum((shift_pos(:)-posXop(:))**2)),0d0)) then
        toadd=.true.
!        ! also remove redundant operations (in terms of rotational axes and planes used for fitting)
!        do j=1,num_atm_sym_ops(atom)
!          op=atm_sym_ops(atom,j)
!          if((rotations(op)%label.eq.'rot'.or.rotations(op)%label.eq.'imp').and. &
!             (rotations(i)%label.eq.'rot'.or.rotations(i)%label.eq.'imp')) then
!            if(equal(sqrt(sum((rotations(op)%axis(:)-rotations(i)%axis(:))**2)) &
!              ,0d0)) then
!              toadd=.false.
!            endif
!          elseif(rotations(op)%label.eq.'ref'.and.rotations(i)%label.eq.'ref')then
!            if(equal(sqrt(sum((rotations(op)%axis(:)-rotations(i)%axis(:))**2)) &
!              ,0d0)) then
!              toadd=.false.
!            endif
!          endif
!        enddo
      endif

      if(toadd) then
        ! include this sym op for charges of this atom
        if(verbose) then
          write(*,'(A,I0,3A,I0)')'include operation ',i,'(',allSymOps(i)%label, &
               ') to fit charges for atom ',atom
        endif
        num_atm_sym_ops(atom)=num_atm_sym_ops(atom)+1
        atm_sym_ops(atom,num_atm_sym_ops(atom))=i
      else
        if(verbose) then
          write(*,'(A,I0,3A,I0)')'reject operation ',i,'(',allSymOps(i)%label, &
               ') while fitting charges for atom ',atom
        endif
      endif
    enddo
    getChargeOps=num_atm_sym_ops(atom)

end function getChargeOps


!-------------------------------------------------------------------------------
! find out how many new charges are spawned by applying consecutive sim ops to 
! an initial charge with position "pos"
function chgsSpawned(atom,pos)
    implicit none
    integer :: atom,chgsSpawned
    integer :: i,j,k,l,op,num_trans_chg
    integer, parameter :: MAX_CHG=9999
    real(rp), dimension(3) :: pos,tpos
    real(rp), dimension(MAX_CHG,3) :: trans_pos ! no idea how many charges we'll have in advance...
    logical :: duplicate

    num_trans_chg=1
    trans_pos(1,:)=pos(:)

    ! apply sim ops to charge and add spawned charge(s) if unique
    do i=1,num_atm_sym_ops(atom)
      op=atm_sym_ops(atom,i)
      do j=1,num_trans_chg
        tpos(:)=trans_pos(j,:)
        do k=1,allSymOps(op)%n
          tpos=matmul(tpos(:),allSymOps(op)%M)
          duplicate=.false.
          do l=1,num_trans_chg
            if(equal(sqrt(sum((trans_pos(l,:)-tpos(:))**2)),0d0))then
              duplicate=.true.
            endif
          enddo !l
          if(.not. duplicate) then
            num_trans_chg=num_trans_chg+1
            if(num_trans_chg.gt.MAX_CHG)then
              call throw_error('ERROR: Too many charges in chgsSpawned()')
            endif
            trans_pos(num_trans_chg,:)=tpos(:)
          endif
        enddo !k
      enddo !j
    enddo !i

!    write(*,'(A,3F6.3)') ' H ',atom_pos(atom,1:3)
!    do i=1,num_trans_chg
!      write(*,'(A,3F6.3)') ' X ',trans_pos(i,1:3)
!    enddo

    chgsSpawned=num_trans_chg

end function chgsSpawned

!-------------------------------------------------------------------------------
! convert a coordinate along a rotation axis relative to
! a nuclear position to global Cartesian coordinates
function rotax_to_cartesian(atom,op,r)
    implicit none
    integer :: atom,op ! op is index in array "allSymOps"
    real(rp), dimension(3) :: rotax_to_cartesian
    real(rp) :: r

    if(allSymOps(op)%typ.ne.'rot'.and.allSymOps(op)%typ.ne.'imp')then
      write(*,'(2A)') 'ERROR, rotax_to_cartesian called for non-rotation operation ',&
         allSymOps(op)%typ
      stop
    endif

    rotax_to_cartesian=atom_pos(atom,:)+r*(allSymOps(op)%axis)

!    write(*,'(A,3F6.3)') ' X ',atom_pos(atom,1:3)
!    write(*,'(A,3F6.3)') ' H ',rotax_to_cartesian(1:3)
end function rotax_to_cartesian

!-------------------------------------------------------------------------------
! convert a coordinate in a mirror plane relative to
! a nuclear position to global Cartesian coordinates
function refplane_to_cartesian(atom,op,rx,ry)
    implicit none
    integer :: i,j,k,atom,op ! op is index in array "allSymOps"
    real(rp), dimension(3) :: xax,yax,refplane_to_cartesian
    real(rp) :: rx,ry
    real(rp) :: small=1.D-6

    if(allSymOps(op)%typ.ne.'ref')then
      call throw_error('ERROR, refplane_to_cartesian called for non-reflection operation!')
    endif

    ! define x-axis as arbitrary vector in mirror plane
    do i=1,3
      if(allSymOps(op)%axis(i).ne.0.d0)then
        do j=1,3
          if(j.ne.i)then
            xax(j)=allSymOps(op)%axis(i)
            xax(i)=-allSymOps(op)%axis(j)
            exit
          endif
        enddo
        exit
      endif
    enddo
    do k=1,3
      if(k.ne.i .and. k.ne.j)then
        if(abs(allSymOps(op)%axis(k)).lt.small)then
          xax(k)=0.d0
        else
          xax(k)=(-allSymOps(op)%axis(i)*xax(i)-allSymOps(op)%axis(j)*xax(j)) &
               /allSymOps(op)%axis(k)
        endif
      endif
    enddo
    ! normalize
    xax=xax/sqrt(sum(xax**2))

    ! define y-axis
    yax=crossProduct(allSymOps(op)%axis,xax)

    ! first along x-axis in mirror plane
    refplane_to_cartesian=atom_pos(atom,:)+rx*(xax)
    ! now add y-axis contribution
    refplane_to_cartesian=refplane_to_cartesian+ry*(yax)

!    write(*,'(A,3F6.3)') ' X ',atom_pos(atom,1:3)
!    write(*,'(A,3F6.3)') ' H ',refplane_to_cartesian(1:3)
end function refplane_to_cartesian

!-------------------------------------------------------------------------------
! convert a random point relative to a nuclear position to global Cartesian coordinates
function point_to_cartesian(atom,rx,ry,rz)
    implicit none
    integer :: atom
    real(rp), dimension(3) :: point_to_cartesian
    real(rp) :: rx,ry,rz

    point_to_cartesian(1) = atom_pos(atom,1)+rx
    point_to_cartesian(2) = atom_pos(atom,2)+ry
    point_to_cartesian(3) = atom_pos(atom,3)+rz

end function point_to_cartesian

!-------------------------------------------------------------------------------
! transform charge coordinates in a fitting solution to another symmetry
! equivalent atom
subroutine sym_map_sea_q_coords(sol,mapsol,atm1,atm2,num_charges)
    implicit none
    integer :: atm1,atm2,num_charges
    real(rp), dimension(:) :: sol,mapsol
    integer :: i,j
    real(rp), dimension(3) :: tpos

    if(atom_sea(sea_ops(atm2,1),1).ne.atm1) call throw_error('sea mismatch in sym_map_sea_q_coords')
    do i=1,num_charges*4-3,4
      tpos=sol(i:i+2)-atom_com !translate back to com coords
      do j=1,sea_ops(atm2,4)
        tpos=matmul(tpos(:),allSymOps(sea_ops(atm2,2))%M) !transform charge coords
      enddo
      if(sea_ops(atm2,3).ne.0)then
        tpos=matmul(tpos(:),allSymOps(sea_ops(atm2,3))%M) !transform charge coords
      endif
      tpos=tpos+atom_com !back to global coords
      mapsol(i:i+2)=tpos
      if(i<num_charges*4-3) mapsol(i+3)=sol(i+3)
    enddo

end subroutine sym_map_sea_q_coords

!-------------------------------------------------------------------------------
! return mapping index of atom operation to full operations array 'allSymOps'
function get_atm_sym_op(atom,i)
    implicit none
    integer :: atom,i,get_atm_sym_op

    get_atm_sym_op=atm_sym_ops(atom,i)

end function get_atm_sym_op

!-------------------------------------------------------------------------------
! check how many charges are spawned by each of the sym ops of each atom
! by new charges placed along rotational axes or in mirror planes
subroutine get_chgs_spawned(a)
    implicit none
    integer :: a
    integer :: tnum_spawned,op,i,j
    real(rp), dimension(3) :: tpos
    real(rp) :: r,rx,ry,rz
    
    if(.not.allocated(num_spawned)) allocate(num_spawned(Natom,num_SymOps))
    do i=1,num_atm_sym_ops(a)
      num_spawned(a,i)=0
      op=get_atm_sym_op(a,i)
      !ToDo, rather than sample randomly N times we should find a single
      !point that lies along the desired axis or in the desired plane and
      !isn't close to an intersection with another axis or plane (handling
      !axes that lie in planes etc.)
      do j=1,50
        if(allSymOps(op)%typ.eq.'rot'.or.allSymOps(op)%typ.eq.'imp') then
          call random_number(r)
          tpos=rotax_to_cartesian(a,op,r*5.D0-2.5D0)
          tnum_spawned=chgsSpawned(a,tpos)
          if(tnum_spawned.gt.num_spawned(a,i))then
            num_spawned(a,i)=tnum_spawned
          endif
        elseif(allSymOps(op)%typ.eq.'ref') then
          call random_number(rx)
          call random_number(ry)
          tpos=refplane_to_cartesian(a,op,rx*5.D0-2.5D0,ry*5.D0-2.5D0)
          tnum_spawned=chgsSpawned(a,tpos)
          if(tnum_spawned.gt.num_spawned(a,i))then
            num_spawned(a,i)=tnum_spawned
          endif
        elseif(allSymOps(op)%typ.eq.'inv') then
          num_spawned(a,i)=1
        elseif(allSymOps(op)%typ.eq.'nuc') then
          num_spawned(a,i)=1
        elseif(allSymOps(op)%typ.eq.'fre') then ! random point not along any rotation axes or in any mirror planes
          call random_number(rx)
          call random_number(ry)
          call random_number(rz)
          tpos=point_to_cartesian(a,rx*5.D0-2.5D0,ry*5.D0-2.5D0,rz*5.D0-2.5D0)
          tnum_spawned=chgsSpawned(a,tpos)
          if(tnum_spawned.gt.num_spawned(a,i))then
            num_spawned(a,i)=tnum_spawned
          endif
        else
          write(*,'(2A)') 'ERROR in get_chgs_spawned: unknown symmetry type: ',&
             allSymOps(op)%typ
          stop
        endif
      enddo !j
      write(*,'(I0,2A)') num_spawned(a,i), &
          ' charges spawned by charges in axis or in plane of operation ', &
          allSymOps(op)%label
    enddo !i

end subroutine get_chgs_spawned

!-------------------------------------------------------------------------------
! selects symmetry axes and mirror planes to place charges for subsequent
! symmetry-constrained atom fitting, using the total allowed number of charges for
! the atom as a constraint
function init_atm_sym_search(num_charges,a)
    implicit none
    integer :: num_charges,a
    integer, parameter :: maxtry=200
    integer :: i,cur_charges
    real(rp) :: ran
    logical :: init_atm_sym_search
    integer, dimension(num_charges) :: num_chgs_per_fit_op ! temporary array

    init_atm_sym_search=.true.

    if(allocated(atm_fit_ops)) deallocate(atm_fit_ops)
    if(allocated(num_atm_fit_ops)) deallocate(num_atm_fit_ops)
    allocate(atm_fit_ops(Natom,num_charges))
    allocate(num_atm_fit_ops(Natom))

    cur_charges=0
    num_atm_fit_ops(a)=0
    do i=1,maxtry
      ! randomly select symmetry operations that satisfy requested total charges for atom
      call random_number(ran)
      ran=ran*(num_atm_sym_ops(a)-1) ! -1 to exclude nuclear charge option for now
      if(num_spawned(a,int(ran)+1).gt.0.and. &
         cur_charges+num_spawned(a,int(ran)+1).le.num_charges)then
        num_atm_fit_ops(a)=num_atm_fit_ops(a)+1
        atm_fit_ops(a,num_atm_fit_ops(a))=atm_sym_ops(a,int(ran)+1)
        cur_charges=cur_charges+num_spawned(a,int(ran)+1)
        num_chgs_per_fit_op(num_atm_fit_ops(a))=num_spawned(a,int(ran)+1)
      endif
      if(cur_charges.eq.num_charges)then !reached the desired number of charges
        exit
      endif
    enddo
    !use the nuclear position if we still need 1 more charge
    if(num_charges-cur_charges.eq.1)then
      num_atm_fit_ops(a)=num_atm_fit_ops(a)+1
      atm_fit_ops(a,num_atm_fit_ops(a))=num_SymOps !nuclear charge stored as last sym op in array
      cur_charges=cur_charges+1
      num_chgs_per_fit_op(num_atm_fit_ops(a))=1
    endif

    !if we found a solution to provide the correct no. of charges then return,
    !otherwise throw error
    if(num_charges.ne.cur_charges)then
      write(*,'(/,A,I0,A,I0,/)') &
          'init_atm_sym_search: no solution found that provides ',&
          num_charges,' charges for atom ',a
      init_atm_sym_search=.false.
    else
      write(*,'(/,2(A,I0),A)') &
        'Fitting sym ops for atom ',a,' with ',num_charges,' charges: '
      do i=1,num_atm_fit_ops(a)
        write(*,'(3A,I0,A,I0,A)') '  ',allSymOps(atm_fit_ops(a,i))%label,&
          ' #',atm_fit_ops(a,i),' (',num_chgs_per_fit_op(i),' charges)'
      enddo
      write(*,'(/)')
    endif

end function init_atm_sym_search


!-------------------------------------------------------------------------------
subroutine init_sym_search_range(search_range,symFitAtms,num_symFitAtms,&
           max_extend,max_charge,sqdim)
    implicit none
    integer :: n,i,j,a,num_symFitAtms,sqdim
    integer, dimension(:) :: symFitAtms
    real(rp) :: max_extend,max_charge
    character(len=3) :: lab
    real(rp), dimension(:,:) :: search_range ! has a minimum and a maximum value for every entry of fitting parameters
    logical :: lasta

    n=1
    lasta=.false. ! is last atom?
    do j=1,num_symFitAtms
      a=symFitAtms(j)
      if(j.eq.num_symFitAtms) lasta=.true. ! last atom
      do i=1,num_atm_fit_ops(a)
        lab=allSymOps(atm_fit_ops(a,i))%typ
        if(lab.eq.'rot'.or.lab.eq.'imp') then
          search_range(1,n  ) = -max_extend
          search_range(2,n  ) = max_extend
          n=n+1
          if(i<num_atm_fit_ops(a).or..not.lasta)then ! constrain last charge
            search_range(1,n+1) = -max_charge
            search_range(2,n+1) = max_charge
            n=n+1
          endif
        elseif(lab.eq.'ref') then
          search_range(1,n  ) = -max_extend
          search_range(2,n  ) = max_extend
          search_range(1,n+1) = -max_extend
          search_range(2,n+1) = max_extend
          n=n+2
          if(i<num_atm_fit_ops(a).or..not.lasta)then ! constrain last charge
            search_range(1,n+2) = -max_charge
            search_range(2,n+2) = max_charge
            n=n+1
          endif
        elseif(lab.eq.'inv') then
          if(i<num_atm_fit_ops(a).or..not.lasta)then ! constrain last charge
            search_range(1,n) = -max_charge
            search_range(2,n) = max_charge
            n=n+1
          endif
        elseif(lab.eq.'fre') then
          search_range(1,n  ) = -max_extend
          search_range(2,n  ) = max_extend
          search_range(1,n+1) = -max_extend
          search_range(2,n+1) = max_extend
          search_range(1,n+2) = -max_extend
          search_range(2,n+2) = max_extend
          n=n+3
          if(i<num_atm_fit_ops(a).or..not.lasta)then ! constrain last charge
            search_range(1,n+3) = -max_charge
            search_range(2,n+3) = max_charge
            n=n+1
          endif
        elseif(lab.eq.'nuc') then
          if(i<num_atm_fit_ops(a).or..not.lasta)then ! constrain last charge
            search_range(1,n  ) = -max_charge
            search_range(2,n  ) = max_charge
            n=n+1
          endif
        else
          write(*,'(/,3A,/)') 'ERROR: unrecognized operation type ',lab,&
             ' in init_sym_search_range'
          stop
        endif
      enddo !i
    enddo !j

    sqdim=n-1

end subroutine init_sym_search_range

!-------------------------------------------------------------------------------
! initializes the population to only feasible solutions
subroutine sym_init_pars(pop,fitAtms,nFitAtms,qatm,scal,radius,atom_num)
    implicit none
    real(rp), dimension(:), intent(out) :: pop
    real(rp), dimension(3) :: ranvec
    real(rp), dimension(:) :: radius, qatm
    real(rp) :: ran, ranx, rany, scal, f
    integer, dimension(:) :: fitAtms, atom_num
    integer  :: i,d,l,op,atm,nFitAtms
    logical :: lasta ! last atom?

    l=1 ! index in pop parameter array
    lasta=.false.

    ! loop over symmetry operations for this atom
    do i=1,nFitAtms
      atm=fitAtms(i)
      if(i.eq.nFitAtms) lasta=.true. !last atom
      do d = 1,num_atm_fit_ops(atm)
        op=atm_fit_ops(atm,d)
        if(allSymOps(op)%typ.eq.'rot'.or.allSymOps(op)%typ.eq.'imp')then
          ! draw a random scaling factor and scale distance along axis
          call random_number(ranx)
          ranx=2.d0*(ranx-0.5d0)
          ranx = ranx*scal*radius(atom_num(atm))
          ! set charge position
          pop(l) = ranx
          l=l+1
          if(d.lt.num_atm_fit_ops(atm).or..not.lasta)then ! constrain last charge
            ! set charge
            call random_number(ran)
            pop(l) = -abs(qatm(atm)) + 2*ran*abs(qatm(atm))
            l=l+1
          endif
        elseif(allSymOps(op)%typ.eq.'ref') then
          ! draw a normalized random scaling vector in mirror plane
          call random_number(ranx)
          ranx=2.d0*(ranx-0.5d0)
          call random_number(rany)
          rany=2.d0*(rany-0.5d0)
          f=sqrt(ranx**2+rany**2)
          ranx=ranx/f
          rany=rany/f
  
          ! draw a random scaling factor and scale vector
          call random_number(ran)
          ranx=ran*scal*radius(atom_num(atm))*ranx
          rany=ran*scal*radius(atom_num(atm))*rany
  
          ! set charge position
          pop(l) = ranx
          pop(l+1) = rany
          l=l+2
  
          if(d.lt.num_atm_fit_ops(atm).or..not.lasta)then ! constrain last charge
            ! set charge
            call random_number(ran)
            pop(l) = -abs(qatm(atm)) + 2*ran*abs(qatm(atm))
            l=l+1
          endif
        elseif(allSymOps(op)%typ.eq.'inv') then
          if(d.lt.num_atm_fit_ops(atm).or..not.lasta)then ! constrain last charge
            call random_number(ran)
            pop(l) = -abs(qatm(atm)) + 2*ran*abs(qatm(atm))
            l=l+1
          endif
        elseif(allSymOps(op)%typ.eq.'nuc')then ! nuclear charge
          if(d.lt.num_atm_fit_ops(atm).or..not.lasta)then ! constrain last charge
            call random_number(ran)
            pop(l) = -abs(qatm(atm)) + 2*ran*abs(qatm(atm))
            l=l+1
          endif
        elseif(allSymOps(op)%typ.eq.'fre')then ! free point
          ! draw a random vector and normalize (random direction)
          call random_number(ranvec)
          ranvec = ranvec - 0.5d0
          ranvec = ranvec/sqrt(sum(ranvec**2))
  
          ! draw a random scaling factor and scale vector
          call random_number(ran)
          ranvec = ran*scal*radius(atom_num(atm))*ranvec
  
          ! set charge position
          pop(l:l+2) = ranvec
          l=l+3
  
          if(d.lt.num_atm_fit_ops(atm).or..not.lasta)then ! constrain last charge
            ! set charge
            call random_number(ran)
            pop(l) = -abs(qatm(atm)) + 2*ran*abs(qatm(atm))
            l=l+1
          endif
        else
          write(*,'(2A)') 'ERROR in sym_init_pars: unknown symmetry type: ',&
             allSymOps(op)%typ
          stop
        endif
      end do
    enddo
end subroutine sym_init_pars

!-------------------------------------------------------------------------------
! initialize atm_fit_ops() and num_atm_fit_ops() for the combination of symmetry
! operations that was found to give the best initial guess for the current total
! number of charges
subroutine sym_init_fit_ops(symFitAtms,num_symFitAtms,best_symatm_combo,num_charges,&
                            sym_solution_ops,num_sym_solution_ops)
    implicit none
    integer, dimension(:) :: symFitAtms,best_symatm_combo
    integer, dimension(:,:) :: num_sym_solution_ops
    integer, dimension(:,:,:) :: sym_solution_ops
    integer :: num_symFitAtms,num_charges,i,atm,nchg

    if(allocated(num_atm_fit_ops)) deallocate(num_atm_fit_ops)
    if(allocated(atm_fit_ops)) deallocate(atm_fit_ops)
    allocate(num_atm_fit_ops(Natom))
    allocate(atm_fit_ops(Natom,num_charges))

    do i=1,num_symFitAtms
      atm=symFitAtms(i)
      nchg=best_symatm_combo(i)
      if(nchg.eq.0)then
        num_atm_fit_ops(atm)=0
      else
        num_atm_fit_ops(atm)=num_sym_solution_ops(nchg,atm)
        atm_fit_ops(atm,1:num_atm_fit_ops(atm))=&
                         sym_solution_ops(nchg,atm,1:num_atm_fit_ops(atm))
      endif
    enddo

end subroutine sym_init_fit_ops


!-------------------------------------------------------------------------------
! applies sym ops to sqin (local sym axis, mirror plane etc. coordinates) to
! populate a new Cartesian array q of all coorindates after applying sym ops
subroutine spawn_sym_chgs(sqin,q,nq,symFitAtms,num_symFitAtms,total_charge,&
               replicate_atoms)
    implicit none
    integer :: atm,atm1,atm2,nq,npts,num_symFitAtms
    integer, dimension(:) :: symFitAtms
    real(rp), dimension(:) :: sqin ! input charges
    integer, dimension(num_symFitAtms) :: num_atm_chgs,n1
    integer, dimension(nq) :: corr_chgs ! hold q() indices of charges to correct in order to maintain total charge
    integer :: a,b,i,j,k,l,n,op,nlast,tlast,nn,t
    logical :: replicate_atoms !whether to replicate charges to other sea's
    real(rp) :: chg,chg_corr,total_charge,local_charge
    real(rp), dimension(size(sqin,dim=1)+1) :: tqin    ! sqin with additional empty charge
    real(rp), dimension(nq*4) :: tq    ! output: spawned charge array in global axis
    real(rp), dimension(:) :: q    ! output: spawned charge array in global axis
    real(rp), dimension(3) :: pt,tpos
    real(rp), dimension(nq,3) :: sym_pts ! hold points generated by sym ops

    tqin(1:size(sqin,dim=1))=sqin(:)
    tqin(size(sqin,dim=1)+1)=0.D0 ! add empty charge to end of fitting parameter array

    n=0 ! no. chgs spawned
    l=1 ! current index in sqin array
    local_charge=0.d0

    do a=1,num_symFitAtms
      atm=symFitAtms(a)
      num_atm_chgs(a)=0
      n1(a)=n*4
      do i=1,num_atm_fit_ops(atm) ! loop over symmetry-constrained charges
        op=atm_fit_ops(atm,i)
        if(allSymOps(op)%typ.eq.'rot'.or.allSymOps(op)%typ.eq.'imp')then
          chg=tqin(l+1)
          pt=rotax_to_cartesian(atm,op,tqin(l))
          call apply_atm_sym_ops(atm,pt,sym_pts,npts)
          nlast=npts
          do j=1,npts
            n=n+1
            num_atm_chgs(a)=num_atm_chgs(a)+1
            tq(n*4-3:n*4-1)=sym_pts(j,:)
            tq(n*4)=chg
            local_charge=local_charge+chg
          enddo
          l=l+2
        elseif(allSymOps(op)%typ.eq.'ref')then
          chg=tqin(l+2)
          pt=refplane_to_cartesian(atm,op,tqin(l),tqin(l+1))
          call apply_atm_sym_ops(atm,pt,sym_pts,npts)
          nlast=npts
          do j=1,npts
            n=n+1
            num_atm_chgs(a)=num_atm_chgs(a)+1
            tq(n*4-3:n*4-1)=sym_pts(j,:)
            tq(n*4)=chg
            local_charge=local_charge+chg
          enddo
          l=l+3
        elseif(allSymOps(op)%typ.eq.'inv')then ! nuclear charge spawns no others
          n=n+1
          num_atm_chgs(a)=num_atm_chgs(a)+1
          tq(n*4-3:n*4-1)=atom_pos(atm,1:3)
          tq(n*4)=tqin(l)
          local_charge=local_charge+tqin(l)
          nlast=1
          l=l+1
        elseif(allSymOps(op)%typ.eq.'nuc')then ! nuclear charge spawns no others
          n=n+1
          num_atm_chgs(a)=num_atm_chgs(a)+1
          tq(n*4-3:n*4-1)=atom_pos(atm,1:3)
          tq(n*4)=tqin(l)
          local_charge=local_charge+tqin(l) ! keep track of total charge of this arrangement
          nlast=1 ! remember how many charges were spawned by the last fitted charge
          l=l+1
        elseif(allSymOps(op)%typ.eq.'fre')then ! no symmetry
          chg=tqin(l+3)
          pt=point_to_cartesian(atm,tqin(l),tqin(l+1),tqin(l+2))
          call apply_atm_sym_ops(atm,pt,sym_pts,npts)
          nlast=npts
          do j=1,npts
            n=n+1
            num_atm_chgs(a)=num_atm_chgs(a)+1
            tq(n*4-3:n*4-1)=sym_pts(j,:)
            tq(n*4)=chg
            local_charge=local_charge+chg
          enddo
          l=l+4
        else
          call throw_error('ERROR: unrecognized operation type in spawn_sym_chgs')
        endif
      enddo !i
    enddo !a

    !store array indices in 'tq()' of charges to correct
    do i=1,nlast
      corr_chgs(i)=n*4-((i-1)*4)
    enddo

    tlast=nlast
    if(replicate_atoms)then !also transform atom charges to other sea's
      do a=1,num_symFitAtms
        atm1=symFitAtms(a)
        do b=2,Natom !loop over sea's
          if(atom_sea(a,b) == 0) exit
          atm2=atom_sea(a,b)
          if(atom_sea(sea_ops(atm2,1),1).ne.atm1) call throw_error('sea mismatch in spawn_sym_chgs')
          do i=1,num_atm_chgs(a)
            nn=n1(a)+4*i-3
            tpos=tq(nn:nn+2)
            do j=1,sea_ops(atm2,4)
              tpos=matmul(tpos(:),allSymOps(sea_ops(atm2,2))%M) !transform charge coords
            enddo
            if(sea_ops(atm2,3).ne.0)then
              tpos=matmul(tpos(:),allSymOps(sea_ops(atm2,3))%M) !transform charge coords
            endif
            n=n+1
            tq(n*4-3:n*4-1)=tpos(:)
            tq(n*4)=tq(nn+3)
            local_charge=local_charge+tq(nn+3)
            do k=1,tlast
              if(nn+3.eq.corr_chgs(k)) then ! if this charge was spawned from the atomic charge used to correct the total charge
                nlast=nlast+1 ! final charge needs correcting
                corr_chgs(nlast)=n*4
                exit
              endif
            enddo
          enddo
        enddo
      enddo
    endif

   ! correct total charge
    chg_corr=(total_charge-local_charge)/nlast
    do i=1,nlast
      tq(corr_chgs(i))=chg_corr
    enddo
    ! transform charges back to global axis:
    do i=1,nq*4-3,4
      tq(i:i+2)=tq(i:i+2)+atom_com
!      write(*,'(A,4F7.3)') ' H ',tq(i:i+2),tq(i+3)
    enddo
!    write(*,'(/)')
    q(1:size(q,dim=1))=tq(1:size(q,dim=1))
end subroutine spawn_sym_chgs

!-------------------------------------------------------------------------------
! applies all symmetry operations for a given atom to a point and returns an
! array of spawned points

subroutine apply_atm_sym_ops(atom,pt,sym_pts,npts)
    implicit none
 
    integer :: atom,npts
    real(rp), dimension(:) :: pt
    real(rp), dimension(:,:) :: sym_pts

    integer :: i,j,k,l,op,t
    logical :: duplicate
    real(rp), dimension(3) :: tpos

    npts=1
    sym_pts(1,1:3)=pt(:)
    do i=1, num_atm_sym_ops(atom)
      op=atm_sym_ops(atom,i)
      j=1
      do while (j.le.npts)
        tpos(:)=sym_pts(j,1:3)
        do k=1, allSymOps(op)%n
          tpos=matmul(tpos(:),allSymOps(op)%M)
          duplicate=.false.
          do l=1,npts
            if(equal(sqrt(sum((sym_pts(l,:)-tpos(:))**2)),0d0))then
              duplicate=.true.
            endif
          enddo !l
          if(.not. duplicate) then
            npts=npts+1
            sym_pts(npts,:)=tpos(:)
          endif
        enddo !k
        j=j+1
      enddo !j
    enddo !i

end subroutine apply_atm_sym_ops

!-------------------------------------------------------------------------------
! read fitted symmetry-constrained parameters from file
function read_sym_file(filename,solutions,sym_solution_ops,num_sym_solution_ops,&
         sqdim,num_charges,atm,maxq)
implicit none
integer :: i,l,ios,atm,sqdim,num_charges,maxq
logical :: read_sym_file
character(len=1024) :: dummy,filename
character(len=4) :: label
character(len=3) :: typ
integer, dimension(:,:) :: num_sym_solution_ops
integer, dimension(:,:,:) :: sym_solution_ops
real(rp), dimension(:) :: solutions

read_sym_file=.false.
sqdim=0
open(30, file=trim(filename), status="old", action="read", iostat = ios)
if(ios == 0) then

  l=1
  if(allocated(num_atm_fit_ops)) deallocate(num_atm_fit_ops)
  if(allocated(atm_fit_ops)) deallocate(atm_fit_ops)
  allocate(num_atm_fit_ops(Natom))  
  allocate(atm_fit_ops(Natom,maxq))
  read(30,*,iostat=ios) num_atm_fit_ops(atm)
  do i=1,num_atm_fit_ops(atm)
    ! read symmetry operation index from allSym array, operation type and label
    read(30,*) dummy,atm_fit_ops(atm,i),typ,label
    ! check the type and label match what we have for this index in allSym array
    if(allSymOps(atm_fit_ops(atm,i))%typ.ne.typ.or. &
       allSymOps(atm_fit_ops(atm,i))%label.ne.label)then
      write(*,'(A,I0,5(X,A))') 'For operation index ',atm_fit_ops(atm,i),&
         ', type and label should be: ',allSymOps(atm_fit_ops(atm,i))%typ,&
         allSymOps(atm_fit_ops(atm,i))%label,', read: ',typ,label
      call throw_error('Mismatch when reading results file '//filename)
    endif
    if(typ.eq.'rot'.or.typ.eq.'imp')then
      if(i.ne.num_atm_fit_ops(atm))then
        read(30,*) dummy,solutions(l:l+1)
        l=l+2
      else
        read(30,*) dummy,solutions(l)
        l=l+1
      endif
    elseif(typ.eq.'ref')then
      if(i.ne.num_atm_fit_ops(atm))then
        read(30,*) dummy,solutions(l:l+2)
        l=l+3
      else
        read(30,*) dummy,solutions(l:l+1)
        l=l+2
      endif
    elseif(typ.eq.'fre')then
      if(i.ne.num_atm_fit_ops(atm))then
        read(30,*) dummy,solutions(l:l+3)
        l=l+4
      else
        read(30,*) dummy,solutions(l:l+2)
        l=l+3
      endif
    elseif(typ.eq.'nuc'.or.typ.eq.'inv')then
      if(i.ne.num_atm_fit_ops(atm))then
        read(30,*) dummy,solutions(l)
        l=l+1
      endif
    else
      call throw_error('Unrecognized sym type in read_sym_file')
    endif
  enddo
  sqdim=l-1
  sym_solution_ops(num_charges,atm,1:num_atm_fit_ops(atm))= &
                            atm_fit_ops(atm,1:num_atm_fit_ops(atm))
  num_sym_solution_ops(num_charges,atm)=num_atm_fit_ops(atm)
  read_sym_file=.true.
!  print*,'previous results: ',solutions(:)

close(30)
 
endif

end function read_sym_file

!-------------------------------------------------------------------------------
! write fitted symmetry-constrained parameters to file
subroutine write_sym_file(charges,a,nchg,filename)
implicit none
integer :: ios,i,l,op,nchg
integer, optional :: a
real(rp), dimension(:), intent(in) :: charges
character(len=*), intent(in), optional :: filename
character(len=1024) :: outfile, dummy, prefix = ''

write(outfile,'(I0)') nchg !number of charges
if(present(a)) then
    write(dummy,'(I0)') a
    outfile = "symatm"//trim(dummy)//"_"//trim(outfile)//"charges.fit"
else
    outfile = "sym_"//trim(outfile)//"charges.fit"
end if
if(trim(prefix) /= '') outfile = trim(prefix)//"/"//trim(outfile)

if(present(filename)) outfile = filename

open(30, file=trim(outfile), status="replace", action="write", iostat = ios)
if(ios /= 0) call throw_error('Could not open "'//trim(outfile)//'" for writing')
write(30,'(I0)') num_atm_fit_ops(a) !number of symmetry-axes, mirror planes etc. used to fit charges

l=1
do i=1,num_atm_fit_ops(a)
  op=atm_fit_ops(a,i)
  write(30,'(A,I0,4A)') 'SymOp ',op,' ',allSymOps(op)%typ,' ',allSymOps(op)%label
  if(allSymOps(op)%typ.eq.'rot'.or.allSymOps(op)%typ.eq.'imp')then
    if(i<num_atm_fit_ops(a))then
      write(30,'(A,2(X,F9.6))') 'r.q: ',charges(l),charges(l+1)
      l=l+2
    else ! final charge (constrained)
      write(30,'(A,2(X,F9.6))') 'r: ',charges(l)
      l=l+1
    endif
  elseif(allSymOps(op)%typ.eq.'ref')then
    if(i<num_atm_fit_ops(a))then
      write(30,'(A,3(X,F9.6))') 'x.y.q: ',charges(l),charges(l+1),charges(l+2)
      l=l+3
    else ! final charge (constrained)
      write(30,'(A,2(X,F9.6))') 'x.y: ',charges(l),charges(l+1)
      l=l+2
    endif
  elseif(allSymOps(op)%typ.eq.'inv'.or.allSymOps(op)%typ.eq.'nuc')then
    if(i<num_atm_fit_ops(a))then
      write(30,'(A,3(X,F9.6))') 'q: ',charges(l)
      l=l+1
    else ! final charge (constrained)
      write(30,'(A,2(X,F9.6))') 'q: '
    endif
  elseif(allSymOps(op)%typ.eq.'fre')then
    if(i<num_atm_fit_ops(a))then
      write(30,'(A,3(X,F9.6))') 'x.y.z.q: ',charges(l),charges(l+1),charges(l+2),&
         charges(l+3)
      l=l+4
    else ! final charge (constrained)
      write(30,'(A,2(X,F9.6))') 'x.y.z: ',charges(l),charges(l+1),charges(l+2)
      l=l+3
    endif
  else
    write(*,'(/,3A,/)') 'ERROR: unrecognized operation type ',allSymOps(op)%typ,&
        ' in write_sym_file'
    stop
  endif
enddo
close(30)

end subroutine write_sym_file
!-------------------------------------------------------------------------------


function getNumSymOps(atom)
    implicit none
    integer :: getNumSymOps,atom

    getNumSymOps=num_atm_sym_ops(atom)
end function getNumSymOps

!===============================================================================
!   BELOW, ALL THE DIFFERENT SYMMETRY OPERATIONS ARE DEFINED AS 3x3 matrices
!===============================================================================

!-------------------------------------------------------------------------------
! returns the inversion matrix
function inversion()
    implicit none
    real(rp), dimension(3,3), parameter  :: inv = reshape((/ -1, 0, 0, 0, -1, 0, 0, 0, -1/), (/3,3/))
    real(rp), dimension(3,3) :: inversion
    inversion = inv
end function inversion
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! returns the inversion matrix
function identity()
    implicit none
    real(rp), dimension(3,3), parameter  :: ide = reshape((/ 1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3,3/))
    real(rp), dimension(3,3) :: identity
    identity = ide
end function identity
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! returns the reflection matrix for a plane defined by the normalized normal vector n
function reflection(n)
    implicit none
    real(rp), parameter :: pi = acos(-1d0)
    real(rp), dimension(3,3)           :: reflection !output rotation matrix
    real(rp), dimension(3), intent(in) :: n     ! normalized normal vector
    integer :: i
    do i = 1,3
        reflection(i,i) = 1d0-2d0*n(i)**2
    end do
    !off diagonal elements
    reflection(1,2) = -2d0*n(1)*n(2)
    reflection(2,1) = reflection(1,2)
    reflection(1,3) = -2d0*n(1)*n(3)
    reflection(3,1) = reflection(1,3)
    reflection(2,3) = -2d0*n(2)*n(3)
    reflection(3,2) = reflection(2,3)
end function reflection
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! returns the improper rotation matrix around the normalized vector "axis" by angle "alpha"
function improper(axis,alpha,n)
    implicit none
    real(rp), dimension(3,3)    :: improper
    real(rp), dimension(3), intent(in) :: axis  ! rotation axis
    real(rp), optional,     intent(in) :: alpha ! rotation angle
    integer, optional,     intent(in) :: n     ! to directly specify a "Cn" rotation 
    if(present(alpha).and..not.present(n)) then
        improper = matmul(rotation(axis, alpha = alpha),reflection(axis))
    else if(present(n).and..not.present(alpha)) then
        improper = matmul(rotation(axis, n = n),reflection(axis))
    else
        call throw_error("improper function must be called with either argument alpha or argument n!")
    end if   

end function improper
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! returns the rotation matrix around the normalized vector "axis" by angle "alpha"
function rotation(axis,alpha,n)
    implicit none
    real(rp), parameter :: pi = acos(-1d0)
    real(rp), dimension(3,3)           :: rotation !output rotation matrix
    real(rp), dimension(3), intent(in) :: axis  ! rotation axis
    real(rp), optional,     intent(in) :: alpha ! rotation angle
    integer, optional,     intent(in) :: n     ! to directly specify a "Cn" rotation
    real(rp) :: tmp, tmp2, cosine, sine
    integer :: i,j
    if(present(alpha).and..not.present(n)) then
        cosine = cos(alpha)
        sine   = sin(alpha)
    else if(present(n).and..not.present(alpha)) then
        cosine = cos(pi/(0.5d0*n))
        sine   = sin(pi/(0.5d0*n))
    else
        call throw_error("rotation function must be called with either argument alpha or argument n!")
    end if
    !calculate diagonal elements
    tmp2 = 0d0
    j = 1
    do i = 1,3
        tmp = axis(i)**2
        if(tmp > tmp2) then
            j = i
            tmp2 = tmp
        end if
        rotation(i,i) = tmp+(1d0-tmp)*cosine
    end do
    !calculate off-diagonal elements,
    !the if is so that all rotations are in the same orientation
    if(axis(j) > 0d0) then
        tmp = axis(1)*axis(2)*(1d0-cosine)
        rotation(1,2) = tmp - axis(3)*sine
        rotation(2,1) = tmp + axis(3)*sine
        tmp = axis(1)*axis(3)*(1d0-cosine)
        rotation(1,3) = tmp + axis(2)*sine
        rotation(3,1) = tmp - axis(2)*sine
        tmp = axis(2)*axis(3)*(1d0-cosine)
        rotation(2,3) = tmp - axis(1)*sine
        rotation(3,2) = tmp + axis(1)*sine
    else
        tmp = (-axis(1))*(-axis(2))*(1d0-cosine)
        rotation(1,2) = tmp + axis(3)*sine
        rotation(2,1) = tmp - axis(3)*sine
        tmp = (-axis(1))*(-axis(3))*(1d0-cosine)
        rotation(1,3) = tmp - axis(2)*sine
        rotation(3,1) = tmp + axis(2)*sine
        tmp = (-axis(2))*(-axis(3))*(1d0-cosine)
        rotation(2,3) = tmp + axis(1)*sine
        rotation(3,2) = tmp - axis(1)*sine    
    end if
end function rotation
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! throws an error message and terminates the code
subroutine throw_error(message)
    implicit none
    character(len=*), intent(in) :: message
    write(*,'(A)') "ERROR: "//message
    stop
end subroutine throw_error


end module symmetry

!=== gfortran punch-to-dcm.f90 -o punch-to-dcm.x -llapack -lblas
! Program to convert multipole exapnsion in global axis (GDMA punch file format)
! to DCM charge distribution in local axis for subsequent use in CHARMM
program solve_eigen
  implicit none

  integer, parameter           :: N = 3
  integer, parameter           :: NATMAX = 100
  integer, parameter           :: MAXCHG = 6
  double precision, parameter  :: d = 0.25d0
  double precision, parameter  :: b2a = 0.52917720859d0
  integer                      :: i, j, k, l

  double precision             :: A(1:N,1:N)
! double precision             :: MAT_CHK(1:N,1:N)
! double precision             :: MAT_EVAL(1:N,1:N)
  double precision             :: EVAL(1:N)
  double precision             :: V(1:N,1:N)
  double precision             :: V3(NATMAX,1:N,1:N)
  double precision             :: SEVAL(1:N) ! for sorted eigenvalues
  double precision             :: SV(1:N,1:N) ! for sorted eigenvectors

  double precision             :: Q00
  double precision             :: Q10, Q11c, Q11s 
  double precision             :: Q20, Q22c, Q22s, Q21c, Q21s
  double precision             :: Q30, Q31c, Q31s, Q32c, Q32s, Q33c, Q33s
  double precision             :: Q10p, Q11cp, Q11sp 
  double precision             :: Q20p, Q22cp, Q22sp, Q21cp, Q21sp
  double precision             :: Qxx, Qyy, Qzz, Qxy, Qxz, Qyz 
  double precision             :: Qxxp, Qyyp, Qzzp
  double precision             :: e(1:NATMAX,1:MAXCHG)
  double precision             :: PTCHGX(1:NATMAX,1:MAXCHG), PTCHGY(1:NATMAX,1:MAXCHG)
  double precision             :: PTCHGZ(1:NATMAX,1:MAXCHG)

  integer                      :: CNT,readstat
  integer                      :: NPTCHG(NATMAX)
  character(len=100)           :: line
  character(len=100)           :: arg
  character(len=2)             :: ATM_SYM(1:NATMAX)
  character(len=2)             :: tmp1
  character(len=4)             :: tmp2
  character(len=100)           :: resname,dcmfile
  double precision             :: X, Y, Z
  double precision             :: RX, RY
  double precision             :: CART_COORD(1:NATMAX,1:3)
  double precision             :: MONO(1:NATMAX)
  double precision             :: DIP(1:NATMAX,1:3)
  double precision             :: QUAD(1:NATMAX,1:5)

  double precision             :: gex(1:3), gey(1:3), gez(1:3)

  double precision             :: EX(1:3,1:3), EY(1:3,1:3), EZ(1:3,1:3)

  double precision             :: px, py, pz
  double precision             :: gmx, gmy, gmz
  double precision             :: lx, ly, lz

  integer                      :: NFRAME, INUM, ANUM(1:NATMAX,1:3)
  integer                      :: RANK
  logical                      :: FLAG_NDC


  if(iargc().ne.2)then
    write(6,*) 'Usage: punch-to-dcm file.punch frame_def.txt'
    return 
  endif
  CALL getarg(1, arg)
  write(*,*) '- Opening file ',TRIM(arg),' for reading'
  open(unit=100,file=arg,IOSTAT=readstat,STATUS='OLD',ACTION='read')
  if(readstat .ne. 0) then
    write(*,*) arg,' cannot be opened'
    stop 1
  end if
  read(100,*) line
  read(100,*) line

  CALL getarg(2, arg)

  CNT = 0
  readloop : do

   read(100,*,iostat=readstat)
   if (readstat .ne. 0) then
     exit readloop
   endif

! First parse the multipoles and nuclear coords from the punch file

   read(100,*,iostat=readstat) tmp1, X, Y, Z, tmp2, RANK
   ! handle 1 blank line at eof
   if (readstat .ne. 0) then
     exit readloop
   endif
   CNT = CNT + 1
   ATM_SYM(CNT) = tmp1
   CART_COORD(CNT,1:3) = (/X, Y, Z/)
   read(100,*) MONO(CNT)
   if(RANK.gt.0)then
     read(100,*) DIP(CNT,1:3)
   else
     DIP(CNT,1:3)=0.d0
   endif
   if(RANK.gt.1)then
     read(100,*) QUAD(CNT,1:5)
   else
     QUAD(CNT,1:5)=0.d0
   endif
  enddo readloop

  ! Get residue name and number of frames from axis file
  NFRAME=0
  write(*,*) '- Opening file ',TRIM(arg),' for reading'
  open(unit=108,file=arg,IOSTAT=readstat,STATUS='OLD',ACTION='read')
  if(readstat .ne. 0) then
    write(*,*) arg,' cannot be opened'
    stop 1
  end if
  read(108,*,iostat=readstat) resname, NFRAME
  close(108)
  write(*,*) '- Residue name is ',TRIM(resname)
  write(*,*) '- Axis file contains ',NFRAME,' frames'

  open(unit=101,file='multipoles_paxis.dat')
  open(unit=102,file='octa_point_charges.dat')
  open(unit=103,file='multipoles_rebuilt.dat')
  open(unit=109,file='coord_point_charges_paxis.xyz')
  open(unit=104,file='coord_point_charges_gaxis.xyz')
  tmp2='.dcm'
  dcmfile=TRIM(resname) // TRIM(tmp2)
  open(unit=105,file=dcmfile)
  open(unit=106,file='coord_mol.xyz')
  open(unit=107,file='octapoles_from_point_charges.dat')
  write(*,*) '- DCM parameter file will be written to ',TRIM(dcmfile)

  write(101,'(a)')'Multipoles in principle axis system'
  write(101,'(a)')'==================================='
  write(101,'(a,9(5x,a4,9x))')'Atoms','Q00','Q10','Q11c','Q11s','Q20','Q21c','Q21s','Q22c','Q22s'
  write(101,'("-----",9(5x,"-------",6x))')

  write(102,'(a)')'Point charges in octahedral model without the central charge'
  write(102,'(a)')'============================================================'

  write(102,'(a,7(8x,a,10x))')'Atoms', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'Total e'
  write(102,'("-----",7(8x,"-------",5x))')

  write(103,'(a)')'Multipoles rebuilt using point charges around each atomic center'
  write(103,'(a)')'================================================================'
  write(103,'(a,9(5x,a4,9x))')'Atoms','Q00','Q10','Q11c','Q11s','Q20','Q21c','Q21s','Q22c','Q22s'
  write(103,'("-----",9(5x,"-------",6x))')

  write(104,'(i3)') 6*CNT
  write(104,'(a,i3,a)')'Cartesian coordinates of point charges in global axis system'

  write(109,'(i3)') 6*CNT
  write(109,'(a,i3,a)')'Cartesian coordinates of point charges in global atomic axis system'

! write(105,'(i3)') NPTCHG*CNT
! write(105,'(a)')'Cartesian coordinates of 6 point charges per atom in local axis system'

! FOR CHARMM INPUT PREPARATION
  write(105,'(i1,10x,a)') 1,'! no. residue types defined here'
  write(105,*)
  write(105,'(a,7x,a)') resname, '! residue name'
  write(105,'(i3,10x,a)') NFRAME,'! no. axis system frames'

  write(106,'(i3)') CNT
  write(106,'(a)')'Cartesian coordinates of the molecule'

  write(107,'(a)')'Octapoles using point charges around each atomic center'
  write(107,'(a)')'======================================================='
  write(107,'(a,8(5x,a4,9x))')'Atoms','Q00','Q30','Q31c','Q31s','Q32c','Q32s','Q33c','Q33s'
  write(107,'("-----",8(5x,"-------",6x))')

  do i = 1, CNT
   write(106,'(a,3f16.10)') ATM_SYM(i), CART_COORD(i,1:3)
  enddo
  
  do i = 1, CNT
   ! Create atomic quadrupole matrix in Cartesian form from spherical moments in
   ! punch file
   Q00  = MONO(i)
   Q10  = DIP(i,1) 
   Q11c = DIP(i,2)
   Q11s = DIP(i,3)
   Q20  = QUAD(i,1)
   Q21c = QUAD(i,2)
   Q21s = QUAD(i,3)
   Q22c = QUAD(i,4)
   Q22s = QUAD(i,5)

   Qxx = (-0.5d0*Q20 + 0.5d0*(sqrt(3.0d0))*Q22c)
   Qyy = (-0.5d0*Q20 - 0.5d0*(sqrt(3.0d0))*Q22c)
   Qzz = Q20
   Qxy = 0.5d0*(sqrt(3.0d0))*Q22s
   Qxz = 0.5d0*(sqrt(3.0d0))*Q21c
   Qyz = 0.5d0*(sqrt(3.0d0))*Q21s
  
   A(1,1) = Qxx
   A(1,2) = Qxy
   A(1,3) = Qxz
   A(2,1) = Qxy
   A(2,2) = Qyy
   A(2,3) = Qyz
   A(3,1) = Qxz
   A(3,2) = Qyz
   A(3,3) = Qzz

!========================================================================================
! PRINT ORIGINAL QUADROPOLE TENSOR
!========================================================================================
!  do j = 1, N
!    write(*,'(3f15.8)')A(j,1:N)
!  enddo
!  print*,'**********'
!========================================================================================
   ! diagonalize quadrupole tensor
   call eigen(N, A, EVAL, V)
!========================================================================================
!  call Jacobi(N, A, EVAL,V)
!========================================================================================
!  print*,'EIGENVALUES'
!  write(*,'(3f15.8)')EVAL
! EIGENVECTOR MATRIX
!  V3(i,:,:) = V
!  print*,'EIGENVECTOR MATRIX'
!  do j = 1, N
!    write(*,'(3f15.8)')V3(i,j,:)
!  enddo
!========================================================================================
! TRANSPOSE OF EIGENVECTOR MATRIX
!========================================================================================
!  print*,'SORT EIGENVALUES AND EIGENVECTORS (MAKE Z LARGEST)'
!  Dropped: risks switching to left-handed axis system
!   SEVAL=EVAL
!   SV=V
!   do k=1,2
!    if(ABS(EVAL(k)).gt.ABS(SEVAL(3)))then
!     SEVAL(k)=SEVAL(3)
!     SEVAL(3)=EVAL(k)
!     SV(k,1:3)=SV(3,1:3)
!     SV(3,1:3)=V(k,1:3)
!    endif
!   enddo
!   EVAL=SEVAL
!   V=SV
   ! store principal axes in V
   V = transpose(V)
!
   V3(i,:,:) = V
!========================================================================================
!  print*,'TRANSPOSE OF EIGENVECTOR MATRIX'
!  do j = 1, N
!    write(*,'(3f15.8)')V3(i,j,:)
!  enddo
!========================================================================================
!TEST MATRIX BACK TRANSFORMATION
!========================================================================================
!  MAT_EVAL(1,1) = EVAL(1)
!  MAT_EVAL(1,2) = 0.0d0
!  MAT_EVAL(1,3) = 0.0d0
!  MAT_EVAL(2,1) = 0.0d0
!  MAT_EVAL(2,2) = EVAL(2)
!  MAT_EVAL(2,3) = 0.0d0
!  MAT_EVAL(3,1) = 0.0d0
!  MAT_EVAL(3,2) = 0.0d0
!  MAT_EVAL(3,3) = EVAL(3)
!  MAT_CHK = matmul( V, matmul(MAT_EVAL,transpose(V) ) )
!  do j = 1, N
!    write(*,'(3f15.8)')MAT_CHK(j,1:N)
!  enddo
!  print*,'**********'
!========================================================================================
   ! Transform dipole components to new principal axes
   Q11cp = dot_product((/Q11c,Q11s,Q10/), V3(i,1,1:3))
   Q11sp = dot_product((/Q11c,Q11s,Q10/), V3(i,2,1:3))
   Q10p  = dot_product((/Q11c,Q11s,Q10/), V3(i,3,1:3))

   ! Convert Cartesian multipoles to spherical multipoles in new principal axes
   ! (off-diagonal terms are zero due to diagonalization)
   Q20p  = EVAL(3)
   Q21cp = (2.0d0/sqrt(3.0d0)) * 0.0d0
   Q21sp = (2.0d0/sqrt(3.0d0)) * 0.0d0
   Q22cp = (1.0d0/sqrt(3.0d0)) * (EVAL(1) - EVAL(2))
   Q22sp = (2.0d0/sqrt(3.0d0)) * 0.0d0

   write(101,'(a,9f18.10)') ATM_SYM(i), Q00, Q10p, Q11cp, Q11sp, Q20p, Q21cp, Q21sp, Q22cp, Q22sp
!=======================================================================================
! Calculate point charges using appropriate model, selected according to which multipole
! moment components (if any) are zeroed
! Define  point charges in principle axis/eigen vector basis                                      |
!==================================================================================================
!              e1     e2    e3    e4     e5     e6
   Qxxp = -0.5d0*Q20p + (sqrt(3.0d0)/2.0d0)*Q22cp
   Qyyp = -0.5d0*Q20p - (sqrt(3.0d0)/2.0d0)*Q22cp
   Qzzp = Q20p

   if(Q10p.eq.0.d0 .and. Q11cp.eq.0.d0 .and. Q11sp.eq.0.d0 .and. Q20p.eq.0.d0 .and. &
      Q22cp.eq.0.d0)then
    ! charge only / no moments
    print*,'ATOM ',I,': Using central chg only'
    NPTCHG(I)=1
    e(i,1)=Q00
    PTCHGX(I,1:NPTCHG(I)) = (/0.0d0/)
    PTCHGY(I,1:NPTCHG(I)) = (/0.0d0/)
    PTCHGZ(I,1:NPTCHG(I)) = (/0.0d0/)
   elseif(Q20p.ne.0.d0.and.Q22cp.ne.0.d0.and.Q11cp.eq.0.d0 .and. Q11sp.eq.0.d0)then
    ! 5-chg model for nonzero Q20, Q22c and no xy dipole component (z-only)
    print*,'ATOM ',I,': Using 5-chg model A'
    NPTCHG(I)=5
    e(i,1) = Q00 - Q20p/d**2 - sqrt(3.d0)*Q22cp/d**2
    e(i,2) = Q10p/(2.d0*d) + Q20p/(2.d0*d**2) + sqrt(3.d0)*Q22cp/(6.d0*d**2)
    e(i,3) = -Q10p/(2.d0*d) + Q20p/(2.d0*d**2) + sqrt(3.d0)*Q22cp/(6.d0*d**2)
    e(i,4) = Q22cp*sqrt(3.d0)/(6.d0*d**2)
    e(i,5) = Q22cp*sqrt(3.d0)/(6.d0*d**2)
    PTCHGX(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0, 1.0d0,-1.0d0/)
    PTCHGY(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
    PTCHGZ(I,1:NPTCHG(I)) = (/0.0d0, 1.0d0,-1.0d0, 0.0d0, 0.0d0/)
   elseif(Q20p.ne.0.d0.and.Q22cp.ne.0.d0.and.(Q11cp.ne.0.d0 .and. Q11sp.eq.0.d0))then
    ! 5-chg model for nonzero Q20, Q22c and xy dipole component along x or y axis
    print*,'ATOM ',I,': Using 5-chg model B1'
    NPTCHG(I)=5
    e(i,1) = Q00 - Q20p/d**2 - (sqrt(3.d0)*Q22cp) / d**2
    e(i,2) = Q10p/(2.d0*d) + Q20p/(2.d0*d**2) + (Q22cp / (2.d0*sqrt(3.d0)*d**2))
    e(i,3) = -Q10p/(2.d0*d) + Q20p/(2.d0*d**2) + (Q22cp / (2.d0*sqrt(3.d0)*d**2))
    e(i,4) = Q11cp/(2.d0*d) + Q22cp / (sqrt(3.d0)*d**2)
    e(i,5) = -Q11cp/(2.d0*d) + Q22cp / (sqrt(3.d0)*d**2)
    PTCHGX(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0, 1.d0,-1.d0 /)
    PTCHGY(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0, 0.d0, 0.d0 /)
    PTCHGZ(I,1:NPTCHG(I)) = (/0.0d0, 1.0d0,-1.0d0, 0.d0, 0.d0 /)
   elseif(Q20p.ne.0.d0.and.Q22cp.ne.0.d0.and.(Q11cp.eq.0.d0 .and. Q11sp.ne.0.d0))then
    ! 5-chg model for nonzero Q20, Q22c and xy dipole component along x or y axis
    print*,'ATOM ',I,': Using 5-chg model B2'
    NPTCHG(I)=5
    e(i,1) = Q00 - Q20p/d**2 + (sqrt(3.d0)*Q22cp) / d**2
    e(i,2) = Q10p/(2.d0*d) + Q20p/(2.d0*d**2) - (Q22cp / (2.d0*sqrt(3.d0)*d**2))
    e(i,3) = -Q10p/(2.d0*d) + Q20p/(2.d0*d**2) - (Q22cp / (2.d0*sqrt(3.d0)*d**2))
    e(i,4) = Q11sp/(2.d0*d) - Q22cp / (sqrt(3.d0)*d**2)
    e(i,5) = -Q11sp/(2.d0*d) - Q22cp / (sqrt(3.d0)*d**2)
    PTCHGX(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0, 0.d0, 0.d0 /)
    PTCHGY(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0, 1.d0,-1.d0 /)
    PTCHGZ(I,1:NPTCHG(I)) = (/0.0d0, 1.0d0,-1.0d0, 0.d0, 0.d0 /)
   elseif(Q11cp.eq.0.d0 .and. Q11sp.eq.0.d0 .and. Q20p.eq.0.d0 .and. Q22cp.eq.0.d0)then
   ! z-dipole only
    print*,'ATOM ',I,': Using 2-chg z-axis model'
    NPTCHG(I)=2
    e(i,1) = Q10p/(2.d0*d)
    e(i,2) = -Q10p/(2.d0*d)
    PTCHGX(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0/)
    PTCHGY(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0/)
    PTCHGZ(I,1:NPTCHG(I)) = (/1.0d0,-1.0d0/)
   elseif(Q11cp.eq.0.d0 .and. Q11sp.eq.0.d0 .and. Q20p.eq.0.d0 .and. Q22cp.eq.0.d0)then
   ! charge and z-dipole only
    print*,'ATOM ',I,': Using 3-chg z-axis model A'
    NPTCHG(I)=3
    e(i,1) = Q00
    e(i,2) = Q10p/(2.d0*d)
    e(i,3) = -Q10p/(2.d0*d)
    PTCHGX(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0/)
    PTCHGY(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0/)
    PTCHGZ(I,1:NPTCHG(I)) = (/0.0d0, 1.0d0,-1.0d0/)
   elseif(Q11cp.eq.0.d0 .and. Q11sp.eq.0.d0 .and. Q20p.ne.0.d0 .and. Q22cp.eq.0.d0)then
   ! charge, z-dipole and Q20 quadrupole only
    print*,'ATOM ',I,': Using 3-chg z-axis model B'
    NPTCHG(I)=3
    e(i,1) = Q00 - Q20p/d**2
    e(i,2) = Q10p/(2.d0*d)+Q20p/(2.d0*d**2)
    e(i,3) = -Q10p/(2.d0*d)+Q20p/(2.d0*d**2)
    PTCHGX(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0/)
    PTCHGY(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0, 0.0d0/)
    PTCHGZ(I,1:NPTCHG(I)) = (/0.0d0, 1.0d0,-1.0d0/)
   else
   ! 6-chg model for general case (equivalent to formula in DCM paper replacing
   ! quadrupole terms with Cartesian form)
    print*,'ATOM ',I,': Using 6-chg general octahedral model'
    NPTCHG(I)=6
    e(i,1) = Q00/6.0d0 + Qxxp/(3.0d0*d**2) + Q11cp/(2.0d0*d)
    e(i,2) = Q00/6.0d0 + Qxxp/(3.0d0*d**2) - Q11cp/(2.0d0*d)
    e(i,3) = Q00/6.0d0 + Qyyp/(3.0d0*d**2) + Q11sp/(2.0d0*d)
    e(i,4) = Q00/6.0d0 + Qyyp/(3.0d0*d**2) - Q11sp/(2.0d0*d)
    e(i,5) = Q00/6.0d0 + Qzzp/(3.0d0*d**2) + Q10p /(2.0d0*d)
    e(i,6) = Q00/6.0d0 + Qzzp/(3.0d0*d**2) - Q10p /(2.0d0*d)
    PTCHGX(I,1:NPTCHG(I)) = (/1.0d0,-1.0d0,0.0d0, 0.0d0,0.0d0, 0.0d0/)
    PTCHGY(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0,1.0d0,-1.0d0,0.0d0, 0.0d0/)
    PTCHGZ(I,1:NPTCHG(I)) = (/0.0d0, 0.0d0,0.0d0, 0.0d0,1.0d0,-1.0d0/)
   endif

   write(102,'(a,2x,7f20.10)') ATM_SYM(i), e(i,1:NPTCHG(I)), sum(e(i,1:NPTCHG(I)))

  enddo

  ! recalculate atomic multipoles from our new charge models for each atom.
  ! Should recover original atomic multipole moments. Use octopole to check d is
  ! small enough
  do i = 1, CNT

   Q00  = sum(e(i,1:NPTCHG(I))) 

   Q10  = d*dot_product (PTCHGZ(I,1:NPTCHG(I)), e(i,1:NPTCHG(I)))
   Q11c = d*dot_product (PTCHGX(I,1:NPTCHG(I)), e(i,1:NPTCHG(I)))
   Q11s = d*dot_product (PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I)))

   Q20 = 0.5d0 * d**2 * ( 3.0d0*dot_product (PTCHGZ(I,1:NPTCHG(I))*PTCHGZ(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                              ( dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) +&
                                dot_product (PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) +&
                                dot_product (PTCHGZ(I,1:NPTCHG(I))*PTCHGZ(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) ) )

   Q21c = sqrt(3.0d0) * d**2 * ( sum ( e(i,1:NPTCHG(I)) * PTCHGX(I,1:NPTCHG(I)) * PTCHGZ(I,1:NPTCHG(I)) )) 
   Q21s = sqrt(3.0d0) * d**2 * ( sum ( e(i,1:NPTCHG(I)) * PTCHGY(I,1:NPTCHG(I)) * PTCHGZ(I,1:NPTCHG(I)) )) 

   Q22c = (sqrt(3.0d0)/2.0d0) * d**2 * ( dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                                         dot_product (PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) )
  
   Q22s = sqrt(3.0d0) * d**2 * ( sum ( e(i,1:NPTCHG(I)) * PTCHGX(I,1:NPTCHG(I)) * PTCHGY(I,1:NPTCHG(I)) ) ) 



   Q30  = 0.5d0 * d**3 * ( 5.0d0*dot_product (PTCHGZ(I,1:NPTCHG(I))*PTCHGZ(I,1:NPTCHG(I)) &
          *PTCHGZ(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                       3.0d0 * ( dot_product (PTCHGZ(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I)) &
          *PTCHGX(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) +&
                                 dot_product (PTCHGZ(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)) &
          *PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) +&
                                 dot_product (PTCHGZ(I,1:NPTCHG(I))*PTCHGZ(I,1:NPTCHG(I)) &
          *PTCHGZ(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) ) )

   Q31c = (1.0d0/4.0d0) * sqrt(6.0d0) * d**3 *&
          ( 4.0d0*dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGZ(I,1:NPTCHG(I))*PTCHGZ(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                  dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                  dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) ) 

   Q31s = (1.0d0/4.0d0) * sqrt(6.0d0) * d**3 *&
          ( 4.0d0*dot_product (PTCHGY(I,1:NPTCHG(I))*PTCHGZ(I,1:NPTCHG(I))*PTCHGZ(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                  dot_product (PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                  dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) ) 

   Q32c = 0.5d0 * sqrt(15.0d0) * d**3 * ( dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                                          dot_product (PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) )
  
   Q32s = sqrt(15.0d0) * d**3 * ( sum ( e(i,1:NPTCHG(I)) * PTCHGX(I,1:NPTCHG(I)) * PTCHGY(I,1:NPTCHG(I)) * PTCHGZ(I,1:NPTCHG(I)) ) ) 


   Q33c = (1.0d0/4.0d0) * sqrt(10.0d0) * d**3 *&
                   (dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
              3.0d0*dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) )

   Q33s = (1.0d0/4.0d0) * sqrt(10.0d0) * d**3 *&
          (3.0d0 * dot_product (PTCHGX(I,1:NPTCHG(I))*PTCHGX(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) -&
                   dot_product (PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I))*PTCHGY(I,1:NPTCHG(I)), e(i,1:NPTCHG(I))) ) 

  write(103,'(a,9f18.10)') ATM_SYM(i), Q00, Q10, Q11c, Q11s, Q20, Q21c, Q21s, Q22c, Q22s

  write(107,'(a,8f18.10)') ATM_SYM(i), Q00, Q30, Q31c, Q31s, Q32c, Q32s, Q33c, Q33s
  enddo

!==================================================================================================
! Axes transformation subroutine called                                                           |
! For C6H5Cl: 6 frames are defined (NFRAME = 6). Each frame contains 3 no. of atoms (INUM = 3).   |
! Double counting of atoms per frame is excluded using logical conditions. FLAG_NDC is true. But  |
! when double counting happens, FLAG_NDC is false, i.e., exclude double counting.                 |
!==================================================================================================
! DEFINE ATOMIC GLOBAL AXIS
  gex = (/1.0d0,0.0d0,0.0d0/)
  gey = (/0.0d0,1.0d0,0.0d0/)
  gez = (/0.0d0,0.0d0,1.0d0/)

  open(unit=108,file=arg,IOSTAT=readstat,STATUS='OLD',ACTION='read')
  if(readstat .ne. 0) then
    write(*,*) arg,' cannot be opened'
    stop 1
  end if

  read(108,*,iostat=readstat) resname
  NFRAME = 0
  readframe : do

   NFRAME = NFRAME +1

   read(108,*,iostat=readstat) ANUM(NFRAME,1:3) 

   if (readstat .ne. 0) then
     exit readframe
   endif

   write(105,'(3(i2,1x),a,2x,a,1x,i2)') ANUM(NFRAME,1:3),'BO','! atom indices involved in frame', NFRAME

   do INUM = 1, 3
    FLAG_NDC = .true.    ! TILL HERE READING OF DUPLICATE IS TRUE 

! PRINT 0 (CHARGE) FOR DUPLICATE ATOM INDEX
     do k = 1, NFRAME-1
      do l = 1, 3
       if ( ANUM(NFRAME,INUM) .eq. ANUM(k,l) ) write(105,'("0 0")')  
      enddo
     enddo

! EXCLUDE DUPLICATES
     do k = 1, NFRAME-1
      do l = 1, 3
       if ( ANUM(NFRAME,INUM) .eq. ANUM(k,l) ) FLAG_NDC = .false.
      enddo
     enddo

    if (FLAG_NDC) then   ! IF EXCLUDING DUPLICATES TRUE  

     write(105,'(i1,1x,i1,10x,a,1x,i2)') NPTCHG(ANUM(NFRAME,INUM)),0,'! no. chgs for atom', ANUM(NFRAME,INUM)

     do j = 1, NPTCHG(ANUM(NFRAME,INUM))  !for NPTCHG point charges
     k=ANUM(NFRAME,INUM)
! IN DIAGONALIZED VECTOR BASIS/PRINCIPLE AXES
     px = d * PTCHGX(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),1,1:3)), gex(1:3) ) +&
          d * PTCHGY(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),2,1:3)), gex(1:3) ) +&
          d * PTCHGZ(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),3,1:3)), gex(1:3) )
                                                                  
     py = d * PTCHGX(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),1,1:3)), gey(1:3) ) +&
          d * PTCHGY(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),2,1:3)), gey(1:3) ) +&
          d * PTCHGZ(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),3,1:3)), gey(1:3) )
                                                                  
     pz = d * PTCHGX(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),1,1:3)), gez(1:3) ) +&
          d * PTCHGY(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),2,1:3)), gez(1:3) ) +&
          d * PTCHGZ(k,j) * dot_product( (V3(ANUM(NFRAME,INUM),3,1:3)), gez(1:3) )

     write(109,'(a,4f20.10)') 'X', px, py, pz, e(ANUM(NFRAME,INUM),j) 

! IN GLOBAL MOLECULAR AXES IN BOHR
     gmx = px + CART_COORD(ANUM(NFRAME,INUM),1)
     gmy = py + CART_COORD(ANUM(NFRAME,INUM),2)
     gmz = pz + CART_COORD(ANUM(NFRAME,INUM),3)

     write(104,'(a,4f20.10)') 'X', gmx, gmy, gmz, e(ANUM(NFRAME,INUM),j) 

! CALL NEW AXIS-TRANSFORMATION, & LOCAL AXES PRINTED IN ANGSTROM FOR CHARMM INPUT
     call axis (NATMAX, ANUM(NFRAME,1:3), CART_COORD, EX, EY, EZ)
     lx = dot_product( EX(INUM,1:3), (/px,py,pz/) )* b2a
     ly = dot_product( EY(INUM,1:3), (/px,py,pz/) )* b2a
     lz = dot_product( EZ(INUM,1:3), (/px,py,pz/) )* b2a

!    write(105,'(a,4f20.10)') 'X', lx, ly, lz, e(ANUM(NFRAME,INUM),j) 

     write(105,'(4f20.10)')lx, ly, lz, e(ANUM(NFRAME,INUM),j)
     
     enddo ! do j = 1, NPTCHG

    endif ! if (FLAG_NDC) then

   enddo ! do INUM = 1, 3

  enddo readframe
  

end program solve_eigen

!==================================================================================================
! Axes transformation frames, each containing 3 atoms                         |
! Returns unit vectors along local x, y and z axes for each atom in arrays
! EX,EY,EZ
!==================================================================================================
subroutine axis (NATMAX, ANUM, CART_COORD, EX, EY, EZ)
  implicit none

  integer, intent(in)           :: NATMAX, ANUM(1:3)
  double precision, intent(in)  :: CART_COORD(1:NATMAX,1:3)
  double precision, intent(out) :: EX(1:3,1:3)
  double precision, intent(out) :: EY(1:3,1:3)
  double precision, intent(out) :: EZ(1:3,1:3)

  double precision              :: ATM1(1:3), ATM2(1:3), ATM3(1:3)
  double precision              :: B1X, B1Y, B1Z, B2X, B2Y, B2Z 
  double precision              :: B1(1:3)
  double precision              :: B2(1:3)
  double precision              :: B3(1:3)
  double precision              :: B4(1:3)
  double precision              :: B5(1:3)
  double precision              :: RB1, RB2


 ATM1 = CART_COORD(ANUM(1),1:3)
 ATM2 = CART_COORD(ANUM(2),1:3)
 ATM3 = CART_COORD(ANUM(3),1:3)

!FIRST DEFINE LOCAL Z AXIS ALONG BOND VECTOR B1
!BOND VECTOR B1
 B1X = ATM1(1) - ATM2(1)
 B1Y = ATM1(2) - ATM2(2)
 B1Z = ATM1(3) - ATM2(3)
 RB1 = SQRT(B1X*B1X + B1Y*B1Y + B1Z*B1Z)
 B1 = (/B1X,B1Y,B1Z/)/RB1

!EZ
 EZ(1,1:3) = B1(1:3) 
 EZ(2,1:3) = EZ(1,1:3)

!BOND VECTOR B2
 B2X = ATM3(1) - ATM2(1)
 B2Y = ATM3(2) - ATM2(2)
 B2Z = ATM3(3) - ATM2(3)
 RB2 =  SQRT(B2X*B2X + B2Y*B2Y + B2Z*B2Z)
 B2 = (/B2X,B2Y,B2Z/)/RB2

 EZ(3,1:3) = B2(1:3) 

!EY
!B3 = B1 cross B2, AND NORMALIZE 
 B3  =(/( B1(2)*B2(3) - B1(3)*B2(2) ), (-( B1(1)*B2(3) - B1(3)*B2(1) ) ), ( B1(1)*B2(2) - B1(2)*B2(1) )/)
 B3 =  B3 / SQRT(SUM(B3*B3))
 EY(1,1:3) = B3(1:3)
 EY(2,1:3) = EY(1,1:3)
 EY(3,1:3) = EY(1,1:3)

!EX = EZ cross EY
!B4 = B1 cross B3, AND NORMALIZE
 B4 = (/( B1(2)*B3(3) - B1(3)*B3(2) ), (-( B1(1)*B3(3) - B1(3)*B3(1) ) ), ( B1(1)*B3(2) - B1(2)*B3(1) )/)
 B4 =  B4 / SQRT(SUM(B4*B4))
 EX(1,1:3) = B4(1:3)
 EX(2,1:3) = EX(1,1:3) 

!B5 = B2 cross B3, AND NORMALIZE
 B5 = (/( B2(2)*B3(3) - B2(3)*B3(2) ), (-( B2(1)*B3(3) - B2(3)*B3(1) ) ), ( B2(1)*B3(2) - B2(2)*B3(1) )/) 
 B5 =  B5 / SQRT(SUM(B5*B5))
 EX(3,1:3) = B5(1:3)

return
end  
  

subroutine eigen(N, MAT, EVAL, EVEC)
  !=======================================================================================
  !
  !  Interface to Lapack 'dsyev'
  !
  !=======================================================================================
  implicit none

  !=== Input, Output data
  integer, intent(in)           :: N
  double precision, intent(in)  :: MAT(N,N)
  double precision, intent(out) :: EVAL(N)
  double precision, intent(out) :: EVEC(N,N)

  !=== Lapack data
  integer                       :: info
  double precision              :: work(3*N)
  character, parameter          :: jobz="V", uplo="U"

  !=== Execution

  EVEC = MAT

  call dsyev(jobz, uplo, N, EVEC, N, EVAL, work, 3*N, info)

  !=======================================================================================
end subroutine eigen


subroutine Jacobi(n,oria,eval,x)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer i, j, k, n
double precision a(n,n),x(n,n)
double precision eval(n), oria(n,n)
double precision abserr, b2, bar
double precision beta, coeff, c, s, cs, sc

abserr = 5.0d-15

a = oria

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.0d0
do i=1,n
  x(i,i) = 1.0d0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0d0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) then
  do i=1,n
    eval(i)=a(i,i)
  end do
  return
endif

! average for off-diagonal elements /2
bar = 0.5d0*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0d0*a(j,i)**2
      bar = 0.5d0*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0d0*a(j,i))
      coeff = 0.5d0*beta/sqrt(1.0d0+beta**2)
      s = sqrt(max(0.5d0+coeff,0.0d0))
      c = sqrt(max(0.5d0-coeff,0.0d0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do

do i = 1, n
  eval(i) = a(i,i)
enddo
return
end


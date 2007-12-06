!==========================================================================
SUBROUTINE get_displ( latt, symm, iprint )
  !==========================================================================
  !
  ! This routine generates a set of displacements, if you calculate
  ! the forces on the cell induced by these displacements than these
  ! forces can be used to construct the force constant matrix.

  USE nrutil
  USE data

  IMPLICIT NONE

  LOGICAL :: found, move, equiva1

  INTEGER :: ndispl, conta, na, naa, nb, isym, natoms, iprint

  REAL(DP) :: dx(3,3), tmp1(3), tmp2(3), at(3,3), bg(3,3)
  REAL(DP), ALLOCATABLE :: x(:,:), xtmp(:,:), lattxtmp(:,:)

  TYPE(lattice)  :: latt
  TYPE(symmetry) :: symm

  open(1,file='SPOSCAR',status='old')
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*) 
  read(1,*) 

  natoms = sum( latt%nions * product( latt%ndim ) )

  allocate( x(3,natoms), xtmp(3,natoms), lattxtmp(3,latt%natoms) )

  do na = 1, natoms
     read(1,*) xtmp(:,na)
  enddo

  at(:,1) = latt%at(:,1) * latt%ndim(1)
  at(:,2) = latt%at(:,2) * latt%ndim(2)
  at(:,3) = latt%at(:,3) * latt%ndim(3)

  call inve( at, bg )
  bg = transpose(bg)

  xtmp = matmul( at, xtmp ) * latt%scale 

  lattxtmp = matmul( latt%at, latt%x ) * latt%scale 

  ndispl = 0
  conta = 0
  
  if( latt%dxstart(2) /= 0._dp .or. latt%dxstart(3) /= 0._dp ) then
     print'(/''------------------------------------------------------------'')'
     print'(''First displacement provided by user: '',3f6.2)',latt%dxstart
  endif

  print'(''The following displacements are needed: (also written in file DISP)'')'
  print'(''(Displacements are in units of lattice vectors)'')'
  open(17,file='DISP',status='unknown')
  !==============================================================================
  ! For each atom in the primitive unit cell check if it is necessary to move it
  !==============================================================================
  do na = 1, latt%natoms
     x = xtmp
     latt%x = lattxtmp
     symm%is = symm%isstart

     !------------------------------------------------------------------------------
     ! is it possible that, with a symmetry operation, I send this atom to an atom
     ! I have already moved?
     !------------------------------------------------------------------------------
     move = .true.
     do_nb: do nb = 1, na-1
        do isym = 1, symm%nsymstart

           !-------------------------------------------
           ! Yes. There is no need to move the atom
           !-------------------------------------------
           if( nb == symm%irt(na,isym) ) then
              move = .false.
              exit do_nb
           endif
        enddo
     enddo do_nb
     !---------------------------------------------
     ! No, move the atom
     !---------------------------------------------
     if( move ) then

        conta = conta + 1

        do nb = 1, natoms
!
! some compilers might not like this, it is safer to pass the vector as in the following line
!
!           if( equiva1( latt%x(:,na), x(:,nb) ) ) exit
           if( equiva1( latt%x(1,na), x(1,nb) ) ) exit
        enddo

        do naa = 1, latt%natoms
           latt%x(:,naa) = latt%x(:,naa) - lattxtmp(:,na)
        enddo
        latt%x = matmul( transpose(latt%bg), latt%x ) / latt%scale
        print'(/''-------------------------------------------------------------'')'
        call sgama1 ( symm, latt, iprint )
        latt%x = matmul( latt%at, latt%x ) * latt%scale

        print'(''Move atom # '',i6,3f12.6,/18x3f12.6/,''Along direction(s):'')', &
             &           nb, matmul( transpose(bg), x(:,nb) ) / latt%scale, x(:,nb)

        dx = 0._dp
        dx(:,1) = latt%dxstart
        dx(:,1) = matmul( transpose(latt%bgs), dx(:,1) )
        print'(''First, Direct: '',3f7.3,10x,''Cartesian:'',3f7.3)', &
             dx(:,1)/vabs(dx(:,1)), matmul( latt%ats, dx(:,1) ) / vabs(matmul( latt%ats, dx(:,1) ))
        write(17,'('' " '',i3,3f12.8,'' " '', ''\ '')') nb,dx(:,1)/vabs(dx(:,1))/latt%disp
        dx(:,1) = matmul( latt%ats, dx(:,1) )
        ndispl = ndispl + 1

        !---------------------------------------------------
        ! second displacement
        !---------------------------------------------------
        do isym = 1, symm%nsym
           dx(:,2) = matmul( symm%is(:,:,isym), dx(:,1) )
           !----------------------------------------------------------------
           ! is dx(:,2) linearly independent from dx(:,1) ?
           !----------------------------------------------------------------
           call vecprod( dx(:,1), dx(:,2), tmp1 )
           if( vabs( tmp1 ) > 1.d-8 ) then
              found = .true.
              exit
           else
              dx(:,2) = 0
              found = .false.
           endif
        enddo
        if(.not. found ) then
           dx(2,2) = 1.0_dp
           dx(:,2) = matmul( transpose(latt%bgs), dx(:,2) )
           print'(''Second, Direct: '',3f7.3,10x,''Cartesian:'',3f7.3)', &
                dx(:,2)/vabs(dx(:,2)), matmul( latt%ats, dx(:,2) )/vabs(matmul( latt%ats, dx(:,2)))
           write(17,'('' " '',i3,3f12.8,'' " '', '' \ '')') nb,dx(:,2)/vabs(dx(:,2))/latt%disp
           dx(:,2) = matmul( latt%ats, dx(:,2) )
           ndispl = ndispl + 1
        endif

        !------------------------------------------------------------------
        ! Third displacement
        !------------------------------------------------------------------
        do isym = 1, symm%nsym
           dx(:,3) = matmul( symm%is(:,:,isym), dx(:,2) )
           !----------------------------------------------------------------
           ! is dx(:,3) linearly independent from dx(:,1) and dx(:,2) ?
           !----------------------------------------------------------------
           call vecprod( dx(:,1), dx(:,2), tmp1 )
           tmp2 = dot_product( tmp1, dx(:,3) )
           if( vabs( tmp2 ) > 1.d-8 ) then
              found=.true.
              exit
           else
              dx(:,3) = 0
              found=.false.
           endif
        enddo
        if( .not. found ) then
           dx(3,3) = 1.0_dp
           call vecprod( dx(:,1), dx(:,2), tmp1 )
           tmp2 = dot_product( tmp1, dx(:,3) )
           if( vabs( tmp2 ) < 1.d-8 ) then
              dx(3,3) = 0._dp
              dx(2,3) = 1._dp
           endif
           call vecprod( dx(:,1), dx(:,2), tmp1 )
           tmp2 = dot_product( tmp1, dx(:,3) )
           if( vabs( tmp2 ) < 1.d-8 ) then
              dx(2,3) = 0._dp
              dx(1,3) = 1._dp
           endif
           call vecprod( dx(:,1), dx(:,2), tmp1 )
           tmp2 = dot_product( tmp1, dx(:,3) )
           if( vabs( tmp2 ) < 1.d-8 ) then
              stop 'error, cannot find third displacement'
           endif
           dx(:,3) = matmul( transpose(latt%bgs), dx(:,3) )
           print'(''Third, Direct: '',3f7.3,10x,''Cartesian:'',3f7.3)', &
                dx(:,3)/vabs(dx(:,3)), matmul( latt%ats, dx(:,3)) / vabs(matmul( latt%ats, dx(:,3)))
           write(17,'('' " '',i3,3f12.8,'' " '', '' \ '')') nb,dx(:,3)/vabs(dx(:,3))/latt%disp
           dx(:,3) = matmul( latt%ats, dx(:,3) )
           ndispl = ndispl + 1
        endif
     endif

  enddo

  close(17)
  write(*,'(/''Total number of displaced atoms = '',i4,/, &
       &           ''Total number of displacements   = '',i4)') conta, ndispl

END SUBROUTINE get_displ







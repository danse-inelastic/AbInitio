!==========================================
! Generate a set of Monkhorst-Pack points 
!==========================================
SUBROUTINE generate_punti( symm, latt, dx, dy, dz )

  USE nrtype
  USE nrutil
  USE data
  IMPLICIT NONE

  LOGICAL :: equiva, lgamma
  LOGICAL, ALLOCATABLE :: mask(:)
  INTEGER ::  npoints, i, j, isym, qa, qb, qc
  REAL(DP) :: bg(3,3), at(3,3), r, p, s, x(3), xm(3), xp(3)
  REAL(DP) :: dx, dy, dz, b
  REAL(DP), PARAMETER :: eps = 1.e-6
  REAL(DP), ALLOCATABLE :: peso(:), kt(:,:), is(:,:,:)


  TYPE(symmetry) :: symm
  TYPE(lattice) :: latt
  
  qa = latt%qa
  qb = latt%qb
  qc = latt%qc
  bg = latt%bg
  lgamma = latt%lgamma

  if( qa < 0 ) return

!---------------------------------------------------------------------------------
! If you don't want to go to the edge of the BZ set dx, and/or dy and/or dz < 1.0
! Defaults are dx=dy=dz=1.0 (normal behaviour)
!---------------------------------------------------------------------------------

  print'(/''Generating IBZ points ....'',3i5)',qa,qb,qc

  npoints = qa * qb * qc
  allocate( mask(npoints), peso(npoints), kt(3,npoints), is(3,3,symm%nsym) )
  kt = 0
!----------------------------
! generate points
!----------------------------
  i = 0
  if(.not.lgamma)then                           ! reticolo generale
     do r = 1, qa
        do p = 1, qb
           do s = 1, qc
              i = i + 1
              kt(1,i) = &
                   ( 2*r - qa - 1 ) /2 / qa * bg(1,1) + &
                   ( 2*p - qa - 1 ) /2 / qa * bg(1,2) + &
                   ( 2*s - qa - 1 ) /2 / qa * bg(1,3)
              if( dx < 0 ) then
                 kt(1,i) = (0.5 - kt(1,i) * dx)
              else
                 kt(1,i) = kt(1,i) * dx
              endif
              kt(2,i) = &
                   ( 2*r - qb - 1 ) /2 / qb * bg(2,1) + &
                   ( 2*p - qb - 1 ) /2 / qb * bg(2,2) + &
                   ( 2*s - qb - 1 ) /2 / qb * bg(2,3)
              if( dy < 0 ) then
                 kt(2,i) = (0.5 - kt(2,i) * dy)
              else
                 kt(2,i) = kt(2,i) * dy
              endif
              kt(3,i) = &
                   ( 2*r - qc - 1 ) /2 / qc * bg(3,1) + &
                   ( 2*p - qc - 1 ) /2 / qc * bg(3,2) + &
                   ( 2*s - qc - 1 ) /2 / qc * bg(3,3)
              if( dz < 0 ) then
                 kt(3,i) = (0.5 - kt(3,i) * dz)
              else
                 kt(3,i) = kt(3,i) * dz
              endif
           enddo
        enddo
     enddo
  else
     do r = 1, qa
        do p = 1, qb
           do s = 1, qc
              i = i + 1
              kt(1,i) = &
                   ( r - 1 ) / qa * bg(1,1) + &
                   ( p - 1 ) / qa * bg(1,2) + &
                   ( s - 1 ) / qa * bg(1,3)
              kt(2,i) = &
                   ( r - 1 ) / qb * bg(2,1) + &
                   ( p - 1 ) / qb * bg(2,2) + &
                   ( s - 1 ) / qb * bg(2,3)
              kt(3,i) = &
                   ( r - 1 ) / qc * bg(3,1) + &
                   ( p - 1 ) / qc * bg(3,2) + &
                   ( s - 1 ) / qc * bg(3,3)
           enddo
        enddo
     enddo
  endif

  call inve( bg, at )
  at = transpose(at)

  kt = matmul( transpose(at), kt )
  kt = kt - aint(2*kt)

  if( .not. symm%lsymm ) then
     mask=.true.
     peso = 1
     goto 10
  endif

  do i=1,symm%nsym
     is(:,:,i) = matmul( transpose(bg), matmul(symm%iscryst(:,:,i), at ) )
     is(:,:,i) = transpose( is(:,:,i) )
  enddo
  
  
 !----------------------------
 ! now reduce them to the IBZ
 !----------------------------
  mask = .true.
  do i = 1, npoints
     peso(i) = 1
     if( mask(i) ) then
        do j = i+1, npoints
           if( mask(j) ) then
              dosym: do isym = 1, symm%nsym
                 x(1) = is(1,1,isym) * kt(1,j) + is(1,2,isym) * kt(2,j) + is(1,3,isym) * kt(3,j) 
                 x(2) = is(2,1,isym) * kt(1,j) + is(2,2,isym) * kt(2,j) + is(2,3,isym) * kt(3,j) 
                 x(3) = is(3,1,isym) * kt(1,j) + is(3,2,isym) * kt(2,j) + is(3,3,isym) * kt(3,j) 
                 equiva = &
                      ( abs( x(1) - kt(1,i)  - nint( x(1) - kt(1,i) ) ) < eps .and. &
                      abs( x(2) - kt(2,i)  - nint( x(2) - kt(2,i) ) ) < eps .and. &
                      abs( x(3) - kt(3,i)  - nint( x(3) - kt(3,i) ) ) < eps )
  !------------------------------------------------------------------------------------
  ! keep only K, eliminate -K (the dynamical matrix has equal eigenvalues at K and -K)
  !------------------------------------------------------------------------------------
                 if( .not. equiva ) equiva = &
                      ( abs( x(1) + kt(1,i)  - nint( x(1) + kt(1,i) ) ) < eps .and. &
                      abs( x(2) + kt(2,i)  - nint( x(2) + kt(2,i) ) ) < eps .and. &
                      abs( x(3) + kt(3,i)  - nint( x(3) + kt(3,i) ) ) < eps ) 
                 if(equiva)then
                    mask(j) = .false.
                    peso(i) =  peso(i) + 1
                    exit dosym
                 endif
              enddo dosym
           endif
        enddo
     endif
  enddo
  
10 continue
  do i = 1, npoints
     if( mask(i) ) then
        do j = i+1, npoints
           if( mask(j) ) then
              equiva = &
                   ( abs( kt(1,j) + kt(1,i)  - nint( kt(1,j) + kt(1,i) ) ) < eps .and. &
                   abs( kt(2,j) + kt(2,i)  - nint( kt(2,j) + kt(2,i) ) ) < eps .and. &
                   abs( kt(3,j) + kt(3,i)  - nint( kt(3,j) + kt(3,i) ) ) < eps )
              if(equiva)then
                 mask(j) = .false.
                 peso(i) =  peso(i) + 1
              endif
           endif
        enddo
     endif
  enddo

  write(*,'(/''Writing in file QPOINTS '',i6,'' MP points, from the grid '',3i4/)') &
       count(mask), qa, qb, qc
  if(lgamma) write(*,'(''The grid passes through GAMMA'')')
  open(1,file='QPOINTS', status='unknown')
  write(1,*) count(mask)
  do i = 1, npoints
     if( mask(i) ) then
        if( abs( peso(i) - nint(peso(i)) ) > eps  ) &
             stop 'something wrong in MP generation, non integer weights, stopping...'
        write(1,'(3f20.16,f10.6)') kt(:,i), peso(i)
     endif
  enddo
  close(1)
  
END SUBROUTINE generate_punti

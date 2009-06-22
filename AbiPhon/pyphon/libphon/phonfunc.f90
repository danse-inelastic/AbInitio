!=========================================================================
!
!  Program PHON: 
!
!
! Force constant matrix calculation program. Uses finite displacements.
!
!
! Dario Alfe`,      November  1998
!
! This program is freely distributed and comes with no warranty. Please
! send any comments to the author (d.alfe@ucl.ac.uk)
! 
! If you use this code to publish scientific results the author would appreciate a
! citation as follows:
! "D. Alfe`, 1998, program available at http://chianti.geol.ucl.ac.uk/~dario"  
!
! 
! 26/2/2003  added partial density of states (USETHIS)   
! 21/8/2003  minor modification in 'numeric_force', which didn't compile
!            properly with the linux-pgf90 compiler
! 30/9/2003  added translational invariance symmetrization
! 16/12/2003 added makefile for linux ifc compiler
! 15/07/2004 removed generation of SPOSCART file (not needed anyway)
! 28/09/2005 corrected DXSTART; increased default of RMAX
! 03/03/2006 increased default of RMAX in the right routine!
! 12/03/2006 Drift in forces is not removed (for non periodic calculations)
! 18/03/2006 Corrected error in generate_punti
! 16/04/2006 Problem with symmetry in numeric_force hopefully fixed
!            (check performed in crystal coordinates rather than carthesian)
!
! subroutinized 25/09/2006 patrickh 
!========================================================================= 

subroutine phonfunc

  USE nrtype
  USE nrutil
  USE diagonalize
  USE data

  IMPLICIT NONE

  CHARACTER(1) :: fl
  CHARACTER(2) :: fl1
  CHARACTER(7) :: dum
  CHARACTER(80), PARAMETER :: version='  PHON, VERSION 1.09 (15/04/2006)'

  LOGICAL :: linverse, lforceout, lrecip, esiste, equiva1, &
       lgamma = .false.

  INTEGER :: nbranches, nd, nq, iprint, ndos, j, k, kmax, ndiff
  INTEGER :: i, i1, i2, na, nb, nrr, nr, ntot, idos, ntyp, conta, nti

  INTEGER, ALLOCATABLE :: list(:,:)

  REAL(DP), ALLOCATABLE :: dos(:)

  REAL(DP) :: temp(3), omegam2, peso, pesotot, pesopart, delta, slope(3), &
       alpha, a, temperature, q(3), dq(3), zeta, tmp1(3), q0, e, dosin, dosend, &
       dosstep, dossmear, eint, tv, tc, tv0, tc0, omega, fakt, free, &
       free_q, rcut, x, cv, zeropoint, fac
  
  REAL(DP) :: dx, dy, dz

  COMPLEX(DP) :: zdotc

  REAL(DP), PARAMETER :: zero = 1.d-6, convert_thz_to_cm = 33.357, convert_thz_to_meV = 4.1357    

  TYPE(lattice)  :: latt
  TYPE(shell)    :: sh
  TYPE(dynmat)   :: dyn
  TYPE(symmetry) :: symm

  call vtime( tv0, tc0 )

  !==========================================================================
  ! read in run time variables
  !==========================================================================

  write(*,*)'----------------------------------------------------------------'
  write(*,*) version
  write(*,*)'----------------------------------------------------------------'
  call reader( latt, dyn, temperature, lforceout, linverse, alpha, a,  lrecip, sh%rmax, &
       sh%cutoff, dosin, dosend, dosstep, dossmear, &
       nd, iprint, symm%lsymm, symm%lfullqgrid, dx, dy, dz, nti )

  fac = evtoj   &        !  D^2 U eV  --> Joule
       /1.d-20  &        !  D^2 R A^2 --> meter 
       /amtokg  &        !  proton mass --> Kg
       /1.d24/twopi**2   !  cicles in THZ^2

  ndos = (dosend - dosin)/dosstep + 1
  allocate( dos(ndos) )
  dos = 0

  !======================================================================
  ! set lattice
  !=====================================================================
  call set_lattice ( latt, symm, iprint, sh%rmax )
  call generate_punti( symm, latt, dx, dy, dz )


  !---------------------------
  ! Number of phonon branches
  !---------------------------
  nbranches = latt%natoms * 3

  !==================================================
  ! For the inverse power problem, ignore otherwise
  !==================================================
  a = 4*a**alpha
  !---------------------------------------
  ! set the adimensional variable zeta
  !---------------------------------------
  if( dyn%lfree ) then
     zeta = a / exp( alpha * log( latt%volume / latt%natoms ) / 3._dp ) / twopi / ( temperature * bolkev )
     if(linverse)print*,'Zeta = ',zeta
  endif
  !=======================================================================

  !===========================
  ! Build a supercell ?
  !===========================
  if(latt%lsuper)then
     call build_supercell( latt )
     call get_displ( latt, symm, iprint )
!     stop
     return
  endif

  !=============================================================================
  ! construct the shell of vectors in real space
  !=============================================================================
  call rshell( sh, latt ) 

  if( linverse )then
!-----------------------------------
! inverse power problem (analytical)
     call force_constant( symm, dyn, sh, latt, a, alpha )
!-----------------------------------
  else
!---------------------------------------------
! numeric forces (small displacement method)
!---------------------------------------------
     call numeric_force( latt, symm, sh, dyn, iprint, nti )
  endif

  !--------------------------
  ! set the list of vectors
  !--------------------------
  ntot = 0
  do nrr=1,sh%nrm
     do na=1,latt%natoms
        do nb=1,latt%natoms
           tmp1 = sh%r(:,nrr) + latt%x(:,nb) - latt%x(:,na)
           if( vabs(tmp1) < sh%rmax ) ntot = ntot + 1
        enddo
     enddo
  enddo
  allocate( list(3,ntot) )
  ntot = 0
  do nrr=1,sh%nrm
     do na=1,latt%natoms
        do nb=1,latt%natoms
           tmp1 = sh%r(:,nrr) + latt%x(:,nb) - latt%x(:,na)
           if( vabs(tmp1) < sh%rmax ) then
              ntot = ntot + 1
              list(1,ntot) = na
              list(2,ntot) = nb
              list(3,ntot) = nrr
           endif
        enddo
     enddo
  enddo

  !==========================================================
  ! Write the force constant matrix in the file HARMONIC ?
  !==========================================================
  if(lforceout)then
     open(1,file='HARMONIC',status='unknown')
     write(1,*) sh%rmax,'   cutoff radius'
     write(1,*) ntot,'   number of vectors'
     do i = 1, ntot
        na = list(1,i)
        nb = list(2,i)
        nrr = list(3,i)
        write(1,'(''vector:'',4i20)') na, nb, nrr, dyn%weight(na,nb,nrr)
        write(1,'(3f20.12)') sh%r(:,nrr) + latt%x(:,nb) - latt%x(:,na)
        write(1,'(''force constant matrix:'')')
        write(1,'(3f20.12)') dyn%dmat(:,:,na,nb,nrr)
     enddo
     write(*,'(/''Force constant matrix file HARMONIC written''/)')
  endif

  !==========================================================================
  ! calculate the free energy as the sum of frequencies, you must provide
  ! the file 'QPOINTS' which contains the special points of the BZ 
  ! IN RECIPROCAL LATTICE COORDINATES
  !=========================================================================

!xx
  return
  if( dyn%lfree ) then

     inquire(file='QPOINTS',exist=esiste)
     if(.not.esiste) stop 'cannot find file QPOINTS, use QA, QB and QC to generate one'
     open(10,file='QPOINTS',status='old')

     read(10,*) nd
     print'(/''Integrating frequencies...'')'
     print'(/''Using '',i6,'' q-points from file QPOINTS (in reciprocal space coordinates)''/)',nd
     pesotot = 0; omegam2 = 0; eint = 0
     omega = 0; free_q = 0; zeropoint = 0
     free = 0; cv = 0
     conta = 0
     do nq = 1, nd
        read(10,*) q, peso
        pesotot = pesotot + peso
     enddo
     rewind(10)
     read(10,*) nd
     do nq = 1, nd
        read(10,*) q, peso
        q = twopi/latt%scale*matmul(latt%bg,q) ! in cartesian coordinates

        !==========================================================================
        ! construct the dynamical matrix
        !==========================================================================
        call dyn_mat( dyn, latt, sh, q, ntot, list )

        if(iprint>2)then
           print'(''Dynamical matrix:'')'
           do j=1,latt%natoms
              do k=1,latt%natoms
                 print'(2i3)',j,k
                 do i1=1,3
                    print'(6f13.8)', ( dyn%cmat( i1 + (j-1)*3, i2 + (k-1)*3 ), i2=1,3) 
                 enddo
              enddo
           enddo
        endif

        !==========================================================================
        ! diagonalize the dynamical matrix
        !==========================================================================
        call diagh( dyn%cmat, dyn%eig, 'L' )

        if(iprint>2)then
           do j = 1, nbranches
              print'(''Eigenvalue '',i3)',j
              print'(f20.12)',dyn%eig(j) 
              print'(''Eigenvector '',i3)',j
              do k=1,latt%natoms
                 print'(''atom '',i3)',k
                 do i1=1,3
                    print'(2f20.12)', dyn%cmat(i1+(k-1)*3,j)
                 enddo
              enddo
           enddo
        endif

        !=====================
        ! density of states
        !=====================
        do i1 = 1, nbranches
           pesopart = 0
           do k = 1, latt%natoms
              if(dyn%usethis(latt%ityp(k))) pesopart = pesopart + &
                   abs(zdotc(3,dyn%cmat(1+(k-1)*3,i1),1, &
                   dyn%cmat(1+(k-1)*3,i1),1))
           enddo
           idos = 0
           kmax = (dosend - dosin)/dosstep
           do k = 1, kmax + 1
              e = dosin + (k-1)*dosstep + dosstep/2
              idos = idos + 1
              if( dyn%eig(i1) >= 0 )then
                 if( sqrt(dyn%eig(i1)) > e - dosstep/2 .and. sqrt(dyn%eig(i1))< e + dosstep/2 )then
                    dos(idos) = dos(idos) + peso * pesopart
                 endif
              endif
           enddo

           !==========================================================================
           ! the phonon frequencies are the squares of the matrix eigenvalues, for the
           ! inverse power problem I need
           ! the adimensional constants ci = m * omega^2 * V^[(alpha+2)/3] / A (V is
           ! the volume per atom) 
           !==========================================================================
           !-----------------------------------------------------------------------
           ! exclude the acustic mode at gamma 
           !-----------------------------------------------------------------------
           if( dyn%eig(i1) < -zero ) then
              print*,'negative frequencies',dyn%eig(i1),i1
              !                  stop 'negative frequencies'
              !----------------------------------------------------------
              ! omegam2 is the average logaritmic square frequency.
              !
              ! nu = sqrt( dyn%eig(i) ) in THZ
              ! omega = 2pi * nu
              !
              ! I want m*omega^2*V^(2/3) to be an energy in eV:
              ! 
              ! omega^2 = omega^2 * 1.d24  and now it is in Rad^2 / sec
              !
              ! m = m * amtokg   in Kg
              !
              ! The volume will be expressed in A^3, so I multiply per 1.d-20,
              ! so that V^(3/2) will be in meter^(3/2).
              !
              ! Now m*omega^2*V^(2/3) is in Joule, I divide by 1.6d-19.
              !
              ! 1.d-20 / 1.6d-19 =  0.06241460122_dp
              !  
              !--------------------------------------------------------------

           elseif( dyn%eig(i1) > zero ) then

              conta = conta + 1
              fakt = 1.d4 * twopi**2 / evtoj 

              omega = omega + peso * pesopart * log( dyn%eig(i1) * fakt )

              free = free + peso * pesopart * bolkev * temperature * &
                   log( hplank * sqrt( dyn%eig(i1) )*1.d12 / temperature / bolk )

              free_q = free_q  + peso * pesopart * ( &
                   0.5_dp * hplank * sqrt( dyn%eig(i1) )*1.d12 / evtoj + &
                   bolkev * temperature * &
                   log( 1 - exp( -hplank * sqrt( dyn%eig(i1) )*1.d12 / temperature / bolk) ) )

              zeropoint = zeropoint + peso * pesopart * &
                   0.5_dp * hplank * sqrt( dyn%eig(i1) )*1.d12 / evtoj

              omegam2 = omegam2 + peso * pesopart * ( log( dyn%eig(i1)*fakt ) + log( latt%mass(1)*amtokg ) )

              eint = eint + peso * pesopart * (hplank/2)*sqrt( dyn%eig(i1) )*1.d12 * &
                   coth( hplank*sqrt( dyn%eig(i1) )*1.d12 / &
                   ( 2*bolk*Temperature ) ) 
              
              x = hplank*sqrt( dyn%eig(i1) )*1.d12 / ( bolk*Temperature ) 
              cv = cv + peso * pesopart * x**2*exp(x)/(exp(x) - 1)**2

           endif
        enddo
     enddo

     omegam2 = omegam2 / pesotot / nbranches
     omega = omega / pesotot / nbranches
     free = free / pesotot
     free_q = free_q / pesotot
     zeropoint = zeropoint / pesotot
     cv = cv / pesotot
     eint = eint / pesotot / evtoj

     if( conta /= nbranches*nd ) &
          write(*,'(/''WARNING, Found '',i8, &
          & '' frequencies different from zero, out of a total of '',i8/)') conta, nbranches*nd

     !---------------------------------------------------------------------------
     ! print out (excess) free energy(classical limit) per atom for the solid
     !---------------------------------------------------------------------------
     write(6,'(/''Your primitive cell contains N ='',i3,'' atoms''/)') latt%natoms
     write(6,'(''Temperature  '',f10.2,'' K'')') temperature
     write(6,'(''Zero point energy           '',f15.5,''  (eV/cell)'')') zeropoint
     write(6,'(''Free energy                 '',f15.5,''  (eV/cell)'', &
          & f12.5,''  (eV/atom)'')') free_q, free_q/latt%natoms 
     write(6,'(''Free energy(classical_limit)'',f15.5,''  (eV/cell)'', &
          & f12.5,''  (eV/atom)'')') free, free/latt%natoms
     write(6,'(''Internal Energy             '',f15.5,''  (eV/cell)'')') eint
     write(6,'(''3 N KT                      '',f15.5,''  (eV/cell)'')') &
          3*temperature*bolkev*latt%natoms
     write(6,'(''Cv                          '',f15.5,''  (K/cell)'')') cv
     write(6,'(''S                           '',f15.5,''  (K/cell)''/)') &
          (eint - free_q) / (temperature*bolkev)

     write(6,'(''Omegam2:'',f15.5)') omegam2
     write(6,'(''Omegabar:'',f15.5)') (omegam2-log(latt%mass(1)*amtokg)-log(1.d0/1.d20/evtoj))/2



     open(4,file='DOS',status='unknown')
     open(44,file='DOS.cm',status='unknown')
     open(45,file='DOS.meV',status='unknown')
     idos = 0
     dos = dos / pesotot / nbranches / dosstep
     call smooth( dos, dossmear, ndos, dosin, dosstep )
     write(6,'(/''Density of states integral = '',f10.5)') sum(dos)*dosstep
     if( abs(sum(dos)*dosstep-1.0_dp) > 1.d-6 ) then
        print*,'renormalizing DOS..'
        dos = dos / ( sum(dos)*dosstep )
     endif        

     kmax = (dosend - dosin)/dosstep
     do k = 1, kmax + 1
        e = dosin + (k-1)*dosstep + dosstep/2
        idos = idos + 1
        write(4,'(2f20.5)') e, dos(idos)
        write(44,'(2f20.5)') e*convert_thz_to_cm, dos(idos)/convert_thz_to_cm
        write(45,'(2f20.5)') e*convert_thz_to_mev, dos(idos)/convert_thz_to_mev
     enddo

     call vtime( tv, tc )
     write(*,'(/''CPU time '',f10.2,10x,''Elapsed time '',f10.2)') tv - tv0, tc - tc0

!    stop     
     return

  endif

  !============================================================================
  ! calculate phonon dispersions
  !===========================================================================

  do i=1,nbranches/3
     if( i <=9 ) then
        write(fl,'(i1)')i
        open(10+i,file='FREQ'//fl,status='unknown')
     else
        write(fl1,'(i2)')i
        open(10+i,file='FREQ'//fl1,status='unknown')
     endif
  enddo
  open(2,file='FREQ',status='unknown')
  open(3,file='FREQ.cm',status='unknown')
  q0=0
  write(*,'(/''Distances:'')')
  do nq = 1, nd
     print'(''point '',i4,f10.6,10x,3f10.6)', nq, q0, dyn%qi(:,nq)
     if( dyn%npoints > 1 ) then
        dq = ( dyn%qf(:,nq) - dyn%qi(:,nq) ) / ( dyn%npoints - 1 )      
     else
        dq = 0
     endif
     if( lrecip ) dq = matmul(latt%bg,dq)
     q0 = q0 - vabs( dq )
     do i = 1, dyn%npoints
        q = dyn%qi(:,nq) 
        if( lrecip ) then
           q = twopi/latt%scale*matmul(latt%bg,q) + &
                (i-1)*dq*twopi/latt%scale
        else
           q = twopi/latt%scale*q + (i-1)*dq*twopi/latt%scale
        endif
        if(iprint > 1)then
           write(*,'(/''Point #'',i3,'' q-vector: '',3f10.5)') i, q*latt%scale/twopi
        endif
        q0 = q0 + vabs(dq)

        !==========================================================================
        ! construct the dynamical matrix
        !==========================================================================
        call dyn_mat( dyn, latt, sh, q, ntot, list )

        if(iprint>2)then
           print'(''Dynamical matrix:'')'
           do j=1,latt%natoms
              do k=1,latt%natoms
                 print'(2i3)',j,k
                 do i1=1,3
                    print'(6f13.8)', &
                         ( dyn%cmat( i1 + (j-1)*3, i2 + (k-1)*3 )/fac, i2=1,3) 
                 enddo
              enddo
           enddo
        endif
        !==========================================================================
        ! diagonalize the dynamical matrix
        !==========================================================================
        call diagh( dyn%cmat, dyn%eig, 'L' )


        do i1=1,nbranches

           !==========================================================================
           ! the phonon frequencies are the square roots of the matrix eigenvalues
           !==========================================================================
           if( dyn%eig(i1) < 0 ) then
              dyn%eig(i1) = -sqrt(-dyn%eig(i1))
           else
              dyn%eig(i1) = sqrt(dyn%eig(i1))
           endif
        enddo
        write(2,'(f12.8,42f12.5)')q0,dyn%eig
        write(3,'(f12.8,42f12.5)')q0,dyn%eig*convert_thz_to_cm
        do i1=1,nbranches/3
           write(10+i1,'(4f12.5)')q0,dyn%eig( 1 + (i1 - 1)*3 : i1*3 )
        enddo

        !--------------------------------------------------
        ! estimates the slopes of the first three branches
        !--------------------------------------------------
        if( iprint > 0 ) then
           if( i == 1 .or. i == dyn%npoints - 1 ) then
              slope = dyn%eig(1:3)
           elseif( i == 2 .or. i == dyn%npoints ) then
              slope = ( dyn%eig(1:3) - slope ) / vabs( dq/latt%scale ) * 100
              write(*,'(''Slopes (m/s) : '',3f10.1)') slope
           endif
        endif

        !=====================
        ! density of states
        !=====================
        do i1 = 1, nbranches
           idos = 0
           do k = 1, ndos
              e = dosin + (k-1)*dosstep
              idos = idos + 1
              if( dyn%eig(i1) >= 0 )then
                 if( dyn%eig(i1) > e - dosstep/2 .and. dyn%eig(i1) < e + dosstep/2 )then
                    dos(idos) = dos(idos) + 1
                 endif
              endif
           enddo
        enddo

        if(iprint>2)then
           do j=1,nbranches
              Print'(''Eigenvalue '',i3,'' : '',f12.5)',j,dyn%eig(j) 
              print'(''Eigenvector '',i3)',j
              do k=1,latt%natoms
                 print'(''atom '',i3)',k
                 do i1=1,3
                    print'(2f10.4)', dyn%cmat(i1+(k-1)*3,j)
                 enddo
              enddo
           enddo
        endif
     enddo
  enddo

  print'(''end point '',f10.6,10x,3f10.6)',  q0, dyn%qf(:,nd)

  open(4,file='DOS',status='unknown')
  open(44,file='DOS.cm',status='unknown')
  open(45,file='DOS.meV',status='unknown')
  idos = 0
  dos = dos / nd / dyn%npoints / nbranches / dosstep
  call smooth( dos, dossmear, ndos, dosin, dosstep )
  write(6,'(/''Density of states integral = '',f10.5)') sum(dos)*dosstep
  if( abs(sum(dos)*dosstep-1.0_dp) > 1.d-6 ) then
     print*,'renormalizing DOS..'
     dos = dos / ( sum(dos)*dosstep )
  endif
  do k = 1, ndos
     e = dosin + (k-1)*dosstep
     idos = idos + 1
     write(4,'(2f20.5)') e, dos(idos)
     write(44,'(2f20.5)') e*convert_thz_to_cm, dos(idos)/convert_thz_to_cm
     write(45,'(2f20.5)') e*convert_thz_to_meV, dos(idos)/convert_thz_to_mev
  enddo

  call vtime( tv, tc )
  write(*,'(/''CPU time '',f10.2,10x,''Elapsed time '',f10.2)') tv - tv0, tc - tc0

END



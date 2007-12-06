!========================================================================
SUBROUTINE force_constant( symm, dyn, sh, latt, a, alpha )
  !========================================================================
  USE nrtype
  USE nrutil
  USE data
  IMPLICIT NONE         

  INTEGER :: i, na, nb, nr, iprint

  REAL(DP) :: one(3,3), temp(3), ftot(3,3)
  REAL(DP), INTENT(IN) :: alpha, a

  TYPE (dynmat), INTENT(OUT) :: dyn
  TYPE (lattice), INTENT(INOUT) :: latt
  TYPE (shell), INTENT(IN)   :: sh
  TYPE (symmetry), INTENT(INOUT) :: symm
  !
  ! allocate space for the dynamical matrix
  !
  allocate(  dyn%dmat(3,3,latt%natoms,latt%natoms,sh%nrm), &
       dyn%weight(latt%natoms,latt%natoms,sh%nrm) )

  !------------------------------------------------
  ! find the symmetry operations of the supercell
  !------------------------------------------------
  call sgama ( symm, latt, iprint )

  !============================================================================
  ! construct the force constant matrix,     
  ! now I need x in cartesian coordinates 
  !============================================================================
  latt%x =  latt%scale * matmul( latt%at, latt%x )

  !============================================================================
  ! put is in carthesian coordinates
  !============================================================================
  do i=1,symm%nsym
     symm%is(:,:,i) = &
          &        matmul(latt%at,matmul(symm%is(:,:,i),transpose(latt%bg)))
     symm%ftau(:,i) = latt%scale * matmul(latt%at,symm%ftau(:,i))
  enddo

  print*,'Memory for dmat: ', &
       &     real(3*3*latt%natoms*latt%natoms*sh%nrm*8)/1000000,&
       &     'MBytes' 

  dyn%dmat = 0;  one=0
  do i=1,3
     one(i,i)=1
  enddo

  do nr = 1, sh%nrm
     do na = 1, latt%natoms
        do nb = 1, latt%natoms
           !
           ! jump the construction of D for the on-site terms
           !
           temp = sh%r(:,nr) + latt%x(:,nb) - latt%x(:,na) 
           if( (nr>1 .or. nb/=na) .and. (vabs(temp) < sh%rmax) )then
              dyn%dmat(:,:,na,nb,nr) = alpha * A * &
                   ( one - ( alpha + 2 ) * outerprod( temp, temp ) / &
                   dot_product( temp, temp ) ) / vabs( temp )**( alpha + 2 )
           endif
        enddo
     enddo
  enddo
  !
  ! on-site terms
  !
  do na = 1, latt%natoms
     do nb = 1, latt%natoms
        do nr=1,sh%nrm
           temp = sh%r(:,nr) + latt%x(:,nb) - latt%x(:,na) 
           if( (nr>1 .or. na/=nb) .and. (vabs(temp) < sh%rmax) ) then
              dyn%dmat(:,:,na,na,1) = dyn%dmat(:,:,na,na,1) - dyn%dmat(:,:,na,nb,nr)
           endif
        enddo
     enddo
  enddo

  !========================================================================
  ! check that the sum of the force constant matrices is zero
  !========================================================================
  ftot = 0
  do na = 1, latt%natoms
     do nb = 1, latt%natoms
        do nr = 1, sh%nrm
           temp = sh%r(:,nr) + latt%x(:,nb) - latt%x(:,na) 
           if(vabs(temp) < sh%rmax) then
              ftot = ftot + dyn%dmat(:,:,na,nb,nr)
           endif
        enddo
     enddo
  enddo

  print'(/''Total sum of the force constants:'')'
  print'(3f16.8)',ftot

!  if( any(abs(ftot)) > 1.d-6 ) print*, ftot

END SUBROUTINE force_constant


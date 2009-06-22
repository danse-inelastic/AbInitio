MODULE data
  USE nrutil
  IMPLICIT NONE

  TYPE lattice  
     CHARACTER*1 :: string
     CHARACTER*50 :: string2
     LOGICAL :: lsuper, lgamma
     INTEGER :: ntypes, natoms, natomss, ndim(3), qa, qb, qc, disp
     INTEGER, POINTER :: nions(:), ityp(:), ityps(:), super_atom(:), prim_atom(:)
     REAL(DP) :: at(3,3), bg(3,3), scale, volume
     REAL(DP) :: ats(3,3), bgs(3,3), volumes
     REAL(DP) :: dxstart(3)
     REAL(DP), POINTER :: x(:,:), xs(:,:), mass(:)
  END TYPE lattice

  TYPE shell
     INTEGER :: nrm, nsh
     REAL(DP) :: rmax, cutoff
     INTEGER, POINTER :: nl(:), ishel(:)
     REAL(DP), POINTER :: r(:,:), rr(:), rl(:)
  END TYPE shell

  TYPE dynmat
     LOGICAL :: lfree
     LOGICAL,      POINTER :: usethis(:)
     INTEGER :: npoints
     INTEGER,      POINTER :: weight(:,:,:)
     REAL(DP),     POINTER :: tmpmat2(:,:,:,:), dmat(:,:,:,:,:), invdmat(:,:,:,:,:), &
          qi(:,:), qf(:,:), eig(:)
     COMPLEX(DPC), POINTER :: cmat(:,:),  invcmat(:,:)
  END TYPE dynmat

  TYPE symmetry
     LOGICAL :: lsymm, invsym, lfullqgrid
     INTEGER :: nsym, nrot, nsymstart
     REAL(DP), POINTER :: is(:,:,:), isstart(:,:,:), iscryst(:,:,:), ftau(:,:)
     INTEGER, POINTER :: irt(:,:), irts(:,:)
  END TYPE symmetry


END MODULE data

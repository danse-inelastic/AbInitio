!==========================================================================
      FUNCTION equiva( a, b )
!==========================================================================
      USE nrtype
      IMPLICIT NONE
      LOGICAL :: equiva
      REAL(DP), PARAMETER :: eps = 1.d-6
      REAL(DP) :: a(3), b(3)

      equiva = abs(  a(1) - b(1)  - nint( a(1) - b(1) ) ) < eps .and. &
     &     abs( a(2) - b(2)  - nint( a(2) - b(2) ) ) < eps .and. &
     &     abs( a(3) - b(3)  - nint( a(3) - b(3) ) ) < eps 
!      equiva = equiva .and. .not. ( any( abs ( nint( a - b) ) > 1 ) )

      END FUNCTION

!==========================================================================
    FUNCTION equiva1( a, b )
!==========================================================================
      USE nrtype
      IMPLICIT NONE
      LOGICAL :: equiva1
      REAL(DP), PARAMETER :: eps = 2.d-6
      REAL(DP) :: a(3), b(3)

      equiva1 = abs(  a(1) - b(1) ) < eps .and. &
           abs( a(2) - b(2) ) < eps .and. &
           abs( a(3) - b(3) ) < eps

    END FUNCTION equiva1

      MODULE nr

      INTERFACE gammaln
         FUNCTION gammaln_d(xx)
         USE nrtype
         REAL(DP), INTENT(IN) :: xx
         REAL(DP) :: gammaln_d
         END FUNCTION gammaln_d
      END INTERFACE

      INTERFACE gammaq
         FUNCTION gammaq_d(a,x)
         USE nrtype
         REAL(DP), INTENT(IN) :: a,x
         REAL(DP) :: gammaq_d
         END FUNCTION gammaq_d
      END INTERFACE

      INTERFACE gcf
         FUNCTION gcf_d(a,x,gln)
         USE nrtype
         REAL(DP), INTENT(IN) :: a,x
         REAL(DP), OPTIONAL, INTENT(OUT) :: gln
         REAL(DP) :: gcf_d
         END FUNCTION gcf_d
      END INTERFACE

      INTERFACE gser
         FUNCTION gser_d(a,x,gln)
         USE nrtype
         REAL(DP), INTENT(IN) :: a,x
         REAL(DP), OPTIONAL, INTENT(OUT) :: gln
         REAL(DP) :: gser_d
         END FUNCTION gser_d
      END INTERFACE

      END MODULE nr

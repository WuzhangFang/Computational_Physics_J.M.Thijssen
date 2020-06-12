MODULE Excorr

USE globals

CONTAINS

  DOUBLE PRECISION FUNCTION Vxc (I, J, K)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K
  INTEGER :: ICnt

  DOUBLE PRECISION :: AParams(4), DAparams(4), &
                      BParams(4), DBParams(4), TA(4), TB(4)

  DATA AParams/0.4581652932831429D0, 2.217058676663745D0, &
                0.7405551735357053D0, 1.968227878617998D-2/
  DATA DAparams/0.119086804055547D0, 0.6157402568883345D0, &
                 0.1574201515892867D0, 3.532336663397157D-3/
  DATA BParams/1.D0, 4.504130959426697D0, 1.110667363742916D0, &
                2.359291751427506D-2/
  DATA DBParams/0.D0, 0.2673612973836267D0, 0.2052004607777787D0, &
                 4.200005045691381D-3/

  DOUBLE PRECISION :: LocDens, Rs, Denom, Numer, DDenom, DNumer

  LocDens = Density_R(I,J,K)
! Hydrogen case, density is up!!  
  Rs = (3/(4.D0*PI*LocDens))**0.33333333333333D0
  TA = AParams !+DAParams
  TB = BParams !+DBParams
  Numer = 0.D0
  DNumer = 0.D0
  Denom = 0.D0
  DDenom = 0.D0
  DO ICnt=1, 4
    Numer = Numer + TA(ICnt)*Rs**(ICnt-1)
    DNumer = DNumer + TA(ICnt)*(ICnt-1)*Rs**(ICnt-2)
    Denom = Denom + TB(ICnt)*Rs**(ICnt+3)
    DDenom = DDenom + TB(ICnt)*(ICnt+3)*Rs**(ICnt+2)
  END DO
  Vxc = (Denom*DNumer-Numer*DDenom)/Denom**2*Rs**4/3.D0
  
  END FUNCTION Vxc


  DOUBLE PRECISION FUNCTION epsilon_xc (I, J, K)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K
  INTEGER :: ICnt

  DOUBLE PRECISION :: AParams(4), DAparams(4), &
                      BParams(4), DBParams(4), TA(4), TB(4)

  DATA AParams/0.4581652932831429D0, 2.217058676663745D0, &
                0.7405551735357053D0, 1.968227878617998D-2/
  DATA DAparams/0.119086804055547D0, 0.6157402568883345D0, &
                 0.1574201515892867D0, 3.532336663397157D-3/
  DATA BParams/1.D0, 4.504130959426697D0, 1.110667363742916D0, &
                2.359291751427506D-2/
  DATA DBParams/0.D0, 0.2673612973836267D0, 0.2052004607777787D0, &
                 4.200005045691381D-3/

  DOUBLE PRECISION :: LocDens, Rs, Denom, Numer

  LocDens = Density_R(I,J,K)
! Hydrogen case, density is up!!  
  Rs = (3/(4.D0*PI*LocDens))**0.33333333333333D0
  TA = AParams!+DAParams
  TB = BParams!+DBParams
  Numer = 0.D0
  Denom = 0.D0
  DO ICnt=1, 4
    Numer = Numer + TA(ICnt)*Rs**(ICnt-1)
    Denom = Denom + TB(ICnt)*Rs**ICnt
  END DO
  epsilon_xc = -Numer/Denom
  
  END FUNCTION epsilon_xc



END MODULE ExCorr

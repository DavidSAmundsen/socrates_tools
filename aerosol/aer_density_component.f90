SUBROUTINE aer_density_component(i_component, density)

  USE realtype_rd, ONLY: RealK
  USE rad_pcf

  IMPLICIT NONE

! Input variables
  INTEGER, INTENT(IN) :: i_component
!   Aerosol component index

! Output variables
  REAL*8, INTENT(OUT) :: density
!   Density of aerosol

  INCLUDE 'aerosol_component.finc'

  density = density_component(i_component)

END SUBROUTINE aer_density_component
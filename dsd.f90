


MODULE dsd

USE defs_eddington


IMPLICIT NONE


!======================================================================

CONTAINS

!--------------------------------------------------------

!--------------------------------------------------
!
!    Gamma lut consists of an array of
!    values for mu, lambda, and the
!    total water content for each value pair.
!    More information can be found in
!    lwc_gamma_LUT.
!
!    Spencer Jones, CSU Atmos., 10/2023
!
!--------------------------------------------------

SUBROUTINE read_gamma_lut

CHARACTER(LEN=100) :: gamma_lut_file = 'LUT/gamma/lwc_gamma_LUT_N08000.tbl'

INTEGER :: nN0, nmu, nlamb

OPEN(UNIT=99, FILE=gamma_lut_file, STATUS='OLD', ACCESS='STREAM', &
        FORM='UNFORMATTED')

READ(99) N0_lut
READ(99) nmu
ALLOCATE(mu_lut(nmu))
READ(99) mu_lut(:)
READ(99) nlamb
ALLOCATE(lamb_lut(nlamb))
READ(99) lamb_lut(:)
ALLOCATE(gamma_lwc_lut(nmu,nlamb))
READ(99) gamma_lwc_lut(:,:)


CLOSE(99)


RETURN

END SUBROUTINE read_gamma_lut


!--------------------------------------------------------

SUBROUTINE gamma_dsd(lwc, N_0, mu, lamb)

!Inputs
REAL :: lwc ![g m^-3]       
REAL :: N_0 ![number mm^-1 m^-3]
REAL :: mu  ![]

!Output
REAL :: lamb ![mm^-1]

REAL    :: mu_min
INTEGER :: mu_indx(1)
INTEGER :: lwc_indx(1)
REAL    :: dmu
INTEGER :: N_0_indx(1)
REAL    :: log_n0


mu_indx  = MINLOC(ABS(mu - mu_lut(:)))

lwc_indx = MINLOC(ABS(lwc - gamma_lwc_lut(mu_indx(1),:)))


lamb = lamb_lut(lwc_indx(1))


RETURN

END SUBROUTINE gamma_dsd


!----------------------------------------------------------



!=========================================================================

END MODULE dsd

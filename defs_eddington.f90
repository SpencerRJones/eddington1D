MODULE defs_eddington

IMPLICIT NONE

SAVE

LOGICAL :: debug = .FALSE.
INTEGER :: gas_abs = 1             ! 0 = MonoRTM, 1 = Eddington Absorption

INTEGER, PARAMETER :: nchans = 13

INTEGER, PARAMETER :: nlevs = 31
INTEGER, PARAMETER :: nlyrs = 30

REAL :: frequency(nchans)          = (/ 10.65, 10.65, 18.7, 18.7, 23.8, 36.64, 36.64, &
                                       89.0, 89.0, 166.0, 166.0, 180.31, 190.31 /)
INTEGER :: polarization(nchans)      = (/ 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0 /)
REAL :: incidence_angle(nchans)   = (/ 52.8, 52.8, 52.8, 52.8, 52.8, 52.8, 52.8, &
                                        52.8, 52.8, 49.2, 49.2, 49.2, 49.2 /)

REAL, PARAMETER :: n0_liq = 8000. !Intercept parameter [mm^-1 m^-3]
REAL, PARAMETER :: n0_sn = 5100.  !Intercept parameter [mm^-1 m^-3]
REAL, PARAMETER :: mu_liq = 0.    !Shape parameter []
REAL, PARAMETER :: rho_sn = 0.4   !Snow particle density [g cm^-3]
REAL, PARAMETER :: rho_ciw = 0.1  !Cloud ice density [g cm^-3]

REAL, PARAMETER :: r_clw = 0.1  !Radius of monodisperse cloud water droplets[mm]
REAL, PARAMETER :: r_ciw = 0.03 !Radius of monodisperse cloud ice particles [mm]




!Vars for MonoRTM:

      ! variables for MonoRTM lookup table

      integer              :: nfreq_lut
      integer              :: nchan_lut
      integer              :: npres_lut
      integer              :: ntemp_lut
      integer              :: nrmix_lut

      integer, allocatable :: ifreq_lut(:)
      integer, allocatable :: ipol_lut(:)
      real, allocatable    :: freq_lut(:)
      real, allocatable    :: pres_lut(:)
      real, allocatable    :: temp_lut(:)
      real, allocatable    :: rmix_lut(:)
      real, allocatable    :: kabs_lut(:,:,:,:)


!Vars for gamma LUT:
        !REAL, ALLOCATABLE :: N0_lut(:)
        REAL              :: N0_lut
        REAL, ALLOCATABLE :: mu_lut(:)
        REAL, ALLOCATABLE :: lamb_lut(:)
        !REAL, ALLOCATABLE :: gamma_lwc_lut(:,:,:)
        REAL, ALLOCATABLE :: gamma_lwc_lut(:,:)

END MODULE defs_eddington

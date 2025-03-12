  MODULE FASTEM

  USE Type_Kinds, ONLY: fp
  USE CRTM_FastemX
  USE CRTM_MWwaterCoeff
  USE Message_Handler

  IMPLICIT NONE

  INTEGER :: FASTEM_MODEL=0

  CONTAINS

!--------------------------------------------------------------------------------
! NAME:
!       Compute_Fastem
!
! PURPOSE:
!       Subroutine to compute the Fastem4, Fastem5, or Fastem6 microwave sea surface
!       emissivity and reflectivity.
!
! CALLING SEQUENCE:
!       CALL Compute_Fastem(
!              Frequency    , &  ! Input
!              Zenith_Angle , &  ! Input
!              Temperature  , &  ! Input
!              Salinity     , &  ! Input
!              Wind_Speed   , &  ! Input
!              Emissivity   , &  ! Output
!              Reflectivity , &  ! Output
!              Azimuth_Angle)
!
!
! INPUTS:
!       Frequency:      Microwave frequency.
!                       UNITS:      GHz
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Zenith_Angle:   Sensor zenith angle at the sea surface
!                       UNITS:      Degrees
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Temperature:    Sea surface temperature
!                       UNITS:      Kelvin, K
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Salinity:       Water salinity
!                       UNITS:      ppt (parts per thousand)
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Wind_Speed:     Sea surface wind speed
!                       UNITS:      m/s
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Azimuth_Angle:  Relative azimuth angle (wind direction - sensor azimuth)
!                       UNITS:      Degrees
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
! OUTPUTS:
!       Emissivity:     The surface emissivity
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Rank-1, 4-elements (n_Stokes)
!                       ATTRIBUTES: INTENT(OUT)
!
!       Reflectivity:   The surface reflectivity.
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Rank-1, 4-elements (n_Stokes)
!                       ATTRIBUTES: INTENT(OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE Compute_Fastem(model,    &
                            offsets,  &
                            diffuse,  &
                            freq,     &
                            eia,      &
                            sst,      &
                            salinity, &
                            wsp,      &
                            azimuth,  &
                            trans,    &
                            emis,     &
                            refl,     &
                            verbose)

  integer            :: model
  integer            :: offsets
  integer            :: diffuse
  real               :: freq
  real               :: eia
  real               :: sst
  real               :: salinity
  real               :: wsp
  real               :: azimuth
  real               :: trans
  real(4)            :: emis(4)
  real(4)            :: refl(4)
  integer, optional  :: verbose

  character(len=200) :: File_Path
  character(len=200) :: MWWaterCoeff_File
  character(len=200) :: msg,pid_msg
  integer            :: err_stat
  logical            :: Quiet
  integer            :: Process_ID
  integer            :: Output_Process_ID
  real(fp)           :: frequency
  real(fp)           :: view_angle
  real(fp)           :: tsfc
  real(fp)           :: salin
  real(fp)           :: wndspd
  real(fp)           :: azimuth_angle
  real(fp)           :: emissivity(4)
  real(fp)           :: reflectivity(4)
  TYPE(iVar_type)    :: iVar
  real(fp)           :: transmittance
  
  ! Added for emissivity corrections
  real :: emis_off(2)
  real :: fq
  integer :: pz

  !if (PRESENT(verbose)) then
  !  Quiet=0
  !else
    Quiet=1
  !endif
  File_Path = 'fastem/'
  if (model .eq. 4) then
    MWwaterCoeff_File = 'FASTEM4.MWwater.EmisCoeff.bin'
  else
    if (model .eq. 5) then
      MWwaterCoeff_File = 'FASTEM5.MWwater.EmisCoeff.bin'
    else
      MWwaterCoeff_File = 'FASTEM6.MWwater.EmisCoeff.bin'
    endif
  endif
  MWwaterCoeff_File  = TRIM(ADJUSTL(File_Path)) // TRIM(MWwaterCoeff_File)



  ! Load coefficient file for MW water

  if (FASTEM_MODEL .eq. 0 .or. model .ne. FASTEM_MODEL) then
    err_stat = CRTM_MWwaterCoeff_Load( MWwaterCoeff_File, &
                 Quiet             = Quiet            , &
                 Process_ID        = Process_ID       , &
                 Output_Process_ID = Output_Process_ID  )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error loading MWwaterCoeff data from '//TRIM(MWwaterCoeff_File)
      CALL Display_Message("FASTEM",TRIM(msg)//TRIM(pid_msg),err_stat )
      STOP
    END IF
    FASTEM_MODEL = model
  endif


  ! FastemX model

  if (verbose) then
    write(6,'("model    = ",i2)') model
    write(6,'("offsets  = ",i2)') offsets
    write(6,'("diffuse  = ",i2)') diffuse
    write(6,'("freq     = ",f7.2)') freq
    write(6,'("eia      = ",f7.2)') eia
    write(6,'("sst      = ",f7.2)') sst
    write(6,'("salinity = ",f7.4)') salinity
    write(6,'("wsp      = ",f7.2)') wsp
    write(6,'("azimuth  = ",f7.2)') azimuth
    write(6,'("trans    = ",f7.4)') trans
  endif
  
  frequency     = freq
  view_angle    = eia
  tsfc          = sst
  salin         = salinity
  wndspd        = wsp
  azimuth_angle = azimuth
  if (diffuse .eq. 1) then
    transmittance = trans
  else
    transmittance = 1.0
  endif

  CALL Compute_FastemX(                 &
         MWwaterC,                      & ! Input model coefficients
         frequency,                     & ! Input frequency (GHz)
         view_angle,                    & ! Input view angle (degrees)
         tsfc,                          & ! Input surface temperature (K)
         salin,                         & ! Input salinity (ppm)
         wndspd,                        & ! Input wind speed (m/sec)
         iVar,                          & ! Internal variable output
         emissivity,                    & ! Output emissivity
         reflectivity,                  & ! Output reflectivity
         Azimuth_Angle = azimuth_angle, & ! Optional input (wind direction - sensor azimuth)
         Transmittance = transmittance)   ! Optional input (atmospheric transmittance)

  emis(:) = emissivity(:)
  refl(:) = reflectivity(:)

  !write(6,'("emis = ",4(f9.4))') emis
  !write(6,'("refl = ",4(f9.4))') refl

  END SUBROUTINE Compute_Fastem
  END MODULE Fastem

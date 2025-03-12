







SUBROUTINE eddington(freqy, polz, eia, height, plevel, tlevel, mix_ratio, cloud_water, rain_water, &
                cloud_ice, snow_water, tsfc, saln, wspd, n0_plw, &
                mu_plw, n0_snow, rho_snow, tbout)

USE MonoRTM
USE defs_eddington
USE fastem
USE mie

IMPLICIT NONE

REAL, PARAMETER :: pi = 3.141592654

!---INPUTS
REAL :: freqy(nchans)        !Vector of channel frequencies [GHz]
INTEGER :: polz(nchans)          !Vector of channel polarizations (0 = vertical, 1 = horizontal)
REAL :: eia(nchans)          !Vector of channel earth incidence angles [degrees]
REAL :: height(nlevs)        !Layer boundary heights from sfc to TOA [km]
REAL :: plevel(nlevs)        !Layer boundary pressures [hPa]
REAL :: tlevel(nlevs)        !Layer boundary temperatures [K]
REAL :: mix_ratio(nlyrs)     !Layer mixing ratio [g/kg]
REAL :: cloud_water(nlyrs)   !Layer cloud water content [g/m^3]
REAL :: rain_water(nlyrs)    !Layer rain water content [g/m^3]
REAL :: cloud_ice(nlyrs)     !Layer cloud ice content [g/m^3]
REAL :: snow_water(nlyrs)    !Layer snow water content [g/m^3]
REAL :: tsfc                 !Surface temperature (or SST) [K]
REAL :: saln                 !Sea surface salinity [-9999.9 for land]
REAL :: wspd                 !Surface Wind speed [m/s]
REAL :: n0_plw               !Gamma distribution intercept parameter [mm^-1 m^-3]
REAL :: mu_plw               !Gamma distribution shape parameter []
REAL :: n0_snow              !Gamma distribution intercept parameter [mm^-1 m^-3]
REAL :: rho_snow             !Snow particle density [g cm^-3]


!---OUTPUTS
REAL :: tbout(nchans)     !Output simulation vector

!---General Vars:
INTEGER :: ichan, ilyr
INTEGER :: pol
REAL    :: freq, view_angle, optdepth
REAL    :: kext(nlyrs), salb(nlyrs), asym(nlyrs)
REAL    :: dhgt(nlyrs)
REAL    :: tavg, pavg
REAL    :: kabs_clr
REAL    :: kext_clw, salb_clw, asym_clw, pbck_clw
REAL    :: kext_plw, salb_plw, asym_plw, pbck_plw
REAL    :: kext_ciw, salb_ciw, asym_ciw, pbck_ciw
REAL    :: kext_piw, salb_piw, asym_piw, pbck_piw
REAL    :: trans
REAL    :: femis(4), frefl(4)
INTEGER :: verbose=0
REAL    :: emis(nchans), refl(nchans), ebar
REAL    :: tb, tbdown

IF (debug) THEN
        WRITE(*,'("Height        = ",31(f7.2))') height
        WRITE(*,'("Pressure      = ",31(f7.2))') plevel
        WRITE(*,'("Temp          = ",31(f7.2))') tlevel
        WRITE(*,'("Mix Ratio     = ",30(f7.2))') mix_ratio
        WRITE(*,'("Cloud Water   = ",30(f7.2))') cloud_water
        WRITE(*,'("Rain Water    = ",30(f7.2))') rain_water
        WRITE(*,'("Cloud Ice     = ",30(f7.2))') cloud_ice
        WRITE(*,'("Snow Water    = ",30(f7.2))') snow_water
        WRITE(*,'("Tsfc          = ",30(f7.2))') tsfc
        WRITE(*,'("Salinity      = ",30(f7.2))') saln
        WRITE(*,'("Wind Speed    = ",30(f7.2))') wspd
        WRITE(*,'("N0_rain       = ",30(f7.2))') n0_plw
        WRITE(*,'("mu_rain       = ",30(f7.2))') mu_plw
        WRITE(*,'("N0_snow       = ",30(f7.2))') n0_snow
        WRITE(*,'("Rho_snow      = ",30(f7.2))') rho_snow
        WRITE(*,'("Frequency:    = ",13(f7.2))') freqy
        WRITE(*,'("Polarization: = ",13(f7.2))') polz
        WRITE(*,'("EIA:          = ",13(f7.2))') eia
ENDIF


!---Loop through channels and compute Tb:

kext(:) = 0.
salb(:) = 0.
asym(:) = 0.

DO ilyr = 1, nlyrs
        dhgt(ilyr) = height(ilyr+1) - height(ilyr)
ENDDO

DO ichan = 1, nchans

        freq       = freqy(ichan)
        view_angle = eia(ichan)
        pol        = polz(ichan) + 1

        optdepth = 0.

        DO ilyr = 1, nlyrs
                
                tavg = 0.5 * (tlevel(ilyr) + tlevel(ilyr+1))
                pavg = (plevel(ilyr) - plevel(ilyr+1)) / LOG(plevel(ilyr)/plevel(ilyr+1))

                !Clear air absorption:
                IF (gas_abs .EQ. 0) THEN  !Use MonoRTM
                        CALL monortm_lut(ichan, pavg, tavg, mix_ratio(ilyr), kabs_clr)
                ELSEIF (gas_abs .EQ. 1) THEN !Use builtin Eddington absorption
                        CALL absorb_clr(freq, tavg, pavg, mix_ratio(ilyr), kabs_clr)
                ELSE
                        WRITE(*,*) 'Error: No gas absorption specified.'
                        STOP
                ENDIF

                !Cloud water scattering (monodisperse):
                CALL mie_clw(freq, tavg, cloud_water(ilyr), r_clw, &
                             kext_clw, salb_clw, asym_clw, pbck_clw)
     
                !write(*,*) kext_clw, salb_clw, asym_clw, pbck_clw
                

                !Rain water scattering (Gamma distribution):
                IF (N0_plw .EQ. 8000. .AND. mu_plw .EQ. 0.) THEN
                        CALL mie_rain_mp(freq, tavg, rain_water(ilyr), &
                                         kext_plw, salb_plw, asym_plw, pbck_plw)
                        !write(*,*) kext_plw, salb_plw, asym_plw, pbck_plw
                ELSE
                        !CALL mie_rain()
                ENDIF

                !Cloud ice scattering (monodisperse):
                CALL mie_ciw(freq, tavg, cloud_ice(ilyr), r_ciw, rho_ciw, &
                             kext_ciw, salb_ciw, asym_ciw, pbck_ciw)
                !write(*,*) kext_ciw, salb_ciw, asym_ciw, pbck_ciw

                !Snow/graupel scattering (Inverse exponential distribution):
                CALL mie_snow(freq, tavg, snow_water(ilyr), n0_snow, rho_snow, &
                              kext_piw, salb_piw, asym_piw, pbck_piw)
                !write(*,*) kext_piw, salb_piw, asym_piw, pbck_piw

      
      
      
                !Add up k_ext from all sources      
                kext(ilyr) = kabs_clr + kext_clw + kext_plw + kext_ciw + kext_piw

                salb(ilyr) = ((salb_clw*kext_clw) + &
                              (salb_plw*kext_plw) + &
                              (salb_ciw*kext_ciw) + &
                              (salb_piw*kext_piw) / kext(ilyr))

                IF (salb(ilyr) .GT. 0) THEN
                        asym(ilyr) = ((asym_clw*salb_clw*kext_clw) +  &
                                      (asym_plw*salb_plw*kext_plw) +  &
                                      (asym_ciw*salb_ciw*kext_ciw) +  &
                                      (asym_piw*salb_piw*kext_piw)) / &
                                      (salb(ilyr) * kext(ilyr))
                ELSE
                        asym(ilyr) = 0.
                ENDIF

                optdepth = optdepth + (kext(ilyr) * dhgt(ilyr))

        ENDDO !Layers

        trans = EXP(-optdepth)

        CALL compute_fastem(6, 0, 1, freq, view_angle, tsfc, saln, &
                            wspd, 0., trans, femis, frefl, verbose)

        IF (pol .EQ. 1) THEN     !Vertical polarization
                emis(ichan) = femis(1)
                refl(ichan) = frefl(1)
        ELSEIF (pol .EQ. 2) THEN !Horizontal polarization
                emis(ichan) = femis(2)
                refl(ichan) = frefl(2)
        ELSE !Should not have any other value for AMSR2
                WRITE(*,*) 'Polarization value out of bounds.'
                STOP
        ENDIF

        ebar = emis(ichan)


        IF (debug) THEN
                WRITE(*,*) 'Into radtran:'
                WRITE(*,*) 'nz:          ', nlyrs
                WRITE(*,*) 'view_angle:  ', view_angle
                WRITE(*,*) 'tsfc:        ', tsfc
                WRITE(*,*) 'hgt:         ', height
                WRITE(*,*) 'tmp:         ', tlevel
                WRITE(*,*) 'kext:        ', kext
                WRITE(*,*) 'salb:        ', salb
                WRITE(*,*) 'emis(ichan): ', emis(ichan)
                WRITE(*,*) 'ebar:        ', ebar
                WRITE(*,*) 'refl(ichan): ', refl(ichan)
        ENDIF

        CALL radtran(nlyrs,view_angle,tsfc,height,tlevel, &
                     kext,salb,asym,emis(ichan),ebar,refl(ichan),tb,tbdown)

        tbout(ichan) = tb


ENDDO !Channel loop

IF (debug) THEN
        WRITE(*,*) 'Tbs Out: ', tbout
ENDIF

RETURN

END SUBROUTINE eddington 

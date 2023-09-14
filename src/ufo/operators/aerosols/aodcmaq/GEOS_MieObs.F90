!
! This module --load the optics tables speficied in a resource file
!             -- destroy the tables
!             -- get AOD at observations location from aer mixing ratio profiles
! ! Interface:
!

    MODULE GEOS_MieObs_mod

! ! Uses:
!
    use Chem_MieMod
    use geos_aero_kinds_mod
    Implicit None
!
    PRIVATE
    PUBLIC get_GEOS_AOD       ! return AOD from aer Mass Mixing Ratio,
                              ! RH and DELP
    PUBLIC get_GEOS_AOD_tl
    PUBLIC get_GEOS_AOD_ad

CONTAINS

!-----------------

 subroutine get_geos_aod ( km, nobs, nch, nq, rcfile, channels, vname,  &
                           qm, rh,  aod_tot, ext, rc )

! Returns AOD at the asking channels at obs location.

  implicit NONE

  integer,               intent(in)  :: km               ! number vertical layers
  integer,               intent(in)  :: nobs             ! number of profiles

  character(len=*),      intent(in)  :: rcfile           ! resource file, e.g., Aod_***.rc
  integer,               intent(in)  :: nch              ! number of wavelengths
  real(kind=kind_real),                  intent(in)  :: channels(nch)    ! Wavelengths

  integer,               intent(in)  :: nq               ! number of tracers
  character(len=*),      intent(in)  :: vname(nq)        ! variable name

  real(kind=kind_real),  intent(in)  :: qm(nq,km,nobs)   ! (mixing ratio) * delp/g -> [kg/m2]
  real(kind=kind_real),  intent(in)  :: rh(km,nobs)      ! relative humidity

  real(kind=kind_real),  optional,  intent(out) :: aod_tot(nch, nobs)   ! total aerosol optical depth (summed over km)
  real(kind=kind_real),  optional,  intent(out) :: ext(km,nch,nq,nobs) ! mass extinction efficiency [m2 (kg dry mass)-1]

  integer,  optional,          intent(out) :: rc

!                               ---

  real                :: idxChannel(nch) ! this should have been integer
  integer             :: idxTable
  integer :: iq, n, m, i, k
  real(kind=kind_real) :: aod(km, nch, nobs)
  real(kind=kind_real) :: bext_, tau_

  type(Chem_Mie) :: mieTables        ! Mie Tables
  character(len=16) :: vname_(nq)        ! variable name

  rc = 0
  ! Create the Mie Tables
 ! ---------------------
  mieTables = Chem_MieCreate(rcfile,rc)
  if ( rc /= 0 ) then
     print *, 'Cannot create Mie tables from '//trim(rcfile)
     return
  end if

! Determine channel indices
! -------------------------
  do n = 1, nch
     idxChannel(n) = -1 ! this is really the channel index
     do m = 1, mieTables%nch
        if ( abs(channels(n) - (1.e9)*mieTables%channels(m)) < 1. ) then
           idxChannel(n) = m
           exit
         end if
     end do
  end do

  if ( any(idxChannel<0) ) then
        print *, 'Mie resource files does not set the required channel'
        print *, 'Channels requested:  ', channels
        print *, 'Channels on RC file: ', 1.e+9 * mieTables%channels
        rc = 99
        return
  end if

! Initialize output arrays to zero
! --------------------------------
  AOD = 0.0_kind_real
  if ( present(aod_tot) ) AOD_tot = 0.0_kind_real
  if ( present(ext) ) ext = 0.0_kind_real
  
! Loop over aerosol species
! -------------------------
  do iq =1, nq

     if(vname(iq) == 'mass_fraction_of_dust001_in_air') then
        vname_(iq) = 'du001'

     else if(vname(iq) == 'mass_fraction_of_dust002_in_air') then
        vname_(iq) = 'du002'
     
     else if(vname(iq) == 'mass_fraction_of_dust003_in_air') then
        vname_(iq) = 'du003'

     else if(vname(iq) == 'mass_fraction_of_dust004_in_air') then
        vname_(iq) = 'du004'

     else if(vname(iq) == 'mass_fraction_of_dust005_in_air') then
        vname_(iq) = 'du005'

     else if(vname(iq) == 'mass_fraction_of_sea_salt001_in_air') then
        vname_(iq) = 'ss001'
     
     else if(vname(iq) == 'mass_fraction_of_sea_salt002_in_air') then
        vname_(iq) = 'ss002'

     else if(vname(iq) == 'mass_fraction_of_sea_salt003_in_air') then
        vname_(iq) = 'ss003'

     else if(vname(iq) == 'mass_fraction_of_sea_salt004_in_air') then
        vname_(iq) = 'ss004'

     else if(vname(iq) == 'mass_fraction_of_sea_salt005_in_air') then
        vname_(iq) = 'ss005'

     else if(vname(iq) == 'mass_fraction_of_hydrophobic_black_carbon_in_air') then
        vname_(iq) = 'bcphobic'
     
     else if(vname(iq) == 'mass_fraction_of_hydrophilic_black_carbon_in_air') then
        vname_(iq) = 'bcphilic'

     else if(vname(iq) == 'mass_fraction_of_hydrophobic_organic_carbon_in_air') then
        vname_(iq) = 'ocphobic'

     else if(vname(iq) == 'mass_fraction_of_hydrophilic_organic_carbon_in_air') then
        vname_(iq) = 'ocphilic'

     else if(vname(iq) == 'mass_fraction_of_sulfate_in_air') then
        vname_(iq) = 'so4'

     else if(vname(iq) == 'mass_fraction_of_nitrate001_in_air') then
        vname_(iq) = 'no3an1'

     else if(vname(iq) == 'mass_fraction_of_nitrate002_in_air') then
        vname_(iq) = 'no3an2'

     else if(vname(iq) == 'mass_fraction_of_nitrate003_in_air') then
        vname_(iq) = 'no3an3'

     endif
  enddo 

  do iq = 1,nq

     idxTable = Chem_MieQueryIdx(mieTables,vname_(iq),rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        print *, 'cannot get Mie index for '//vname_(iq)
        return
     end if

!    Loop over nobs, km, nch
!    --------------------------
     do i = 1, nobs
        do n = 1, nch
           do k =1, km


                call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                               qm(iq,k,i), rh(k,i), bext=bext_, tau=tau_, rc= rc)
                            AOD(k,n,i) = AOD(k,n,i) + bext_*qm(iq,k,i)  ! sum over the tracers
                            if ( present(ext))  EXT(k,n,iq,i) = bext_
                

           end do  ! end km
          end do ! end nch
      end do  ! end nobs
  end do ! end tracers

 ! sum over the layers to get total AOD for all obs and wavelengths
  if (present(aod_tot)) then
    do k = 1, km
       AOD_tot = AOD_tot + AOD(k,:,:)
    enddo
  endif
 
 ! Clean up
  call Chem_MieDestroy(mieTables, rc)
  if ( rc /= 0 ) then
     print *, 'Cannot destroy MieTables'
     return
  end if

  end subroutine get_geos_aod

!-------------------------------------

  subroutine get_geos_aod_tl(km, nobs, nch, nq, bext, qm_tl, aod_tot_tl)

  implicit none
  integer, intent(in)    :: km                      ! number of layers
  integer, intent(in)    :: nobs                    ! number of profiles
  integer, intent(in)    :: nch                     ! number of wavelengths
  integer, intent(in)    :: nq                      ! number of tracers
  real(kind=kind_real),    intent(in)    :: bext(km, nch, nq, nobs) ! mass extinction efficiency [m2 (kg dry mass)-1]

  real(kind=kind_real),    intent(in)    :: qm_tl( nq, km, nobs)    ! (mixing ratio) * delp/g -> [kg/m2]
  real(kind=kind_real),    intent(inout) :: aod_tot_tl(nch, nobs)   ! aerosol optical depth tangent linear

  integer :: ob, ch, tr, lv

  aod_tot_tl = 0.0_kind_real
  do ob = 1, nobs
    do ch = 1, nch
      do lv = 1,km
        do tr = 1,nq
          aod_tot_tl(ch,ob) = aod_tot_tl(ch,ob) + bext(lv,ch,tr,ob) * qm_tl(tr,lv,ob)
        enddo
      enddo
    enddo
  enddo      
  end subroutine get_geos_aod_tl

! -----------------------------------
  subroutine get_geos_aod_ad(km, nobs, nch, nq, bext,  aod_tot_ad, qm_ad)

  implicit none
  integer, intent(in)    :: km                       ! number of layers
  integer, intent(in)    :: nobs                     ! number of profiles
  integer, intent(in)    :: nch                      ! number of wavelengths
  integer, intent(in)    :: nq                       ! number of tracers
  real(kind=kind_real),    intent(in)    :: bext(km, nch, nq, nobs)  ! mass extinction efficiency
  real(kind=kind_real),    intent(in)    :: aod_tot_ad(nch, nobs)    ! aerosol optical depth adjoint
  real(kind=kind_real),    intent(out)   :: qm_ad( nq, km, nobs)      ! aerosol concentration

  integer :: ob, ch, tr, lv

  qm_ad = 0.0_kind_real

  do ob=nobs,1,-1
    do ch=nch,1,-1
      do lv=km,1,-1
        do tr=nq,1,-1
          qm_ad(tr, lv, ob) = qm_ad(tr, lv, ob) + bext(lv, ch, tr, ob)*aod_tot_ad(ch, ob)
        end do
      end do
    end do
  end do

  end subroutine get_geos_aod_ad

end module GEOS_MieObs_mod

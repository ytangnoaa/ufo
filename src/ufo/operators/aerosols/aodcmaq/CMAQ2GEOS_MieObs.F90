!
! This module --load the optics tables specified in a resource file
!             -- destroy the tables
!             -- get AOD at observations location from aerosol mixing ratio profiles
! ! Interface:
!

    MODULE CMAQ2GEOS_MieObs_mod

! ! Uses:
!
    use Chem_MieMod
    use aero_kinds_mod
    Implicit None
!
    PRIVATE
    PUBLIC get_CMAQ_AOD       ! return AOD from aerosol Mass Mixing Ratio, RH and DELP
    PUBLIC get_CMAQ_AOD_tl
    PUBLIC get_CMAQ_AOD_ad

CONTAINS

!-----------------

 subroutine get_CMAQ_AOD ( km, nobs, nch, nq, rcfile, channels, cmaq2geos,  &
                           qm, rh,  aod_tot, ext, rc )

! Returns AOD at the asking channels at obs location.

  implicit NONE

  integer,               intent(in)  :: km               ! number vertical layers
  integer,               intent(in)  :: nobs             ! number of profiles

  character(len=*),      intent(in)  :: rcfile           ! resource file, e.g., Aod_***.rc
  integer,               intent(in)  :: nch              ! number of wavelengths
  real(kind=kind_real),  intent(in)  :: channels(nch)    ! Wavelengths

  integer,               intent(in)  :: nq               ! number of tracers
  integer,               intent(in)  :: cmaq2geos(nq)    ! Map cmaq to geos

  real(kind=kind_real),  intent(in)  :: qm(nq,km,nobs)   ! (mixing ratio) * delp/g -> [kg/m2]
  real(kind=kind_real),  intent(in)  :: rh(km,nobs)      ! relative humidity

  real(kind=kind_real),  optional,  intent(out) :: aod_tot(nch, nobs)   ! total aerosol optical depth (summed over km)
  real(kind=kind_real),  optional,  intent(out) :: ext(km,nch,nq,nobs)  ! mass extinction efficiency [m2 (kg dry mass)-1]

  integer,  optional,          intent(out) :: rc

!                               ---

  real                :: idxChannel(nch) ! this should have been integer
  integer             :: idxTable
  integer :: iq, n, m, i, k
  real(kind=kind_real) :: aod(km, nch, nobs)
  real(kind=kind_real) :: bext_, tau_

  type(Chem_Mie) :: mieTables        ! Mie Tables
  character(len=16) :: cmaq2geos_(nq)        ! variable name

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

     if(cmaq2geos(iq) == 1) then
        cmaq2geos_(iq) = 'du001'

     else if(cmaq2geos(iq) == 2) then
        cmaq2geos_(iq) = 'ss001'
     
     else if(cmaq2geos(iq) == 3) then
        cmaq2geos_(iq) = 'bcphilic'

     else if(cmaq2geos(iq) == 4) then
        cmaq2geos_(iq) = 'ocphilic'

     else if(cmaq2geos(iq) == 5) then
        cmaq2geos_(iq) = 'so4'

     else if(cmaq2geos(iq) == 6) then
        cmaq2geos_(iq) = 'no3an1'

     else if(cmaq2geos(iq) == 7) then
        cmaq2geos_(iq) = 'bcphobic'

     else if(cmaq2geos(iq) == 8) then
        cmaq2geos_(iq) = 'ocphobic'

     else if(cmaq2geos(iq) == 9) then
        cmaq2geos_(iq) = 'no3an2'

     else if(cmaq2geos(iq) == 10) then
        cmaq2geos_(iq) = 'ss002'

     else if(cmaq2geos(iq) == 11) then
        cmaq2geos_(iq) = 'ss003'

     else if(cmaq2geos(iq) == 12) then
        cmaq2geos_(iq) = 'ss004'

     else if(cmaq2geos(iq) == 13) then
        cmaq2geos_(iq) = 'du003'
     endif
  enddo 

  do iq = 1,nq

     idxTable = Chem_MieQueryIdx(mieTables,cmaq2geos_(iq),rc)
     if(idxTable == -1) cycle
     if ( rc/=0 ) then
        print *, 'cannot get Mie index for '//cmaq2geos_(iq)
        return
     end if

!    Loop over nobs, km, nch
!    --------------------------
     do i = 1, nobs
        do n = 1, nch
           do k =1, km


                call Chem_MieQuery(mieTables, idxTable, idxChannel(n), &
                               qm(iq,k,i), rh(k,i), tau=tau_, bext=bext_, rc=rc)
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

  end subroutine get_CMAQ_AOD

!-------------------------------------

  subroutine get_CMAQ_AOD_tl(km, nobs, nch, nq, bext, qm_tl, aod_tot_tl)

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
  end subroutine get_CMAQ_AOD_tl

! -----------------------------------
  subroutine get_CMAQ_AOD_ad(km, nobs, nch, nq, bext,  aod_tot_ad, qm_ad)

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

  end subroutine get_CMAQ_AOD_ad

end module CMAQ2GEOS_MieObs_mod

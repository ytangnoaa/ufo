! (C) Copyright 2021.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for aodcmaq observation operator

module ufo_aodcmaq_mod

 use iso_c_binding
 use kinds
 use CMAQ2GEOS_MieObs_mod
 use oops_variables_mod
 use ufo_vars_mod
 use ufo_basis_mod, only: ufo_basis

 implicit none
 private

!> Fortran derived type for the observation type
 type, public :: ufo_aodcmaq
 private
   type(oops_variables), public :: geovars
   type(oops_variables), public :: obsvars
   integer, public              :: ntracers
   integer, public, allocatable :: cmaq2geosrc(:)
   real(kind_real), public, allocatable :: wavelength(:)
   character(len=maxvarlen),public :: rcfile
 contains
   procedure :: setup  => ufo_aodcmaq_setup
   procedure :: simobs => ufo_aodcmaq_simobs
   final :: destructor
 end type ufo_aodcmaq

!> Default variables required from model
 character(len=maxvarlen), dimension(2), parameter :: varindefault = (/var_delp, var_rh/)

contains

! ------------------------------------------------------------------------------
subroutine ufo_aodcmaq_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_aodcmaq), intent(inout)     :: self
type(fckit_configuration), intent(in) :: f_conf

!Locals
integer :: iq, nvars
character(kind=c_char,len=:), allocatable :: tracer_variables(:)
character(kind=c_char,len=:), allocatable :: str
character(len=maxvarlen) :: err_msg

  ! Fill in geovars: variables we need from the model
  ! Let users choose specific aerosol species (defined in the yaml file) needed in the AOD calculation.
  ! Followed by slots for delp and RH.

  call f_conf%get_or_die("tracer_geovals",tracer_variables)
  self%ntracers = f_conf%get_size("tracer_geovals")

  ! CMAQ aerosol species 
  do iq = 1, self%ntracers
     call self%geovars%push_back(tracer_variables(iq))      
  enddo

  ! delp and RH 
  call self%geovars%push_back(varindefault)                      

  ! Size of variables (number of obs type (wavelength for AOD))
  nvars = self%obsvars%nvars()

  ! List of wavelengths for the calculation of AOD
  allocate(self%wavelength(nvars))
  call f_conf%get_or_die("wavelengths", self%wavelength)
  
  ! RC File needed for ChemBase 
  call f_conf%get_or_die("RCFile",str)
  self%rcfile = str

  ! Species mapping info, CMAQ to GEOS RC File
  call f_conf%get_or_die("cmaq2geosrc",self%cmaq2geosrc)  
  if(f_conf%get_size("cmaq2geosrc") .ne. self%ntracers) then
  write(err_msg, *) "species mapping information missing for some CMAQ aerosol species"
  call abor1_ftn(err_msg)
  end if

  deallocate(tracer_variables, str)

end subroutine ufo_aodcmaq_setup

! ------------------------------------------------------------------------------
subroutine destructor(self)
implicit none
type(ufo_aodcmaq), intent(inout) :: self

  if (allocated(self%wavelength)) deallocate(self%wavelength)

end subroutine destructor

! ------------------------------------------------------------------------------
subroutine ufo_aodcmaq_simobs(self, geovals, obss, nvars, nlocs, hofx)
use kinds
use ufo_constants_mod, only: grav
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use iso_c_binding

implicit none
class(ufo_aodcmaq), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: aer_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: rh_profile
integer :: nlayers, rc, iq

real(kind_real), dimension(:,:,:), allocatable :: qm    ! aerosol mass mix ratio(ug/kg) that becomes concentration (*delp/g) profiles at obs loc
real(kind_real), dimension(:,:),   allocatable :: rh    ! relative humidity profile interp at obs loc
real(kind_real), dimension(:,:),   allocatable :: delp  ! air pressure thickness profiles at obs loc

character(len=maxvarlen) :: geovar
character(len=maxvarlen), dimension(:), allocatable:: tracer_name

  ! Get delp and rh from model interp at obs loc (from geovals)
  call ufo_geovals_get_var(geovals, var_delp, delp_profile)
  nlayers = delp_profile%nval                            ! number of model layers

  allocate(delp(nlayers,nlocs))
  delp = delp_profile%vals

  ! Get RH from geovals
  allocate(rh(nlayers,nlocs))
  call ufo_geovals_get_var(geovals, var_RH, rh_profile)
  rh = rh_profile%vals

  ! Get aerosol profiles interpolated at obs loc
  allocate(qm(self%ntracers, nlayers, nlocs))
  allocate(tracer_name(self%ntracers))
  do iq = 1, self%ntracers
     geovar = self%geovars%variable(iq)                   !self%geovars contains tracers 
     tracer_name(iq) = geovar
     call ufo_geovals_get_var(geovals, geovar, aer_profile)
     qm(iq,:,:) = aer_profile%vals                   ! aerosol mass mixing ratio
     qm(iq,:,:) = qm(iq,:,:) * delp / grav / 1e9     ! aerosol concentration (kg/m2) 
  enddo
 
  ! Call observation operator code
  ! -----------------------------
  hofx(:,:) = 0.0
  call get_CMAQ_AOD(nlayers, nlocs, nvars, self%ntracers, self%rcfile,  &
                    self%wavelength, self%cmaq2geosrc, qm, rh,       &
                    aod_tot = hofx, rc = rc)  

  ! cleanup memory
  ! --------
  deallocate(qm, rh, delp, tracer_name)

end subroutine ufo_aodcmaq_simobs


! ------------------------------------------------------------------------------

end module ufo_aodcmaq_mod

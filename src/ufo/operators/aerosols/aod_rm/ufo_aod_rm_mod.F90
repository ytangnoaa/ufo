! (C) Copyright 2022 NOAA, EPA, UCAR and George Mason Univeristy.
! This method uses user-assigned aerosol extinciton coefficent per mass 
! in unit m^2/g (ug/m3 -> /Mm ) to construct AOD operator for CMAQ aerosols,
! or the reconstruction method.
! https://doi.org/10.1029/2001JD001409
! https://doi.org/10.1029/2018JD029009 
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for aod_rm observation operator

module ufo_aod_rm_mod

 use iso_c_binding
 use kinds
 use oops_variables_mod
 use ufo_vars_mod
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use missing_values_mod
 use ufo_basis_mod, only: ufo_basis

 implicit none
 private
 integer, parameter :: max_string=800

 type, public :: ufo_aod_rm
 private
   type(oops_variables), public :: obsvars
   type(oops_variables), public :: geovars
   integer, public              :: ntracers
   real(kind_real), public, allocatable :: extincpm(:) ! in unit m^2/g for ug/m3 -> /Mm
 contains
   procedure :: setup  => ufo_aod_rm_setup
   procedure :: simobs => ufo_aod_rm_simobs
   final :: destructor
 end type ufo_aod_rm
!> Default variables required from model
 character(len=maxvarlen), dimension(2), parameter :: varindefault = (/var_delp, var_rh/)

! air humidity scaling factors at 1% RH intervals
         real, public, parameter :: humfac( 99 ) = (/                                &
     	   1.0000e+00,  1.0000e+00,  1.0000e+00,  1.0000e+00,  1.0000e+00,   &
     	   1.0000e+00,  1.0000e+00,  1.0000e+00,  1.0000e+00,  1.0000e+00,   &
     	   1.0000e+00,  1.0000e+00,  1.0000e+00,  1.0001e+00,  1.0001e+00,   &
     	   1.0004e+00,  1.0006e+00,  1.0024e+00,  1.0056e+00,  1.0089e+00,   &
     	   1.0097e+00,  1.0105e+00,  1.0111e+00,  1.0115e+00,  1.0118e+00,   &
     	   1.0122e+00,  1.0126e+00,  1.0130e+00,  1.0135e+00,  1.0139e+00,   &
     	   1.0173e+00,  1.0206e+00,  1.0254e+00,  1.0315e+00,  1.0377e+00,   &
     	   1.0486e+00,  1.0596e+00,  1.0751e+00,  1.0951e+00,  1.1151e+00,   &
     	   1.1247e+00,  1.1343e+00,  1.1436e+00,  1.1525e+00,  1.1615e+00,   &
     	   1.1724e+00,  1.1833e+00,  1.1955e+00,  1.2090e+00,  1.2224e+00,   &
     	   1.2368e+00,  1.2512e+00,  1.2671e+00,  1.2844e+00,  1.3018e+00,   &
     	   1.3234e+00,  1.3450e+00,  1.3695e+00,  1.3969e+00,  1.4246e+00,   &
     	   1.4628e+00,  1.5014e+00,  1.5468e+00,  1.5992e+00,  1.6516e+00,   &
     	   1.6991e+00,  1.7466e+00,  1.7985e+00,  1.8549e+00,  1.9113e+00,   &
     	   1.9596e+00,  2.0080e+00,  2.0596e+00,  2.1146e+00,  2.1695e+00,   &
     	   2.2630e+00,  2.3565e+00,  2.4692e+00,  2.6011e+00,  2.7330e+00,   &
     	   2.8461e+00,  2.9592e+00,  3.0853e+00,  3.2245e+00,  3.3637e+00,   &
     	   3.5743e+00,  3.7849e+00,  4.0466e+00,  4.3594e+00,  4.6721e+00,   &
     	   5.3067e+00,  5.9412e+00,  6.9627e+00,  8.3710e+00,  9.7793e+00,   &
     	   1.2429e+01,  1.5078e+01,  1.8059e+01,  2.1371e+01/)

! Update from Christian Hogrefe (U.S. EPA)
! sea salt humidity scaling factors at 1% RH values
! based on "fss(Rh)" in Table 1 of
! "revised IMPROVE algorithm for estimating light extinction
! from particle speciation data", available at
! http://vista.cira.colostate.edu/improve/Publications/GrayLit/gray_literature.htm
! values for RH96%-99% were linearly extrapolated by using the 94% and 95% values:
!     f(96%)=f(95%)+(f(95%)-f(94%))
!     f(97%)=f(96%)+(f(96%)-f(95%))
!     f(98%)=f(97%)+(f(97%)-f(96%))
!     f(99%)=f(98%)+(f(98%)-f(97%))

         real, public, parameter :: humfac_ss( 99 ) = (/  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 1.00, 1.00, 1.00, 1.00,	  &
     	     1.00, 2.36, 2.38, 2.42, 2.45,	  &
     	     2.48, 2.50, 2.51, 2.53, 2.56,	  &
     	     2.58, 2.59, 2.62, 2.66, 2.69,	  &
     	     2.73, 2.78, 2.83, 2.83, 2.86,	  &
     	     2.89, 2.91, 2.95, 3.01, 3.05,	  &
     	     3.13, 3.17, 3.21, 3.25, 3.27,	  &
     	     3.35, 3.42, 3.52, 3.57, 3.63,	  &
     	     3.69, 3.81, 3.95, 4.04, 4.11,	  &
     	     4.28, 4.49, 4.61, 4.86, 5.12,	  &
     	     5.38, 5.75, 6.17, 6.72, 7.35,	  &
     	     7.98, 8.61, 9.24, 9.87/)
contains

subroutine ufo_aod_rm_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_aod_rm), intent(inout)     :: self
type(fckit_configuration), intent(in) :: f_conf
character(len=MAXVARLEN) :: err_msg

!Locals
integer :: iq, nvars
character(kind=c_char,len=:), allocatable :: tracer_variables(:)
  ! Fill in geovars: variables we need from the model
  ! Need slots for RH and delp 
  ! Let user choose specific aerosols needed in aod calculation.

  call f_conf%get_or_die("tracer_geovals",tracer_variables)
  self%ntracers = f_conf%get_size("tracer_geovals")

  do iq = 1, self%ntracers
     call self%geovars%push_back(tracer_variables(iq))      ! aerosols
  enddo
  call self%geovars%push_back(varindefault)                 ! delp and rh (for concentration)

  ! Species mapping info, CMAQ to GEOS RC File
  call f_conf%get_or_die('extincpm',self%extincpm)  
  if(f_conf%get_size('extincpm') .ne. self%ntracers) then
   write(err_msg, *)'species mapping information missing for some CMAQ aerosol species'
   call abor1_ftn(err_msg)
  endif
  deallocate(tracer_variables)
end subroutine ufo_aod_rm_setup

! ------------------------------------------------------------------------------
! TODO: add cleanup of your observation operator (optional)
subroutine destructor(self)
implicit none
type(ufo_aod_rm), intent(inout) :: self

end subroutine destructor

subroutine ufo_aod_rm_simobs(self, geovals, obss, nvars, nlocs, hofx)
use kinds
use ufo_constants_mod, only: grav
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use iso_c_binding
use obsspace_mod
implicit none
class(ufo_aod_rm), intent(in)    :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: aer_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: rh_profile
integer :: nlayers, rc, iq, k, n

real(kind_real), dimension(:,:,:), allocatable :: qm    ! aer mass mix ratio(g/kg) that becomes concentration (*delp/g) profiles at obs loc
real(kind_real), dimension(:,:),   allocatable :: frh,frh_ss, rh    ! relative humidity profile interp at obs loc
real(kind_real), dimension(:,:),   allocatable :: delp  ! air pressure thickness profiles at obs loc

character(len=MAXVARLEN) :: geovar
character(len=MAXVARLEN), dimension(:), allocatable:: tracer_name

! get some metadata from obsspace
  call ufo_geovals_get_var(geovals, var_delp, delp_profile)
  nlayers = delp_profile%nval                            ! number of model layers

  allocate(delp(nlayers,nlocs))
  delp = delp_profile%vals

  ! Get RH from geovals
  allocate(rh(nlayers,nlocs),frh(nlayers,nlocs),frh_ss(nlayers,nlocs))
  call ufo_geovals_get_var(geovals, var_RH, rh_profile)
  rh = rh_profile%vals
  if (maxval(rh).le.1.1) rh(:,:)=rh(:,:)*100
  rh(:,:)=amin1(amax1(1.,rh(:,:)),99.)
  do n=1, nlocs
   do k=1,nlayers
     frh(k,n)=humfac(int(rh(k,n)))
     frh_ss(k,n)=humfac_ss(int(rh(k,n)))
   enddo
  enddo
     
  ! Get Aer profiles interpolated at obs loc
  allocate(qm(self%ntracers, nlayers, nlocs))
  allocate(tracer_name(self%ntracers))
  hofx(:,:) = 0.0
  do iq = 1, self%ntracers
     geovar = self%geovars%variable(iq)                   !self%geovars contains tracers 
     tracer_name(iq) = geovar
     call ufo_geovals_get_var(geovals, geovar, aer_profile)
     qm(iq,:,:) = aer_profile%vals             ! mass mixing ratio (ug/kg)
     qm(iq,:,:) = qm(iq,:,:) * abs(delp) / grav /1e6     ! dp=-dens*g*dz aer concentration (g/m2) 
     do k=1,nlayers
       if(any(tracer_name(iq)(1:4).eq.['aso4','ano3','anh4'])) then
        hofx(1,:)=hofx(1,:)+qm(iq,k,:)*frh(k,:)*self%extincpm(iq)
       else if(any(tracer_name(iq)(1:3).eq.['acl','ana','ase'])) then
        hofx(1,:)=hofx(1,:)+qm(iq,k,:)*frh_ss(k,:)*self%extincpm(iq)
       else
        hofx(1,:)=hofx(1,:)+qm(iq,k,:)*self%extincpm(iq)
       endif
     enddo  
  enddo
 
  ! cleanup memory
  ! --------
  deallocate(qm, rh, frh, frh_ss, delp, tracer_name)

end subroutine ufo_aod_rm_simobs


! ------------------------------------------------------------------------------

end module ufo_aod_rm_mod

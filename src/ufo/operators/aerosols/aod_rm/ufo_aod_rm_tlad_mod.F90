! (C) Copyright 2022 NOAA, EPA, UCAR and George Mason Univeristy.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for aod_rm tl/ad observation operator

module ufo_aod_rm_tlad_mod

 use iso_c_binding

 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_constants_mod, only: grav

 use ufo_aod_rm_mod, only: humfac,humfac_ss
 use oops_variables_mod
 use ufo_vars_mod

 implicit none
 private
 integer, parameter :: max_string=800
 
 type, public :: ufo_aod_rm_tlad
  integer :: nlocs, nlayers, ntracers, nvars
  type(oops_variables), public :: obsvars
  type(oops_variables), public :: geovars
  real(kind=kind_real), allocatable :: bext(:,:,:,:)
  real(kind=kind_real), public, allocatable :: extincpm(:)
  real(kind=kind_real), dimension(:,:), allocatable :: delp(:,:)
 contains
  procedure :: setup  => ufo_aod_rm_tlad_setup
  procedure :: settraj => ufo_aod_rm_tlad_settraj
  procedure :: simobs_tl  => ufo_aod_rm_simobs_tl
  procedure :: simobs_ad  => ufo_aod_rm_simobs_ad
  procedure :: cleanup => ufo_aod_rm_tlad_cleanup
  final :: destructor
 end type ufo_aod_rm_tlad

contains


subroutine ufo_aod_rm_tlad_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_aod_rm_tlad), intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf

!Locals
integer :: iq
character(len=:), allocatable :: tracer_variables(:)
character(len=maxvarlen) :: err_msg

 call f_conf%get_or_die("tracer_geovals",tracer_variables)
 self%ntracers = f_conf%get_size("tracer_geovals")
 ! CMAQ aerosol species
 do iq = 1, self%ntracers
    call self%geovars%push_back(tracer_variables(iq))	  
 enddo

 self%nvars = self%obsvars%nvars()

 call f_conf%get_or_die("extincpm",self%extincpm)
 if(f_conf%get_size("extincpm") .ne. self%ntracers) then
 write(err_msg, *) "species mapping information missing for some CMAQ aerosol species"
 call abor1_ftn(err_msg)
 end if

 deallocate(tracer_variables)

end subroutine ufo_aod_rm_tlad_setup

subroutine destructor(self)
 implicit none
 type(ufo_aod_rm_tlad), intent(inout) :: self
 call self%cleanup()
end subroutine destructor

subroutine ufo_aod_rm_tlad_cleanup(self)
implicit none
class(ufo_aod_rm_tlad), intent(inout) :: self

 if (allocated(self%bext))    deallocate(self%bext)
 if (allocated(self%delp))    deallocate(self%delp)

end subroutine ufo_aod_rm_tlad_cleanup

subroutine ufo_aod_rm_tlad_settraj(self, geovals, obss)
use iso_c_binding
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
implicit none

class(ufo_aod_rm_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss
!!type(ufo_geovals),       intent(inout) :: hofxdiags    !non-h(x) diagnostics
! Local variables
type(ufo_geoval), pointer :: aer_profile
type(ufo_geoval), pointer :: delp_profile
type(ufo_geoval), pointer :: rh_profile
integer :: rc, iq, n, k
character(len=MAXVARLEN) :: geovar
character(len=MAXVARLEN), dimension(:), allocatable:: tracer_name

real(kind=kind_real), dimension(:,:),   allocatable :: frh,frh_ss, rh

 ! Get number of locations
 self%nlocs = obsspace_get_nlocs(obss)

 ! Get delp and rh from model interp at obs loc (from geovals)
 call ufo_geovals_get_var(geovals, var_delp, delp_profile)
 self%nlayers = delp_profile%nval                                          ! number of model layers
 allocate(self%delp(self%nlayers,self%nlocs))
 self%delp = delp_profile%vals

 ! Get RH from geovals
 allocate(rh(self%nlayers,self%nlocs),frh(self%nlayers,self%nlocs), &
    frh_ss(self%nlayers,self%nlocs))
 call ufo_geovals_get_var(geovals, var_RH, rh_profile)
 rh = rh_profile%vals
 if (maxval(rh).le.1.1) rh(:,:)=rh(:,:)*100
 rh(:,:)=amin1(amax1(1.,rh(:,:)),99.)
  do n=1,self%nlocs
   do k=1,self%nlayers
     frh(k,n)=humfac(int(rh(k,n)))
     frh_ss(k,n)=humfac_ss(int(rh(k,n)))
   enddo
  enddo

! Get Aer profiles interpolated at obs loc
 allocate(tracer_name(self%ntracers))
 allocate(self%bext(self%nlayers, self%nvars, self%ntracers, self%nlocs)) !mass extinction efficiency 
 do iq = 1, self%ntracers
    geovar = self%geovars%variable(iq)                   !self%geovars contains tracers 
    tracer_name(iq) = geovar
     do k=1,self%nlayers
       if(any(tracer_name(iq)(1:4).eq.['aso4','ano3','anh4'])) then
        self%bext(k,1,iq,:)=frh(k,:)*self%extincpm(iq)
       else if(any(tracer_name(iq)(1:3).eq.['acl','ana','ase'])) then
        self%bext(k,1,iq,:)=frh_ss(k,:)*self%extincpm(iq)
       else 	
        self%bext(k,1,iq,:)=self%extincpm(iq)
       endif
     enddo
 enddo

 deallocate(rh,frh,frh_ss)
 deallocate(tracer_name)

end subroutine ufo_aod_rm_tlad_settraj


subroutine ufo_aod_rm_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)
use iso_c_binding
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
implicit none
class(ufo_aod_rm_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss
integer :: iq,k,n
real(kind_real), dimension(:,:,:), allocatable :: qm_tl
type(ufo_geoval), pointer :: aer_profile

character(len=MAXVARLEN) :: geovar

 ! Get Aer profiles interpolated at obs loc
 allocate(qm_tl(self%ntracers, self%nlayers, nlocs))
 hofx=0.
 do iq = 1, self%ntracers
    geovar = self%geovars%variable(iq)                      !self%geovars contains tracers 
    call ufo_geovals_get_var(geovals, geovar, aer_profile)
    qm_tl(iq,:,:) = aer_profile%vals                         ! aer mass mixing ratio
    qm_tl(iq,:,:) = qm_tl(iq,:,:) * self%delp / grav /1e6         ! aer concentration
    do k=1,self%nlayers
     hofx(1,:)=hofx(1,:)+qm_tl(iq,k,:)*self%bext(k,1,iq,:)
    enddo 
 enddo

 deallocate(qm_tl)


end subroutine ufo_aod_rm_simobs_tl


subroutine ufo_aod_rm_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)
use iso_c_binding
use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use obsspace_mod
implicit none
class(ufo_aod_rm_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(in)    :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss

integer :: iq,nv,k,loc
real(kind_real), dimension(:,:,:), allocatable :: qm_ad
type(ufo_geoval), pointer :: aer_profile
character(len=MAXVARLEN) :: geovar

 allocate(qm_ad(self%ntracers, self%nlayers, nlocs))   
 qm_ad=0.
 
 do loc=nlocs,1,-1
  do k=self%nlayers,1,-1
   do iq = self%ntracers,1,-1
    qm_ad(iq,k,loc)=qm_ad(iq,k,loc)+self%bext(k,1,iq,loc)*hofx(1,loc)
   enddo
  enddo
 enddo   
     
 do iq = self%ntracers,1,-1

   geovar = self%geovars%variable(iq)                   !self%geovars contains tracers 
   call ufo_geovals_get_var(geovals, geovar, aer_profile)  ! assign pointer:  aer_profile => geovals(ivar)
      
   qm_ad(iq,:,:) = qm_ad(iq,:,:) * self%delp / grav /1e6
   aer_profile%vals = qm_ad(iq,:,:)

 enddo
deallocate(qm_ad)
end subroutine ufo_aod_rm_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_aod_rm_tlad_mod

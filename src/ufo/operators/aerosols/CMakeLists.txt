# (C) Copyright 2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


add_subdirectory( AOP )
add_subdirectory( AODMetOffice )
add_subdirectory( aod_rm )

if( ${geos-aero_FOUND} )
  add_subdirectory( MR )
  add_subdirectory( aodcmaq )
endif( ${geos-aero_FOUND} )

PREPEND( _p_aodext_files  "operators/aerosols/AOP"  ${aodext_src_files} )
PREPEND( _p_aodmetoffice_files  "operators/aerosols/AODMetOffice"  ${aodmetoffice_src_files} )
PREPEND( _p_aod_rm_files  "operators/aerosols/aod_rm" ${aod_rm_src_files} )

if( ${geos-aero_FOUND} )
  PREPEND( _p_aodgeos_files   "operators/aerosols/MR"   ${aodgeos_src_files} )
  PREPEND( _p_aodcmaq_files   "operators/aerosols/aodcmaq"  ${aodcmaq_src_files} )
endif( ${geos-aero_FOUND} )

set ( aerosols_src_files
    ${_p_aodext_files}
    ${_p_aodgeos_files}
    ${_p_aodmetoffice_files}
    ${_p_aod_rm_files}
    ${_p_aodcmaq_files}
    PARENT_SCOPE
)

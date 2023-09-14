/*
 * (C) Copyright 2021.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_AEROSOLS_AODCMAQ_OBSAODCMAQPARAMETERS_H_
#define UFO_AEROSOLS_AODCMAQ_OBSAODCMAQPARAMETERS_H_

#include <string>
#include <vector>

// List of Parameter classes
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// Configuration options recognized by the AodCmaq operator.
class ObsAodCmaqParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAodCmaqParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter <std::vector<std::string>> tracerGeovals
    {"tracer_geovals",
     "Names of model tracer variables",
     this};

  oops::RequiredParameter <std::vector<int>> Cmaq2Geosrc
    {"cmaq2geosrc",
     "CMAQ to GEOS aero: "
     "must be within 1-13",
     this};

  oops::Parameter <std::string> Rcfile
    {"RCFile",
     "GEOS resource file",
     "geosaod.rc", this};

  oops::RequiredParameter <std::vector<float>> Wavelengths
    {"wavelengths",
     "wavelengths for AOD",
     this};
};

}  // namespace ufo
#endif  // UFO_AEROSOLS_AODCMAQ_OBSAODCMAQPARAMETERS_H_

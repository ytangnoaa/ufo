/*
 * (C) Copyright 2022 NOAA, EPA, UCAR and George Mason Univeristy.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_AEROSOLS_AOD_RM_OBSAOD_RMPARAMETERS_H_
#define UFO_AEROSOLS_AOD_RM_OBSAOD_RMPARAMETERS_H_

#include <string>
#include <vector>

// TODO: modify the list of Parameter classes to include
// depending on what is used below.
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// Configuration options recognized by the AOD_RM operator.
class ObsAOD_RMParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsAOD_RMParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter <std::vector<std::string>> tracerGeovals
    {"tracer_geovals",
     "Names of model tracer variables",
     this};

  oops::RequiredParameter <std::vector<float>> extincpm
    {"extincpm",
     "Extinction per mass (m^2/g):",
     this};
};

}  // namespace ufo
#endif  // UFO_AEROSOLS_AOD_RM_OBSAOD_RMPARAMETERS_H_

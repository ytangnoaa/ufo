/*
 * (C) Copyright 2022 NOAA, EPA, UCAR and George Mason Univeristy.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_AEROSOLS_AOD_RM_OBSAOD_RM_H_
#define UFO_AEROSOLS_AOD_RM_OBSAOD_RM_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"

#include "ufo/operators/aerosols/aod_rm/ObsAOD_RM.interface.h"
#include "ufo/operators/aerosols/aod_rm/ObsAOD_RMParameters.h"
#include "ufo/ObsOperatorBase.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// AOD_RM observation operator class
class ObsAOD_RM : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAOD_RM> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsAOD_RMParameters Parameters_;

  static const std::string classname() {return "ufo::ObsAOD_RM";}

  ObsAOD_RM(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAOD_RM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_AEROSOLS_AOD_RM_OBSAOD_RM_H_

/*
 * (C) Copyright 2022 NOAA, EPA, UCAR and George Mason Univeristy.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_AEROSOLS_AOD_RM_OBSAOD_RMTLAD_H_
#define UFO_AEROSOLS_AOD_RM_OBSAOD_RMTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/operators/aerosols/aod_rm/ObsAOD_RMParameters.h"
#include "ufo/operators/aerosols/aod_rm/ObsAOD_RMTLAD.interface.h"
#include "ufo/LinearObsOperatorBase.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// AOD_RM TL/AD observation operator class
class ObsAOD_RMTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsAOD_RMTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsAOD_RMParameters Parameters_;

  static const std::string classname() {return "ufo::ObsAOD_RMTLAD";}

  ObsAOD_RMTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAOD_RMTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_AEROSOLS_AOD_RM_OBSAOD_RMTLAD_H_

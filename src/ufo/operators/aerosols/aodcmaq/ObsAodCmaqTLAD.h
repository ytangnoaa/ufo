/*
 * (C) Copyright 2021.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_AEROSOLS_AODCMAQ_OBSAODCMAQTLAD_H_
#define UFO_AEROSOLS_AODCMAQ_OBSAODCMAQTLAD_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/operators/aerosols/aodcmaq/ObsAodCmaqParameters.h"
#include "ufo/operators/aerosols/aodcmaq/ObsAodCmaqTLAD.interface.h"
#include "ufo/LinearObsOperatorBase.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// AodCmaq TL/AD observation operator class
class ObsAodCmaqTLAD : public LinearObsOperatorBase,
                       private util::ObjectCounter<ObsAodCmaqTLAD> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the LinearObsOperatorFactory.
  typedef ObsAodCmaqParameters Parameters_;

  static const std::string classname() {return "ufo::ObsAodCmaqTLAD";}

  ObsAodCmaqTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodCmaqTLAD();

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
#endif  // UFO_AEROSOLS_AODCMAQ_OBSAODCMAQTLAD_H_

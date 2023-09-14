/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_AEROSOLS_AODCMAQ_OBSAODCMAQ_H_
#define UFO_AEROSOLS_AODCMAQ_OBSAODCMAQ_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/operators/aerosols/aodcmaq/ObsAodCmaq.interface.h"
#include "ufo/operators/aerosols/aodcmaq/ObsAodCmaqParameters.h"
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
/// AodCmaq observation operator class
class ObsAodCmaq : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAodCmaq> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsAodCmaqParameters Parameters_;

  static const std::string classname() {return "ufo::ObsAodCmaq";}

  ObsAodCmaq(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsAodCmaq();

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
#endif  // UFO_AEROSOLS_AODCMAQ_OBSAODCMAQ_H_

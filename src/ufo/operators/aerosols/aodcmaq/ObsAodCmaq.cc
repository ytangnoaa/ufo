/*
 * (C) Copyright 2021.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/aerosols/aodcmaq/ObsAodCmaq.h"

#include <ostream>
#include <set>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsAodCmaq> makerAodCmaq_("AodCmaq");
// -----------------------------------------------------------------------------

ObsAodCmaq::ObsAodCmaq(const ioda::ObsSpace & odb,
                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  ufo_aodcmaq_setup_f90(keyOper_, parameters.toConfiguration(), odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsAodCmaq created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodCmaq::~ObsAodCmaq() {
  ufo_aodcmaq_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAodCmaq destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCmaq::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                             ObsDiagnostics &) const {
  ufo_aodcmaq_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsAodCmaq: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCmaq::print(std::ostream & os) const {
  os << "ObsAodCmaq::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo

/*
 * (C) Copyright 2022 NOAA, EPA, UCAR and George Mason Univeristy.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/aerosols/aod_rm/ObsAOD_RM.h"

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
static ObsOperatorMaker<ObsAOD_RM> makerAOD_RM_("AOD_RM");
// -----------------------------------------------------------------------------

ObsAOD_RM::ObsAOD_RM(const ioda::ObsSpace & odb,
                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOper_(0), odb_(odb), varin_()
{
  ufo_aod_rm_setup_f90(keyOper_, parameters.toConfiguration(), odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsAOD_RM created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsAOD_RM::~ObsAOD_RM() {
  ufo_aod_rm_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAOD_RM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAOD_RM::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                             ObsDiagnostics &) const {
  ufo_aod_rm_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsAOD_RM: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAOD_RM::print(std::ostream & os) const {
  os << "ObsAOD_RM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo

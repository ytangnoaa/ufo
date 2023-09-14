/*
 * (C) Copyright 2022 NOAA, EPA, UCAR and George Mason Univeristy.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/aerosols/aod_rm/ObsAOD_RMTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAOD_RMTLAD> makerAOD_RMTL_("AOD_RM");
// -----------------------------------------------------------------------------

ObsAOD_RMTLAD::ObsAOD_RMTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  ufo_aod_rm_tlad_setup_f90(keyOper_, parameters.toConfiguration(), odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsAOD_RMTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAOD_RMTLAD::~ObsAOD_RMTLAD() {
  ufo_aod_rm_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAOD_RMTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAOD_RMTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics & ) {
  ufo_aod_rm_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsAOD_RMTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAOD_RMTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_aod_rm_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAOD_RMTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAOD_RMTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_aod_rm_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAOD_RMTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAOD_RMTLAD::print(std::ostream & os) const {
  os << "ObsAOD_RMTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
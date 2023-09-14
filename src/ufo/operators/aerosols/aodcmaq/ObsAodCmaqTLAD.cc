/*
 * (C) Copyright 2021.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/aerosols/aodcmaq/ObsAodCmaqTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsAodCmaqTLAD> makerAodCmaqTL_("AodCmaq");
// -----------------------------------------------------------------------------

ObsAodCmaqTLAD::ObsAodCmaqTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), keyOper_(0), varin_()
{
  ufo_aodcmaq_tlad_setup_f90(keyOper_, parameters.toConfiguration(), odb.obsvariables(), varin_);
  oops::Log::trace() << "ObsAodCmaqTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsAodCmaqTLAD::~ObsAodCmaqTLAD() {
  ufo_aodcmaq_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsAodCmaqTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCmaqTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &) {
  ufo_aodcmaq_tlad_settraj_f90(keyOper_, geovals.toFortran(), obsspace());
  oops::Log::trace() << "ObsAodCmaqTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCmaqTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec) const {
  ufo_aodcmaq_simobs_tl_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodCmaqTLAD: TL observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCmaqTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec) const {
  ufo_aodcmaq_simobs_ad_f90(keyOper_, geovals.toFortran(), obsspace(),
                            ovec.nvars(), ovec.nlocs(), ovec.toFortran());
  oops::Log::trace() << "ObsAodCmaqTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsAodCmaqTLAD::print(std::ostream & os) const {
  os << "ObsAodCmaqTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

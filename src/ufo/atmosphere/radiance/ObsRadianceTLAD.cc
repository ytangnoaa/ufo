/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/atmosphere/radiance/ObsRadianceTLAD.h"

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/utils/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadianceTLAD> makerAmsuaTL_("AMSU-A");
static LinearObsOperatorMaker<ObsRadianceTLAD> makerAvhrrTL_("AVHRR");
// -----------------------------------------------------------------------------

ObsRadianceTLAD::ObsRadianceTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperRadiance_(0), varin_(), odb_(odb)
{
  const std::vector<std::string> vv{"air_temperature"};
  varin_.reset(new oops::Variables(vv));

  // parse channels from the config and create variable names
  std::string chlist = config.getString("channels");
  std::set<int> channels = parseIntSet(chlist);
  channels_.reserve(channels.size());
  for (const int jj : channels) {
    channels_.push_back(jj);
  }

  const eckit::LocalConfiguration obsOptions(config, "ObsOptions");
  const eckit::Configuration * configc = &obsOptions;
  ufo_radiance_tlad_setup_f90(keyOperRadiance_, &configc);
  oops::Log::trace() << "ObsRadianceTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceTLAD::~ObsRadianceTLAD() {
  ufo_radiance_tlad_delete_f90(keyOperRadiance_);
  oops::Log::trace() << "ObsRadianceTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_radiance_tlad_settraj_f90(keyOperRadiance_, geovals.toFortran(), odb_,
                                channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsRadianceTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                    const ObsBiasIncrement & bias) const {
  ufo_radiance_simobs_tl_f90(keyOperRadiance_, geovals.toFortran(), odb_,
                             ovec.size(), ovec.toFortran(),
                             channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsRadianceTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                    ObsBiasIncrement & bias) const {
  ufo_radiance_simobs_ad_f90(keyOperRadiance_, geovals.toFortran(), odb_,
                             ovec.size(), ovec.toFortran(),
                             channels_.size(), channels_[0]);
}

// -----------------------------------------------------------------------------

void ObsRadianceTLAD::print(std::ostream & os) const {
  os << "ObsRadianceTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

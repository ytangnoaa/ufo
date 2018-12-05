/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/marine/insitutemperature/ObsInsituTemperature.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"


namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsInsituTemperature> makerInsituTemperature_("InsituTemperature");
// -----------------------------------------------------------------------------

ObsInsituTemperature::ObsInsituTemperature(const ioda::ObsSpace & odb,
                                           const eckit::Configuration & config)
  : keyOper_(0), odb_(odb), varin_(), varout_()
{
  const std::vector<std::string> vvin{"ocean_potential_temperature",
                                      "ocean_salinity", "ocean_layer_thickness"};
  varin_.reset(new oops::Variables(vvin));
  const std::vector<std::string> vvout{"insitu_temperature"};
  varout_.reset(new oops::Variables(vvout));
  const eckit::Configuration * configc = &config;
  ufo_insitutemperature_setup_f90(keyOper_, &configc);
  oops::Log::trace() << "ObsInsituTemperature created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsInsituTemperature::~ObsInsituTemperature() {
  ufo_insitutemperature_delete_f90(keyOper_);
  oops::Log::trace() << "ObsInsituTemperature destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsInsituTemperature::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                              const ObsBias & bias) const {
  ufo_insitutemperature_simobs_f90(keyOper_, gv.toFortran(), odb_, ovec.size(), ovec.toFortran(),
                      bias.toFortran());
  oops::Log::trace() << "ObsInsituTemperature: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

Locations * ObsInsituTemperature::locateObs(const util::DateTime & t1,
                                            const util::DateTime & t2) const {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  int keylocs;
  ufo_insitutemperature_locateobs_f90(keyOper_, odb_, &p1, &p2, keylocs);

  return new Locations(keylocs);
}

// -----------------------------------------------------------------------------

void ObsInsituTemperature::print(std::ostream & os) const {
  os << "ObsInsituTemperature::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo

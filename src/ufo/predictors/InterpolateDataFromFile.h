/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_INTERPOLATEDATAFROMFILE_H_
#define UFO_PREDICTORS_INTERPOLATEDATAFROMFILE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/obsfunctions/DrawValueFromFile.h"
#include "ufo/predictors/PredictorBase.h"

namespace ufo {

class VariableCorrectionParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableCorrectionParameters, Parameters)

 public:
  /// Name of the variable to be bias-corrected.
  oops::RequiredParameter<std::string> name{"name", this};

  /// Options controlling bias correction of this variable.
  DrawValueFromFileParametersWithoutGroup details{this};
};

/// \brief Parameters recognized in the `options` section of the Configuration passed to
/// the InterpolateDataFromFile constructor.
class InterpolateDataFromFileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InterpolateDataFromFileParameters, Parameters)

 public:
  oops::Parameter<std::vector<VariableCorrectionParameters>> correctedVariables{
      "corrected variables", {}, this};
};

/// \brief A predictor returning values produced by the DrawValueFromFile ObsFunction.
///
/// Let's start with a simple example. Suppose this predictor is configured with the following
/// YAML snippet:
///
/// \code{.yaml}
/// name: interpolate_data_from_file
/// options:
///   corrected variables:
///   - name: air_temperature
///     file: Data/ufo/testinput_tier_1/air_temperature_bias.nc4
///     interpolation:
///     - name: station_id@MetaData
///       method: exact
/// \endcode
///
/// and the `air_temperature_bias.nc4` file contains a 1D array `air_temperature@ObsBias` indexed
/// by a `station_id@MetaData` coordinate. The DrawValueFromFile ObsFunction will then load this
/// file and for each location produce the element of the `air_temperature@ObsBias` array
/// corresponding to the element of the `station_id@MetaData` array matching the value of the
/// `station_id@MetaData` ObsSpace variable at that location. This will also be the value produced
/// by this predictor.
///
/// The predictor will produce zeros for all bias-corrected variables missing from the `corrected
/// variables` list.
///
/// It is possible to make the bias correction dependent on more than one ObsSpace variable and to
/// use a different matching method than `exact`. For more details, see the documentation of
/// DrawValueFromFile and DataExtractor.
class InterpolateDataFromFile : public PredictorBase {
 public:
  InterpolateDataFromFile(const eckit::Configuration &, const oops::Variables &);

  void compute(const ioda::ObsSpace &, const GeoVaLs &,
               const ObsDiagnostics &, ioda::ObsVector &) const override;

 private:
  /// `obsFunctions_[varName]` is the ObsFunction that will calculate the predictions for variable
  /// `varName`.
  // The map is storing unique_ptrs to make it possible to compile this code with GCC 4.8.5,
  // whose STL implementation is affected by LWG 2397
  // (http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-defects.html#2397), later resolved by
  // amending the C++11 standard as described in N4387
  // (http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/n4387.html).
  std::map<std::string, std::unique_ptr<DrawValueFromFile>> obsFunctions_;
};

}  // namespace ufo

#endif  // UFO_PREDICTORS_INTERPOLATEDATAFROMFILE_H_

/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORCSVBACKEND_H_
#define UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORCSVBACKEND_H_

#include <string>

#include "ufo/utils/dataextractor/DataExtractorBackend.h"

namespace ufo
{

/// \brief Produces input for a DataExtractor by loading data from a CSV file.
///
/// The file should have the following structure (described in more detail in the example below):
/// * First line: comma-separated column names in ioda-v1 style (`var@Group`) or ioda-v2 style
///   (`Group/var`)
/// * Second line: comma-separated column data types (datetime, float, int or string)
/// * Further lines: comma-separated data entries
/// The number of entries in each line should be the same.
///
/// Here's an example of a file that could be read by this backend and used for bias correction:
///
///     station_id@MetaData,air_pressure@MetaData,air_temperature@ObsBias
///     string,float,float
///     ABC,30000,0.1
///     ABC,60000,0.2
///     ABC,90000,0.3
///     XYZ,40000,0.4
///     XYZ,80000,0.5
///
/// One of the columns (above, air_temperature@ObsBias) contains the values to be extracted (also
/// known as the _payload_). The payload column is identified by the group it belongs to, i.e. the
/// part of its name following the `@` sign (ioda-v1 style) or preceding the last `/` sign (ioda-v2
/// style); this group is specified in the call to the loadData() member function. The values from
/// the other columns (_coordinates_) are compared against ObsSpace variables with the same names
/// to determine the row or rows from which the payload value should be extracted for each
/// observation. The details of this comparison (e.g. whether an exact match is required, the
/// nearest match is used, or piecewise linear interpolation is performed) depend on how the class
/// using the extracted data (e.g. the DrawValueFromFile ObsFunction) is configured. The data type
/// of each column must match the data type of the corresponding ObsSpace variable. The payload
/// column must be of type `float` or `int`. The column order does not matter.
///
/// Special case: a column containing channel numbers (which aren't stored in a separate ObsSpace
/// variable) should be labelled `channel_number@MetaData` or `MetaData/channel_number`.
///
/// To continue the example above, suppose the file shown earlier is passed to the
/// DrawValueFromFile ObsFunction configured in the following way:
///
///     name: DrawValueFromFile@ObsFunction
///     options:
///       file: ...       # path to the CSV file
///       group: ObsBias  # group with the payload variable
///       interpolation:
///       - name: station_id@MetaData
///         method: exact
///       - name: air_pressure@MetaData
///         method: linear
///
/// For an observation taken by station XYZ at pressure 60000 the function would be evaluated in the
/// following way:
/// * First, find all rows in the CSV file with a value of `XYZ` in the `station_id@MetaData`
///   column.
/// * Then take the values of the `air_pressure@MetaData` and `air_temperature@ObsBias` columns
///   in these rows and use them to construct a piecewise linear interpolant. Evaluate this
///   interpolant at pressure 60000. This produces the value of 0.45.
///
/// Refer to the documentation of the DrawValueFromFile ObsFunction for more information about the
/// available extraction methods.
class DataExtractorCSVBackend : public DataExtractorBackend {
 public:
  /// \brief Create a new instance.
  ///
  /// \param filepath Path to the CSV file that will be read by loadData().
  explicit DataExtractorCSVBackend(const std::string &filepath);

  DataExtractorInput loadData(const std::string &payloadGroup) const override;

 private:
  std::string filepath_;
};

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_DATAEXTRACTORCSVBACKEND_H_

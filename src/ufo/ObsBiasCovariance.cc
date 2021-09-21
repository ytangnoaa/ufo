/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <Eigen/Core>
#include <cmath>
#include <fstream>
#include <memory>
#include <random>
#include <set>

#include "ufo/ObsBiasCovariance.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/Engines/Factory.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/predictors/PredictorBase.h"
#include "ufo/utils/IodaGroupIndices.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasCovariance::ObsBiasCovariance(ioda::ObsSpace & odb,
                                     const Parameters_ & params)
  : odb_(odb), prednames_(0), vars_(odb.obsvariables()), variances_(),
    preconditioner_(0),
    ht_rinv_h_(0), obs_num_(0), analysis_variances_(0), minimal_required_obs_number_(0) {
  oops::Log::trace() << "ObsBiasCovariance::Constructor starting" << std::endl;

  // Predictor factory
  for (const PredictorParametersWrapper &wrapper :
       params.variationalBC.value().predictors.value()) {
    std::shared_ptr<PredictorBase> pred(PredictorFactory::create(wrapper.predictorParameters,
                                                                 vars_));
    prednames_.push_back(pred->name());
  }

  if (prednames_.size()*vars_.size() > 0) {
    if (params.covariance.value() == boost::none)
      throw eckit::UserError("obs bias.covariance section missing from the YAML file");
    const ObsBiasCovarianceParameters &biasCovParams = *params.covariance.value();

    // Get the minimal required filtered obs number
    minimal_required_obs_number_ = biasCovParams.minimalRequiredObsNumber;

    // Override the variance range if provided
    {
      const std::vector<double> &range = biasCovParams.varianceRange.value();
      ASSERT(range.size() == 2);
      smallest_variance_ = range[0];
      largest_variance_ = range[1];
    }

    // Override the preconditioning step size if provided
    step_size_ = biasCovParams.stepSize;

    // Override the largest analysis variance if provided
    largest_analysis_variance_ = biasCovParams.largestAnalysisVariance;

    // Initialize the variances to upper limit
    variances_ = Eigen::VectorXd::Constant(prednames_.size()*vars_.size(), largest_variance_);

    // Initialize the hessian contribution to zero
    ht_rinv_h_.resize(prednames_.size() * vars_.size());
    std::fill(ht_rinv_h_.begin(), ht_rinv_h_.end(), 0.0);

    // Initialize the preconditioner to default step size
    preconditioner_.resize(prednames_.size() * vars_.size());
    std::fill(preconditioner_.begin(), preconditioner_.end(), step_size_);

    // Initialize obs_num_ to ZERO
    obs_num_.resize(vars_.size());
    std::fill(obs_num_.begin(), obs_num_.end(), 0);

    // Initialize analysis error variances to the upper limit
    analysis_variances_.resize(prednames_.size() * vars_.size());
    std::fill(analysis_variances_.begin(), analysis_variances_.end(), largest_variance_);

    // Initializes from given prior
    if (biasCovParams.prior.value() != boost::none) {
      const ObsBiasCovariancePriorParameters &priorParams = *biasCovParams.prior.value();

      // Get default inflation ratio
      const double inflation_ratio = priorParams.inflation.value().ratio;

      // Check the large inflation ratio when obs number < minimal_required_obs_number
      const double large_inflation_ratio = priorParams.inflation.value().ratioForSmallDataset;

      // read in Variances prior (analysis_variances_) and number of obs. (obs_num_)
      // from previous cycle
      this->read(priorParams);

      // set variances for bias predictor coeff. based on diagonal info
      // of previous analysis error variance
      std::size_t ii;
      for (std::size_t j = 0; j < vars_.size(); ++j) {
        const double inflation = (obs_num_[j] <= minimal_required_obs_number_) ?
                                 large_inflation_ratio : inflation_ratio;
        for (std::size_t p = 0; p < prednames_.size(); ++p) {
          ii = j*prednames_.size() + p;
          if (inflation > inflation_ratio)
            analysis_variances_[ii] = inflation * analysis_variances_[ii] + smallest_variance_;
          variances_[ii] = inflation * analysis_variances_[ii] + smallest_variance_;
          if (variances_[ii] > largest_variance_) variances_[ii] = largest_variance_;
          if (analysis_variances_[ii] > largest_analysis_variance_)
            analysis_variances_[ii] = largest_analysis_variance_;
        }
      }
    }
  }

  oops::Log::trace() << "ObsBiasCovariance::Constructor is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::read(const ObsBiasCovariancePriorParameters & params) {
  oops::Log::trace() << "ObsBiasCovariance::read from file " << std::endl;

  if (params.inputFile.value() != boost::none) {
    // Open an hdf5 file, read only
    ioda::Engines::BackendNames  backendName = ioda::Engines::BackendNames::Hdf5File;
    ioda::Engines::BackendCreationParameters backendParams;
    backendParams.fileName = *params.inputFile.value();
    backendParams.action   = ioda::Engines::BackendFileActions::Open;
    backendParams.openMode = ioda::Engines::BackendOpenModes::Read_Only;

    // Create the backend and attach it to an ObsGroup
    // Use the None DataLyoutPolicy for now to accommodate the current file format
    ioda::Group backend = constructBackend(backendName, backendParams);
    ioda::ObsGroup obsgroup = ioda::ObsGroup(backend,
                   ioda::detail::DataLayoutPolicy::generate(
                         ioda::detail::DataLayoutPolicy::Policies::None));

    // Read coefficients error variances into the Eigen array
    ioda::Variable bcerrvar = obsgroup.vars["bias_coeff_errors"];
    Eigen::ArrayXXf allbcerrors;
    bcerrvar.readWithEigenRegular(allbcerrors);

    // Read nobs into Eigen array
    ioda::Variable nobsvar = obsgroup.vars["number_obs_assimilated"];
    Eigen::ArrayXf nobsassim;
    nobsvar.readWithEigenRegular(nobsassim);

    // Find indices of predictors and variables/channels that we need in the data read from the file
    const std::vector<int> pred_idx = getRequiredVariableIndices(obsgroup, "predictors",
                                              prednames_.begin(), prednames_.end());
    const std::vector<int> var_idx = getRequiredVarOrChannelIndices(obsgroup, vars_);

    // Filter predictors and channels that we need
    // FIXME: may be possible by indexing allbcerrors(pred_idx, chan_idx) when Eigen 3.4
    // is available
    for (size_t jvar = 0; jvar < var_idx.size(); ++jvar) {
      obs_num_[jvar] = nobsassim(var_idx[jvar]);
      for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
        analysis_variances_[jvar*pred_idx.size()+jpred] =
             allbcerrors(pred_idx[jpred], var_idx[jvar]);
      }
    }
  }

  oops::Log::trace() << "ObsBiasCovariance::read is done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::write(const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsBiasCovariance::write to file " << std::endl;
  oops::Log::trace() << "ObsBiasCovariance::write is not implemented " << std::endl;
  oops::Log::trace() << "ObsBiasCovariance::write is done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::linearize(const ObsBias & bias, const eckit::Configuration & innerConf) {
  oops::Log::trace() << "ObsBiasCovariance::linearize starts" << std::endl;
  if (bias) {
    // Retrieve the QC flags and do statistics from second outer loop
    const int jouter = innerConf.getInt("iteration");
    if (jouter >= 1) {
      std::unique_ptr<ioda::Accumulator<std::vector<size_t>>> obs_num_accumulator =
          odb_.distribution()->createAccumulator<size_t>(obs_num_.size());

      // Retrieve the QC flags of previous outer loop and recalculate the number of effective obs.
      const std::string qc_group_name = "EffectiveQC" + std::to_string(jouter-1);
      const std::vector<std::string> vars = odb_.obsvariables().variables();
      std::vector<int> qc_flags(odb_.nlocs(), 999);
      for (std::size_t jvar = 0; jvar < vars.size(); ++jvar) {
        if (odb_.has(qc_group_name, vars[jvar])) {
          odb_.get_db(qc_group_name, vars[jvar], qc_flags);
          for (std::size_t jloc = 0; jloc < qc_flags.size(); ++jloc)
            if (qc_flags[jloc] == 0)
              obs_num_accumulator->addTerm(jloc, jvar, 1);
        } else {
          throw eckit::UserError("Unable to find QC flags : " + vars[jvar] + "@" + qc_group_name);
        }
      }

      // Sum across the processors
      obs_num_ = obs_num_accumulator->computeResult();

      const float missing = util::missingValue(missing);

      // compute the hessian contribution from Jo bias terms channel by channel
      // retrieve the effective error (after QC) for this channel
      const std::string err_group_name = "EffectiveError" + std::to_string(jouter-1);
      ioda::ObsVector r_inv(odb_, err_group_name);

      // compute \mathrm{R}^{-1}
      std::size_t nvars = r_inv.nvars();
      for (size_t vv = 0; vv < nvars; ++vv) {
        for (size_t ii = 0; ii < r_inv.nlocs(); ++ii) {
          if (r_inv[ii*nvars + vv] != missing) {
            r_inv[ii*nvars + vv] = 1.0f / pow(r_inv[ii*nvars + vv], 2);
          } else {
            r_inv[ii*nvars + vv] = 0.0f;
          }
        }
      }

      // compute \mathrm{H}_\beta^\intercal \mathrm{R}^{-1} \mathrm{H}_\beta
      // -----------------------------------------
      std::unique_ptr<ioda::Accumulator<std::vector<double>>> ht_rinv_h_accumulator =
          odb_.distribution()->createAccumulator<double>(ht_rinv_h_.size());
      for (std::size_t p = 0; p < prednames_.size(); ++p) {
        // retrieve the predictors
        const ioda::ObsVector predx(odb_, prednames_[p] + "Predictor");

        // for each variable
        ASSERT(r_inv.nlocs() == predx.nlocs());
        std::size_t nvars = predx.nvars();
        // only keep the diagnoal
        for (size_t vv = 0; vv < nvars; ++vv) {
          for (size_t ii = 0; ii < predx.nlocs(); ++ii)
            ht_rinv_h_accumulator->addTerm(ii,
                                           vv*prednames_.size() + p,
                                           pow(predx[ii*nvars + vv], 2) * r_inv[ii*nvars + vv]);
        }
      }

      // Sum the hessian contributions across the tasks
      ht_rinv_h_ = ht_rinv_h_accumulator->computeResult();
    }

    // reset variances for bias predictor coeff. based on current data count
    for (std::size_t j = 0; j < obs_num_.size(); ++j) {
      if (obs_num_[j] <= minimal_required_obs_number_) {
        for (std::size_t p = 0; p < prednames_.size(); ++p)
          variances_[j*prednames_.size() + p] = smallest_variance_;
      }
    }

    // set a coeff. factor for variances of control variables
    for (std::size_t j = 0; j < vars_.size(); ++j) {
      for (std::size_t p = 0; p < prednames_.size(); ++p) {
        const std::size_t index = j*prednames_.size() + p;
        preconditioner_[index] = step_size_;
        // L = \mathrm{A}^{-1}
        if (obs_num_[j] > 0)
          preconditioner_[index] = 1.0 / (1.0 / variances_[index] + ht_rinv_h_[index]);
        if (obs_num_[j] > minimal_required_obs_number_) {
          if (ht_rinv_h_[index] > 0.0) {
            analysis_variances_[index] = 1.0 / (1.0 / variances_[index] + ht_rinv_h_[index]);
          } else {
            analysis_variances_[index] = largest_analysis_variance_;
          }
        }
      }
    }
  }
  oops::Log::trace() << "ObsBiasCovariance::linearize is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::multiply(const ObsBiasIncrement & dx1,
                                 ObsBiasIncrement & dx2) const {
  oops::Log::trace() << "ObsBiasCovariance::multiply starts" << std::endl;

  dx2.data().array() = dx1.data().array() * variances_.array();

  oops::Log::trace() << "ObsBiasCovariance::multiply is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::inverseMultiply(const ObsBiasIncrement & dx1,
                                        ObsBiasIncrement & dx2) const {
  oops::Log::trace() << "ObsBiasCovariance::inverseMultiply starts" << std::endl;

  dx2.data().array() = dx1.data().array() / variances_.array();

  oops::Log::trace() << "ObsBiasCovariance::inverseMultiply is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::randomize(ObsBiasIncrement & dx) const {
  oops::Log::trace() << "ObsBiasCovariance::randomize starts" << std::endl;
  if (dx) {
    static util::NormalDistribution<double> dist(variances_.size());
    for (std::size_t jj = 0; jj < variances_.size(); ++jj) {
      dx.data()[jj] = dist[jj] * std::sqrt(variances_[jj]);
    }
  }
  oops::Log::trace() << "ObsBiasCovariance::randomize is done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

#include <atomic>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "asa103.hpp"

// C++ string formatting library
//#include "spdlog/fmt/fmt.h"

#include "Optimizer.hpp"
#include "eigen3/Eigen/Dense"

// intelligently chosen value adopted from
// https://github.com/pachterlab/kallisto/blob/master/src/EMAlgorithm.h#L18
constexpr double minEQClassWeight = std::numeric_limits<double>::denorm_min();
constexpr double minWeight = std::numeric_limits<double>::denorm_min();
// A bit more conservative of a minimum as an argument to the digamma function.
constexpr double digammaMin = 1e-2;

template <typename VecT>
double truncateCountVector(VecT& alphas, double cutoff) {
  // Truncate tiny expression values
  double alphaSum = 0.0;

  for (size_t i = 0; i < alphas.size(); ++i) {
    if (alphas[i] <= cutoff) {
      alphas[i] = 0.0;
    }
    alphaSum += alphas[i];
  }
  return alphaSum;
}

template <typename VecT>
double truncateCountVector(VecT& alphas, std::vector<double>& cutoff) {
  // Truncate tiny expression values
  double alphaSum = 0.0;

  for (size_t i = 0; i < alphas.size(); ++i) {
    if (alphas[i] <= cutoff[i]) {
      alphas[i] = 0.0;
    }
    alphaSum += alphas[i];
  }
  return alphaSum;
}

/**
 *  Populate the prior parameters for the VBEM
 *  Note: effLens *must* be valid before calling this function.
 */
std::vector<double> populatePriorAlphas_(
    Eigen::VectorXd& effLens, // current effective length estimate
    double priorValue,        // the per-nucleotide prior value to use
    bool perTranscriptPrior   // true if prior is per-txp, else per-nucleotide
    ) {
  // start out with the per-txp prior
  std::vector<double> priorAlphas(effLens.size(), priorValue);

  // If the prior is per-nucleotide (default, then we need a potentially
  // different
  // value for each transcript based on its length).
  if (!perTranscriptPrior) {
    for (size_t i = 0; i < effLens.size(); ++i) {
      priorAlphas[i] = priorValue * effLens(i);
    }
  }
  return priorAlphas;
}

/**
 * Single-threaded EM-update routine for use in bootstrapping
 */
template <typename VecT>
void EMUpdate_(std::vector<std::vector<uint32_t>>& txpGroupLabels,
               std::vector<std::vector<double>>& txpGroupCombinedWeights,
               std::vector<size_t>& txpGroupCounts, const VecT& alphaIn,
               VecT& alphaOut) {

  size_t N = alphaIn.size();
  assert(alphaIn.size() == alphaOut.size());

  size_t numEqClasses = txpGroupLabels.size();
  for (size_t eqID = 0; eqID < numEqClasses; ++eqID) {
    uint64_t count = txpGroupCounts[eqID];
    // for each transcript in this class
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    double denom = 0.0;
    size_t groupSize = txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (groupSize > 1) {
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        double v = alphaIn[tid] * aux;
        denom += v;
      }

      if (denom <= ::minEQClassWeight) {
        // tgroup.setValid(false);
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          double v = alphaIn[tid] * aux;
          if (!std::isnan(v)) {
            alphaOut[tid] += v * invDenom;
          }
        }
      }
    } else {
      alphaOut[txps.front()] += count;
    }
  }
}

/*
void FixCounts_(
                 std::vector<std::vector<uint32_t>>& txpGroupLabels,
                 std::vector<std::vector<double>>& txpGroupCombinedWeights,
                 std::vector<uint64_t>& txpGroupCounts,
                 std::vector<bool>& fix,
                 const VecT& alphaIn,
                 VecT& fixedCounts) {
    for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
        if (fix[eqId]) { continue; }
    uint64_t count = txpGroupCounts[eqID];
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    double denom = 0.0;
    size_t groupSize = txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (groupSize > 1) {
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        if (alphaIn[tid] > 0.0) {
          double v = alphaIn[tid] * aux;
          denom += v;
        }
      }
      if (denom <= ::minEQClassWeight) {
        // tgroup.setValid(false);
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          if (alphaIn[tid] > 0.0) {
            double v = alphaIn[tid] * aux;
            fixedCounts[tid] += v * invDenom;
          }
        }
      }

    } else {
      fixedCounts[txps.front()] += count;
    }
  }
}
*/

/**
 * Single-threaded VBEM-update routine for use in bootstrapping
 */
template <typename VecT>
void VBEMUpdate_(
    const std::vector<std::vector<uint32_t>>& txpGroupLabels,
    const std::vector<std::vector<double>>& txpGroupCombinedWeights,
    const std::vector<size_t>& txpGroupCounts, const VecT& priorAlphas,
    const VecT& flatPriorAlphas, const VecT& alphaIn, VecT& alphaOut,
    VecT& expTheta, const VecT& factors) {

  size_t N = alphaIn.size();
  assert(alphaIn.size() == alphaOut.size());

  size_t numEQClasses = txpGroupLabels.size();

  double alphaSumInformative = {0.0};
  alphaSumInformative = (alphaIn + priorAlphas).sum();

  int ifault;
  double logNorm = digamma(alphaSumInformative, &ifault);
  for (size_t i = 0; i < N; ++i) {
    double apInfo = 0.0;
    apInfo = alphaIn[i] + priorAlphas[i];

    if (apInfo > ::digammaMin) {
      expTheta[i] = std::exp(digamma(apInfo, &ifault) - logNorm);
      if (ifault or std::isinf(expTheta[i]) or std::isnan(expTheta[i])) {
        std::cerr << "alpha[" << i << "] = " << alphaIn[i] << ", expTheat[" << i
                  << "] = " << expTheta[i] << "\n";
      }
    } else {
      expTheta[i] = 0.0;
    }
    alphaOut[i] = 0.0;
  }

  for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
    uint64_t count = txpGroupCounts[eqID];
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    double denom = 0.0;
    size_t groupSize = txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (groupSize > 1) {
      double minFactor{1.0};
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        if (expTheta[tid] > 0.0) {
          double v = expTheta[tid] * aux;
          denom += v;
        }
      }
      if (denom <= ::minEQClassWeight) {
        // tgroup.setValid(false);
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          if (expTheta[tid] > 0.0) {
            double v =  expTheta[tid] * aux;
            alphaOut[tid] += v * invDenom;
          }
        }
      }

    } else {
      alphaOut[txps.front()] += count;
    }
  }
}

/**
 * Single-threaded VBEM-update routine for use in bootstrapping
 */
template <typename VecT>
void EMUpdateAdaptive_(
                         const std::vector<std::vector<uint32_t>>& txpGroupLabels,
                         const std::vector<std::vector<double>>& txpGroupCombinedWeights,
                         const std::vector<size_t>& txpGroupCounts,
                         const VecT& priorAlphas,
                         const VecT& flatPriorAlphas,
                         const VecT& alphaIn,
                         const VecT& alphaFlat,
			 //const VecT& effLens,
                         VecT& alphaOut,
                         VecT& expTheta,
                         const VecT& factors,
                         const size_t itNum) {
    size_t N = alphaIn.size();
    assert(alphaIn.size() == alphaOut.size());
    VecT expThetaFlat(N); expThetaFlat.setZero();
    size_t numEQClasses = txpGroupLabels.size();

    double decay = std::exp(-static_cast<double>(itNum));
    //auto thetas = alphaIn / effLens;
    //thetas /= thetas.sum();

for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
    uint64_t count = txpGroupCounts[eqID];
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    double denom = 0.0;
    size_t groupSize = txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (groupSize > 1) {
      double minFactor{1.0};
      //double minFactor{0.0};
      double localWeightFlat{0.0};
      double localWeightInfo{0.0};
      double atot{0.0};
      for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          localWeightFlat += flatPriorAlphas[tid];
          localWeightInfo += priorAlphas[tid];
          atot += alphaIn[tid];
          if (factors[tid] < minFactor) {
              minFactor = factors[tid];
          }
          //minFactor += factors[tid];
      }
      //minFactor /= groupSize;

      //std::cerr << "local weight flat = " << localWeightFlat << ", local weight info = " << localWeightInfo << "\n";
      double lwFlat = decay;//0.05 * (atot / localWeightFlat);
      double localWeight = (decay*localWeightFlat) / localWeightInfo;

      for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto apFlat = alphaIn[tid] + lwFlat * flatPriorAlphas[tid];
          auto apInfo = alphaIn[tid] + localWeight * priorAlphas[tid];
          //std::cerr << "apFlat = " << apFlat << ", apInfo = " << apInfo << "\n";
          expThetaFlat[tid] = apFlat;//(invFlat > 0.0) ? (apFlat * invFlat) : 0.0;
          expTheta[tid] = apInfo;//(invInfo > 0.0) ? (apInfo * invInfo) : 0.0;
      }

      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        if ((expTheta[tid] + expThetaFlat[tid]) > 0.0) {
            double v = ((1.0 - minFactor) * expThetaFlat[tid] +
                      minFactor * expTheta[tid]) *
                     aux;
          denom += v;
        }
      }
      if (denom <= ::minEQClassWeight) {
        // tgroup.setValid(false);
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          if ((expTheta[tid] + expThetaFlat[tid]) > 0.0) {
              double v = ((1.0 - minFactor) * expThetaFlat[tid] +
                        minFactor * expTheta[tid]) *
                       aux;
            alphaOut[tid] += v * invDenom;
          }
        }
      }

    } else {
      alphaOut[txps.front()] += count;
    }
  }
}

/**
 * Single-threaded VBEM-update routine for use in bootstrapping
 */
template <typename VecT>
void VBEMUpdateAdaptive_(
    const std::vector<std::vector<uint32_t>>& txpGroupLabels,
    const std::vector<std::vector<double>>& txpGroupCombinedWeights,
    const std::vector<size_t>& txpGroupCounts,
    const VecT& priorAlphas,
    const VecT& flatPriorAlphas,
    const VecT& alphaIn,
    const VecT& alphaFlat,
    VecT& alphaOut,
    VecT& expTheta,
    const VecT& factors,
    const size_t itNum) {

    double decay = 1.0;//std::exp(-static_cast<double>(itNum + 1));

  size_t N = alphaIn.size();
  assert(alphaIn.size() == alphaOut.size());

  size_t numEQClasses = txpGroupLabels.size();

  double alphaSumInformative = {0.0};
  double alphaSumFlat = {0.0};

  alphaSumInformative = (alphaIn + decay * priorAlphas).sum();
  alphaSumFlat = (alphaIn + decay * flatPriorAlphas).sum();

  int ifault;
  double logNorm = digamma(alphaSumInformative, &ifault);
  double logNormFlat = digamma(alphaSumFlat, &ifault);
  VecT expThetaFlat(N); expThetaFlat.setZero();

  for (size_t i = 0; i < N; ++i) {
    //if (active.find(i) == active.end()) {
    //  expTheta[i] = 0.0;
    //  expThetaFlat[i] = 0.0;
    //} else {
    double apFlat = 0.0;
    double apInfo = 0.0;
    apInfo = alphaIn[i] + decay * priorAlphas[i];
    apFlat = alphaIn[i] + decay * flatPriorAlphas[i];

    if (apInfo > ::digammaMin) {
      expTheta[i] = std::exp(digamma(apInfo, &ifault) - logNorm);
      if (ifault or std::isinf(expTheta[i]) or std::isnan(expTheta[i])) {
        std::cerr << "alpha[" << i << "] = " << alphaIn[i] << ", expTheat[" << i
                  << "] = " << expTheta[i] << "\n";
      }
    } else {
      expTheta[i] = 0.0;
    }

    if (apFlat > ::digammaMin) {
      expThetaFlat[i] = std::exp(digamma(apFlat, &ifault) - logNormFlat);
      if (ifault or std::isinf(expThetaFlat[i]) or
          std::isnan(expThetaFlat[i])) {
        std::cerr << "alpha[" << i << "] = " << alphaIn[i] << ", expTheat[" << i
                  << "] = " << expThetaFlat[i] << "\n";
      }
    } else {
      expThetaFlat[i] = 0.0;
    }
    //}

    alphaOut[i] = 0.0;
  }

  for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
    uint64_t count = txpGroupCounts[eqID];
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    double denom = 0.0;
    size_t groupSize = txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (groupSize > 1) {
      double minFactor{1.0};
      for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];

          if (factors[tid] < minFactor) {
              minFactor = factors[tid];
	  }

      }
      double wflat = 1.0 - minFactor;
      double winfo = minFactor;

      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        if (expTheta[tid] > 0.0) {
	  double v = (wflat * expThetaFlat[tid] +
		     (winfo * expTheta[tid])) * aux;
          denom += v;
        }
      }
      if (denom > ::minEQClassWeight) {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          if (expTheta[tid] > 0.0) {
              double v = (wflat * expThetaFlat[tid] +
			 (winfo * expTheta[tid])) * aux;
            alphaOut[tid] += v * invDenom;
          }
        }
      }

    } else {
      alphaOut[txps.front()] += count;
    }
  }
}

Optimizer::Optimizer() {}

Eigen::VectorXd
Optimizer::optimize(EquivCollection& eqc, Eigen::VectorXd& alphas,
                    Eigen::VectorXd& lengths,
                    Eigen::VectorXd& effLens,
                    Eigen::VectorXd& priorAlphas,
                    Eigen::VectorXd& flatPriorAlphas,
                    Eigen::VectorXd& factors, Eigen::VectorXd& alphasFlat, OptimizationType ot,
                    double relDiffTolerance, uint32_t maxIter) {

  uint32_t minIter = 50;
  if (minIter > maxIter) { minIter = maxIter; }
  size_t M = alphas.size();
  using VecT = Eigen::VectorXd;

  //VecType alphasFlat(M); alphasFlat.setOnes();
  //alphasFlat *= 1.0 / M;
  //VecType alphasPrimeFlat(M); alphasPrimeFlat.setZero();

  VecType alphasPrime(M); alphasPrime.setZero();
  VecType expTheta(M); expTheta.setZero();

  bool useVBEM{true};
  bool perTranscriptPrior{false};

  // auto jointLog = sopt.jointLog;
  // auto& fragStartDists = readExp.fragmentStartPositionDistributions();
  double totalNumFrags = alphas.sum();
  double totalLen{0.0};

  bool useEffectiveLengths = true;
  bool converged{false};
  double maxRelDiff = -std::numeric_limits<double>::max();
  size_t itNum = 0;

  double minAlpha = 1e-8;
  double alphaCheckCutoff = 1e-2;
  double cutoff = minAlpha;
  while (itNum < minIter or (itNum < maxIter and !converged)) {

      switch (ot) {
          case OptimizationType::EM:
              EMUpdate_(eqc.labels_, eqc.auxProbs_, eqc.counts_, alphas, alphasPrime);
              break;
          case OptimizationType::VBEM:
              VBEMUpdate_(eqc.labels_, eqc.auxProbs_, eqc.counts_, priorAlphas,
                          flatPriorAlphas, alphas, alphasPrime, expTheta, factors);
              break;
          case OptimizationType::VBEM_ADAPTIVE:
              VBEMUpdateAdaptive_(eqc.labels_, eqc.auxProbs_, eqc.counts_, priorAlphas,
                                flatPriorAlphas, alphas, alphasFlat, alphasPrime, expTheta, factors, itNum);
              break;
      }


    converged = true;
    maxRelDiff = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < M; ++i) {
      if (alphasPrime[i] > alphaCheckCutoff) {
        double relDiff = std::abs(alphas[i] - alphasPrime[i]) / alphasPrime[i];
        maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
        if (relDiff > relDiffTolerance) {
          converged = false;
        }
      }
      alphas[i] = alphasPrime[i];
      alphasPrime[i] = 0.0;
    }

    ++itNum;
    if (itNum % 10 == 0) {
      std::cerr << "\r\riteration = " << itNum
                << " | max rel diff. = " << maxRelDiff;
    }
  }
  std::cerr << "\n\n";
  // Truncate tiny expression values
  double alphaSum = 0.0;
  if (useVBEM and !perTranscriptPrior) {
    std::vector<double> cutoffs(M, 0.0);
    for (size_t i = 0; i < M; ++i) {
      // cutoffs[i] = priorAlphas[i] + minAlpha;
      cutoffs[i] = minAlpha;
    }
    alphaSum = truncateCountVector(alphas, cutoffs);
  } else {
    // Truncate tiny expression values
    alphaSum = truncateCountVector(alphas, cutoff);
  }

  if (alphaSum < minWeight) {
    std::cerr << "Total alpha weight was too small! "
              << "Make sure you ran salmon correclty.\n";
    return alphas;
  }

  bool useScaledCounts{false};
  if (useScaledCounts) {
    /*
    double mappedFragsDouble = static_cast<double>(numMappedFrags);
    double alphaSum = 0.0;
    for (auto a : alphas) {
      alphaSum += a;
    }
    if (alphaSum > ::minWeight) {
      double scaleFrac = 1.0 / alphaSum;
      // scaleFrac converts alpha to nucleotide fraction,
      // and multiplying by numMappedFrags scales by the total
      // number of mapped fragments to provide an estimated count.
      for (auto& a : alphas) {
        a = mappedFragsDouble * (a * scaleFrac);
      }
    */
  } else { // This shouldn't happen!
           /*
               sopt.jointLog->error(
                   "Bootstrap had insufficient number of fragments!"
                   "Something is probably wrong; please check that you "
                   "have run salmon correctly and report this to GitHub.");
           */
  }
  return alphas;
}


template
void EMUpdate_<Eigen::VectorXd>(std::vector<std::vector<uint32_t>>& txpGroupLabels,
                                std::vector<std::vector<double>>& txpGroupCombinedWeights,
                                std::vector<size_t>& txpGroupCounts, const Eigen::VectorXd& alphaIn,
                                Eigen::VectorXd& alphaOut);
template
double truncateCountVector<Eigen::VectorXd>(Eigen::VectorXd& alphas, double cutoff);

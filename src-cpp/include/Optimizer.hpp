#ifndef _OPTIMIZER_HPP
#define _OPTIMIZER_HPP

#include "eigen3/Eigen/Dense"
#include "EquivCollection.hpp"
#include <vector> 

enum class OptimizationType : uint8_t {
    EM = 0,
    VBEM = 1,
    VBEM_ADAPTIVE = 2,
    EM_ADAPTIVE = 3
};

class Optimizer {
public:
  using VecType = Eigen::VectorXd;
  using SerialVecType = std::vector<double>;
  Optimizer();

  Eigen::VectorXd optimize(
			   EquivCollection& eqc,
			   Eigen::VectorXd& alphas,
			   Eigen::VectorXd& lengths,
			   Eigen::VectorXd& effLens,
			   Eigen::VectorXd& priorAlphas,
			   Eigen::VectorXd& flatPriorAlphas,
               Eigen::VectorXd& factors,
               Eigen::VectorXd& estCounts,
               OptimizationType ot,
			   double tolerance = 0.01,
			   uint32_t maxIter = 1000);
};

#endif // _OPTIMIZER_HPP

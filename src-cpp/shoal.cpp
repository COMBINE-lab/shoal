#include <fstream>
#include <iostream>
#include <stack>
#include <unordered_set>

#include "EquivCollection.hpp"
#include "Optimizer.hpp"
#include "csv.h"
#include "eigen3/Eigen/Dense"
#include "filesystem/path.h"
#include "popl.hpp"
#include "sparsepp.h"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"

struct PriorEntry {
  double factor;
  double mean;
  double stddev;
};
struct QuantEntry {
  uint32_t len;
  double efflen;
  double tpm;
  double numReads;
};

void writeSFFile(std::string fname, std::vector<std::string>& names,
                 Eigen::VectorXd& counts, Eigen::VectorXd& lengths,
                 Eigen::VectorXd& effLens) {
  double million = 1000000.0;
  size_t M = counts.size();
  Eigen::VectorXd tpm(M); tpm.setZero();
  double denom{0.0};
  for (size_t i = 0; i < M; ++i) {
    denom += counts[i] / effLens[i];
  }
  for (size_t i = 0; i < M; ++i) {
    tpm[i] += million * (counts[i] / effLens[i]) / denom;
  }

  std::ofstream ofile(fname);
  ofile << "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";
  for (size_t i = 0; i < M; ++i) {
    ofile << names[i] << '\t' << lengths[i] << '\t' << effLens[i] << '\t'
          << tpm[i] << '\t' << counts[i] << '\n';
  }
  ofile.close();
}

template <typename VecT>
void EMUpdate_(std::vector<std::vector<uint32_t>>& txpGroupLabels,
               std::vector<std::vector<double>>& txpGroupCombinedWeights,
               std::vector<size_t>& txpGroupCounts, const VecT& alphaIn,
               VecT& alphaOut);

template <typename VecT>
double truncateCountVector(VecT& alphas, double cutoff);

Eigen::VectorXd optAdaptEst(EquivCollection& ec, spp::sparse_hash_map<std::string, QuantEntry>& quantMap,
                              spp::sparse_hash_map<std::string, PriorEntry>& priorMap, double weight) {
    auto console = spdlog::get("console");
    size_t N = quantMap.size();
    Eigen::VectorXd prior(N);
    Eigen::VectorXd priorInform(N);
    Eigen::VectorXd alphas(N);
    Eigen::VectorXd effLens(N);
    Eigen::VectorXd lengths(N);
    Eigen::VectorXd tpms(N);
    Eigen::VectorXd estCounts(N);
    Eigen::VectorXd factors(N);
    Eigen::VectorXd factorsInform(N);

    size_t i{0};
    for (auto& tname : ec.tnames_) {
        auto& quantEnt = quantMap[tname];
        lengths[i] = static_cast<double>(quantEnt.len);
        effLens[i] = quantEnt.efflen;
        tpms[i] = quantEnt.tpm;
        estCounts[i] = quantEnt.numReads;
        alphas[i] = 1.0 / N;
        auto& priorEnt = priorMap[tname];
        prior[i] = 1e-3 * effLens[i];
        priorInform[i] = (weight * priorEnt.mean);
        factors[i] = 1.0;
        factorsInform[i] = priorEnt.factor;
        ++i;
    }

    console->info("num txps = {}", priorMap.size());
    Optimizer opt;
    Eigen::VectorXd alphaOptInform;

	alphaOptInform = opt.optimize(ec, alphas, lengths, effLens, priorInform, prior, factorsInform, estCounts, OptimizationType::VBEM);
    auto alphaOpt = opt.optimize(ec, alphas, lengths, effLens, prior, prior, factors, estCounts, OptimizationType::VBEM);



    Eigen::VectorXd merged(N);
    for (size_t i = 0; i < N; ++i) {
        auto& priorEnt = priorMap[ec.tnames_[i]];
        merged[i] = (1.0 - priorEnt.factor) * alphaOpt[i] +  (priorEnt.factor * alphaOptInform[i]);
    }
    Eigen::VectorXd mergedOut(N); mergedOut.setZero();
    EMUpdate_(ec.labels_, ec.auxProbs_, ec.counts_, merged, mergedOut);
    truncateCountVector(merged, 1e-8);
    return mergedOut;
}



Eigen::VectorXd optAdaptPrior(EquivCollection& ec, spp::sparse_hash_map<std::string, QuantEntry>& quantMap,
                              spp::sparse_hash_map<std::string, PriorEntry>& priorMap, double weight) {

    auto console = spdlog::get("console");
    size_t N = quantMap.size();
    Eigen::VectorXd prior(N); prior.setZero();
    Eigen::VectorXd flatPrior(N); flatPrior.setZero();
    Eigen::VectorXd alphas(N);
    Eigen::VectorXd effLens(N);
    Eigen::VectorXd lengths(N);
    Eigen::VectorXd estCounts(N); estCounts.setZero();
    Eigen::VectorXd factors(N);

    auto numEQClasses = ec.labels_.size();
    for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
        uint64_t count = ec.counts_[eqID];
        const std::vector<uint32_t>& txps = ec.labels_[eqID];
        size_t groupSize = txps.size();
        // If this is a single-transcript group,
        // then it gets the full count.  Otherwise,
        // update according to our VBEM rule.
        if (groupSize == 1) {
            estCounts[txps.front()] = count;
        }
    }

    size_t i{0};
    double flatSum = 0.0;
    double infoSum = 0.0;
    for (auto& tname : ec.tnames_) {
        auto& quantEnt = quantMap[tname];
        lengths[i] = static_cast<double>(quantEnt.len);
        effLens[i] = quantEnt.efflen;
        estCounts[i] += 1e-3 * effLens[i];//quantEnt.numReads;
        alphas[i] = 1.0/N;//estCounts[i];
        auto& priorEnt = priorMap[tname];
        factors[i] = priorEnt.factor;
        ++i;
    }
    /*
    alphas = Eigen::VectorXd::Random(N);
    double mc = alphas.minCoeff();
    alphas = alphas.array() + MC;
    */

    std::unordered_set<uint32_t> active;
    for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
      const std::vector<uint32_t>& txps = ec.labels_[eqID];
      size_t groupSize = txps.size();
      for (size_t i = 0; i < groupSize; ++i) {
	active.insert(txps[i]);
      }
    }

    //std::ofstream fsr("ratios.txt");
    //std::ofstream fsf("flat.txt");
    //std::ofstream fsi("info.txt");
    double mf{0.0};
    double mi{0.0};
    std::vector<double> ratios;
    for (auto t : active) {
      auto& priorEnt = priorMap[ec.tnames_[t]];
      auto f = 1e-3 * effLens[t];
      auto info = priorEnt.mean;
      flatSum += f;
      infoSum += info;
      if (f >= mf) { mf = f; }
      if (info >= mi) { mi = info; }
      double r = (info > 0.0) ? (f / info) : ((f > 0.0) ? 1.0 : 0.0);
      //fsr << r << '\t';
      //fsf << 1e-3 * effLens[t] << '\t';
      //fsi << info << '\t';
      ratios.push_back(r);
    }
    //fsr.close();
    //fsf.close();
    //fsi.close();

    console->info("flatSum = {}, infoSum = {}, ratio = {}", flatSum, infoSum, flatSum / infoSum);
    console->info("num txps = {}", priorMap.size());
    auto autoWeight = 0.15 * (flatSum / infoSum);
    //autoWeight = (flatSum / infoSum);
    //autoWeight = 1.0;//mf / mi;


    //double dfact = 0.1;
    //std::nth_element(ratios.begin(), ratios.begin() + (dfact * ratios.size()), ratios.end());
    //autoWeight = ratios[(dfact * ratios.size())];

    for (auto t : active ) {
      auto& tname = ec.tnames_[t];
      auto& priorEnt = priorMap[tname];
      auto& quantEnt = quantMap[tname];
      prior[t] = autoWeight * (priorEnt.mean);
      flatPrior[t] = 1e-3 * effLens[t];
    }

    Optimizer opt;
    auto alphaOpt = opt.optimize(ec, alphas, lengths, effLens, prior, flatPrior, factors, estCounts, OptimizationType::VBEM_ADAPTIVE);
    /*
    i = 0;
    for (auto& tname : ec.tnames_) {
        auto& quantEnt = quantMap[tname];
        if (alphaOpt[i] <= 1e-8 and quantEnt.numReads >= 1e-8) {
            alphaOpt[i] = 1e-3 * quantEnt.efflen;
        }
        ++i;
    }
    */
    return alphaOpt;
}



Eigen::VectorXd optNoPrior(EquivCollection& ec, spp::sparse_hash_map<std::string, QuantEntry>& quantMap, double weight, std::string optBaseStr) {

    auto console = spdlog::get("console");
    size_t N = quantMap.size();
    Eigen::VectorXd prior(N);
    Eigen::VectorXd alphas(N);
    Eigen::VectorXd effLens(N);
    Eigen::VectorXd lengths(N);
    Eigen::VectorXd tpms(N);
    Eigen::VectorXd estCounts(N);
    Eigen::VectorXd factors(N);

    size_t i{0};
    for (auto& tname : ec.tnames_) {
        auto& quantEnt = quantMap[tname];
        lengths[i] = static_cast<double>(quantEnt.len);
        effLens[i] = quantEnt.efflen;
        tpms[i] = quantEnt.tpm;
        estCounts[i] = quantEnt.numReads;
        alphas[i] = 1.0 / N;
        prior[i] = weight * effLens[i];
        factors[i] = 1.0;
        ++i;
    }

    console->info("num txps = {}", quantMap.size());
    Optimizer opt;
    if (optBaseStr == "em"){
	    return opt.optimize(ec, alphas, lengths, effLens, prior, prior, factors, estCounts, OptimizationType::EM);
    }
    else{
	    return opt.optimize(ec, alphas, lengths, effLens, prior, prior, factors, estCounts, OptimizationType::VBEM);
    }

}


namespace spd = spdlog;
int main(int argc, char* argv[]) {
  using namespace popl;

  std::string priorFile;
  std::string sampleDir;
  std::string outFile;
  std::string optStr;
  std::string optBaseStr;
  double weight{0.001};

  Switch helpOption("h", "help", "produce help message");
  Value<std::string> optBaseType("b", "optBaseType", "The type of base level optimization to perform one of {vbem, or em}",
                             "vbem", &optBaseStr);
  Value<std::string> optType("t", "optType", "The type of optimization to perform one of {adapt-prior, or adapt-est}",
                             "adapt-prior", &optStr);
  Value<std::string> priorOpt("p", "prior", "file containing prior", "",
                              &priorFile);
  Value<std::string> sampleOpt("s", "sample", "file containing quantifications",
                               "", &sampleDir);
  Value<std::string> outOpt("o", "output", "output file", "", &outFile);
  Value<double> weightOpt("w", "weight", "informative prior weight", 0.001,
                          &weight);

  OptionParser op("Options");
  op.add(helpOption)
      .add(optType)
      .add(optBaseType)
      .add(priorOpt)
      .add(sampleOpt)
      .add(outOpt)
      .add(weightOpt);

  try {
    // Console logger with color
    auto console = spd::stdout_color_mt("console");
    op.parse(argc, argv);

    if (helpOption.isSet()) {
      std::cout << op << "\n";
      std::exit(0);
    }

    std::unordered_set<std::string> validOptTypes = {"adapt-prior", "adapt-est"};
    if (validOptTypes.find(optStr) == validOptTypes.end()) {
        console->critical("Do not recognize optType {}!", optStr);
        std::exit(1);
    }

    std::unordered_set<std::string> validOptBaseTypes = {"em", "vbem"};
    if (validOptBaseTypes.find(optBaseStr) == validOptBaseTypes.end()) {
        console->critical("Do not recognize optBaseType {}!", optBaseStr);
        std::exit(1);
    }

    std::string sampleFile = sampleDir + "/quant.sf";
    auto sampleFilePath = filesystem::path(sampleFile);
    if (!sampleOpt.isSet() || !sampleFilePath.exists()) {
      console->critical("Could not find the \"quant.sf\" file. Give the path to "
                       "parent directory containing the \"quant.sf\" file");
      std::exit(0);
    }

    auto outpath = filesystem::path(outFile).parent_path();
    std::stack<filesystem::path> createStack;
    while (!outpath.exists() and !outpath.empty()) {
      createStack.push(outpath);
      outpath = outpath.parent_path();
    }
    while (!createStack.empty()) {
      filesystem::create_directory(createStack.top());
      console->info("creating directory: {}", createStack.top());
      createStack.pop();
    }

    std::string eqFile = sampleDir + "/aux_info/eq_classes.txt";
    if (!filesystem::path(eqFile).exists()){
        console->critical("Could not find the \"eq_classes.txt\" file."
                        "Make sure you ran right version of salmon");
        std::exit(0);
    }
    EquivCollection ec;
    ec.fromFile(eqFile);
    console->info("num txp names = {}", ec.tnames_.size());

    spp::sparse_hash_map<std::string, PriorEntry> priorMap;
    spp::sparse_hash_map<std::string, QuantEntry> quantMap;
    std::string name;

    bool havePrior{priorOpt.isSet()};
    if (havePrior) {
      io::CSVReader<4, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(
          priorFile);
      in.read_header(io::ignore_missing_column, "Name", "factor", "mean",
                     "stddev");
      double factor, mean, stddev;
      while (in.read_row(name, factor, mean, stddev)) {
        priorMap[name] = {factor, mean, stddev};
      }
    }

    io::CSVReader<5, io::trim_chars<' '>, io::no_quote_escape<'\t'>> quant(
        sampleFile);
    quant.read_header(io::ignore_missing_column, "Name", "Length",
                      "EffectiveLength", "TPM", "NumReads");
    uint32_t len;
    double efflen, tpm, numReads;
    while (quant.read_row(name, len, efflen, tpm, numReads)) {
      quantMap[name] = {len, efflen, tpm, numReads};
    }

    size_t N = quantMap.size();
    Eigen::VectorXd effLens(N);
    Eigen::VectorXd lengths(N);
    size_t i{0};
    for (auto& tname : ec.tnames_) {
        auto& quantEnt = quantMap[tname];
        lengths[i] = static_cast<double>(quantEnt.len);
        effLens[i] = quantEnt.efflen;
        ++i;
    }

    Eigen::VectorXd merged;
    if (havePrior) {
        if (optStr == "adapt-prior") {
            merged = optAdaptPrior(ec, quantMap, priorMap, weight);
        } else if (optStr == "adapt-est") {
            merged = optAdaptEst(ec, quantMap, priorMap, weight);
        }
    } else {
        merged = optNoPrior(ec, quantMap, weight, optBaseStr);
    }

    writeSFFile(outFile, ec.tnames_, merged, lengths, effLens);
    // Release and close all loggers
    spdlog::drop_all();
  } catch (const spd::spdlog_ex& ex) {
    std::cout << "Log init failed: " << ex.what() << std::endl;
    return 1;
  }
}

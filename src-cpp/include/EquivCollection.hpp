#ifndef _EQUIV_COLLECTION_H_
#define _EQUIV_COLLECTION_H_

#include <vector>
#include <fstream>

class EquivCollection {
public:
  EquivCollection() {}
  
  void fromFile(const std::string& fn) {
    std::ifstream ifile(fn);
    size_t numT;
    size_t numEC;

    ifile >> numT >> numEC;
    for (size_t i = 0; i < numT; ++i) {
      std::string tname;
      ifile >> tname;
      tnames_.push_back(tname);
    }

    for (size_t i = 0; i < numEC; ++i) {
      size_t nlabels;
      ifile >> nlabels;

      std::vector<uint32_t> labels;
      std::vector<double> auxs;
      for (size_t j = 0; j < nlabels; ++j) {
	uint32_t label;
	ifile >> label;
	labels.push_back(label);
      }
      labels_.push_back(labels);
      
      for (size_t j = 0; j < nlabels; ++j) {
	double p;
	ifile >> p;
	auxs.push_back(p);
      }
      auxProbs_.push_back(auxs);

      size_t count;
      ifile >> count;
      counts_.push_back(count);
    }
    ifile.close();
  }

  std::vector<std::string> tnames_; 
  std::vector<std::vector<uint32_t>> labels_; 
  std::vector<std::vector<double>> auxProbs_; 
  std::vector<size_t> counts_; 
};

#endif // _EQUIV_COLLECTION_H_

// =============================================================================
//
//   Copyright (C) 2020-2024 Vasileios Vavourakis (vasvav@gmail.com)
//   All Rights Reserved.
//
//   Licensed under the GNU General Public License v3.0 (the "License").
//   See the LICENSE file provided in this project details the License.
//   You cannot use this file except in compliance with the License.
//
// =============================================================================

// =============================================================================
#ifndef _IDA_H_
#define _IDA_H_
// =============================================================================
#include "./global.h"
// =============================================================================
class InverseDistanceAlgorithm {
//
public:
  explicit InverseDistanceAlgorithm(double c =1.0) : is_init_(false), coeff_(c), pnt_(0) {}
  InverseDistanceAlgorithm(double c, const std::vector<bdm::Double3>& p) : coeff_(c) {
    this->init(p);
  }
  ~InverseDistanceAlgorithm() {}
  //
  inline
  void init(const std::vector<bdm::Double3>& p) {
    this->is_init_ = true;
    this->pnt_ = &p;
  }
  //
  inline
  double interpolate(const bdm::Double3& a, const std::vector<double>& input) const {
    const int npnt = this->pnt_->size();
    // sanity check
    ASSERT_(npnt==input.size(),"size of input data is wrong");
    // calculate the weights
    std::vector<double> phi;
    this->calculate_(a, phi);
    // ...and now interpolate the input data
    double out;
    out = 0.0;
    for (int i=0; i<npnt; i++)
      out += phi[i] * input[i];
    return out;
  }
//
private:
  inline
  void calculate_(const bdm::Double3& a, std::vector<double>& phi) const {
      const int npnt = this->pnt_->size();
      std::vector<double> weights(npnt);
      for (int i=0; i<npnt; i++)
        {
          const bdm::Double3 delta = (*this->pnt_)[i] - a;
          double d = sqrt(pow(delta[0],2)+pow(delta[1],2)+pow(delta[2],2));
          weights[i] = 1.0 / pow(d, this->coeff_);
        }
      const double weights_tot = accumulate(weights.begin(), weights.end(), 0.0);
      // calculate the normalized weights
      phi.resize(npnt);
      for (int i=0; i<npnt; i++)
        phi[i] = weights[i] / weights_tot;
  }
//
private:
    bool is_init_;
    // (reciprocal) exponent of the inverse distance weight function
    double coeff_;
    // input point (cloud) used to interpolate the data
    const std::vector<bdm::Double3>* pnt_;
};
// =============================================================================
#endif // _IDA_H_
// =============================================================================

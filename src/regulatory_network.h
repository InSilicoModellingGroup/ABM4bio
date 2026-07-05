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
#ifndef _REGULATORY_NETWORK_H_
#define _REGULATORY_NETWORK_H_
// =============================================================================
#include "./global.h"
#include "boost/numeric/odeint.hpp"
#include "boost/phoenix/core.hpp"
#include "boost/phoenix/operator.hpp"
// =============================================================================
typedef boost::numeric::ublas::vector<double>  bvector_t;
typedef boost::numeric::ublas::matrix<double>  bmatrix_t;
// =============================================================================
struct RegulatoryNetworkData {
  // pseudo-time for ODE(s) time integration
  mutable
  double time = 0.0;
  // time-step for ODE(s) time integration
  double time_step = 1.0;
  int time_step_subdivision = 100;
  // current solution of the species concentration
  mutable
  bvector_t current_species = {};
  // previous solution of the species concentration
  mutable
  bvector_t previous_species = {};
  // parameters of the regulatory network
  std::vector<double> params, params_a, params_i;
};
// =============================================================================
inline
int read_regulatory_network_data(const std::string& fname, RegulatoryNetworkData& rn)
{
  std::ifstream fin(fname);
  ASSERT_(fin.good(),"file \""+fname+"\" cannot be accessed");
  //
  int number_of_species = 0;
  //
  fin >> rn.time_step >> rn.time_step_subdivision >> number_of_species;
  //
  rn.time = 0.0;
  // now read the regulatory network parameters respectively
  //
  rn.params.resize(number_of_species+1);
  for (int p=0; p<number_of_species+1; p++)
    fin >> rn.params[p];
  //
  rn.params_a.resize(number_of_species);
  for (int p=0; p<number_of_species; p++)
    fin >> rn.params_a[p];
  //
  rn.params_i.resize(number_of_species);
  for (int p=0; p<number_of_species; p++)
    fin >> rn.params_i[p];
  //
  // also read the initial conditions of the species for the regulatory network
  rn.current_species.resize(number_of_species);
  rn.previous_species.resize(number_of_species);
  for (int l=0; l<number_of_species; l++)
    {
      fin >> rn.current_species[l];
      rn.previous_species[l] = rn.current_species[l];
    }
  //
  return number_of_species;
}
// =============================================================================
#endif // _REGULATORY_NETWORK_H_
// =============================================================================

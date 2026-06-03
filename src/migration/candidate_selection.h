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

#ifndef _MIGRATION_CANDIDATE_SELECTION_H_
#define _MIGRATION_CANDIDATE_SELECTION_H_
// =============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>
#include "core/container/math_array.h"
#include "global.h"
// =============================================================================
namespace bdm {
class Random;
namespace migration {
// -----------------------------------------------------------------------------
/// Unit directions for local candidate moves (26 in 3D, 8 in 2D), shuffled.
inline
std::vector<Double3> BuildShuffledLocalDirections(const bool simulation_domain_is_2d,
                                                  Random* rg)
{
  std::vector<Double3> local_directions;
  for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
      for (int dz = -1; dz <= 1; dz++)
        {
          if (simulation_domain_is_2d && dz) continue;
          if (!dx && !dy && !dz) continue;
          Double3 direction = {
            static_cast<double>(dx),
            static_cast<double>(dy),
            static_cast<double>(dz)
          };
          if (!normalize(direction, direction)) continue;
          local_directions.push_back(direction);
        }
  for (size_t i = local_directions.size() - 1; i > 0; --i)
    {
      const size_t j = static_cast<size_t>(rg->Uniform(0, i + 1));
      std::swap(local_directions[i], local_directions[j]);
    }
  return local_directions;
}
// -----------------------------------------------------------------------------
/// Chemotaxis score: sign(χ)[C(x_j)-C(x)] plus optional crowding adjustment.
inline
double ComputeChemotaxisScore(const double chemotaxis_sign,
                              const double concentration_at_origin,
                              const double concentration_at_candidate,
                              const double crowding_adjustment)
{
  return chemotaxis_sign * (concentration_at_candidate - concentration_at_origin)
         + crowding_adjustment;
}
// -----------------------------------------------------------------------------
/// Greedy local-best selection over feasible candidates.
inline
bool SelectLocalBestCandidate(
    const std::vector<Double3>& local_directions,
    const Double3& current_position,
    const double step_length,
    const bool simulation_domain_is_2d,
    const double chemotaxis_sign,
    const double concentration,
    const double tie_epsilon,
    const double migration_tolerance,
    const std::function<bool(const Double3&, double&)>& evaluate_candidate,
    const std::function<double(const Double3&)>& sample_concentration,
    Double3& chosen_dvec)
{
  double best_score = tie_epsilon;
  chosen_dvec = {0.0, 0.0, 0.0};
  bool found_better_candidate = false;
  for (const auto& direction : local_directions)
    {
      Double3 point = current_position + direction * step_length;
      if (simulation_domain_is_2d) point[2] = current_position[2];
      double crowding_adj = 0.0;
      if (!evaluate_candidate(point, crowding_adj)) continue;
      const double C_cand = sample_concentration(point);
      const double score =
        ComputeChemotaxisScore(chemotaxis_sign, concentration, C_cand, crowding_adj);
      if (score > best_score)
        {
          best_score = score;
          chosen_dvec = point - current_position;
          found_better_candidate = true;
        }
    }
  return found_better_candidate && L2norm(chosen_dvec) > migration_tolerance;
}
// -----------------------------------------------------------------------------
/// Boltzmann/softmax stochastic selection over feasible candidates.
inline
bool SelectStochasticLocalCandidate(
    const std::vector<Double3>& local_directions,
    const Double3& current_position,
    const double step_length,
    const bool simulation_domain_is_2d,
    const double chemotaxis_sign,
    const double concentration,
    const double beta,
    const double migration_tolerance,
    Random* rg,
    const std::function<bool(const Double3&, double&)>& evaluate_candidate,
    const std::function<double(const Double3&)>& sample_concentration,
    Double3& chosen_dvec)
{
  std::vector<double> scores_vec;
  std::vector<Double3> dvecs_vec;
  for (const auto& direction : local_directions)
    {
      Double3 point = current_position + direction * step_length;
      if (simulation_domain_is_2d) point[2] = current_position[2];
      double crowding_adj = 0.0;
      if (!evaluate_candidate(point, crowding_adj)) continue;
      const double C_cand = sample_concentration(point);
      scores_vec.push_back(
        ComputeChemotaxisScore(chemotaxis_sign, concentration, C_cand, crowding_adj));
      dvecs_vec.push_back(point - current_position);
    }
  if (dvecs_vec.empty()) return false;
  double max_s = scores_vec[0];
  for (double s : scores_vec)
    if (s > max_s) max_s = s;
  double sum_w = 0.0;
  std::vector<double> weights;
  weights.reserve(scores_vec.size());
  for (double s : scores_vec)
    {
      const double w = std::exp(beta * (s - max_s));
      weights.push_back(w);
      sum_w += w;
    }
  double r = rg->Uniform(0.0, sum_w);
  size_t chosen = dvecs_vec.size() - 1;
  double cumulative = 0.0;
  for (size_t k = 0; k < weights.size(); ++k)
    {
      cumulative += weights[k];
      if (r < cumulative) { chosen = k; break; }
    }
  chosen_dvec = dvecs_vec[chosen];
  return L2norm(chosen_dvec) > migration_tolerance;
}
// -----------------------------------------------------------------------------
/// Whether a substance-specific chemotaxis response is active at `concentration`.
inline
bool IsChemotaxisActive(const double concentration, const double threshold)
{
  return (threshold == 0.0) ||
         (threshold > 0.0 && concentration >= threshold) ||
         (threshold < 0.0 && concentration <= -threshold);
}
// =============================================================================
}  // namespace migration
}  // namespace bdm
// =============================================================================
#endif  // _MIGRATION_CANDIDATE_SELECTION_H_

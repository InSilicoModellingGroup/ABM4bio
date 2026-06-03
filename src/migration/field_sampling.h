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

#ifndef _MIGRATION_FIELD_SAMPLING_H_
#define _MIGRATION_FIELD_SAMPLING_H_
// =============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include "core/diffusion/diffusion_grid.h"
#include "neuroscience/param.h"
// =============================================================================
namespace bdm {
// -----------------------------------------------------------------------------
/// How a scalar field on a structured grid is sampled at a continuous position.
/// Independent of simulation parameters (CAP, cell phenotypes, etc.).
enum class SamplingMode {
  kNearestCell,           ///< Voxel containing the position (BioDynaMo GetValue).
  kLinearInterpolation,   ///< Trilinear weights from grid-corner-aligned voxels.
  kTrilinearInterpolation,///< Trilinear weights from cell-centre-aligned voxels.
  kSafeBoundaryClamp      ///< Cell-centre trilinear with position clamped to the grid.
};
// -----------------------------------------------------------------------------
/// Sample a diffusion-grid scalar field at `position` using `mode`.
/// When `sample_in_2d` is true, interpolation modes use bilinear weights in x–y
/// and hold the z voxel fixed (no through-plane interpolation).
inline
double SampleFieldAtPosition(const DiffusionGrid* grid,
                             const Real3& position,
                             const SamplingMode mode,
                             const bool sample_in_2d = false)
{
  if (!grid) return 0.0;
  const size_t res = grid->GetResolution();
  if (res <= 1 || mode == SamplingMode::kNearestCell)
    return static_cast<double>(grid->GetValue(position));
  const auto dims = grid->GetDimensions();
  const double grid_min = static_cast<double>(dims[0]),
               grid_max = static_cast<double>(dims[1]),
               box_len = grid->GetBoxLength();
  if (box_len <= 0.0)
    return static_cast<double>(grid->GetValue(position));
  const bool use_cell_centers =
    (mode == SamplingMode::kTrilinearInterpolation)
    || (mode == SamplingMode::kSafeBoundaryClamp);
  const double origin = use_cell_centers ? (grid_min + 0.5 * box_len) : grid_min,
               extent_max = use_cell_centers ? (grid_max - 0.5 * box_len) : grid_max;
  auto sample_axis = [&](const double coord, uint32_t& i0, uint32_t& i1, double& t) {
    double mapped = coord;
    if (mode == SamplingMode::kSafeBoundaryClamp)
      mapped = std::clamp(mapped, origin, extent_max);
    const double g = (mapped - origin) / box_len;
    int lower = static_cast<int>(std::floor(g));
    lower = std::clamp(lower, 0, static_cast<int>(res) - 2);
    i0 = static_cast<uint32_t>(lower);
    i1 = static_cast<uint32_t>(lower + 1);
    t = std::clamp(g - static_cast<double>(lower), 0.0, 1.0);
  };
  auto sample_axis_nearest = [&](const double coord, uint32_t& iz) {
    double mapped = coord;
    if (mode == SamplingMode::kSafeBoundaryClamp)
      mapped = std::clamp(mapped, origin, extent_max);
    const double g = (mapped - origin) / box_len;
    int index = static_cast<int>(std::floor(g + 0.5));
    index = std::clamp(index, 0, static_cast<int>(res) - 1);
    iz = static_cast<uint32_t>(index);
  };
  auto concentration = [&](const uint32_t x, const uint32_t y, const uint32_t z) {
    const std::array<uint32_t, 3> box = {x, y, z};
    return static_cast<double>(grid->GetConcentration(grid->GetBoxIndex(box)));
  };
  uint32_t x0 = 0, x1 = 0, y0 = 0, y1 = 0, z0 = 0, z1 = 0;
  double tx = 0.0, ty = 0.0, tz = 0.0;
  sample_axis(static_cast<double>(position[0]), x0, x1, tx);
  sample_axis(static_cast<double>(position[1]), y0, y1, ty);
  if (sample_in_2d)
    {
      uint32_t z_fixed = 0;
      sample_axis_nearest(static_cast<double>(position[2]), z_fixed);
      const double c00 = concentration(x0, y0, z_fixed),
                   c10 = concentration(x1, y0, z_fixed),
                   c01 = concentration(x0, y1, z_fixed),
                   c11 = concentration(x1, y1, z_fixed);
      const double c0 = c00 * (1.0 - tx) + c10 * tx,
                   c1 = c01 * (1.0 - tx) + c11 * tx;
      return c0 * (1.0 - ty) + c1 * ty;
    }
  sample_axis(static_cast<double>(position[2]), z0, z1, tz);
  const double c000 = concentration(x0, y0, z0),
               c100 = concentration(x1, y0, z0),
               c010 = concentration(x0, y1, z0),
               c110 = concentration(x1, y1, z0),
               c001 = concentration(x0, y0, z1),
               c101 = concentration(x1, y0, z1),
               c011 = concentration(x0, y1, z1),
               c111 = concentration(x1, y1, z1);
  const double c00 = c000 * (1.0 - tx) + c100 * tx,
               c10 = c010 * (1.0 - tx) + c110 * tx,
               c01 = c001 * (1.0 - tx) + c101 * tx,
               c11 = c011 * (1.0 - tx) + c111 * tx;
  const double c0 = c00 * (1.0 - ty) + c10 * ty,
               c1 = c01 * (1.0 - ty) + c11 * ty;
  return c0 * (1.0 - tz) + c1 * tz;
}
// -----------------------------------------------------------------------------
inline
Real3 SampleFieldGradientAtPosition(const DiffusionGrid* grid,
                                  const Real3& position,
                                  const SamplingMode mode,
                                  const bool sample_in_2d = false)
{
  Real3 grad = {0.0, 0.0, 0.0};
  if (!grid) return grad;
  const size_t res = grid->GetResolution();
  const bool interpolate =
    (mode != SamplingMode::kNearestCell) && res > 1;
  if (!interpolate)
    {
      grid->GetGradient(position, &grad);
      if (sample_in_2d) grad[2] = 0.0;
      return grad;
    }
  const double h = grid->GetBoxLength();
  if (h <= 0.0)
    {
      grid->GetGradient(position, &grad);
      if (sample_in_2d) grad[2] = 0.0;
      return grad;
    }
  const real_t half_h = static_cast<real_t>(0.5 * h);
  const Real3 px0({position[0] - half_h, position[1], position[2]}),
            px1({position[0] + half_h, position[1], position[2]}),
            py0({position[0], position[1] - half_h, position[2]}),
            py1({position[0], position[1] + half_h, position[2]}),
            pz0({position[0], position[1], position[2] - half_h}),
            pz1({position[0], position[1], position[2] + half_h});
  const double gx0 = SampleFieldAtPosition(grid, px0, mode, sample_in_2d),
               gx1 = SampleFieldAtPosition(grid, px1, mode, sample_in_2d),
               gy0 = SampleFieldAtPosition(grid, py0, mode, sample_in_2d),
               gy1 = SampleFieldAtPosition(grid, py1, mode, sample_in_2d);
  grad[0] = static_cast<real_t>((gx1 - gx0) / h);
  grad[1] = static_cast<real_t>((gy1 - gy0) / h);
  if (!sample_in_2d)
    {
      const double gz0 = SampleFieldAtPosition(grid, pz0, mode, false),
                   gz1 = SampleFieldAtPosition(grid, pz1, mode, false);
      grad[2] = static_cast<real_t>((gz1 - gz0) / h);
    }
  return grad;
}
// -----------------------------------------------------------------------------
inline
SamplingMode ParseSamplingModeString(const std::string& mode_name,
                                     const SamplingMode default_mode)
{
  if (mode_name == "nearest_cell") return SamplingMode::kNearestCell;
  if (mode_name == "linear_interpolation"
      || mode_name == "bilinear_interpolation")
    return SamplingMode::kLinearInterpolation;
  if (mode_name == "trilinear_interpolation") return SamplingMode::kTrilinearInterpolation;
  if (mode_name == "safe_boundary_clamp") return SamplingMode::kSafeBoundaryClamp;
  return default_mode;
}
// -----------------------------------------------------------------------------
// Parameter-aware field queries (CAP / diffusion_grid config).
inline
SamplingMode ResolveFieldSamplingMode(const ::Parameters* params)
{
  SamplingMode mode = SamplingMode::kNearestCell;
  if (nullptr == params)
    return mode;
  if (params->have_parameter<std::string>("diffusion_grid/sampling_mode"))
    return ParseSamplingModeString(
      params->get<std::string>("diffusion_grid/sampling_mode"), mode);
  const bool cap_enabled =
    params->have_parameter<bool>("CAP/enabled")
    && params->get<bool>("CAP/enabled");
  if (cap_enabled)
    return SamplingMode::kNearestCell;
  if (params->have_parameter<bool>("diffusion_grid/trilinear_interpolation")
      && params->get<bool>("diffusion_grid/trilinear_interpolation"))
    return SamplingMode::kSafeBoundaryClamp;
  return mode;
}
// -----------------------------------------------------------------------------
inline
bool SimulationDomainIs2D(const ::Parameters* params)
{
  return nullptr != params
    && params->have_parameter<bool>("simulation_domain_is_2D")
    && params->get<bool>("simulation_domain_is_2D");
}
// -----------------------------------------------------------------------------
inline
double GetInterpolatedValue(const DiffusionGrid* dg,
                            const Double3& position,
                            const ::Parameters* params)
{
  const Real3 pos({
    static_cast<real_t>(position[0]),
    static_cast<real_t>(position[1]),
    static_cast<real_t>(position[2])});
  return SampleFieldAtPosition(
    dg, pos, ResolveFieldSamplingMode(params), SimulationDomainIs2D(params));
}
// -----------------------------------------------------------------------------
inline
Double3 GetInterpolatedGradient(const DiffusionGrid* dg,
                                const Double3& position,
                                const ::Parameters* params)
{
  const Real3 pos({
    static_cast<real_t>(position[0]),
    static_cast<real_t>(position[1]),
    static_cast<real_t>(position[2])});
  const Real3 grad = SampleFieldGradientAtPosition(
    dg, pos, ResolveFieldSamplingMode(params), SimulationDomainIs2D(params));
  return {static_cast<double>(grad[0]),
          static_cast<double>(grad[1]),
          static_cast<double>(grad[2])};
}
// =============================================================================
}  // namespace bdm
// =============================================================================
#endif  // _MIGRATION_FIELD_SAMPLING_H_

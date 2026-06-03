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
#ifndef _ABM4bio_H_
#define _ABM4bio_H_
// =============================================================================
#include "./global.h"
#include "./csv_io.h"
#include "./biology.h"
#include "./biological_cell.h"
#include "./cell_protrusion.h"
#include "./vessel.h"
namespace bdm {
inline
double GetInterpolatedValue(const bdm::DiffusionGrid* dg,
                            const bdm::Double3& position,
                            const ::Parameters* params)
{
  if (!dg) return 0.0;
  const bool use_trilinear =
    nullptr != params
    && (!params->have_parameter<bool>("CAP/enabled")
        || !params->get<bool>("CAP/enabled"))
    && params->have_parameter<bool>("diffusion_grid/trilinear_interpolation")
    && params->get<bool>("diffusion_grid/trilinear_interpolation");
  if (!use_trilinear)
    return dg->GetValue(position);
  const size_t res = dg->GetResolution();
  if (res<=1)
    return dg->GetValue(position);
  const auto dims = dg->GetDimensions();
  const double grid_min = static_cast<double>(dims[0]),
               grid_max = static_cast<double>(dims[1]),
               box_len = dg->GetBoxLength(),
               center_min = grid_min + 0.5*box_len,
               center_max = grid_max - 0.5*box_len;
  auto sample_axis = [=] (const double coord, uint32_t& i0, uint32_t& i1, double& t) {
    const double clamped = std::clamp(coord, center_min, center_max);
    const double g = (clamped-center_min)/box_len;
    int lower = static_cast<int>(std::floor(g));
    lower = std::clamp(lower, 0, static_cast<int>(res)-2);
    i0 = static_cast<uint32_t>(lower);
    i1 = static_cast<uint32_t>(lower+1);
    t = std::clamp(g-static_cast<double>(lower), 0.0, 1.0);
  };
  uint32_t x0=0, x1=0, y0=0, y1=0, z0=0, z1=0;
  double tx=0.0, ty=0.0, tz=0.0;
  sample_axis(position[0], x0, x1, tx);
  sample_axis(position[1], y0, y1, ty);
  sample_axis(position[2], z0, z1, tz);
  auto concentration = [&] (const uint32_t x, const uint32_t y, const uint32_t z) {
    const std::array<uint32_t,3> box = {x, y, z};
    return dg->GetConcentration(dg->GetBoxIndex(box));
  };
  const double c000 = concentration(x0, y0, z0),
               c100 = concentration(x1, y0, z0),
               c010 = concentration(x0, y1, z0),
               c110 = concentration(x1, y1, z0),
               c001 = concentration(x0, y0, z1),
               c101 = concentration(x1, y0, z1),
               c011 = concentration(x0, y1, z1),
               c111 = concentration(x1, y1, z1);
  const double c00 = c000*(1.0-tx) + c100*tx,
               c10 = c010*(1.0-tx) + c110*tx,
               c01 = c001*(1.0-tx) + c101*tx,
               c11 = c011*(1.0-tx) + c111*tx;
  const double c0 = c00*(1.0-ty) + c10*ty,
               c1 = c01*(1.0-ty) + c11*ty;
  return c0*(1.0-tz) + c1*tz;
}
inline
double GetNumericParameterFlexible(const ::Parameters& params,
                                   const std::string& parameter_name,
                                   const double fallback_value);

inline
double GetCapInitialFallbackValue(const ::Parameters& params,
                                  const std::string& biochemical_name)
{
  const bool cap_enabled =
    params.have_parameter<bool>("CAP/enabled")
    && params.get<bool>("CAP/enabled");
  if (!cap_enabled)
    return 0.0;

  const std::string cap_application_mode =
    params.have_parameter<std::string>("CAP/application_mode")
    ? params.get<std::string>("CAP/application_mode")
    : std::string("boundary_faces");
  const bool uniform_mode =
    (cap_application_mode == "uniform_grid")
    || (cap_application_mode == "uniform_domain");
  if (!uniform_mode)
    return 0.0;

  const double dt_h =
    params.have_parameter<double>("time_step")
    ? params.get<double>("time_step")
    : 1.0;
  const int start_step =
    params.have_parameter<int>("CAP/start_step")
    ? params.get<int>("CAP/start_step")
    : 0;
  const double start_time_h =
    params.have_parameter<double>("CAP/start_time_h")
    ? params.get<double>("CAP/start_time_h")
    : static_cast<double>(start_step) * dt_h;
  if (start_time_h > 1.0e-12)
    return 0.0;

  const std::string cap_parameter_name =
    biochemical_name == "H2O2"
    ? std::string("CAP/H2O2/concentration")
    : std::string("CAP/NO2_/concentration");
  return GetNumericParameterFlexible(params, cap_parameter_name, 0.0);
}

inline
double GetInitialValueMin(const ::Parameters& params,
                          const std::string& biochemical_name)
{
  const std::string scalar_name = biochemical_name+"/initial_value";
  const std::string min_name = biochemical_name+"/initial_value/min";
  if (params.have_parameter<double>(min_name))
    return params.get<double>(min_name);
  if (params.have_parameter<double>(scalar_name))
    return params.get<double>(scalar_name);
  if (biochemical_name == "H2O2" || biochemical_name == "NO2_")
    return GetCapInitialFallbackValue(params, biochemical_name);
  ABORT_("biochemical \""+biochemical_name+"\" has no initial value configured");
  return 0.0;
}
inline
double GetInitialValueMax(const ::Parameters& params,
                          const std::string& biochemical_name)
{
  const std::string scalar_name = biochemical_name+"/initial_value";
  const std::string max_name = biochemical_name+"/initial_value/max";
  if (params.have_parameter<double>(max_name))
    return params.get<double>(max_name);
  if (params.have_parameter<double>(scalar_name))
    return params.get<double>(scalar_name);
  if (biochemical_name == "H2O2" || biochemical_name == "NO2_")
    return GetCapInitialFallbackValue(params, biochemical_name);
  ABORT_("biochemical \""+biochemical_name+"\" has no initial value configured");
  return 0.0;
}
inline
double GetNumericParameterFlexible(const ::Parameters& params,
                                   const std::string& parameter_name,
                                   const double fallback_value)
{
  if (params.have_parameter<double>(parameter_name))
    return params.get<double>(parameter_name);
  if (params.have_parameter<std::string>(parameter_name))
    {
      const std::string& raw_value = params.get<std::string>(parameter_name);
      try
        {
          size_t parsed = 0;
          const double value = std::stod(raw_value, &parsed);
          ASSERT_(parsed == raw_value.size(),
                  "parameter \""+parameter_name+"\" must be numeric");
          return value;
        }
      catch (...)
        {
          ABORT_("parameter \""+parameter_name+"\" must be numeric");
        }
    }
  return fallback_value;
}
inline
void ValidateCapConcentrations(const double cap_h2o2,
                               const double cap_no2)
{
  ASSERT_(cap_h2o2 >= 0.0 && cap_h2o2 <= 1.0,
          "parameter \"CAP/H2O2/concentration\" must be within [0, 1]");
  ASSERT_(cap_no2 >= 0.0 && cap_no2 <= 1.0,
          "parameter \"CAP/NO2_/concentration\" must be within [0, 1]");
  ASSERT_(cap_h2o2 + cap_no2 <= 1.0 + 1.0e-12,
          "CAP species concentrations must sum to 1 or less");
}
inline
void ValidateNormalizedConcentrationValue(const double value,
                                          const std::string& parameter_name)
{
  ASSERT_(value >= 0.0 && value <= 1.0,
          "parameter \""+parameter_name+"\" must be within [0, 1]");
}
inline
void ValidateNormalizedConcentrationRange(const double min_value,
                                          const double max_value,
                                          const std::string& parameter_prefix)
{
  ASSERT_(min_value >= 0.0 && min_value <= 1.0,
          "parameter \""+parameter_prefix+"\" minimum must be within [0, 1]");
  ASSERT_(max_value >= 0.0 && max_value <= 1.0,
          "parameter \""+parameter_prefix+"\" maximum must be within [0, 1]");
  ASSERT_(min_value <= max_value,
          "parameter \""+parameter_prefix+"\" minimum must not exceed maximum");
}
inline
bool IsCapTransportSpecies(const std::string& biochemical_name)
{
  return biochemical_name == "H2O2" || biochemical_name == "NO2_";
}
inline
std::string GetCapTransportSolver(const ::Parameters& params)
{
  const std::string solver =
    params.have_parameter<std::string>("CAP/transport_solver")
    ? params.get<std::string>("CAP/transport_solver")
    : std::string("explicit");
  ASSERT_(solver == "explicit" || solver == "implicit",
          "parameter \"CAP/transport_solver\" must be \"explicit\" or \"implicit\"");
  return solver;
}
inline
bool UseImplicitCapTransport(const ::Parameters& params,
                             const std::string& biochemical_name)
{
  return params.have_parameter<bool>("CAP/enabled")
         && params.get<bool>("CAP/enabled")
         && IsCapTransportSpecies(biochemical_name)
         && GetCapTransportSolver(params) == "implicit";
}
inline
size_t FlattenGridIndex(const int x, const int y, const int z, const int res)
{
  return static_cast<size_t>((z * res + y) * res + x);
}
inline
void ApplyImplicitTransportStep(bdm::DiffusionGrid* dg,
                                const ::Parameters& params,
                                const std::string& biochemical_name,
                                const double dt_h)
{
  if (!dg) return;
  if (dt_h <= 0.0) return;
  const int res = static_cast<int>(dg->GetResolution());
  if (res <= 0) return;
  const double box_length = dg->GetBoxLength();
  if (box_length <= 0.0) return;
  const double diffusion_coefficient =
    params.have_parameter<double>(biochemical_name+"/diffusion_coefficient")
    ? params.get<double>(biochemical_name+"/diffusion_coefficient")
    : 0.0;
  const double dissipation_coefficient =
    params.have_parameter<double>(biochemical_name+"/dissipation_coefficient")
    ? params.get<double>(biochemical_name+"/dissipation_coefficient")
    : 0.0;
  ASSERT_(diffusion_coefficient >= 0.0,
          "diffusion coefficient must be non-negative for implicit CAP transport");
  ASSERT_(dissipation_coefficient >= 0.0,
          "dissipation coefficient must be non-negative for implicit CAP transport");
  if (diffusion_coefficient == 0.0 && dissipation_coefficient == 0.0) return;
  const int iterations =
    params.have_parameter<int>("CAP/implicit_transport/iterations")
    ? params.get<int>("CAP/implicit_transport/iterations")
    : 50;
  const double tolerance =
    params.have_parameter<double>("CAP/implicit_transport/tolerance")
    ? params.get<double>("CAP/implicit_transport/tolerance")
    : 1.0e-8;
  const double relaxation =
    params.have_parameter<double>("CAP/implicit_transport/relaxation")
    ? params.get<double>("CAP/implicit_transport/relaxation")
    : 1.0;
  ASSERT_(iterations > 0,
          "parameter \"CAP/implicit_transport/iterations\" must be positive");
  ASSERT_(tolerance >= 0.0,
          "parameter \"CAP/implicit_transport/tolerance\" must be non-negative");
  ASSERT_(relaxation > 0.0 && relaxation <= 1.0,
          "parameter \"CAP/implicit_transport/relaxation\" must be in (0, 1]");
  const size_t n = static_cast<size_t>(res) * static_cast<size_t>(res) * static_cast<size_t>(res);
  std::vector<double> previous(n, 0.0);
  std::vector<double> current(n, 0.0);
  for (int z=0; z<res; ++z)
    for (int y=0; y<res; ++y)
      for (int x=0; x<res; ++x)
        {
          const size_t flat = FlattenGridIndex(x, y, z, res);
          const std::array<uint32_t,3> box = {static_cast<uint32_t>(x),
                                              static_cast<uint32_t>(y),
                                              static_cast<uint32_t>(z)};
          previous[flat] = dg->GetConcentration(dg->GetBoxIndex(box));
          current[flat] = previous[flat];
        }
  const double alpha = diffusion_coefficient * dt_h / (box_length * box_length);
  const double beta = dissipation_coefficient * dt_h;
  for (int it=0; it<iterations; ++it)
    {
      double max_delta = 0.0;
      for (int z=0; z<res; ++z)
        for (int y=0; y<res; ++y)
          for (int x=0; x<res; ++x)
            {
              const size_t flat = FlattenGridIndex(x, y, z, res);
              double sum = 0.0;
              int n_neighbors = 0;
              auto add_neighbor = [&] (const int nx, const int ny, const int nz) {
                if (nx < 0 || nx >= res || ny < 0 || ny >= res || nz < 0 || nz >= res)
                  return;
                sum += current[FlattenGridIndex(nx, ny, nz, res)];
                ++n_neighbors;
              };
              add_neighbor(x-1, y, z);
              add_neighbor(x+1, y, z);
              add_neighbor(x, y-1, z);
              add_neighbor(x, y+1, z);
              add_neighbor(x, y, z-1);
              add_neighbor(x, y, z+1);
              const double denominator = 1.0 + beta + alpha * static_cast<double>(n_neighbors);
              const double candidate = (previous[flat] + alpha * sum) / denominator;
              const double updated =
                std::clamp(current[flat] + relaxation * (candidate - current[flat]), 0.0, 1.0);
              max_delta = std::max(max_delta, std::abs(updated - current[flat]));
              current[flat] = updated;
            }
      if (max_delta <= tolerance) break;
    }
  for (int z=0; z<res; ++z)
    for (int y=0; y<res; ++y)
      for (int x=0; x<res; ++x)
        {
          const size_t flat = FlattenGridIndex(x, y, z, res);
          const std::array<uint32_t,3> box = {static_cast<uint32_t>(x),
                                              static_cast<uint32_t>(y),
                                              static_cast<uint32_t>(z)};
          const size_t b = dg->GetBoxIndex(box);
          const double concentration = dg->GetConcentration(b);
          dg->ChangeConcentrationBy(b, current[flat] - concentration);
        }
}
}
#include "./biological_cell-inline.h"
#include "./cell_protrusion-inline.h"
#include "./vessel-inline.h"
#include "./biology4biologicalcell-inline.h"
#include "./biology4cellprotrusion-inline.h"
#include "./biology4vessel-inline.h"
#include "./obstacles.h"
#include "./io_flux.h"
// =============================================================================
// BioDynaMo model parameters
static
Parameters params;
static
std::vector<bdm::Double3> all_agents;
static
SimulationObstacles obstacles;
static
SimulationIOFlux io_flux;
static
std::vector<bdm::Double3> dg_vec;
static
std::map<int, int> vessel_map__ID_age;
// =============================================================================
inline
void save_stats(bdm::Simulation& sim,
                const std::map<int, std::string>& cells,
                std::ofstream& fout )
{
  if (!fout.is_open()) return;
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  //
  unsigned int n_cell = 0;
  unsigned int n_vessel = 0;
  unsigned int n_soma = 0;
  //
  std::map<int, unsigned int> n_cell_per_phenotype;
  // ...per phase of the cell circle
  std::map<int, unsigned int> n_cell_per_phenotype__Ap;
  std::map<int, unsigned int> n_cell_per_phenotype__G1;
  std::map<int, unsigned int> n_cell_per_phenotype__Sy;
  std::map<int, unsigned int> n_cell_per_phenotype__G2;
  std::map<int, unsigned int> n_cell_per_phenotype__Di;
  std::map<int, unsigned int> n_cell_per_phenotype__Tr;
  // iterate for all cell phenotypes
  for ( std::map<int, std::string>::const_iterator
        ci=cells.begin(); ci!=cells.end(); ci++ )
    {
      // obtain the ID & name for this cell phenotype
      const int CP_ID = ci->first;
      //const std::string& CP_name = ci->second;
      // initialize the map (counter) for this cell phenotype
      n_cell_per_phenotype[CP_ID] = 0;
      n_cell_per_phenotype__Ap[CP_ID] = 0;
      n_cell_per_phenotype__G1[CP_ID] = 0;
      n_cell_per_phenotype__Sy[CP_ID] = 0;
      n_cell_per_phenotype__G2[CP_ID] = 0;
      n_cell_per_phenotype__Di[CP_ID] = 0;
      n_cell_per_phenotype__Tr[CP_ID] = 0;
    }
  unsigned int n_protrusion = 0;
  //
  // Variables for tumor volume calculation (length * width^2 / 2)
  std::vector<bdm::Double3> cancer_positions;
  double tumor_volume = 0.0;
  //
  rm->ForEachAgent([&] (bdm::Agent* a) {
    if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
      {
        if (nullptr!=protrusion->GetCell())
          ++n_protrusion;
      }
    else if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
      {
        ++n_cell;
        // ID for this cell phenotype
        const int CP_ID = cell->GetPhenotype();
        // phase of the cell circle
        const bdm::BiologicalCell::Phase CP_Ph =
          static_cast<bdm::BiologicalCell::Phase>( cell->GetPhase() );
        // update this component on the map..
        n_cell_per_phenotype[CP_ID] += 1;
        //
        // Store cancer cell positions for volume calculation (exclude necrotic cells with ID=0)
        if (CP_ID > 0) {
          cancer_positions.push_back(cell->GetPosition());
        }
        //
        if      ( bdm::BiologicalCell::Phase::Ap==CP_Ph )
          n_cell_per_phenotype__Ap[CP_ID] += 1;
        else if ( bdm::BiologicalCell::Phase::I0==CP_Ph ||
                  bdm::BiologicalCell::Phase::G1==CP_Ph )
          n_cell_per_phenotype__G1[CP_ID] += 1;
        else if ( bdm::BiologicalCell::Phase::Sy==CP_Ph )
          n_cell_per_phenotype__Sy[CP_ID] += 1;
        else if ( bdm::BiologicalCell::Phase::G2==CP_Ph )
          n_cell_per_phenotype__G2[CP_ID] += 1;
        else if ( bdm::BiologicalCell::Phase::Di==CP_Ph )
          n_cell_per_phenotype__Di[CP_ID] += 1;
        else if ( bdm::BiologicalCell::Phase::Tr==CP_Ph )
          n_cell_per_phenotype__Tr[CP_ID] += 1;
      }
    else if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
      {
        ++n_vessel;
      }
    else if (auto* soma = dynamic_cast<bdm::neuroscience::NeuronSoma*>(a))
      {
        ++n_soma;
      }
  });
  
  // Calculate tumor volume using formula: length * width^2 / 2
  if (cancer_positions.size() > 1) {
    // Find bounding box of cancer cells
    double min_x = cancer_positions[0][0], max_x = cancer_positions[0][0];
    double min_y = cancer_positions[0][1], max_y = cancer_positions[0][1];
    double min_z = cancer_positions[0][2], max_z = cancer_positions[0][2];
    
    for (const auto& pos : cancer_positions) {
      min_x = std::min(min_x, pos[0]);
      max_x = std::max(max_x, pos[0]);
      min_y = std::min(min_y, pos[1]);
      max_y = std::max(max_y, pos[1]);
      min_z = std::min(min_z, pos[2]);
      max_z = std::max(max_z, pos[2]);
    }
    
    // Calculate dimensions
    double length = max_x - min_x;  // longest dimension
    double width_y = max_y - min_y;
    double width_z = max_z - min_z;
    
    // Use the largest width dimension
    double width = std::max(width_y, width_z);
    
    // Calculate volume: length * width^2 / 2 (ellipsoidal approximation)
    tumor_volume = length * width * width / 2.0;
  }
  //
  // simulation statistics output data - write the CSV file header
  // only once (@zero time step)
  if ( 0.0 == params.get<double>("current time") )
    {
      fout << "current_time, N_vessels, N_cells, tumor_volume_mm3";
      for ( std::map<int, std::string>::const_iterator
            ci=cells.begin(); ci!=cells.end(); ci++ )
        {
          // obtain the ID & name for this cell phenotype
          const int CP_ID = ci->first;
          //+const std::string& CP_name = ci->second;
          fout << ", N_cells_pheno_" << CP_ID << "";
          // by default the necrotic cells have ID equal to zero
          if (CP_ID) // hence, we ignore taking the stats below
            {
              fout << ", N_cells_pheno_" << CP_ID << "_Ap";
              fout << ", N_cells_pheno_" << CP_ID << "_G1";
              fout << ", N_cells_pheno_" << CP_ID << "_Sy";
              fout << ", N_cells_pheno_" << CP_ID << "_G2";
              fout << ", N_cells_pheno_" << CP_ID << "_Di";
              fout << ", N_cells_pheno_" << CP_ID << "_Tr";
            }
        }
      fout << ", N_cell_protrusions";
      fout << std::endl;
    }
  //
  fout << params.get<double>("current time");
  fout //<< ',' << n_soma
       << ',' << n_vessel;
  fout << ',' << n_cell;
  fout << ',' << tumor_volume * 1.0e-9; // Add tumor volume output
  // iterate for all cell phenotypes
  for ( std::map<int, std::string>::const_iterator
        ci=cells.begin(); ci!=cells.end(); ci++ )
    {
      // obtain the ID & name for this cell phenotype
      const int CP_ID = ci->first;
      //+const std::string& CP_name = ci->second;
      fout << ',' << n_cell_per_phenotype[CP_ID];
      // by default the necrotic cells have ID equal to zero
      if (CP_ID) // hence, we ignore taking the stats below
        {
          fout << ',' << n_cell_per_phenotype__Ap[CP_ID];
          fout << ',' << n_cell_per_phenotype__G1[CP_ID];
          fout << ',' << n_cell_per_phenotype__Sy[CP_ID];
          fout << ',' << n_cell_per_phenotype__G2[CP_ID];
          fout << ',' << n_cell_per_phenotype__Di[CP_ID];
          fout << ',' << n_cell_per_phenotype__Tr[CP_ID];
        }
    }
  fout << ',' << n_protrusion;
  fout << std::endl;
}
// =============================================================================
inline
void read_csv_file(const std::string& fn,
                   std::map<int, std::string>& cells,
                   std::vector<std::string>& biochem)
{
  csv_io::CSVReader<3> in(fn);
  //
  in.read_header(csv_io::ignore_extra_column, "parameter_name", "type", "parameter_value");
  // after reading the header, now read the rest of the lines in the file
  std::string parameter_name, type, parameter_value;
  while (in.read_row(parameter_name, type, parameter_value))
    {
      // ignore this line if...
      if ( "#" == parameter_name && "#" == type && "#" == parameter_value )
        continue;
      if ( "#" == parameter_name || '#' == parameter_name[0] )
        continue;
      // distinguish amongst the different input types
      if ( "string" == type )
        {
          params.set<std::string>(parameter_name) = parameter_value;
        }
      else if ( "bool" == type )
        {
          if      ("true" ==parameter_value) params.set<bool>(parameter_name) = true;
          else if ("false"==parameter_value) params.set<bool>(parameter_name) = false;
          else ABORT_("boolean type of data unable to read");
        }
      else if ( "int" == type )
        {
          params.set<int>(parameter_name) = std::stoi(parameter_value);
        }
      else if ( "float" == type )
        {
          params.set<double>(parameter_name) = std::stod(parameter_value);
        }
      else
        {
          ABORT_("unrecognized type of data to read");
        }
      // ...end of line to process loop
    }
  // ...at this point, the CSV file is complete!!!
  // now, proceed with some important initializations as well as
  // with some important checks prior to starting the simulation
  //
  // some crucial sanity checks
  if (!params.have_parameter<bool>("clean_output_directory"))
    params.set<bool>("clean_output_directory") = false;
  if (!params.have_parameter<std::string>("output_directory"))
    ABORT_("fatal error; parameter \"output_directory\" not inserted");
  if (!params.have_parameter<std::string>("simulation_title"))
    ABORT_("fatal error; parameter \"simulation_title\" not inserted");
  // create the simulation directory and sub-directories
  {
    std::string cmd;
    //
    if (params.get<bool>("clean_output_directory"))
      {
        cmd = "rm -rf " + params.get<std::string>("output_directory");
        std::system(cmd.c_str());
      }
    //
    cmd = "mkdir " + params.get<std::string>("output_directory");
    if (params.get<bool>("clean_output_directory"))
      {
        ASSERT_(0==std::system(cmd.c_str()),
                "could not create directory: "+params.get<std::string>("output_directory"));
      }
    else
      std::system(cmd.c_str());
    //
    cmd = "mkdir " + params.get<std::string>("output_directory") + "/in";
    if (params.get<bool>("clean_output_directory"))
      {
        ASSERT_(0==std::system(cmd.c_str()),
                "could not create directory: "+params.get<std::string>("output_directory"));
      }
    else
      std::system(cmd.c_str());
  }
  //
  // print-out all parameters read
  if ( false )
    print_helper(std::cout, &params);
  // preset initializations
  {
    params.set<int>("N_biochemicals") = 0;
    params.set<int>("N_phenotypes") = 0;
    //
    params.set<double>("current time") = 0.0;
    //
    if (!params.have_parameter<bool>("simplify_output"))
      params.set<bool>("simplify_output") = true;
    //
    if (!params.have_parameter<int>("number_of_time_steps"))
      ABORT_("model parameter \"number_of_time_steps\" is not provided");
    if (params.get<int>("number_of_time_steps")<1)
      ABORT_("model parameter \"number_of_time_steps\" is initialized wrong");
    if (!params.have_parameter<int>("statistics_interval"))
      params.set<int>("statistics_interval") = 1;
    if (params.get<int>("statistics_interval")<1)
      ABORT_("model parameter \"statistics_interval\" is initialized wrong");
    if (!params.have_parameter<int>("visualization_interval"))
      params.set<int>("visualization_interval") = 1;
    if (params.get<int>("visualization_interval")<1)
      ABORT_("model parameter \"visualization_interval\" is initialized wrong");
    //
    if (params.get<double>("max_boundary")<=params.get<double>("min_boundary"))
      ABORT_("model parameters \"min_boundary\", \"max_boundary\" are initialized wrong");
    //
    if (!params.have_parameter<int>("diffusion_grid/spatial_resolution"))
      ABORT_("model parameter \"diffusion_grid/spatial_resolution\" is not initialized");
    if (params.get<int>("diffusion_grid/spatial_resolution")<=2)
      ABORT_("model parameter \"diffusion_grid/spatial_resolution\" is initialized wrong");
    //
    if (!params.have_parameter<bool>("diffusion_grid/save_gradients"))
      params.set<bool>("diffusion_grid/save_gradients") = true;
    //
    if (!params.have_parameter<double>("domain_tolerance"))
      params.set<double>("domain_tolerance") = 0.01*(params.get<double>("max_boundary")
                                                    -params.get<double>("min_boundary"));
    if (params.get<double>("domain_tolerance")<=0.0)
      ABORT_("model parameter \"domain_tolerance\" is initialized wrong");
    //
    const double tol = (params.get<double>("max_boundary")-params.get<double>("min_boundary"))
                     / (params.get<int>("diffusion_grid/spatial_resolution")-1);
    if (params.get<double>("domain_tolerance")>tol)
      std::cout << "Warning: model parameter \"domain_tolerance\" needs inspection (wrt value: " << tol << ")." << std::endl;
    // {
    //   params.set<double>("domain_tolerance") = tol;
    //   std::cout << "Warning: model parameter \"domain_tolerance\" has been updated!" << std::endl;
    // }
    //
    if (!params.have_parameter<double>("migration_tolerance"))
      params.set<double>("migration_tolerance") = 0.01*params.get<double>("domain_tolerance");
    if (params.get<double>("migration_tolerance")<=0.0)
      ABORT_("model parameter \"migration_tolerance\" is initialized wrong");
    //
    if (!params.have_parameter<double>("cell/max_displacement"))
      params.set<double>("cell/max_displacement") = 0.0;
    if (params.get<double>("cell/max_displacement")<0.0)
      ABORT_("model parameter \"cell/max_displacement\" is initialized wrong");
    //
    if (!params.have_parameter<bool>("simulation_models_vessels"))
      params.set<bool>("simulation_models_vessels") = false;
    //
    if (!params.have_parameter<int>("simulation_obstacles"))
      params.set<int>("simulation_obstacles") = 0;
    if (params.get<int>("simulation_obstacles")<0)
      ABORT_("model parameter \"simulation_obstacles\" is initialized wrong");
    if (!params.have_parameter<bool>("simulation_obstacles/update"))
      params.set<bool>("simulation_obstacles/update") = false;
    // produce a pointer to parameter of the simulation obstacles object
    params.set<SimulationObstacles*>("simulation_obstacles_data") = &obstacles;
    //
    if (!params.have_parameter<double>("safe_distance_ratio"))
      params.set<double>("safe_distance_ratio") = 0.5;
    else
      {
        if (params.get<double>("safe_distance_ratio")==0.0e+0)
          { ; }
        else if (params.get<double>("safe_distance_ratio")<1.0e-2)
          { ABORT_("model parameter \"safe_distance_ratio\" is initialized wrong"); }
        else if (params.get<double>("safe_distance_ratio")>2.0e+0)
          { ABORT_("model parameter \"safe_distance_ratio\" is initialized wrong"); }
        else if (params.get<double>("safe_distance_ratio")<0.0e+0)
          { ABORT_("model parameter \"safe_distance_ratio\" is initialized wrong"); }
        else
          { ; }
      }
    //
    if (params.get<bool>("simulation_domain_is_periodic"))
      {
        if (!params.have_parameter<double>("simulation_domain_is_periodic/antisymmetry"))
          params.set<double>("simulation_domain_is_periodic/antisymmetry") = true;
      }
  }
  //
  if ( params.get<bool>("simulation_models_vessels") )
    {
      if (!params.have_parameter<double>("default_vessel_stiffness"))
        ABORT_("model parameter \"default_vessel_stiffness\" is initialized wrong");
      if (!params.have_parameter<double>("default_vessel_adherence"))
        ABORT_("model parameter \"default_vessel_adherence\" is initialized wrong");
      if (!params.have_parameter<double>("default_vessel_spring_constant"))
        ABORT_("model parameter \"default_vessel_spring_constant\" is initialized wrong");
      if (!params.have_parameter<double>("min_vessel_length"))
        ABORT_("model parameter \"min_vessel_length\" is initialized wrong");
      if (!params.have_parameter<double>("max_vessel_length"))
        ABORT_("model parameter \"max_vessel_length\" is initialized wrong");
    }
  else
    {
      params.set<double>("default_vessel_stiffness") = 1.0e+2;
      params.set<double>("default_vessel_adherence") = 0.1;
      params.set<double>("default_vessel_spring_constant") = 1.0e+1;
      if (!params.have_parameter<double>("min_vessel_length"))
        params.set<double>("min_vessel_length") = 1.0;
      if (!params.have_parameter<double>("max_vessel_length"))
        params.set<double>("max_vessel_length") = 100.0;
    }
  // biochemical (cues) and cell phenotypes identifiction
  if ( true )
    {
      std::stringstream ss(params.get<std::string>("diffusion_grid/biochemicals"));
      //
      biochem.clear();
      while (ss.good())
        {
          std::string name;
          getline(ss, name, ' ');
          biochem.push_back( name );
        }
    }
  //
  if ( true )
    {
      std::stringstream ss(params.get<std::string>("cell/phenotypes"));
      //
      cells.clear();
      while (ss.good())
        {
          std::string name;
          getline(ss, name, ' ');
          const int ID = params.get<int>(name+"/phenotype_ID");
          cells.insert( std::make_pair(ID, name) );
        }
    }
  // iterate for all cell phenotypes
  for ( std::map<int, std::string>::const_iterator
        ci=cells.begin(); ci!=cells.end(); ci++ )
    {
      // obtain the ID & name for this cell phenotype
      const int CP_ID = ci->first;
      const std::string& CP_name = ci->second;
      // ...by default necrotic cells (ID:0) are ignored!
      if (CP_ID<1) continue;
      //
      if (!params.have_parameter<int>(CP_name+"/io_flux")) continue;
      if (0==params.get<int>(CP_name+"/io_flux")) continue;
      //
      for (int s=0; s<params.get<int>(CP_name+"/io_flux"); s++)
        {
          const std::string fid = std::to_string(s+1);
          //
          const int dt =
            params.have_parameter<int>(CP_name+"/io_flux/"+fid+"/time_step") ?
            params.get<int>(CP_name+"/io_flux/"+fid+"/time_step") : 1;
          ASSERT_(dt>0,
                  "erroneous parameter value for: \""+CP_name+"/io_flux/"+fid+"/time_step\"");
          //
          const int time_function = (int)
            params.get<double>(CP_name+"/io_flux/"+fid+"/time_function/0");
          //
          int n_F_params = 0;
          {
            if      (time_function== 0) n_F_params = 2; // flat function
            else if (time_function==10) n_F_params = 2; // linear function
            else if (time_function==20) n_F_params = 4; // step function
            else if (time_function==30) n_F_params = 5; // ramp function
            else if (time_function==40) n_F_params = 4; // gaussian function
            else if (time_function==50) n_F_params = 5; // logistic function
            else if (time_function==11) n_F_params = 3; // linear periodic function
            else if (time_function==21) n_F_params = 5; // step periodic function
            else if (time_function==31) n_F_params = 6; // ramp periodic function
            else if (time_function==41) n_F_params = 5; // gaussian periodic function
            else if (time_function==51) n_F_params = 6; // logistic periodic function
            else
              ABORT_("erroneous parameter value for: \""+CP_name+"/io_flux/"+fid+"/time_function/0\"");
          }
          //
          std::vector<double> F_params(n_F_params);
          for (int p=0; p<n_F_params; p++)
            {
              F_params[p] =
                params.get<double>(CP_name+"/io_flux/"+fid+"/time_function/"+std::to_string(p));
            }
          //
          std::string stl;
          stl = params.get<std::string>(CP_name+"/io_flux/"+fid+"/STL");
          //
          Surface surf;
          surf.init(dt, F_params, CP_ID, stl);
          //
          io_flux.surface.push_back(surf);
          //
          // create a copy of the file just processed
          const std::string cmd = "cp " + stl + "  "
                                + params.get<std::string>("output_directory")
                                + "/in/" + CP_name + ".io_flux." + fid + ".stl";
          ASSERT_(0==std::system(cmd.c_str()),
                  "could not save a copy of a data file");
          //
        }
      // ...end of cell phenotypes loop
    }
  // produce a pointer to parameter of the simulation io-flux surfaces object
  params.set<SimulationIOFlux*>("simulation_io_flux") = &io_flux;
}
// =============================================================================
inline
void init_obstacles()
{
  const unsigned int n_obstacles = params.get<int>("simulation_obstacles");
  // if there are no obstacles in the simulation, then exit normally
  if (0==n_obstacles) return;
  // iterate for all simulation obstacles
  for (unsigned int l=0; l<n_obstacles; l++)
    {
      const std::string oid = std::to_string(l+1);
      // define the pattern of the simulation obstacle (use template, or load from STL file)
      const std::string pattern = params.get<std::string>("simulation_obstacle/"+oid+"/pattern");
      //
      if ("scaffold"==pattern)
        {
          //
          std::string fn;
          fn = params.get<std::string>("simulation_obstacle/"+oid+"/pattern/scaffold");
          //
          ObstacleScaffold obs;
          obs.init(pattern, fn);
          //
          obstacles.scaffold.push_back(obs);
          //
          // create a copy of the file just processed
          const std::string cmd = "cp " + fn + "  "
                                + params.get<std::string>("output_directory")
                                + "/in/simulation_obstacle." + oid + ".scaffold";
          ASSERT_(0==std::system(cmd.c_str()),
                  "could not save a copy of a data file");
          //
        }
      else if ("box/inside"==pattern || "box/outside"==pattern)
        {
          //
          std::vector<bdm::Double3> vertex(8);
          vertex[0] = { params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_A/0") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_A/1") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_A/2") };
          vertex[1] = { params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_B/0") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_B/1") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_B/2") };
          vertex[2] = { params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_C/0") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_C/1") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_C/2") };
          vertex[3] = { params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_D/0") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_D/1") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_D/2") };
          vertex[4] = { params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_E/0") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_E/1") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_E/2") };
          vertex[5] = { params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_F/0") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_F/1") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_F/2") };
          vertex[6] = { params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_G/0") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_G/1") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_G/2") };
          vertex[7] = { params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_H/0") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_H/1") ,
                        params.get<double>("simulation_obstacle/"+oid+"/pattern/box/point_H/2") };
          //
          ObstacleBox obs;
          obs.init(pattern, vertex);
          //
          obstacles.box.push_back(obs);
          //
        }
      else if ("sphere/inside"==pattern || "sphere/outside"==pattern)
        {
          //
          bdm::Double3 center;
          center = { params.get<double>("simulation_obstacle/"+oid+"/pattern/sphere/center/0") ,
                     params.get<double>("simulation_obstacle/"+oid+"/pattern/sphere/center/1") ,
                     params.get<double>("simulation_obstacle/"+oid+"/pattern/sphere/center/2") };
          //
          double radius;
          radius = params.get<double>("simulation_obstacle/"+oid+"/pattern/sphere/radius");
          //
          ObstacleSphere obs;
          obs.init(pattern, center, radius);
          //
          obstacles.sphere.push_back(obs);
          //
        }
      else if ("STL"==pattern)
        {
          //
          std::string fn;
          fn = params.get<std::string>("simulation_obstacle/"+oid+"/pattern/STL");
          //
          ObstacleSTL obs;
          obs.init(pattern, fn);
          //
          obstacles.surface.push_back(obs);
          //
          // create a copy of the file just processed
          const std::string cmd = "cp " + fn + "  "
                                + params.get<std::string>("output_directory")
                                + "/in/simulation_obstacle." + oid + ".stl";
          ASSERT_(0==std::system(cmd.c_str()),
                  "could not save a copy of the data file");
          //
        }
      else
        ABORT_("model parameter \""+pattern+"\" is initialized wrong");
      // ...end of simulation obstacles loop
    }
}
// =============================================================================
inline
void reinit_obstacles(const int time)
{
  const unsigned int n_obstacles = params.get<int>("simulation_obstacles");
  // if there are no obstacles in the simulation, then exit normally
  if (0==n_obstacles) return;
  // if there is no need to update the obstacles in the simulation, then again
  // exit normally
  if (!params.get<bool>("simulation_obstacles/update")) return;
  // iterate for all simulation obstacles
  for (unsigned int l=0; l<n_obstacles; l++)
    {
      const std::string oid = std::to_string(l+1);
      const std::string T = std::to_string(time);
      // define the pattern of the simulation obstacle (use template, or load from STL file)
      const std::string pattern = params.get<std::string>("simulation_obstacle/"+oid+"/pattern");
      //
      if ("scaffold"==pattern)
        {
          if (params.have_parameter<std::string>("simulation_obstacle/"+oid+"/pattern/scaffold/"+T))
            {
              //
              std::string fn;
              fn = params.get<std::string>("simulation_obstacle/"+oid+"/pattern/scaffold/"+T);
              //
              ObstacleScaffold obs;
              obs.init(pattern, fn);
              //
              obstacles.scaffold.push_back(obs);
              //
              // create a copy of the file just processed
              const std::string cmd = "cp " + fn + "  "
                                    + params.get<std::string>("output_directory")
                                    + "/in/simulation_obstacle." + oid + ".scaffold." + T;
              ASSERT_(0==std::system(cmd.c_str()),
                      "could not save a copy of a data file");
              //
            }
        }
      else if ("box/inside"==pattern || "box/outside"==pattern)
        {
          ;
        }
      else if ("sphere/inside"==pattern || "sphere/outside"==pattern)
        {
          ;
        }
      else if ("STL"==pattern)
        {
          if (params.have_parameter<std::string>("simulation_obstacle/"+oid+"/pattern/STL/"+T))
            {
              //
              std::string fn;
              fn = params.get<std::string>("simulation_obstacle/"+oid+"/pattern/STL/"+T);
              //
              ObstacleSTL obs;
              obs.init(pattern, fn);
              //
              obstacles.surface.push_back(obs);
              //
              // create a copy of the file just processed
              const std::string cmd = "cp " + fn + "  "
                                    + params.get<std::string>("output_directory")
                                    + "/in/simulation_obstacle." + oid + ".stl." + T;
              ASSERT_(0==std::system(cmd.c_str()),
                      "could not save a copy of the data file");
              //
            }
        }
      else
        ABORT_("model parameter \""+pattern+"\" is initialized wrong");
    }
}
// =============================================================================
inline
void set_bdm_params(bdm::Param* p)
{
  // simulation data to output
  const std::vector<std::string> out_obj{"GenericCell"};
  const std::vector<std::string> out_vars{"phenotype_","age_","diameter_","trail_","volume_",
                                          "n_divisions_","n_trasformations_"};
  const std::vector<std::string> out_type{"Int32","Int32","Float64","Float64","Float64",
                                          "Int32","Int32"};
  // BioDynaMo parameters that will be used for all simulations
  p->simulation_time_step = params.get<double>("time_step");
  p->bound_space = bdm::Param::BoundSpaceMode::kClosed;
  //p->bound_space = params.get<bool>("simulation_domain_is_bounded")
  //               ? bdm::Param::BoundSpaceMode::kClosed : bdm::Param::BoundSpaceMode::kOpen;
  p->min_bound = params.get<double>("min_boundary");
  p->max_bound = params.get<double>("max_boundary");
  p->simulation_max_displacement = params.get<double>("cell/max_displacement");
  p->detect_static_agents = params.get<double>("cell/max_displacement") ? true : false;
  p->calculate_gradients = params.get<bool>("diffusion_grid/save_gradients");
  p->environment = "uniform_grid";
  p->diffusion_method = "euler";
  // Selection of boundaries:
  // Dirichlet - Fixed concentration at boundaries (continuous O2 supply)
  // Neumann - Fixed gradient at boundaries
  // open - Substances can leak out (default)
  // closed - No flux through boundaries
  // Periodic - Wrap-around boundaries
  p->diffusion_boundary_condition = params.have_parameter<std::string>("diffusion_boundary_condition")
                                   ? params.get<std::string>("diffusion_boundary_condition")
                                   : "open";
  p->show_simulation_step = false;
  p->export_visualization = false;
  p->output_dir = params.get<std::string>("output_directory");
  p->remove_output_dir_contents = false;
  p->visualization_interval = params.get<int>("visualization_interval");
  for (size_t o=0; o<out_obj.size(); o++)
    {
      p->visualize_agents[out_obj[o]] =
        std::set<std::string>(out_vars.begin(), out_vars.end());
    }
  // BioDynaMo parameters associated with the vessels in the simulation
  {
    auto* p_ns = p->Get<bdm::neuroscience::Param>();
    //
    p_ns->neurite_default_spring_constant = params.get<double>("default_vessel_stiffness");
    p_ns->neurite_default_adherence       = params.get<double>("default_vessel_adherence");
    p_ns->neurite_default_spring_constant = params.get<double>("default_vessel_spring_constant");
    p_ns->neurite_min_length = params.get<double>("min_vessel_length");
    p_ns->neurite_max_length = params.get<double>("max_vessel_length");
    // complete list below (default values shown):
    //double  neurite_default_actual_length = 1.0
    //double  neurite_default_density = 1.0
    //double  neurite_default_diameter = 1.0
    //double  neurite_default_adherence = 0.1
    //double  neurite_default_tension = 0.0
    //double  neurite_minimial_bifurcation_length = 0
  }
}
// =============================================================================
inline
void set_convection(bdm::Simulation& sim, const int time = 0)
{
  //+if ( false )
  if (time==1)
    {
        const int N = params.get<int>("diffusion_grid/spatial_resolution");
        const double S_min = params.get<double>("min_boundary"),
                     S_max = params.get<double>("max_boundary"),
                     DS = (S_max-S_min) / N;
        // points of the Cartesian grid
        std::vector<bdm::Double3> points;
        for (int K=1; K<N; K++) {
          const double Z = S_min + DS * K;
          for (int J=1; J<N; J++) {
            const double Y = S_min + DS * J;
            for (int I=1; I<N; I++) {
              const double X = S_min + DS * I;
              // upload this point into the container
              points.push_back( {X, Y, Z} );
            }
          }
        }
        /*
        // check first for type of simulation domain
        if (params.get<bool>("simulation_domain_is_2D"))
        {
          // iterate for all points of the Cartesian grid
          for (int K=N/2-1; K<=N/2+1; K++) {
            const double Z = S_min + DS * K;
            for (int J=1; J<N; J++) {
              const double Y = S_min + DS * J;
              for (int I=1; I<N; I++) {
                const double X = S_min + DS * I;
                // upload this point into container
                C_vec.push_back( {X, Y, Z} );
              }
            }
          }
        }
        else
        {
          // iterate for all points of the Cartesian grid
          for (int K=1; K<N; K++) {
            const double Z = S_min + DS * K;
            for (int J=1; J<N; J++) {
              const double Y = S_min + DS * J;
              for (int I=1; I<N; I++) {
                const double X = S_min + DS * I;
                // upload this point into container
                C_vec.push_back( {X, Y, Z} );
              }
            }
          }
        }
        */
        // ...initialization is complete
        const std::string name = params.get<std::string>("output_directory")
                               + "/diffusion_grid.dat";
        std::ofstream fout(name);
        fout.precision(6);
        fout.setf(std::ios::scientific);
        //
        const int n_points = points.size();
        fout << n_points << std::endl;
        for (int P=0; P<n_points; P++)
          {
            fout << ' ' << points[P][0]
                 << ' ' << points[P][1]
                 << ' ' << points[P][2] << std::endl;
          }
        //
    }
  //
  if ( ! params.have_parameter<std::string>("convection/dynamic/from_file") )
    return;
  //
  const std::string Time = std::to_string(time);
  const int SpaceDimension = params.get<bool>("simulation_domain_is_2D") ? 2 : 3;
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  //
  std::string fn = params.get<std::string>("convection/dynamic/from_file");
  //
  const std::string s2f = "*.",
                    s2r = Time+".";
  size_t found = fn.find(s2f);
  if (std::string::npos != found)
    fn.replace(found, s2f.length(), s2r);
  else if (0==time)
    ; // simply use this filename...
  else if (0<time)
    return; // no need to reset the convection field...
  else
    ABORT_("convection has erroneous dynamic filename");
  //
  std::ifstream fin(fn);
  ASSERT_(fin.good(),"file \""+fn+"\" cannot be accessed");
  //
  std::set<size_t> scanned_boxes[3];
  //
  int npnt = 0;
  fin >> npnt;
  for (int p=0; p<npnt; p++)
    {
      bdm::Double3 xyz;
      fin >> xyz[0] >> xyz[1] >> xyz[2];
      bdm::Double3 v;
      fin >> v[0] >> v[1] >> v[2];
      //
      for (int ispdm=0; ispdm<SpaceDimension; ispdm++)
        {
          const std::string name = "convection_" + std::to_string(ispdm);
          //
          // access the BioDynaMo diffusion grid
          auto* dg = rm->GetDiffusionGrid(name);
          //
          dg->SetLowerThreshold(-1.0e+20);
          dg->SetUpperThreshold(+1.0e+20);
          //
          const size_t index = dg->GetBoxIndex(xyz);
          // WARNING: check identical piece of code writen in "init_biochemicals"
          if (scanned_boxes[ispdm].end()!=scanned_boxes[ispdm].find(index))
            continue;
          else
            scanned_boxes[ispdm].insert(index);
          //
          const double velocity = v[ispdm];
          const double original_velocity = dg->GetValue(xyz);
          //
          dg->ChangeConcentrationBy(xyz, velocity-original_velocity);
          // ...end of the simulation space dimension loop
        }
      // ...end of data points loop
    }
  // create a copy of the file just processed
  const std::string cmd = "cp " + fn + "  "
                        + params.get<std::string>("output_directory")
                        + "/in/convection." + Time + ".dat";
  ASSERT_(0==std::system(cmd.c_str()),
          "could not save a copy of a data file");
}
// =============================================================================
inline
void init_biochemicals(bdm::Simulation& sim,
                       const std::vector<std::string>& biochem)
{
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = sim.GetRandom();
  //
  // iterate for all biochemicals (substances)
  for (unsigned int icue=0; icue<biochem.size(); icue++)
    {
      const std::string& BC_name = biochem[icue];
      //
      Biochemical bc(Biochemical::N_A);
      // distinguish below for the kind of biochemical cue involved...
      {
        // radiation (not considered a biochemical cue!)
        if      ( BC_name == "RAD" ) bc = Biochemical::RAD;
        // fundamental biochemical cues
        else if ( BC_name == "OH_"  ) bc = Biochemical::OH_;
        else if ( BC_name == "O2"   ) bc = Biochemical::O2;
        else if ( BC_name == "O3"   ) bc = Biochemical::O3;
        else if ( BC_name == "H2O"  ) bc = Biochemical::H2O;
        else if ( BC_name == "H2O2" ) bc = Biochemical::H2O2;
        else if ( BC_name == "N2"   ) bc = Biochemical::N2;
        else if ( BC_name == "NO_"  ) bc = Biochemical::NO_;
        else if ( BC_name == "NO2"  ) bc = Biochemical::NO2;
        else if ( BC_name == "NO3"  ) bc = Biochemical::NO3;
        else if ( BC_name == "NO2_" ) bc = Biochemical::NO2_; // nitrite ion (NO2-)
        else if ( BC_name == "Gluc" ) bc = Biochemical::Gluc;
        // epidermal-/vessel-related biochemical cues
        else if ( BC_name == "VEGF" ) bc = Biochemical::VEGF;
        else if ( BC_name == "PDGF" ) bc = Biochemical::PDGF;
        else if ( BC_name == "PlGF" ) bc = Biochemical::PlGF;
        else if ( BC_name == "Ang1" ) bc = Biochemical::Ang1;
        else if ( BC_name == "Ang2" ) bc = Biochemical::Ang2;
        else if ( BC_name == "EGF"  ) bc = Biochemical::EGF;
        else if ( BC_name == "TGFa" ) bc = Biochemical::TGFa;
        else if ( BC_name == "TGFb" ) bc = Biochemical::TGFb;
        else if ( BC_name == "bFGF" ) bc = Biochemical::bFGF;
        // cancer-related biochemical cues
        else if ( BC_name == "TNF"  ) bc = Biochemical::TNF;
        // CAP-specific ICD markers and inflammatory cytokines
        else if ( BC_name == "CRT"  ) bc = Biochemical::CRT;   // calreticulin
        else if ( BC_name == "HMGB1") bc = Biochemical::HMGB1; // high mobility group box 1
        else if ( BC_name == "HSP70") bc = Biochemical::HSP70; // heat shock protein 70
        else if ( BC_name == "IL1b" ) bc = Biochemical::IL1b;  // interleukin-1 beta
        else if ( BC_name == "IL6"  ) bc = Biochemical::IL6;   // interleukin-6
        else if ( BC_name == "IL12" ) bc = Biochemical::IL12;  // interleukin-12
        else if ( BC_name == "CCL2" ) bc = Biochemical::CCL2;  // C-C motif chemokine ligand 2
        else if ( BC_name == "CCL4" ) bc = Biochemical::CCL4;  // C-C motif chemokine ligand 4
        // neuron-related biochemical cues
        else if ( BC_name == "NGF"  ) bc = Biochemical::NGF;
        else if ( BC_name == "BDNF" ) bc = Biochemical::BDNF;
        // drugs
        else if ( BC_name == "Drug_1" ) bc = Biochemical::Drug_1;
        else if ( BC_name == "Drug_2" ) bc = Biochemical::Drug_2;
        else if ( BC_name == "Drug_3" ) bc = Biochemical::Drug_3;
        // extracellular matrix (volume ratio / density)
        else if ( BC_name == "ECM" ) bc = Biochemical::ECM;
        // ...and exception is caught
        else
          ABORT_("unrecognized type \""+BC_name+"\" of a biochemical");
      }
      //
      const int sr = params.get<int>("diffusion_grid/spatial_resolution");
      // sanity checks
      ASSERT_(sr>=10,"spatial resolution for diffusion grid is erroneous");
      ASSERT_(sr%2==0,"spatial resolution for diffusion grid must be an even number");
      //
      // check which parameters haven't been loaded and set default values
      if (!params.have_parameter<double>(BC_name+"/diffusion_coefficient"))
        params.set<double>(BC_name+"/diffusion_coefficient") = 0.0;
      if (!params.have_parameter<double>(BC_name+"/dissipation_coefficient"))
        params.set<double>(BC_name+"/dissipation_coefficient") = 0.0;
      //
      const bool is_cap_species = bdm::IsCapTransportSpecies(BC_name);
      const bool implicit_cap_transport = bdm::UseImplicitCapTransport(params, BC_name);
      // diffusion (rate) and dissipation (rate) coefficient
      double dc =0.0, mu =0.0;
      if ( Biochemical::RAD != bc )
        {
          dc = params.get<double>(BC_name+"/diffusion_coefficient");
          mu = params.get<double>(BC_name+"/dissipation_coefficient");
          if (implicit_cap_transport)
            {
              dc = 0.0;
              mu = 0.0;
            }
        }
      // insert this biochemical (cue) in the BioDynaMo simulation
      bdm::ModelInitializer::DefineSubstance(bc, BC_name, dc, mu, sr);
      //
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(BC_name);
      dg->Initialize();
      //
      // set boundary condition type: closed for CAP species, Dirichlet for others
      const bool inlet_enabled = params.have_parameter<bool>(BC_name+"/inlet/enabled")
                               ? params.get<bool>(BC_name+"/inlet/enabled")
                               : false;

      const bool linear_gradient_enabled =
        params.have_parameter<bool>(BC_name+"/initial_value/linear_gradient/enabled")
        ? params.get<bool>(BC_name+"/initial_value/linear_gradient/enabled")
        : false;

      const std::string cap_application_mode =
        params.have_parameter<std::string>("CAP/application_mode")
        ? params.get<std::string>("CAP/application_mode")
        : std::string("boundary_faces");
      ASSERT_(cap_application_mode == "boundary_faces"
              || cap_application_mode == "boundary_dirichlet"
              || cap_application_mode == "boundary_neumann"
              || cap_application_mode == "uniform_grid"
              || cap_application_mode == "uniform_domain",
              "parameter \"CAP/application_mode\" must be one of: boundary_faces, boundary_dirichlet, boundary_neumann, uniform_grid, uniform_domain");

      if (is_cap_species)
        {
          if (cap_application_mode == "boundary_dirichlet")
            {
              dg->SetBoundaryConditionType(bdm::BoundaryConditionType::kDirichlet);
              dg->SetBoundaryCondition(
                std::make_unique<bdm::ConstantBoundaryCondition>(0.0));
            }
          else if (cap_application_mode == "boundary_neumann")
            {
              dg->SetBoundaryConditionType(bdm::BoundaryConditionType::kNeumann);
              dg->SetBoundaryCondition(
                std::make_unique<bdm::ConstantBoundaryCondition>(0.0));
            }
          else
            {
              dg->SetBoundaryConditionType(bdm::BoundaryConditionType::kClosedBoundaries);
            }
        }
      else if (inlet_enabled || linear_gradient_enabled)
        {
          dg->SetBoundaryConditionType(bdm::BoundaryConditionType::kClosedBoundaries);
        }
      else
        {
          // Dirichlet boundary: use a single ambient value. If min==max (usual), that
          // value applies; if a [min,max] range is used for interior initialization,
          // pin boundaries to the upper end (e.g. normoxic O2 at max of the 0..1 scale).
          const double min_iv = bdm::GetInitialValueMin(params, BC_name),
                       max_iv = bdm::GetInitialValueMax(params, BC_name);
          bdm::ValidateNormalizedConcentrationRange(min_iv, max_iv, BC_name+"/initial_value");
          const double boundary_value = (min_iv == max_iv) ? min_iv : max_iv;
          dg->SetBoundaryCondition(
            std::make_unique<bdm::ConstantBoundaryCondition>(boundary_value));
        }
      //
      // set the lower and upper threshold for the diffusion grid value range
      dg->SetLowerThreshold(0.0);
      dg->SetUpperThreshold(1.0);
      if (params.have_parameter<double>(BC_name+"/threshold/min"))
        {
          const double threshold_min = params.get<double>(BC_name+"/threshold/min");
          bdm::ValidateNormalizedConcentrationValue(threshold_min, BC_name+"/threshold/min");
          dg->SetLowerThreshold(threshold_min);
        }
      if (params.have_parameter<double>(BC_name+"/threshold/max"))
        {
          const double threshold_max = params.get<double>(BC_name+"/threshold/max");
          bdm::ValidateNormalizedConcentrationValue(threshold_max, BC_name+"/threshold/max");
          dg->SetUpperThreshold(threshold_max);
        }
      ASSERT_(dg->GetLowerThreshold() <= dg->GetUpperThreshold(),
              "biochemical \""+BC_name+"\" has erroneous threshold bounds");
      //
      // print-out the diffusion grid properties
      if ( false )
      if (0==icue)
        {
          std::cout << "Simulation reaction-diffusion domain..." << std::endl;
          std::cout << "Number_of_Boxes=" << dg->GetNumBoxes() << "; ";
          std::cout << "Grid_Size[0]=" << dg->GetGridSize()[0] << "; ";
          std::cout << "Grid_Size[1]=" << dg->GetGridSize()[1] << "; ";
          std::cout << "Grid_Size[2]=" << dg->GetGridSize()[2] << "; ";
          std::cout << "Resolution=" << dg->GetResolution() << "; ";
          std::cout << "Box_Length=" << dg->GetBoxLength() << "; ";
          std::cout << "Box_Volume=" << dg->GetBoxVolume() << "; ";
          auto bounds = sim.GetEnvironment()->GetDimensionThresholds();
          std::cout << "Bounds[0]=" << bounds[0] << "; ";
          std::cout << "Bounds[1]=" << bounds[1] << "; ";
          std::cout << "Delta=" << (bounds[1]-bounds[0])/(dg->GetGridSize()[0]-1.0) << "; ";
          std::cout << "Lower_Threshold=" << dg->GetLowerThreshold() << "; ";
          std::cout << "Upper_Threshold=" << dg->GetUpperThreshold() << "; ";
          std::cout << std::endl;
        }
      //
      const double minBC = bdm::GetInitialValueMin(params, BC_name),
                   maxBC = bdm::GetInitialValueMax(params, BC_name);
      // sanity check
      bdm::ValidateNormalizedConcentrationRange(minBC, maxBC, BC_name+"/initial_value");
      //
      if ( params.have_parameter<std::string>(BC_name+"/dynamic/from_file") )
        {
          std::string fn = params.get<std::string>(BC_name+"/dynamic/from_file");
          //
          const std::string s2f = "*.",
                            s2r = "0.";
          size_t found = fn.find(s2f);
          if (std::string::npos != found)
            fn.replace(found, s2f.length(), s2r);
          else
            ABORT_("biochemical \""+BC_name+"\" has erroneous dynamic filename");
          //
          params.set<std::string>(BC_name+"/initial_value/from_file") = fn;
        }
      //
      // initialize the concentration before simulation starts
      // Use AddInitializer + RunInitializers so values populate BOTH c1_ and c2_
      if ( ! params.have_parameter<std::string>(BC_name+"/initial_value/from_file") )
        {
          // Determine inlet parameters for the initializer lambda
          const std::string inlet_face = inlet_enabled && params.have_parameter<std::string>(BC_name+"/inlet/face")
                                       ? params.get<std::string>(BC_name+"/inlet/face")
                                       : std::string("");
          const double inlet_value = inlet_enabled && params.have_parameter<double>(BC_name+"/inlet/value")
                                   ? params.get<double>(BC_name+"/inlet/value")
                                   : maxBC;
          if (inlet_enabled)
            bdm::ValidateNormalizedConcentrationValue(inlet_value, BC_name+"/inlet/value");
          const double box_len = dg->GetBoxLength();
          const double grid_max = params.get<double>("max_boundary");
          const double grid_min = params.get<double>("min_boundary");
          const double center_min = grid_min + 0.5 * box_len;
          const double center_max = grid_max - 0.5 * box_len;
          
          // Check for linear gradient initialization along an axis
          const bool linear_gradient_enabled =
            params.have_parameter<bool>(BC_name+"/initial_value/linear_gradient/enabled")
            ? params.get<bool>(BC_name+"/initial_value/linear_gradient/enabled")
            : false;
          
          if ( linear_gradient_enabled )
            {
              const std::string grad_axis =
                params.get<std::string>(BC_name+"/initial_value/linear_gradient/axis");
              const double grad_min =
                params.get<double>(BC_name+"/initial_value/linear_gradient/min_value1");
              const double grad_max =
                params.get<double>(BC_name+"/initial_value/linear_gradient/max_value1");
              bdm::ValidateNormalizedConcentrationRange(
                grad_min, grad_max, BC_name+"/initial_value/linear_gradient");
              //
              ASSERT_(grad_axis=="x" || grad_axis=="y" || grad_axis=="z",
                      "linear_gradient axis must be \"x\", \"y\", or \"z\"");
              //
              // check for optional secondary / tertiary gradient axes
              const bool secondary_enabled =
                params.have_parameter<std::string>(BC_name+"/initial_value/linear_gradient/axis2");
              const std::string grad_axis2 = secondary_enabled
                ? params.get<std::string>(BC_name+"/initial_value/linear_gradient/axis2")
                : std::string("");
              const double grad2_min = secondary_enabled
                ? params.get<double>(BC_name+"/initial_value/linear_gradient/min_value2")
                : 0.0;
              const double grad2_max = secondary_enabled
                ? params.get<double>(BC_name+"/initial_value/linear_gradient/max_value2")
                : 0.0;
              if (secondary_enabled)
                bdm::ValidateNormalizedConcentrationRange(
                  grad2_min, grad2_max, BC_name+"/initial_value/linear_gradient/axis2");
              if (secondary_enabled)
                ASSERT_(grad_axis2=="x" || grad_axis2=="y" || grad_axis2=="z",
                        "linear_gradient axis2 must be \"x\", \"y\", or \"z\"");
              const bool tertiary_enabled =
                params.have_parameter<std::string>(BC_name+"/initial_value/linear_gradient/axis3");
              const std::string grad_axis3 = tertiary_enabled
                ? params.get<std::string>(BC_name+"/initial_value/linear_gradient/axis3")
                : std::string("");
              const double grad3_min = tertiary_enabled
                ? params.get<double>(BC_name+"/initial_value/linear_gradient/min_value3")
                : 0.0;
              const double grad3_max = tertiary_enabled
                ? params.get<double>(BC_name+"/initial_value/linear_gradient/max_value3")
                : 0.0;
              if (tertiary_enabled)
                bdm::ValidateNormalizedConcentrationRange(
                  grad3_min, grad3_max, BC_name+"/initial_value/linear_gradient/axis3");
              if (tertiary_enabled)
                ASSERT_(grad_axis3=="x" || grad_axis3=="y" || grad_axis3=="z",
                        "linear_gradient axis3 must be \"x\", \"y\", or \"z\"");
              //
              dg->AddInitializer([=](double x, double y, double z) -> double {
                auto axis_coord = [&] (const std::string& axis) -> double {
                  if      (axis == "x") return x;
                  else if (axis == "y") return y;
                  else                   return z;
                };
                auto axis_value = [&] (const std::string& axis, double vmin, double vmax) -> double {
                  const double coord = axis_coord(axis);
                  double t = (coord - center_min) / (center_max - center_min);
                  t = std::max(0.0, std::min(1.0, t));
                  return vmin + t * (vmax - vmin);
                };

                double value = axis_value(grad_axis, grad_min, grad_max);
                int n_axes = 1;
                if (secondary_enabled) {
                  value += axis_value(grad_axis2, grad2_min, grad2_max);
                  n_axes += 1;
                }
                if (tertiary_enabled) {
                  value += axis_value(grad_axis3, grad3_min, grad3_max);
                  n_axes += 1;
                }
                // Keep the configured value scale independent of number of axes.
                return value / static_cast<double>(n_axes);
              });
            }
          // Add initializer that sets uniform value (or random if min!=max) plus inlet on selected face
          else if ( minBC == maxBC )
            {
              dg->AddInitializer([=](double x, double y, double z) -> double {
                // Check if on inlet face
                if (!inlet_face.empty()) {
                  if (inlet_face == "z_max" && z >= grid_max - box_len) return inlet_value;
                  if (inlet_face == "z_min" && z <= grid_min + box_len) return inlet_value;
                  if (inlet_face == "x_max" && x >= grid_max - box_len) return inlet_value;
                  if (inlet_face == "x_min" && x <= grid_min + box_len) return inlet_value;
                  if (inlet_face == "y_max" && y >= grid_max - box_len) return inlet_value;
                  if (inlet_face == "y_min" && y <= grid_min + box_len) return inlet_value;
                }
                return maxBC;
              });
            }
          else
            {
              const int res = static_cast<int>(dg->GetResolution());
              const std::string inlet_face_local = inlet_enabled && params.have_parameter<std::string>(BC_name+"/inlet/face")
                                                 ? params.get<std::string>(BC_name+"/inlet/face")
                                                 : std::string("z_max");
              const double inlet_val_local = inlet_enabled && params.have_parameter<double>(BC_name+"/inlet/value")
                                           ? params.get<double>(BC_name+"/inlet/value")
                                           : maxBC;

              dg->AddInitializer([=](double x, double y, double z) -> double {
                const int ix = std::max(0, std::min(res - 1,
                    static_cast<int>(std::llround((x - center_min) / box_len))));
                const int iy = std::max(0, std::min(res - 1,
                    static_cast<int>(std::llround((y - center_min) / box_len))));
                const int iz = std::max(0, std::min(res - 1,
                    static_cast<int>(std::llround((z - center_min) / box_len))));

                // Deterministic voxel-wise pseudo-random initialization so c1_/c2_ match.
                uint64_t h = 1469598103934665603ULL;
                h ^= static_cast<uint64_t>(ix + 1); h *= 1099511628211ULL;
                h ^= static_cast<uint64_t>(iy + 1); h *= 1099511628211ULL;
                h ^= static_cast<uint64_t>(iz + 1); h *= 1099511628211ULL;
                h ^= static_cast<uint64_t>(icue + 1);
                const double u01 = static_cast<double>(h & 0xFFFFFFFFULL) / 4294967295.0;
                double value = minBC + u01 * (maxBC - minBC);

                if (inlet_enabled) {
                  if (inlet_face_local == "z_max" && z >= center_max) value = inlet_val_local;
                  else if (inlet_face_local == "z_min" && z <= center_min) value = inlet_val_local;
                  else if (inlet_face_local == "x_max" && x >= center_max) value = inlet_val_local;
                  else if (inlet_face_local == "x_min" && x <= center_min) value = inlet_val_local;
                  else if (inlet_face_local == "y_max" && y >= center_max) value = inlet_val_local;
                  else if (inlet_face_local == "y_min" && y <= center_min) value = inlet_val_local;
                }
                return value;
              });
            }
          // Note: RunInitializers() is called by the scheduler, which will populate both c1_ and c2_
        }
      else
        {
          const std::string fn = params.get<std::string>(BC_name+"/initial_value/from_file");
          //
          std::ifstream fin(fn);
          ASSERT_(fin.good(),"file \""+fn+"\" cannot be accessed");
          //
          std::set<size_t> scanned_boxes;
          //
          int npnt = 0;
          fin >> npnt;
          for (int p=0; p<npnt; p++)
            {
              bdm::Double3 xyz;
              fin >> xyz[0] >> xyz[1] >> xyz[2];
              double concentration;
              fin >> concentration;
              // sanity check
              if ( minBC>concentration || maxBC<concentration )
                ABORT_("biochemical \""+BC_name+"\" has been fed erroneous initial values");
              //
              const size_t index = dg->GetBoxIndex(xyz);
              // WARNING: check identical piece of code writen in "reinit_biochemicals"
              if (scanned_boxes.end()!=scanned_boxes.find(index))
                // ABORT_("biochemical \""+BC_name+"\" has multiple insertion for a DG-box");
                continue;
              else
                scanned_boxes.insert(index);
              //
              dg->ChangeConcentrationBy(xyz, concentration);
            }
          // create a copy of the file just processed
          const std::string cmd = "cp " + fn + "  "
                                + params.get<std::string>("output_directory")
                                + "/in/" + BC_name + ".0.dat";
          ASSERT_(0==std::system(cmd.c_str()),
                  "could not save a copy of a data file");
          // ...end of this if-case
        }
      // Note: inlet values are now set via AddInitializer (for uniform init) or in reinit_biochemicals (per step)
      // ...end of biochemicals (substances) loop
    }
}
// =============================================================================
inline
void init_convection(bdm::Simulation& sim)
{
  if (!params.have_parameter<std::string>("convection/dynamic/from_file"))
    return;
  const int SpaceDimension = params.get<bool>("simulation_domain_is_2D") ? 2 : 3;
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  //
  for (int ispdm=0; ispdm<SpaceDimension; ispdm++)
    {
      const std::string name = "convection_" + std::to_string(ispdm);
      const int id = 1000 + ispdm;
      //
      const int sr = params.get<int>("diffusion_grid/spatial_resolution");
      //
      // insert this biochemical (cue) in the BioDynaMo simulation
      bdm::ModelInitializer::DefineSubstance(id, name, 0.0, 0.0, sr);
      //
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(name);
      dg->Initialize();
      // ...end of the simulation space dimension loop
    }
  //
  // set the initial values for the convection field
  set_convection(sim, 0);
}
// =============================================================================
inline
void init_cells(bdm::Simulation& sim,
                const std::map<int, std::string>& cells,
                const std::vector<std::string>& biochem)
{
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = sim.GetRandom();
  // initialize the random number generator
  srand(static_cast<int>(time(0)));
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = params.get<double>("min_boundary"),
               maxCOORD = params.get<double>("max_boundary"),
               meanCOORD = (maxCOORD+minCOORD)/2.0,
               deltaCOORD = maxCOORD - meanCOORD,
               tol = params.get<double>("domain_tolerance");
  // calculate the min/max radius of all cells
  std::vector<double> cell_Dmin, cell_Dmax;
  // setup the spherical coordinates for all escaping cells:
  // radius, inclination (theta) and azimuth (phi)
  std::vector<double> radius, theta, phi;
  // cell phase for all escaping cells:
  std::vector<int> phase;
  // iterate for all cell phenotypes
  for ( std::map<int, std::string>::const_iterator
        ci=cells.begin(); ci!=cells.end(); ci++ )
    {
      // obtain the ID & name for this cell phenotype
      const int CP_ID = ci->first;
      const std::string& CP_name = ci->second;
      //
      if (CP_ID>0) // ignore necrotic cell whose default phenotype ID = 0
        {
          cell_Dmin.push_back(params.get<double>(CP_name+"/diameter/min"));
          cell_Dmax.push_back(params.get<double>(CP_name+"/diameter/max"));
        }
      //
      params.set<std::vector<double>>(CP_name+"/escaped_cells/radius") = radius;
      params.set<std::vector<double>>(CP_name+"/escaped_cells/theta")  = theta;
      params.set<std::vector<double>>(CP_name+"/escaped_cells/phi")    = phi;
      //
      params.set<std::vector<int>>(CP_name+"/escaped_cells/phase") = phase;
      //
      if (! params.have_parameter<double>(CP_name+"/initial_population/phase_Ap/percentage"))
        params.set<double>(CP_name+"/initial_population/phase_Ap/percentage") =
          0.0;
      if (! params.have_parameter<double>(CP_name+"/initial_population/phase_G1/percentage"))
        params.set<double>(CP_name+"/initial_population/phase_G1/percentage") =
          100.0 -
          params.get<double>(CP_name+"/initial_population/phase_Ap/percentage");
      if (! params.have_parameter<double>(CP_name+"/initial_population/phase_Sy/percentage"))
        params.set<double>(CP_name+"/initial_population/phase_Sy/percentage") =
          100.0 -
          params.get<double>(CP_name+"/initial_population/phase_G1/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_Ap/percentage");
      if (! params.have_parameter<double>(CP_name+"/initial_population/phase_G2/percentage"))
        params.set<double>(CP_name+"/initial_population/phase_G2/percentage") =
          100.0 -
          params.get<double>(CP_name+"/initial_population/phase_G1/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_Ap/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_Sy/percentage");
      if (! params.have_parameter<double>(CP_name+"/initial_population/phase_Di/percentage"))
        params.set<double>(CP_name+"/initial_population/phase_Di/percentage") =
          100.0 -
          params.get<double>(CP_name+"/initial_population/phase_G1/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_Ap/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_Sy/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_G2/percentage");
      if (! params.have_parameter<double>(CP_name+"/initial_population/phase_Tr/percentage"))
        params.set<double>(CP_name+"/initial_population/phase_Tr/percentage") =
          100.0 -
          params.get<double>(CP_name+"/initial_population/phase_G1/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_Ap/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_Sy/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_G2/percentage") -
          params.get<double>(CP_name+"/initial_population/phase_Di/percentage");
      //
      // ...end of cell phenotypes loop
    }
  const double min_radius = 0.5*(*std::min_element(cell_Dmin.begin(),cell_Dmin.end())),
               max_radius = 0.5*(*std::max_element(cell_Dmax.begin(),cell_Dmax.end()));
  const double safe_distance = (min_radius + max_radius)
                             * params.get<double>("safe_distance_ratio");
  //
  // iterate for all cell phenotypes
  for ( std::map<int, std::string>::const_iterator
        ci=cells.begin(); ci!=cells.end(); ci++ )
    {
      // obtain the ID & name for this cell phenotype
      const int CP_ID = ci->first;
      const std::string& CP_name = ci->second;
      //
      // check which parameters haven't been loaded and set default values
      if (! params.have_parameter<bool>(CP_name+"/can_apoptose"))
        params.set<bool>(CP_name+"/can_apoptose") = false;
      if (! params.have_parameter<bool>(CP_name+"/can_grow"))
        params.set<bool>(CP_name+"/can_grow") = false;
      if (! params.have_parameter<bool>(CP_name+"/can_divide"))
        params.set<bool>(CP_name+"/can_divide") = false;
      if (! params.have_parameter<bool>(CP_name+"/can_migrate"))
        params.set<bool>(CP_name+"/can_migrate") = false;
      if (! params.have_parameter<bool>(CP_name+"/can_transform"))
        params.set<bool>(CP_name+"/can_transform") = false;
      if (! params.have_parameter<bool>(CP_name+"/can_polarize"))
        params.set<bool>(CP_name+"/can_polarize") = false;
      if (! params.have_parameter<bool>(CP_name+"/can_protrude"))
        params.set<bool>(CP_name+"/can_protrude") = false;
      // default parameter(s) value
      if (params.get<bool>(CP_name+"/can_apoptose"))
        {
          if (! params.have_parameter<double>(CP_name+"/can_apoptose/probability"))
            params.set<double>(CP_name+"/can_apoptose/probability") = 0.0;
          if (! params.have_parameter<double>(CP_name+"/can_apoptose/probability_increment_with_age"))
            params.set<double>(CP_name+"/can_apoptose/probability_increment_with_age") = 0.0;
        }
      // default parameter(s) value
      if (params.get<bool>(CP_name+"/can_grow"))
        {
          if (! params.have_parameter<double>(CP_name+"/can_grow/probability"))
            params.set<double>(CP_name+"/can_grow/probability") = 0.0;
        }
      // default parameter(s) value
      if (params.get<bool>(CP_name+"/can_divide"))
        {
          if (! params.have_parameter<double>(CP_name+"/can_divide/probability"))
            params.set<double>(CP_name+"/can_divide/probability") = 0.0;
          if (! params.have_parameter<double>(CP_name+"/can_divide/probability_increment_with_age"))
            params.set<double>(CP_name+"/can_divide/probability_increment_with_age") = 0.0;
          if (! params.have_parameter<double>(CP_name+"/can_divide/influence_ratio"))
            params.set<double>(CP_name+"/can_divide/influence_ratio") = 0.0;
          // sanity check...
          ASSERT_(params.get<double>(CP_name+"/can_divide/influence_ratio")>=0.0,
                  "\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has erroneous value for \"can_divide/influence_ratio\"");
        
          // CAP-induced division probability modulation (default: no effect)
          if (! params.have_parameter<double>(CP_name+"/can_divide/CAP_sensitivity"))
            params.set<double>(CP_name+"/can_divide/CAP_sensitivity") = 0.0;
          if (! params.have_parameter<double>(CP_name+"/can_divide/CAP_saturation_time"))
            params.set<double>(CP_name+"/can_divide/CAP_saturation_time") = 60.0; // seconds
        }
      // default parameter(s) value
      if (params.get<bool>(CP_name+"/can_transform"))
        {
          if (! params.have_parameter<double>(CP_name+"/can_transform/probability"))
            params.set<double>(CP_name+"/can_transform/probability") = 0.0;
        }
      // default parameter(s) value
      if (params.get<bool>(CP_name+"/can_polarize"))
        {
          if (! params.have_parameter<double>(CP_name+"/can_polarize/probability"))
            params.set<double>(CP_name+"/can_polarize/probability") = 0.0;
        }
      // default parameter(s) value
      if (params.get<bool>(CP_name+"/can_migrate"))
        {
          if (! params.have_parameter<double>(CP_name+"/can_migrate/probability"))
            params.set<double>(CP_name+"/can_migrate/probability") = 0.0;
          if (! params.have_parameter<bool>(CP_name+"/can_migrate/accumulate_path"))
            params.set<bool>(CP_name+"/can_migrate/accumulate_path") = true;
          // encompass the effect of adhesion for the migratory cell to the ECM
          // by modulating the convection field respectively
          // ...default value for non-migratory cells
          params.set<double>(CP_name+"/can_migrate/max_adhesion/convection") = 0.0;
          //
          if (! params.have_parameter<std::string>("convection/dynamic/from_file"))
            params.set<double>(CP_name+"/can_migrate/max_adhesion/convection") = 0.0;
          // sanity check...
          ASSERT_(params.get<double>(CP_name+"/can_migrate/max_adhesion/convection")>=0.0,
                  "\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has erroneous value for \"max_adhesion/convection\"");
          //
          params.set<double>(CP_name+"/can_migrate/max_adhesion/displacement") =
              params.get<double>("time_step") *
              params.get<double>(CP_name+"/can_migrate/max_adhesion/convection");
        }
      // default parameter(s) value
      if (params.get<bool>(CP_name+"/can_protrude"))
        {
          if (! params.have_parameter<double>(CP_name+"/can_protrude/probability"))
            params.set<double>(CP_name+"/can_protrude/probability") = 0.0;
          if (! params.have_parameter<int>(CP_name+"/can_protrude/time_repeats"))
            params.set<int>(CP_name+"/can_protrude/time_repeats") = 1;
          if (! params.have_parameter<double>(CP_name+"/can_protrude/sprout/probability"))
            params.set<double>(CP_name+"/can_protrude/sprout/probability") = 0.0;
          if (! params.have_parameter<double>(CP_name+"/can_protrude/dissect/probability"))
            params.set<double>(CP_name+"/can_protrude/dissect/probability") = 0.0;
          if (! params.have_parameter<double>(CP_name+"/can_protrude/branch/probability"))
            params.set<double>(CP_name+"/can_protrude/branch/probability") = 0.0;
        }
      //
      if (CP_ID>0) // ignore necrotic cell whose default phenotype ID = 0
        {
          // identify the principal directions' order
          if (! params.have_parameter<double>(CP_name+"/principal/0") &&
              ! params.have_parameter<double>(CP_name+"/principal/1") &&
              ! params.have_parameter<double>(CP_name+"/principal/2") )
            params.set<double>(CP_name+"/principal/0") =
            params.set<double>(CP_name+"/principal/1") =
            params.set<double>(CP_name+"/principal/2") = 1.0;
          //
          // principal directions of the cell polarization matrix
          const double pd0 = (CP_ID==0 ? 1.0 : params.get<double>(CP_name+"/principal/0")),
                       pd1 = (CP_ID==0 ? 1.0 : params.get<double>(CP_name+"/principal/1")),
                       pd2 = (CP_ID==0 ? 1.0 : params.get<double>(CP_name+"/principal/2"));
          //
          std::vector<int> permutation(3);
          if      ( pd0 >= pd1 && pd1 >= pd2 ) permutation = { 0, 1, 2 };
          else if ( pd0 >= pd2 && pd2 >= pd1 ) permutation = { 0, 2, 1 };
          else if ( pd1 >= pd2 && pd2 >= pd0 ) permutation = { 1, 2, 0 };
          else if ( pd1 >= pd0 && pd0 >= pd2 ) permutation = { 1, 0, 2 };
          else if ( pd2 >= pd0 && pd0 >= pd1 ) permutation = { 2, 0, 1 };
          else if ( pd2 >= pd1 && pd1 >= pd0 ) permutation = { 2, 1, 0 };
          // ...and exception is caught
          else
            ABORT_("unrecognized order to put cell principal directions in order");
          //
          params.set<std::vector<int>>(CP_name+"/principal/permutation") = permutation;
        }
      // iterate for all biochemicals (substances)
      for (unsigned int icue=0; icue<biochem.size(); icue++)
        {
          const std::string& BC_name = biochem[icue];
          //
          if (! params.have_parameter<double>(CP_name+"/"+BC_name+"/secretion/net_balance"))
            {
              params.set<double>(CP_name+"/"+BC_name+"/secretion/net_balance") = 0.0;
              params.set<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std") = 0.0;
              params.set<double>(CP_name+"/"+BC_name+"/secretion/saturation") = 0.0;
            }
          else
            {
              if (! params.have_parameter<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std"))
                params.set<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std") = 0.0;
              if (! params.have_parameter<bool>(CP_name+"/"+BC_name+"/secretion/dependency"))
                params.set<bool>(CP_name+"/"+BC_name+"/secretion/dependency") = true;
              if (! params.have_parameter<double>(CP_name+"/"+BC_name+"/secretion/saturation"))
                params.set<double>(CP_name+"/"+BC_name+"/secretion/saturation") = 0.0;
            }
          //
          // --- Michaelis-Menten kinetics model (optional, per-substance) ---
          // When michaelis_menten_model = true, the standard net_balance pathway
          // is REPLACED by a self-contained concentration-dependent model:
          //   R = Vmax * C / (Km + C)
          // The sign of Vmax determines the direction:
          //   Vmax < 0 → consumption,  Vmax > 0 → production.
          // No net_balance, saturation, or dependency parameters are needed.
          if (! params.have_parameter<bool>(CP_name+"/"+BC_name+"/secretion/michaelis_menten_model"))
            params.set<bool>(CP_name+"/"+BC_name+"/secretion/michaelis_menten_model") = false;
          if (params.get<bool>(CP_name+"/"+BC_name+"/secretion/michaelis_menten_model"))
            {
              // Vmax: required — warn if missing, default to -0.001 (consumption)
              if (! params.have_parameter<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/Vmax"))
                {
                  params.set<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/Vmax") = -0.001;
                  bdm::Log::Warning("ABM4bio",
                    "Michaelis-Menten Vmax not provided for "+CP_name+"/"+BC_name
                    +"; defaulting to -0.001 (consumption)");
                }
              // Km: default to 1.0 if not provided
              if (! params.have_parameter<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/Km"))
                {
                  params.set<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/Km") = 1.0;
                  bdm::Log::Warning("ABM4bio",
                    "Michaelis-Menten Km not provided for "+CP_name+"/"+BC_name
                    +"; defaulting to 1.0");
                }
              //
              bdm::Log::Info("ABM4bio",
                "Michaelis-Menten model ENABLED for "+CP_name+"/"+BC_name
                +" (Vmax="+std::to_string(params.get<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/Vmax"))
                +", Km="+std::to_string(params.get<double>(CP_name+"/"+BC_name+"/secretion/michaelis_menten/Km"))+")");
            }
          // ...end of biochemicals (substances) loop
        }
      //
      // --- Default-substance chemotaxis shorthand ---
      // If user defines "CP/can_migrate/chemotaxis" as a double (without substance name),
      // automatically expand it to the first (default) substance with sensible defaults.
      // This avoids the need to explicitly specify substance-specific chemotaxis parameters.
      if (params.get<bool>(CP_name+"/can_migrate"))
        if (params.have_parameter<double>(CP_name+"/can_migrate/chemotaxis"))
          {
            const std::string& default_substance = biochem[0];
            const double chemo_weight = params.get<double>(CP_name+"/can_migrate/chemotaxis");
            // map to the first substance if not already explicitly defined
            if (! params.have_parameter<double>(CP_name+"/can_migrate/chemotaxis/"+default_substance))
              params.set<double>(CP_name+"/can_migrate/chemotaxis/"+default_substance) = chemo_weight;
            // set sensible defaults for sub-parameters if not explicitly provided
            if (! params.have_parameter<double>(CP_name+"/can_migrate/chemotaxis/"+default_substance+"/threshold"))
              params.set<double>(CP_name+"/can_migrate/chemotaxis/"+default_substance+"/threshold") = 1.0e-12;
            if (! params.have_parameter<double>(CP_name+"/can_migrate/chemotaxis/"+default_substance+"/probability"))
              params.set<double>(CP_name+"/can_migrate/chemotaxis/"+default_substance+"/probability") = 1.0;
            if (! params.have_parameter<bool>(CP_name+"/can_migrate/chemotaxis/"+default_substance+"/normalize_gradient"))
              params.set<bool>(CP_name+"/can_migrate/chemotaxis/"+default_substance+"/normalize_gradient") = false;
            //
            bdm::Log::Info("ABM4bio",
              "Chemotaxis shorthand expanded for "+CP_name+": using default substance \""+default_substance
              +"\" (weight="+std::to_string(chemo_weight)+")");
          }
      //
      // --- Default-substance transformation shorthand ---
      // If user defines "CP/can_transform/..." without a substance name, automatically
      // map it to the first (default) substance so single-substance models do not need
      // to repeat the biochemical name explicitly.
      if (params.get<bool>(CP_name+"/can_transform"))
        if (params.have_parameter<int>(CP_name+"/can_transform/new_phenotype"))
          {
            const std::string& default_substance = biochem[0];
            if (! params.have_parameter<int>(CP_name+"/can_transform/"+default_substance+"/new_phenotype"))
              params.set<int>(CP_name+"/can_transform/"+default_substance+"/new_phenotype") =
                params.get<int>(CP_name+"/can_transform/new_phenotype");
            if (params.have_parameter<double>(CP_name+"/can_transform/threshold") &&
                ! params.have_parameter<double>(CP_name+"/can_transform/"+default_substance+"/threshold"))
              params.set<double>(CP_name+"/can_transform/"+default_substance+"/threshold") =
                params.get<double>(CP_name+"/can_transform/threshold");
            if (params.have_parameter<double>(CP_name+"/can_transform/probability") &&
                ! params.have_parameter<double>(CP_name+"/can_transform/"+default_substance+"/probability"))
              params.set<double>(CP_name+"/can_transform/"+default_substance+"/probability") =
                params.get<double>(CP_name+"/can_transform/probability");
            if (params.have_parameter<int>(CP_name+"/can_transform/time_window_open") &&
                ! params.have_parameter<int>(CP_name+"/can_transform/"+default_substance+"/time_window_open"))
              params.set<int>(CP_name+"/can_transform/"+default_substance+"/time_window_open") =
                params.get<int>(CP_name+"/can_transform/time_window_open");
            if (params.have_parameter<int>(CP_name+"/can_transform/time_window_close") &&
                ! params.have_parameter<int>(CP_name+"/can_transform/"+default_substance+"/time_window_close"))
              params.set<int>(CP_name+"/can_transform/"+default_substance+"/time_window_close") =
                params.get<int>(CP_name+"/can_transform/time_window_close");
            bdm::Log::Info("ABM4bio",
              "Transformation shorthand expanded for "+CP_name+": using default substance \""+
              default_substance+"\"");
          }
      //
      // check if to adapt the cell initial population based on some pattern,
      // based on some random distribution, or based on some used-defined pattern
      // that is provided through a (ASCII) file
      params.set<bool>(CP_name+"/initial_population/pattern/box/inside")  = false;
      params.set<bool>(CP_name+"/initial_population/pattern/box/outside") = false;
      params.set<bool>(CP_name+"/initial_population/pattern/sphere/inside")  = false;
      params.set<bool>(CP_name+"/initial_population/pattern/sphere/outside") = false;
      //
      if (! params.have_parameter<std::string>(CP_name+"/initial_population/from_file"))
        {
          int number_of_cells = params.get<int>(CP_name+"/initial_population");
          // sanity check...
          ASSERT_(number_of_cells>=0,
                  "\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has erroneous initial population");
          // check if to consider simulation varying initial cell population
          if ( number_of_cells > 0 )
            if ( params.have_parameter<int>(CP_name+"/initial_population/std") )
              {
                const int noc_std = params.get<int>(CP_name+"/initial_population/std");
                // sanity check...
                ASSERT_(noc_std<=number_of_cells,
                        "\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has erroneous deviation from mean initial population");
                number_of_cells = uniform_distro(number_of_cells-noc_std, number_of_cells+noc_std);
              }
          //
          if ( params.have_parameter<std::string>(CP_name+"/initial_population/pattern") )
            {
              const std::string& pattern = params.get<std::string>(CP_name+"/initial_population/pattern");
              if      ("box/inside" ==pattern)
                params.set<bool>(CP_name+"/initial_population/pattern/box/inside")  = true;
              else if ("box/outside"==pattern)
                params.set<bool>(CP_name+"/initial_population/pattern/box/outside") = true;
              else if ("sphere/inside" ==pattern)
                params.set<bool>(CP_name+"/initial_population/pattern/sphere/inside")  = true;
              else if ("sphere/outside"==pattern)
                params.set<bool>(CP_name+"/initial_population/pattern/sphere/outside") = true;
              else
                ABORT_("\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has unrecognized pattern");
            }
          bdm::Double3 pntA, pntB;
          bdm::Double3 cntr; bdm::real_t rad = 0.0;
          if ( params.get<bool>(CP_name+"/initial_population/pattern/box/inside")  ||
               params.get<bool>(CP_name+"/initial_population/pattern/box/outside") )
            {
              pntA = { params.get<double>(CP_name+"/initial_population/pattern/box/point_A/0") ,
                       params.get<double>(CP_name+"/initial_population/pattern/box/point_A/1") ,
                       params.get<double>(CP_name+"/initial_population/pattern/box/point_A/2") };
              pntB = { params.get<double>(CP_name+"/initial_population/pattern/box/point_B/0") ,
                       params.get<double>(CP_name+"/initial_population/pattern/box/point_B/1") ,
                       params.get<double>(CP_name+"/initial_population/pattern/box/point_B/2") };
            }
          else if ( params.get<bool>(CP_name+"/initial_population/pattern/sphere/inside")  ||
                    params.get<bool>(CP_name+"/initial_population/pattern/sphere/outside") )
            {
              cntr = { params.get<double>(CP_name+"/initial_population/pattern/sphere/center/0") ,
                       params.get<double>(CP_name+"/initial_population/pattern/sphere/center/1") ,
                       params.get<double>(CP_name+"/initial_population/pattern/sphere/center/2") };
              rad = params.get<double>(CP_name+"/initial_population/pattern/sphere/radius");
            }
          //
          // generate the initial population of cells and initialize them
          std::vector<bdm::Double3> XYZ_;
          for (int i=0; i<number_of_cells; i++)
            {
              // cell coordinates (must fall within simulation domain)
              bdm::Double3 xyz{0.0, 0.0, 0.0};
              // check if a cell should follow any user-defined (spatial)
              // restrictions (wrt primitive shapes: box, sphere, cylinder)
              if      ( params.get<bool>(CP_name+"/initial_population/pattern/box/inside") )
                {
                  xyz = { rg->Uniform(pntA[0], pntB[0]) ,
                          rg->Uniform(pntA[1], pntB[1]) ,
                          rg->Uniform(pntA[2], pntB[2]) };
                }
              else if ( params.get<bool>(CP_name+"/initial_population/pattern/box/outside") )
                {
                  xyz = { coin_flip() ? rg->Uniform(minCOORD+tol, pntA[0]) : rg->Uniform(pntB[0], maxCOORD-tol) ,
                          coin_flip() ? rg->Uniform(minCOORD+tol, pntA[1]) : rg->Uniform(pntB[1], maxCOORD-tol) ,
                          coin_flip() ? rg->Uniform(minCOORD+tol, pntA[2]) : rg->Uniform(pntB[2], maxCOORD-tol) };
                }
              else if ( params.get<bool>(CP_name+"/initial_population/pattern/sphere/inside") )
                {
                  bdm::real_t r = rg->Uniform(0.0, rad);
                  bdm::real_t theta = rg->Uniform(0.0, 2.0*bdm::Math::kPi);
                  bdm::real_t u = rg->Uniform(-1.0, 1.0);
                  // https://mathworld.wolfram.com/SpherePointPicking.html
                  xyz = { cntr[0]+r*sqrt(1-u*u)*cos(theta) ,
                          cntr[1]+r*sqrt(1-u*u)*sin(theta) ,
                          cntr[2]+r*u                      };
                }
              else if ( params.get<bool>(CP_name+"/initial_population/pattern/sphere/outside") )
                {
                  bdm::real_t r = rg->Uniform(rad, maxCOORD-tol);
                  bdm::real_t theta = rg->Uniform(0.0, 2.0*bdm::Math::kPi);
                  bdm::real_t u = rg->Uniform(-1.0, 1.0);
                  // https://mathworld.wolfram.com/SpherePointPicking.html
                  xyz = { cntr[0]+r*sqrt(1-u*u)*cos(theta) ,
                          cntr[1]+r*sqrt(1-u*u)*sin(theta) ,
                          cntr[2]+r*u                      };
                }
              else
                {
                  xyz = { rg->Uniform(minCOORD+tol, maxCOORD-tol) ,
                          rg->Uniform(minCOORD+tol, maxCOORD-tol) ,
                          rg->Uniform(minCOORD+tol, maxCOORD-tol) };
                }
              // in case the simulation is in 2D then ensure that all
              // cells are positioned on a X-Y plane
              if ( params.get<bool>("simulation_domain_is_2D") )
                xyz[2] = meanCOORD;
              // in case of polar domain simulation, then check if
              // cell position is within a 3D sphere or a 2D circle
              if ( params.get<bool>("simulation_domain_is_polar") )
                {
                  const bdm::Double3 vec = { xyz[0]-meanCOORD ,
                                             xyz[1]-meanCOORD ,
                                             xyz[2]-meanCOORD };
                  // check radial distance (from the domain center)
                  if ( magnitude(vec) > deltaCOORD )
                    {
                      i -= 1;
                      continue;
                    }
                }
              // check if cell is positioned withing simulations bounds
              if ( xyz[0]<(minCOORD+tol) || xyz[0]>(maxCOORD-tol) ||
                   xyz[1]<(minCOORD+tol) || xyz[1]>(maxCOORD-tol) ||
                   xyz[2]<(minCOORD+tol) || xyz[2]>(maxCOORD-tol) )
                {
                  i -= 1;
                  continue;
                }
              // check also if cell position is suffiently apart to the rest
              // (recently created) of the cells
              bool is_valid = true;
              for (unsigned int c=0; c<all_agents.size(); c++)
                {
                  const bdm::Double3 vec = { xyz[0]-all_agents[c][0] ,
                                             xyz[1]-all_agents[c][1] ,
                                             xyz[2]-all_agents[c][2] };
                  // check distance with respect to other cells
                  if ( L2norm(vec) < safe_distance )
                    {
                      is_valid = false;
                      break;
                    }
                }
              if ( ! is_valid ) // now do what you must
                {
                  i -= 1;
                  continue;
                }
              //
              // save coordinates in temporary container
              XYZ_.push_back( xyz );
              // save spatial coordinates in the local container
              all_agents.push_back( xyz );
            }
          //
          const int number_of_cells__phase_Ap = CP_ID ?
                    number_of_cells*params.get<double>(CP_name+"/initial_population/phase_Ap/percentage")/100.0 : 0,
                    number_of_cells__phase_G1 = CP_ID ?
                    number_of_cells*params.get<double>(CP_name+"/initial_population/phase_G1/percentage")/100.0 : 0,
                    number_of_cells__phase_Sy = CP_ID ?
                    number_of_cells*params.get<double>(CP_name+"/initial_population/phase_Sy/percentage")/100.0 : 0,
                    number_of_cells__phase_G2 = CP_ID ?
                    number_of_cells*params.get<double>(CP_name+"/initial_population/phase_G2/percentage")/100.0 : 0,
                    number_of_cells__phase_Di = CP_ID ?
                    number_of_cells*params.get<double>(CP_name+"/initial_population/phase_Di/percentage")/100.0 : 0,
                    number_of_cells__phase_Tr = CP_ID ?
                    number_of_cells*params.get<double>(CP_name+"/initial_population/phase_Tr/percentage")/100.0 : 0;
          int ncells__phase_Ap = 0,
              ncells__phase_G1 = 0,
              ncells__phase_Sy = 0,
              ncells__phase_G2 = 0,
              ncells__phase_Di = 0,
              ncells__phase_Tr = 0;
          // order of mechanisms that define behaviour
          const int mo = params.get<int>(CP_name+"/mechanism_order");
          // create all cells for this phenotype ID
          for (int i=0; i<number_of_cells; i++)
            {
              // cell position
              const bdm::Double3& xyz = XYZ_[i];
              // cell cycle phase
              int ccp = 0;
              // by default the necrotic cells have ID equal to zero
              if (CP_ID) // hence, we ignore taking the stats below
                if (11==mo)
                  {
                    // check if the total number of cells for each
                    // cell cycle phase has been initiated, otherwise
                    // check for the next cell cycle phase and so on...
                    if (ncells__phase_Ap<number_of_cells__phase_Ap)
                      {
                        ++ncells__phase_Ap;
                        ccp = -1;
                      }
                    else if (ncells__phase_G1<number_of_cells__phase_G1)
                      {
                        ++ncells__phase_G1;
                        ccp = 1;
                      }
                    else if (ncells__phase_Sy<number_of_cells__phase_Sy)
                      {
                        ++ncells__phase_Sy;
                        ccp = 2;
                      }
                    else if (ncells__phase_G2<number_of_cells__phase_G2)
                      {
                        ++ncells__phase_G2;
                        ccp = 3;
                      }
                    else if (ncells__phase_Di<number_of_cells__phase_Di)
                      {
                        ++ncells__phase_Di;
                        ccp = 4;
                      }
                    else if (ncells__phase_Tr<number_of_cells__phase_Tr)
                      {
                        ++ncells__phase_Tr;
                        ccp = 5;
                      }
                  }
              // cell diameter
              const double dia = ( CP_ID==0 ? (min_radius+max_radius)
                                 : rg->Uniform(params.get<double>(CP_name+"/diameter/min"),
                                               params.get<double>(CP_name+"/diameter/max")));
              // principal directions of the cell polarization matrix
              const double pd0 = (CP_ID==0 ? 1.0 : params.get<double>(CP_name+"/principal/0")),
                           pd1 = (CP_ID==0 ? 1.0 : params.get<double>(CP_name+"/principal/1")),
                           pd2 = (CP_ID==0 ? 1.0 : params.get<double>(CP_name+"/principal/2"));
              // cell polarization matrix
              bdm::Double3x3 pl_3x3 = diag(pd0, pd1, pd2);
              //
              bdm::BiologicalCell* cell = new bdm::BiologicalCell(CP_ID, xyz);
              cell->SetPhase(ccp);
              cell->SetDiameter(dia);
              cell->SetAdherence(params.have_parameter<double>(CP_name+"/adherence")
                               ? params.get<double>(CP_name+"/adherence") : 0.0);
              if (params.have_parameter<double>(CP_name+"/density"))
                cell->SetDensity(params.get<double>(CP_name+"/density"));
              cell->SetParametersPointer(&params);
              cell->SetCanApoptose(params.get<bool>(CP_name+"/can_apoptose"));
              cell->SetCanGrow(params.get<bool>(CP_name+"/can_grow"));
              cell->SetCanDivide(params.get<bool>(CP_name+"/can_divide"));
              cell->SetCanMigrate(params.get<bool>(CP_name+"/can_migrate"));
              cell->SetCanTransform(params.get<bool>(CP_name+"/can_transform"));
              cell->SetCanPolarize(params.get<bool>(CP_name+"/can_polarize"));
              cell->SetCanProtrude(params.get<bool>(CP_name+"/can_protrude"));
              cell->SetPolarization(pl_3x3);
              if      (10==mo)
                cell->AddBehavior(new bdm::Biology4BiologicalCell_10());
              else if (11==mo)
                cell->AddBehavior(new bdm::Biology4BiologicalCell_11());
              else
                ABORT_("\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has unrecognized behavior");
              // store this cell into BioDynaMo's resource manager
              rm->AddAgent(cell);
            }
        }
      else
        {
          const std::string fn = params.get<std::string>(CP_name+"/initial_population/from_file");
          //
          std::ifstream fin(fn);
          ASSERT_(fin.good(),"file \""+fn+"\" cannot be accessed");
          //
          if ( params.have_parameter<std::string>(CP_name+"/initial_population/pattern") )
            {
              const std::string& pattern = params.get<std::string>(CP_name+"/initial_population/pattern");
              if      ("box/inside" ==pattern)
                params.set<bool>(CP_name+"/initial_population/pattern/box/inside")  = true;
              else if ("box/outside"==pattern)
                params.set<bool>(CP_name+"/initial_population/pattern/box/outside") = true;
              else if ("sphere/inside" ==pattern)
                params.set<bool>(CP_name+"/initial_population/pattern/sphere/inside")  = true;
              else if ("sphere/outside"==pattern)
                params.set<bool>(CP_name+"/initial_population/pattern/sphere/outside") = true;
              else
                ABORT_("\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has unrecognized pattern");
            }
          bdm::Double3 pntA, pntB;
          bdm::Double3 cntr; double rad = 0.0;
          if ( params.get<bool>(CP_name+"/initial_population/pattern/box/inside")  ||
               params.get<bool>(CP_name+"/initial_population/pattern/box/outside") )
            {
              pntA = { params.get<double>(CP_name+"/initial_population/pattern/box/point_A/0") ,
                       params.get<double>(CP_name+"/initial_population/pattern/box/point_A/1") ,
                       params.get<double>(CP_name+"/initial_population/pattern/box/point_A/2") };
              pntB = { params.get<double>(CP_name+"/initial_population/pattern/box/point_B/0") ,
                       params.get<double>(CP_name+"/initial_population/pattern/box/point_B/1") ,
                       params.get<double>(CP_name+"/initial_population/pattern/box/point_B/2") };
            }
          else if ( params.get<bool>(CP_name+"/initial_population/pattern/sphere/inside")  ||
                    params.get<bool>(CP_name+"/initial_population/pattern/sphere/outside") )
            {
              cntr = { params.get<double>(CP_name+"/initial_population/pattern/sphere/center/0") ,
                       params.get<double>(CP_name+"/initial_population/pattern/sphere/center/1") ,
                       params.get<double>(CP_name+"/initial_population/pattern/sphere/center/2") };
              rad = params.get<double>(CP_name+"/initial_population/pattern/sphere/radius");
            }
          // order of mechanisms that define behaviour
          const int mo = params.get<int>(CP_name+"/mechanism_order");
          //
          int number_of_cells = 0;
          fin >> number_of_cells;
          ////int number_of_cells_loaded = 0;
          for (int c=0; c<number_of_cells; c++)
            {
              // cell position
              bdm::Double3 xyz;
              fin >> xyz[0] >> xyz[1] >> xyz[2];
              // cell diameter
              double dia;
              fin >> dia;
              if ( CP_ID==0 )
                {
                  if ( dia<2.0*min_radius ||
                       dia>2.0*max_radius )
                    ABORT_("\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has erroneous diameter");
                }
              else
                {
                  if ( dia<params.get<double>(CP_name+"/diameter/min") ||
                       dia>params.get<double>(CP_name+"/diameter/max") )
                    ABORT_("\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has erroneous diameter");
                }
              // principal directions of the cell polarization matrix
              bdm::Double3 p0, p1, p2;
              fin >> p0[0] >> p0[1] >> p0[2];
              fin >> p1[0] >> p1[1] >> p1[2];
              fin >> p2[0] >> p2[1] >> p2[2];
              // cell polarization matrix
              bdm::Double3x3 pl_3x3 = tensor(p0,p0)+tensor(p1,p1)+tensor(p2,p2);
              // cell cycle phase
              int ccp;
              fin >> ccp;
              // important checks before uploading this cell
              if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol )
                continue;
              if ( xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol )
                continue;
              if ( xyz[2] < minCOORD+tol || xyz[2] > maxCOORD-tol )
                continue;
              // check if a cell should follow any user-defined (spatial)
              // restrictions (wrt primitive shapes: box, sphere, cylinder)
              bool is_valid = true;
              if      ( params.get<bool>(CP_name+"/initial_population/pattern/box/inside") )
                {
                  is_valid = true;
                  if ( xyz[0]<pntA[0] || xyz[0]>pntB[0] ) is_valid = false;
                  if ( xyz[1]<pntA[1] || xyz[1]>pntB[1] ) is_valid = false;
                  if ( xyz[2]<pntA[2] || xyz[2]>pntB[2] ) is_valid = false;
                }
              else if ( params.get<bool>(CP_name+"/initial_population/pattern/box/outside") )
                {
                  is_valid = false;
                  if ( xyz[0]<pntA[0] || xyz[0]>pntB[0] ) is_valid = true;
                  if ( xyz[1]<pntA[1] || xyz[1]>pntB[1] ) is_valid = true;
                  if ( xyz[2]<pntA[2] || xyz[2]>pntB[2] ) is_valid = true;
                }
              else if ( params.get<bool>(CP_name+"/initial_population/pattern/sphere/inside") )
                {
                  const bdm::Double3 diff = xyz - cntr;
                  is_valid = true;
                  if ( L2norm(diff)>rad ) is_valid = false;
                }
              else if ( params.get<bool>(CP_name+"/initial_population/pattern/sphere/outside") )
                {
                  const bdm::Double3 diff = xyz - cntr;
                  is_valid = false;
                  if ( L2norm(diff)>rad ) is_valid = true;
                }
              //
              if ( ! is_valid ) continue;
              ////++number_of_cells_loaded;
              //
              bdm::BiologicalCell* cell = new bdm::BiologicalCell(CP_ID, xyz);
              cell->SetPhase(ccp);
              cell->SetDiameter(dia);
              cell->SetAdherence(params.have_parameter<double>(CP_name+"/adherence")
                               ? params.get<double>(CP_name+"/adherence") : 0.0);
              if (params.have_parameter<double>(CP_name+"/density"))
                cell->SetDensity(params.get<double>(CP_name+"/density"));
              cell->SetParametersPointer(&params);
              cell->SetCanApoptose(params.get<bool>(CP_name+"/can_apoptose"));
              cell->SetCanGrow(params.get<bool>(CP_name+"/can_grow"));
              cell->SetCanDivide(params.get<bool>(CP_name+"/can_divide"));
              cell->SetCanMigrate(params.get<bool>(CP_name+"/can_migrate"));
              cell->SetCanTransform(params.get<bool>(CP_name+"/can_transform"));
              cell->SetCanPolarize(params.get<bool>(CP_name+"/can_polarize"));
              cell->SetCanProtrude(params.get<bool>(CP_name+"/can_protrude"));
              cell->SetPolarization(pl_3x3);
              if      (10==mo)
                cell->AddBehavior(new bdm::Biology4BiologicalCell_10());
              else if (11==mo)
                cell->AddBehavior(new bdm::Biology4BiologicalCell_11());
              else
                ABORT_("\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has unrecognized behavior");
              // store this cell into BioDynaMo's resource manager
              rm->AddAgent(cell);
            }
          // create a copy of the file just processed
          const std::string cmd = "cp " + fn + "  "
                                + params.get<std::string>("output_directory")
                                + "/in/" + CP_name + ".0.dat";
          ASSERT_(0==std::system(cmd.c_str()),
                  "could not save a copy of a data file");
          // ...end of this if-case
        }
      // update the map containing the ID for this phenotype
      params.set<std::string>("phenotype_ID/"+std::to_string(CP_ID)) = CP_name;
      // ...end of cell phenotypes loop
    }
}
// =============================================================================
inline
void init_vessels(bdm::Simulation& sim,
                  const std::vector<std::string>& biochem)
{
  if (!params.get<bool>("simulation_models_vessels")) return;
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  //
  // default parameter(s) value
  if (! params.have_parameter<double>("vessel/can_grow/probability"))
    params.set<double>("vessel/can_grow/probability") = 0.0;
  //
  bdm::Vessel v;
  //
  // iterate for all biochemicals (substances)
  for (unsigned int icue=0; icue<biochem.size(); icue++)
    {
      const std::string& BC_name = biochem[icue];
      //
      // default parameter(s) value
      if (! params.have_parameter<double>("vessel/can_branch/"+BC_name+"/probability"))
        params.set<double>("vessel/can_branch/"+BC_name+"/probability") = 0.0;
      if (! params.have_parameter<double>("vessel/can_sprout/"+BC_name+"/probability"))
        params.set<double>("vessel/can_sprout/"+BC_name+"/probability") = 0.0;
      //
      if (! params.have_parameter<double>("vessel/"+BC_name+"/secretion/net_balance"))
        {
          params.set<double>("vessel/"+BC_name+"/secretion/net_balance") = 0.0;
          params.set<double>("vessel/"+BC_name+"/secretion/net_balance/std") = 0.0;
          params.set<double>("vessel/"+BC_name+"/secretion/saturation") = 0.0;
        }
      else
        {
          if (! params.have_parameter<double>("vessel/"+BC_name+"/secretion/net_balance/std"))
            params.set<double>("vessel/"+BC_name+"/secretion/net_balance/std") = 0.0;
          if (! params.have_parameter<bool>("vessel/"+BC_name+"/secretion/dependency"))
            params.set<bool>("vessel/"+BC_name+"/secretion/dependency") = true;
          if (! params.have_parameter<double>("vessel/"+BC_name+"/secretion/saturation"))
            params.set<double>("vessel/"+BC_name+"/secretion/saturation") = 0.0;
        }
      //
      // --- Michaelis-Menten kinetics model (optional for vessels) ---
      // Self-contained model: R = Vmax * C / (Km + C)
      // Sign of Vmax determines direction: Vmax < 0 → consumption, Vmax > 0 → production.
      if (! params.have_parameter<bool>("vessel/"+BC_name+"/secretion/michaelis_menten_model"))
        params.set<bool>("vessel/"+BC_name+"/secretion/michaelis_menten_model") = false;
      if (params.get<bool>("vessel/"+BC_name+"/secretion/michaelis_menten_model"))
        {
          if (! params.have_parameter<double>("vessel/"+BC_name+"/secretion/michaelis_menten/Vmax"))
            {
              params.set<double>("vessel/"+BC_name+"/secretion/michaelis_menten/Vmax") = -0.001;
              bdm::Log::Warning("ABM4bio",
                "Michaelis-Menten Vmax not provided for vessel/"+BC_name
                +"; defaulting to -0.001 (consumption)");
            }
          if (! params.have_parameter<double>("vessel/"+BC_name+"/secretion/michaelis_menten/Km"))
            {
              params.set<double>("vessel/"+BC_name+"/secretion/michaelis_menten/Km") = 1.0;
              bdm::Log::Warning("ABM4bio",
                "Michaelis-Menten Km not provided for vessel/"+BC_name
                +"; defaulting to 1.0");
            }
          //
          bdm::Log::Info("ABM4bio",
            "Michaelis-Menten model ENABLED for vessel/"+BC_name
            +" (Vmax="+std::to_string(params.get<double>("vessel/"+BC_name+"/secretion/michaelis_menten/Vmax"))
            +", Km="+std::to_string(params.get<double>("vessel/"+BC_name+"/secretion/michaelis_menten/Km"))+")");
        }
      // ...end of biochemicals (substances) loop
    }
  //
  if (! params.have_parameter<std::string>("vessel/initial_configuration/from_file"))
    {
      ABORT_("initial configuration of vessels should be read from file");
    }
  else
    {
      const std::string fn = params.get<std::string>("vessel/initial_configuration/from_file");
      //
      std::ifstream fin(fn);
      ASSERT_(fin.good(),"file \""+fn+"\" cannot be accessed");
      //
      int n_vessel = 0;
      fin >> n_vessel;
      //
      for (int ivessel=0; ivessel<n_vessel; ivessel++)
        {
          // set the vessel ID
          const int ID = ivessel;
          // read the spatial coordinates (starting point) of the vessel
          bdm::Double3 start;
          fin >> start[0] >> start[1] >> start[2];
          // read the age of this vessel
          int age;
          fin >> age;
          // save spatial coordinates in the local container
          all_agents.push_back( start );
          //
          vessel_map__ID_age.insert( std::make_pair(ID, age) );
          //
          auto* ns = new bdm::neuroscience::NeuronSoma(start);
          ns->SetAdherence(params.get<double>("default_vessel_adherence"));
          rm->AddAgent(ns);
          //
          // read the number of segments for this vessel
          int n_segm;
          fin >> n_segm;
          // sanity check
          ASSERT_(n_segm>=2,
                  "vessel ID loaded from from file is invalid");
          //
          double dia;
          bdm::Double3 end, axis;
          bool can_grow, can_branch, can_sprout;
          //
          // read the diameter for this vessel - for the first vessel!
          fin >> dia;
          // read the spatial coordinates (mass / end point) for this vessel
          fin >> end[0] >> end[1] >> end[2];
          // calculate the axis of the vessel
          axis = end - start;
          // read vessel features
          fin >> can_grow >> can_branch >> can_sprout;
          // save spatial coordinates in the local container
          all_agents.push_back( (start+end)*0.5 );
          all_agents.push_back( end );
          //
          // sanity check
          if (can_branch && can_sprout)
            ABORT_("vessel cannot branch and sprout simultaneously");
          //
          auto* vessel = bdm::bdm_static_cast<bdm::Vessel*>(ns->ExtendNewNeurite(axis, &v));
          vessel->SetVesselID(ID);
          vessel->SetDiameter(dia);
          vessel->SetAge(age);
          vessel->SetCanGrow(can_grow);
          vessel->SetCanBranch(can_branch);
          vessel->SetCanSprout(can_sprout);
          vessel->SetParametersPointer(&params);
          vessel->AddBehavior(new bdm::Biology4Vessel());
          // ...reset this for the next vessel segment
          start = end;
          // iterate for the rest of the vessel segments
          for (int isegm=1; isegm<n_segm; isegm++)
            {
              // read the diameter for this vessel
              fin >> dia;
              // read the spatial coordinates (mass / end point) for this vessel
              fin >> end[0] >> end[1] >> end[2];
              // calculate the axis of the vessel
              axis = end - start;
              // read vessel features
              fin >> can_grow >> can_branch >> can_sprout;
              // save spatial coordinates in the local container
              all_agents.push_back( (start+end)*0.5 );
              all_agents.push_back( end );
              //
              const double extend_rate = L2norm(axis)
                                       / params.get<double>("time_step");
              //
              vessel->ElongateTerminalEnd(extend_rate, normalize(axis));
              vessel->RunDiscretization();
              // ...reset this for the next vessel segment
              start = end;
            }
          // ...end of vessel loop
        }
      // initialize this counter
      params.set<int>("N_vessel_IDs") = n_vessel;
      // create a copy of the file just processed
      const std::string cmd = "cp " + fn + "  "
                            + params.get<std::string>("output_directory")
                            + "/in/vessels.0.dat";
      ASSERT_(0==std::system(cmd.c_str()),
              "could not save a copy of a data file");
      // ...end of this if-case
    }
}
// =============================================================================
inline
void reinit_cells(bdm::Simulation& sim,
                  const std::map<int,std::string>& cells)
{
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = sim.GetRandom();
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = params.get<double>("min_boundary"),
               maxCOORD = params.get<double>("max_boundary"),
               meanCOORD = (maxCOORD+minCOORD)/2.0;
  const double tol = params.get<double>("domain_tolerance");
  //
  // iterate for all cell phenotypes
  for ( std::map<int, std::string>::const_iterator
        ci=cells.begin(); ci!=cells.end(); ci++ )
    {
      // obtain the ID & name for this cell phenotype
      const int CP_ID = ci->first;
      const std::string& CP_name = ci->second;
      //
      // reset the spherical coordinates for all escaping cells:
      // radius, inclination (theta) and azimuth (phi)
      std::vector<double>& escaped_cells_radius =
        params.set<std::vector<double>>(CP_name+"/escaped_cells/radius");
      std::vector<double>& escaped_cells_theta  =
        params.set<std::vector<double>>(CP_name+"/escaped_cells/theta");
      std::vector<double>& escaped_cells_phi    =
        params.set<std::vector<double>>(CP_name+"/escaped_cells/phi");
      //
      std::vector<int>& escaped_cells_phase =
        params.set<std::vector<int>>(CP_name+"/escaped_cells/phase");
      //
      const unsigned int n_escaped_cells = escaped_cells_phase.size();
      //
      // by design we ignore any escaping cells that are necrotic!!!
      if ( 0 == CP_ID )
        {
          // enforce to clear memory
          escaped_cells_radius.clear();
          escaped_cells_theta.clear();
          escaped_cells_phi.clear();
          escaped_cells_phase.clear();
          // ...and then ignore the following computations
          continue;
        }
      //
      if ( n_escaped_cells!=escaped_cells_radius.size() ||
           n_escaped_cells!=escaped_cells_theta.size() ||
           n_escaped_cells!=escaped_cells_phi.size() ) continue;
      //
      if ( 0 == n_escaped_cells )
        {
          // enforce to clear memory
          escaped_cells_radius.clear();
          escaped_cells_theta.clear();
          escaped_cells_phi.clear();
          // ...and then ignore the following computations
          continue;
        }
      // order of mechanisms that define behaviour
      const int mo = params.get<int>(CP_name+"/mechanism_order");
      //
      if ( ! params.get<bool>("simulation_domain_is_bounded") &&
             params.get<bool>("simulation_domain_is_periodic") )
        {
          const bool antisymmetry =
            params.get<double>("simulation_domain_is_periodic/antisymmetry");
          //
          for (unsigned int icell=0; icell<n_escaped_cells; icell++)
            {
              double radius, theta, phi;
              if (antisymmetry)
                {
                  double angle_theta = 0.0;
                  if (! params.get<bool>("simulation_domain_is_2D"))
                    angle_theta = (uniform_distro(-2,4)>0.0 ? 0.5*bdm::Math::kPi : 0.0);
                  double angle_phi = (uniform_distro(-2,4)>0.0 ? bdm::Math::kPi : 0.0);
                  // spherical and cartesian coordinates respectively
                  radius = escaped_cells_radius[icell] - rg->Uniform(0.0, tol),
                  theta  = escaped_cells_theta[icell] + angle_theta,
                  phi    = escaped_cells_phi[icell] + angle_phi;
                }
              else
                {
                  // spherical and cartesian coordinates respectively
                  radius = escaped_cells_radius[icell] - rg->Uniform(0.0, tol),
                  theta  = escaped_cells_theta[icell],
                  phi    = escaped_cells_phi[icell];
                }
              const double x = meanCOORD + radius * sin(theta) * cos(phi),
                           y = meanCOORD + radius * sin(theta) * sin(phi),
                           z = meanCOORD + radius * cos(theta);
              // cell position
              const bdm::Double3 xyz = {x, y, z};
              // cell cycle phase
              const int ccp = escaped_cells_phase[icell];
              // cell diameter
              const double dia = rg->Uniform(params.get<double>(CP_name+"/diameter/min"),
                                             params.get<double>(CP_name+"/diameter/max"));
              // principal directions of the cell polarization matrix
              const double pd0 = params.get<double>(CP_name+"/principal/0"),
                           pd1 = params.get<double>(CP_name+"/principal/1"),
                           pd2 = params.get<double>(CP_name+"/principal/2");
              // cell polarization matrix
              bdm::Double3x3 pl_3x3 = diag(pd0, pd1, pd2);
              //
              bdm::BiologicalCell* cell = new bdm::BiologicalCell(CP_ID, xyz);
              cell->SetPhase(ccp);
              cell->SetDiameter(dia);
              cell->SetAdherence(params.have_parameter<double>(CP_name+"/adherence")
                               ? params.get<double>(CP_name+"/adherence") : 0.0);
              if (params.have_parameter<double>(CP_name+"/density"))
                cell->SetDensity(params.get<double>(CP_name+"/density"));
              cell->SetParametersPointer(&params);
              cell->SetCanApoptose(params.get<bool>(CP_name+"/can_apoptose"));
              cell->SetCanGrow(params.get<bool>(CP_name+"/can_grow"));
              cell->SetCanDivide(params.get<bool>(CP_name+"/can_divide"));
              cell->SetCanMigrate(params.get<bool>(CP_name+"/can_migrate"));
              cell->SetCanTransform(params.get<bool>(CP_name+"/can_transform"));
              cell->SetCanPolarize(params.get<bool>(CP_name+"/can_polarize"));
              cell->SetCanProtrude(params.get<bool>(CP_name+"/can_protrude"));
              cell->SetPolarization(pl_3x3);
              if      (10==mo)
                cell->AddBehavior(new bdm::Biology4BiologicalCell_10());
              else if (11==mo)
                cell->AddBehavior(new bdm::Biology4BiologicalCell_11());
              else
                ABORT_("\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has unrecognized behavior");
              // store this cell into BioDynaMo's resource manager
              rm->AddAgent(cell);
            }
          // ...finished inserting incoming cells due to periodic boundary
        }
      //
      escaped_cells_radius.clear();
      escaped_cells_theta.clear();
      escaped_cells_phi.clear();
      escaped_cells_phase.clear();
      // ...end of cell phenotypes loop
    }
}
// =============================================================================
inline
void save_snapshot(bdm::Simulation& sim, const int time = 0)
{
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  //
  const std::string name = params.get<std::string>("output_directory")
                         + "/" + params.get<std::string>("simulation_title");
  const int n_time = params.get<int>("number_of_time_steps"),
            viz_step = params.get<int>("visualization_interval");
  const double time_step = params.get<double>("time_step");
  const std::vector<std::string>& substances =
    params.get<std::vector<std::string>>("substances");
  const bool simplify_output = params.get<bool>("simplify_output");
  // sanity check
  if (substances.empty())
    ABORT_("list of biochemicals is empty or uninitialiazed");
  //
  // create a PVD file to encompass all time-steps
  // for the cells
  if ( 0 == time )
    {
      std::ofstream fout(name+".C.pvd");
      //
      fout << "<?xml version=\"1.0\"?>" << std::endl;
      fout << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fout << "<Collection>" << std::endl;
      for (int t=0; t<=n_time; t++)
        {
          if (0!=t%viz_step) continue;
          const double timestep = t * time_step;
          const std::string file = params.get<std::string>("simulation_title")
                                 + "/cells." + std::to_string(t) + ".vtu";
          fout << "  <DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\" file=\"" << file << "\"/>" << std::endl;
        }
      fout << "</Collection>" << std::endl;
      fout << "</VTKFile>" << std::endl;
      //... end of PVD file
    }
  // create a VTU file for this time-step
  // for the cells
  ////if ( 0 != time )
    {
      unsigned int n_VTK_points = 0;
      unsigned int n_VTK_cells = 1;
      //
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            ++n_VTK_points;
          }
      });
      //
      std::ofstream fout(name+"/cells."+std::to_string(time)+".vtu");
      //
      fout << "<?xml version=\"1.0\"?>" << std::endl;
      fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fout << "  <UnstructuredGrid>" << std::endl;
      fout << "    <Piece NumberOfPoints=\"" << n_VTK_points << "\" NumberOfCells=\"" << n_VTK_cells << "\">" << std::endl;
      fout << "      <Points>" << std::endl;
      fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            const bdm::Double3& vec = cell->GetPosition();
            fout << ' ' << vec[0] << ' ' << vec[1] << ' ' << vec[2];
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </Points>" << std::endl;
      fout << "      <PointData>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"phenotype\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetPhenotype();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"phase\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetPhase();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"age\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetAge();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"diameter\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetDiameter();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"volume\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetVolume();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"trail\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetTrail();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetDisplacement(0) << ' ' << cell->GetDisplacement(1) << ' ' << cell->GetDisplacement(2);
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      //++++++++++++++++++++++++
      if ( ! simplify_output ) {
      //++++++++++++++++++++++++
      fout << "        <DataArray type=\"Int32\" Name=\"can_apoptose\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << static_cast<int>(cell->GetCanApoptose());
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_grow\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << static_cast<int>(cell->GetCanGrow());
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_divide\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << static_cast<int>(cell->GetCanDivide());
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_migrate\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << static_cast<int>(cell->GetCanMigrate());
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_transform\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << static_cast<int>(cell->GetCanTransform());
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_polarize\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << static_cast<int>(cell->GetCanPolarize());
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_protrude\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << static_cast<int>(cell->GetCanProtrude());
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      //++++++++++++++++++++++++
      } // ...end if-statement
      //++++++++++++++++++++++++
      fout << "        <DataArray type=\"Int32\" Name=\"n_divisions\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetNumberOfDivisions();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"n_trasformations\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetNumberOfTrasformations();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"n_protrusions\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetNumberOfProtrusions();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      //++++++++++++++++++++++++
      if ( ! simplify_output ) {
      //++++++++++++++++++++++++
      fout << "        <DataArray type=\"Float64\" Name=\"polarize\" NumberOfComponents=\"9\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            fout << ' ' << cell->GetPolarization(0,0) << ' ' << cell->GetPolarization(0,1) << ' ' << cell->GetPolarization(0,2);
            fout << ' ' << cell->GetPolarization(1,0) << ' ' << cell->GetPolarization(1,1) << ' ' << cell->GetPolarization(1,2);
            fout << ' ' << cell->GetPolarization(2,0) << ' ' << cell->GetPolarization(2,1) << ' ' << cell->GetPolarization(2,2);
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      //++++++++++++++++++++++++
      } // ...end if-statement
      //++++++++++++++++++++++++
      fout << "        <DataArray type=\"Float64\" Name=\"shape\" NumberOfComponents=\"9\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
          {
            const double& r = 0.5*cell->GetDiameter();
            const bdm::Double3x3& pm = cell->GetPolarization();
            fout << ' ' << r*pm[0][0] << ' ' << r*pm[0][1] << ' ' << r*pm[0][2];
            fout << ' ' << r*pm[1][0] << ' ' << r*pm[1][1] << ' ' << r*pm[1][2];
            fout << ' ' << r*pm[2][0] << ' ' << r*pm[2][1] << ' ' << r*pm[2][2];
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </PointData>" << std::endl;
      fout << "      <Cells>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      int offset = 0;
      for (unsigned int ipnt=0; ipnt<n_VTK_points; ipnt++)
        {
          // VTK_POLY_VERTEX
          ++offset;
        }
      fout << ' ' << offset;
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      for (unsigned int ipnt=0; ipnt<n_VTK_points; ipnt++)
        {
          // VTK_POLY_VERTEX
          fout << ' ' << ipnt;
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
        {
          // VTK_POLY_VERTEX
          fout << ' ' << 2;
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </Cells>" << std::endl;
      fout << "    </Piece>" << std::endl;
      fout << "  </UnstructuredGrid>" << std::endl;
      fout << "</VTKFile>" << std::endl;
      //... end of VTU file
    }
  //
  // create a PVD file to encompass all time-steps
  // for the cell protrusions
  if ( 0 == time )
    {
      std::ofstream fout(name+".CP.pvd");
      //
      fout << "<?xml version=\"1.0\"?>" << std::endl;
      fout << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fout << "<Collection>" << std::endl;
      for (int t=0; t<=n_time; t++)
        {
          if (0!=t%viz_step) continue;
          const double timestep = t * time_step;
          const std::string file = params.get<std::string>("simulation_title")
                                 + "/cell_protrusions." + std::to_string(t) + ".vtu";
          fout << "  <DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\" file=\"" << file << "\"/>" << std::endl;
        }
      fout << "</Collection>" << std::endl;
      fout << "</VTKFile>" << std::endl;
      //... end of PVD file
    }
  // create a VTU file for this time-step
  // for the cell protrusions
  ////if ( 0 != time )
    {
      unsigned int n_VTK_points = 0;
      unsigned int n_VTK_cells = 0;
      //
      // connectivity for all cell protrusions (line segments)
      std::vector<std::pair<int,int>> conn;
      // counter of protrusion nodes
      int node = 0;
      //
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
          if (nullptr!=protrusion->GetCell())
            {
              n_VTK_points += 2;
              n_VTK_cells += 1;
              //
              std::pair<int,int> nodes(node+0, node+1);
              conn.push_back(nodes);
              //
              node += 2;
            }
      });
      // sanity check
      if ( conn.size() != n_VTK_cells )
        ABORT_("unexpected initialiazation error occurred");
      //
      bool set_a_dummy_protrusion = false;
      if ( ! n_VTK_points )
        { // set a dummy protrusion
          n_VTK_points += 2;
          n_VTK_cells += 1;
          //
          std::pair<int,int> nodes(0, 1);
          conn.push_back(nodes);
          //
          set_a_dummy_protrusion = true;
        }
      //
      std::ofstream fout(name+"/cell_protrusions."+std::to_string(time)+".vtu");
      //
      fout << "<?xml version=\"1.0\"?>" << std::endl;
      fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fout << "  <UnstructuredGrid>" << std::endl;
      fout << "    <Piece NumberOfPoints=\"" << n_VTK_points << "\" NumberOfCells=\"" << n_VTK_cells << "\">" << std::endl;
      fout << "      <Points>" << std::endl;
      fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
      if ( set_a_dummy_protrusion )
        {
          const double a = params.get<double>("min_boundary") / 1.0e+4,
                       b = params.get<double>("max_boundary") / 1.0e+4;
          fout << ' ' << a << ' ' << a << ' ' << a;
          fout << ' ' << b << ' ' << b << ' ' << b;
        }
      else
      rm->ForEachAgent([&] (bdm::Agent* a) {
        bdm::Double3 node;
        if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
          if (nullptr!=protrusion->GetCell())
            {
              const bdm::Double3 p = protrusion->GetPosition(),
                                 n = protrusion->GetSpringAxis();
              node = p - n * 0.5;
              fout << ' ' << node[0] << ' ' << node[1] << ' ' << node[2];
              node = p + n * 0.5;
              fout << ' ' << node[0] << ' ' << node[1] << ' ' << node[2];
            }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </Points>" << std::endl;
      fout << "      <Cells>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      int offset = 0;
      for (int s=0; s<n_VTK_cells; s++)
        {
          // VTK_LINE
          offset += 2;
          fout << ' ' << offset;
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      for (int s=0; s<n_VTK_cells; s++)
        {
          const std::pair<int,int>& line = conn[s];
          // VTK_LINE
          fout << ' ' << line.first << ' ' << line.second;
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      for (int s=0; s<n_VTK_cells; s++)
        {
          // VTK_LINE
          fout << ' ' << 3;
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </Cells>" << std::endl;
      fout << "      <CellData>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"age\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      if ( set_a_dummy_protrusion )
        {
          fout << ' ' << 0;
        }
      else
        rm->ForEachAgent([&] (bdm::Agent* a) {
          if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
            if (nullptr!=protrusion->GetCell())
              {
                fout << ' ' << protrusion->GetAge();
              }
        });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      if ( set_a_dummy_protrusion )
        {
          fout << ' ' << 0.0;
        }
      else
        rm->ForEachAgent([&] (bdm::Agent* a) {
          if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
            if (nullptr!=protrusion->GetCell())
            {
              fout << ' ' << 0.5*protrusion->GetDiameter();
            }
        });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"length\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      if ( set_a_dummy_protrusion )
        {
          fout << ' ' << 0.0;
        }
      else
        rm->ForEachAgent([&] (bdm::Agent* a) {
          if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
            if (nullptr!=protrusion->GetCell())
            {
              fout << ' ' << protrusion->GetActualLength();
            }
        });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"length2branch\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      if ( set_a_dummy_protrusion )
        {
          fout << ' ' << 0.0;
        }
      else
        rm->ForEachAgent([&] (bdm::Agent* a) {
          if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
            if (nullptr!=protrusion->GetCell())
            {
              fout << ' ' << protrusion->LengthToProximalBranchingPoint();
            }
        });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"phenotype\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      if ( set_a_dummy_protrusion )
        {
          fout << ' ' << 0;
        }
      else
        rm->ForEachAgent([&] (bdm::Agent* a) {
          if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
            if (nullptr!=protrusion->GetCell())
            {
              fout << ' ' << protrusion->GetCell()->GetPhenotype();
            }
        });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_sprout\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      if ( set_a_dummy_protrusion )
        {
          fout << ' ' << 0;
        }
      else
        rm->ForEachAgent([&] (bdm::Agent* a) {
          if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
            if (nullptr!=protrusion->GetCell())
            {
              fout << ' ' << protrusion->GetCanSprout();
            }
        });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_branch\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      if ( set_a_dummy_protrusion )
        {
          fout << ' ' << 0;
        }
      else
        rm->ForEachAgent([&] (bdm::Agent* a) {
          if (auto* protrusion = dynamic_cast<bdm::CellProtrusion*>(a))
            if (nullptr!=protrusion->GetCell())
            {
              fout << ' ' << protrusion->GetCanBranch();
            }
        });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </CellData>" << std::endl;
      fout << "    </Piece>" << std::endl;
      fout << "  </UnstructuredGrid>" << std::endl;
      fout << "</VTKFile>" << std::endl;
    }
  //
  // create a PVD file to encompass all time-steps
  // for the reaction-diffusion simulator
  if ( 0 == time )
    {
      std::ofstream fout(name+".DG.pvd");
      //
      fout << "<?xml version=\"1.0\"?>" << std::endl;
      fout << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fout << "<Collection>" << std::endl;
      for (int t=0; t<=n_time; t++)
        {
          if (0!=t%viz_step) continue;
          const double timestep = t * time_step;
          const std::string file = params.get<std::string>("simulation_title")
                                 + "/diffusion_grid." + std::to_string(t) + ".vtu";
          fout << "  <DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\" file=\"" << file << "\"/>" << std::endl;
        }
      fout << "</Collection>" << std::endl;
      fout << "</VTKFile>" << std::endl;
      //... end of PVD file
    }
  // create a VTU file for this time-step
  // for the reaction-diffusion simulator
  ////if ( 0 != time )
    {
      const int N = params.get<int>("diffusion_grid/spatial_resolution");
      const double S_min = params.get<double>("min_boundary"),
                   S_max = params.get<double>("max_boundary"),
                   DS = (S_max-S_min) / N;
      //
      // space vectors corresponding to the diffusion grid
      // the size of which in principle is the same for all substances
      // NOTE: query points are placed at voxel CENTERS (edge + DS/2)
      //       to match BioDynaMo's floor-based GetBoxIndex mapping
      if ( dg_vec.empty() )
        {
          // check first for type of simulation domain
          if (params.get<bool>("simulation_domain_is_2D"))
            {
              // iterate for all points of the Cartesian grid
              for (int K=N/2-1; K<=N/2+1; K++)
                {
                  const double Z = S_min + DS * K + DS * 0.5;
                  for (int J=0; J<N; J++)
                    {
                      const double Y = S_min + DS * J + DS * 0.5;
                      for (int I=0; I<N; I++)
                        {
                          const double X = S_min + DS * I + DS * 0.5;
                          // upload this point into container
                          dg_vec.push_back( {X, Y, Z} );
                        }
                    }
                }
            }
          else
            {
              // iterate for all points of the Cartesian grid
              for (int K=0; K<N; K++)
                {
                  const double Z = S_min + DS * K + DS * 0.5;
                  for (int J=0; J<N; J++)
                    {
                      const double Y = S_min + DS * J + DS * 0.5;
                      for (int I=0; I<N; I++)
                        {
                          const double X = S_min + DS * I + DS * 0.5;
                          // upload this point into container
                          dg_vec.push_back( {X, Y, Z} );
                        }
                    }
                }
            }
          // ...initialization is complete
        }
      //
      const int n_points = dg_vec.size();
      //
      std::ofstream fout(name+"/diffusion_grid."+std::to_string(time)+".vtu");
      //
      fout << "<?xml version=\"1.0\"?>" << std::endl;
      fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fout << "  <UnstructuredGrid>" << std::endl;
      fout << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells=\"" << 1 << "\">" << std::endl;
      fout << "      <Points>" << std::endl;
      fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
      for (int P=0; P<n_points; P++)
        {
          const bdm::Double3& vec = dg_vec[P];
          fout << ' ' << vec[0] << ' ' << vec[1] << ' ' << vec[2];
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </Points>" << std::endl;
      fout << "      <PointData>" << std::endl;
      for ( std::vector<std::string>::const_iterator
            ci=substances.begin(); ci!=substances.end(); ci++ )
        {
          // access the BioDynaMo diffusion grid
          auto* dg = rm->GetDiffusionGrid(*ci);
          //
          fout << "        <DataArray type=\"Float64\" Name=\""+dg->GetContinuumName()+"\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
          for (int P=0; P<n_points; P++)
            {
              const double& biochem = dg->GetValue(dg_vec[P]);
              fout << ' ' << biochem;
            }
          fout << std::endl
               << "        </DataArray>" << std::endl;
          //
          fout << "        <DataArray type=\"Float64\" Name=\"Grad_"+dg->GetContinuumName()+"\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
          for (int P=0; P<n_points; P++)
            {
              bdm::Double3 biochem_grad;
              dg->GetGradient(dg_vec[P], &biochem_grad);
              fout << ' ' << biochem_grad[0] << ' ' << biochem_grad[1] << ' ' << biochem_grad[2];
            }
          fout << std::endl
               << "        </DataArray>" << std::endl;
          //...end loop for this substance
        }
      if ( params.have_parameter<std::string>("convection/dynamic/from_file") )
      {
        const int SpaceDimension = params.get<bool>("simulation_domain_is_2D") ? 2 : 3;
        //
        fout << "        <DataArray type=\"Float64\" Name=\"convection\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
        if (2==SpaceDimension)
          {
            for (int P=0; P<n_points; P++)
              {
                for (int ispdm=0; ispdm<SpaceDimension; ispdm++)
                  {
                    const std::string name = "convection_" + std::to_string(ispdm);
                    // access the BioDynaMo diffusion grid
                    auto* dg = rm->GetDiffusionGrid(name);
                    // print out the convection component
                    fout << ' ' << dg->GetValue(dg_vec[P]);
                  }
                {
                  // print out the convection component
                  fout << ' ' << 0.0;
                }
              }
          }
        else
          {
            for (int P=0; P<n_points; P++)
              {
                for (int ispdm=0; ispdm<SpaceDimension; ispdm++)
                  {
                    const std::string name = "convection_" + std::to_string(ispdm);
                    // access the BioDynaMo diffusion grid
                    auto* dg = rm->GetDiffusionGrid(name);
                    // print out the convection component
                    fout << ' ' << dg->GetValue(dg_vec[P]);
                  }
              }
          }
        fout << std::endl
             << "        </DataArray>" << std::endl;
      }
      fout << "      </PointData>" << std::endl;
      fout << "      <Cells>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      int offset = 0;
      // VTK_POLY_VERTEX
      offset += n_points;
      fout << ' ' << offset;
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      for (int P=0; P<n_points; P++) {
        // VTK_POLY_VERTEX
        fout << ' ' << P;
      }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      // VTK_POLY_VERTEX
      fout << ' ' << 2;
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </Cells>" << std::endl;
      fout << "    </Piece>" << std::endl;
      fout << "  </UnstructuredGrid>" << std::endl;
      fout << "</VTKFile>" << std::endl;
      //... end of VTU file
    }
  //
  // do subsequent printout if vessels are simulated
  if ( params.get<bool>("simulation_models_vessels") )
  // create a PVD file to encompass all time-steps
  // for the vessels
  if ( 0 == time )
    {
      std::ofstream fout(name+".V.pvd");
      //
      fout << "<?xml version=\"1.0\"?>" << std::endl;
      fout << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fout << "<Collection>" << std::endl;
      for (int t=0; t<=n_time; t++)
        {
          if (0!=t%viz_step) continue;
          const double timestep = t * time_step;
          const std::string file = params.get<std::string>("simulation_title")
                                 + "/vessels." + std::to_string(t) + ".vtu";
          fout << "  <DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\" file=\"" << file << "\"/>" << std::endl;
        }
      fout << "</Collection>" << std::endl;
      fout << "</VTKFile>" << std::endl;
    }
  //
  // do subsequent printout if vessels are simulated
  if ( params.get<bool>("simulation_models_vessels") )
  // create a VTU file for this time-step
  // for the vessels
  ////if ( 0 != time )
    {
      unsigned int n_VTK_points = 0;
      unsigned int n_VTK_cells = 0;
      //
      // connectivity for all vessels (line segments)
      std::vector<std::pair<int,int>> conn;
      // counter of vascular nodes
      int node = 0;
      //
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            n_VTK_points += 2;
            n_VTK_cells += 1;
            //
            std::pair<int,int> nodes(node+0, node+1);
            conn.push_back(nodes);
            //
            node += 2;
          }
      });
      // sanity check
      if ( conn.size() != n_VTK_cells )
        ABORT_("unexpected initialiazation error occurred");
      //
      std::ofstream fout(name+"/vessels."+std::to_string(time)+".vtu");
      //
      fout << "<?xml version=\"1.0\"?>" << std::endl;
      fout << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fout << "  <UnstructuredGrid>" << std::endl;
      fout << "    <Piece NumberOfPoints=\"" << n_VTK_points << "\" NumberOfCells=\"" << n_VTK_cells << "\">" << std::endl;
      fout << "      <Points>" << std::endl;
      fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        bdm::Double3 node;
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            const bdm::Double3 p = vessel->GetPosition(),
                               n = vessel->GetSpringAxis();
            node = p - n * 0.5;
            fout << ' ' << node[0] << ' ' << node[1] << ' ' << node[2];
            node = p + n * 0.5;
            fout << ' ' << node[0] << ' ' << node[1] << ' ' << node[2];
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </Points>" << std::endl;
      fout << "      <Cells>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      int offset = 0;
      for (int s=0; s<n_VTK_cells; s++)
        {
          // VTK_LINE
          offset += 2;
          fout << ' ' << offset;
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      for (int s=0; s<n_VTK_cells; s++)
        {
          const std::pair<int,int>& line = conn[s];
          // VTK_LINE
          fout << ' ' << line.first << ' ' << line.second;
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      for (int s=0; s<n_VTK_cells; s++)
        {
          // VTK_LINE
          fout << ' ' << 3;
        }
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </Cells>" << std::endl;
      fout << "      <CellData>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"age\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->GetAge();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << 0.5*vessel->GetDiameter();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"length\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->GetActualLength();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"length2branch\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->LengthToProximalBranchingPoint();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Float64\" Name=\"volume\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->GetVolume();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"vessel_ID\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->GetVesselID();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_grow\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->GetCanGrow();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_branch\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->GetCanBranch();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"can_sprout\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->GetCanSprout();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "        <DataArray type=\"Int32\" Name=\"is_terminal\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
      rm->ForEachAgent([&] (bdm::Agent* a) {
        if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
          {
            fout << ' ' << vessel->IsTerminal();
          }
      });
      fout << std::endl
           << "        </DataArray>" << std::endl;
      fout << "      </CellData>" << std::endl;
      fout << "    </Piece>" << std::endl;
      fout << "  </UnstructuredGrid>" << std::endl;
      fout << "</VTKFile>" << std::endl;
    }
}
// =============================================================================
inline
void reinit_biochemicals(bdm::Simulation& sim,
                         const std::vector<std::string>& biochem, const int time = 0)
{
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  //
  // iterate for all biochemicals (cues)
  for (unsigned int icue=0; icue<biochem.size(); icue++)
    {
      const std::string& BC_name = biochem[icue];
      const std::string Time = std::to_string(time);
      //
      // skip following steps if no dynamic initialization for this biochemical
      if (! params.have_parameter<std::string>(BC_name+"/dynamic/from_file") )
        continue;
      //
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(BC_name);
      //
      const double minBC = bdm::GetInitialValueMin(params, BC_name),
                   maxBC = bdm::GetInitialValueMax(params, BC_name);
      //
      std::string fn = params.get<std::string>(BC_name+"/dynamic/from_file");
      //
      const std::string s2f = "*.",
                        s2r = Time+".";
      size_t found = fn.find(s2f);
      if (std::string::npos != found)
        fn.replace(found, s2f.length(), s2r);
      else
        ABORT_("biochemical \""+BC_name+"\" has erroneous dynamic filename");
      //
      std::ifstream fin(fn);
      ASSERT_(fin.good(), "file \""+BC_name+"\" cannot be accessed");
      //
      std::set<size_t> scanned_boxes;
      //
      int npnt = 0;
      fin >> npnt;
      for (int p=0; p<npnt; p++)
        {
          bdm::Double3 xyz;
          fin >> xyz[0] >> xyz[1] >> xyz[2];
          double concentration;
          fin >> concentration;
          // sanity check
          if ( minBC>concentration || maxBC<concentration )
            ABORT_("biochemical \""+BC_name+"\" has been fed erroneous initial values");
          //
          const size_t index = dg->GetBoxIndex(xyz);
          // WARNING: check identical piece of code writen in "init_biochemicals"
          if (scanned_boxes.end()!=scanned_boxes.find(index))
            // ABORT_("biochemical \""+BC_name+"\" has multiple insertion for a DG-box");
            continue;
          else
            scanned_boxes.insert(index);
          //
          const double original_concentration = dg->GetValue(xyz);
          //
          dg->ChangeConcentrationBy(xyz, concentration-original_concentration);
        }
      // create a copy of the file just processed
      const std::string cmd = "cp " + fn + "  "
                            + params.get<std::string>("output_directory")
                            + "/in/" + BC_name + "." + Time + ".dat";
      ASSERT_(0==std::system(cmd.c_str()),
              "could not save a copy of the convection field data file");
      // ...end of biochemical (cue) loop
    }
  // ---------------------------------------------------------------------------
  // CAP scheduler: apply CAP boundary values over an exact physical-time window
  if ( params.have_parameter<bool>("CAP/enabled") && params.get<bool>("CAP/enabled") )
    {
      const int start_step    = params.have_parameter<int>("CAP/start_step")
                              ? params.get<int>("CAP/start_step") : 0;
      const int duration_steps = params.have_parameter<int>("CAP/duration_steps")
                               ? params.get<int>("CAP/duration_steps") : 0;
      const bool all_faces = params.have_parameter<bool>("CAP/all_faces")
                           ? params.get<bool>("CAP/all_faces")
                           : false;
      const std::string cap_face = params.have_parameter<std::string>("CAP/face")
                                 ? params.get<std::string>("CAP/face")
                                 : std::string("z_max");
      const bool has_start_time_h = params.have_parameter<double>("CAP/start_time_h");
      const bool has_duration_h = params.have_parameter<double>("CAP/duration_h");
      const double dt_h = params.get<double>("time_step");
      ASSERT_(dt_h > 0.0, "parameter \"time_step\" must be positive");

      const bool valid_cap_face = (cap_face == "x_min" || cap_face == "x_max" ||
                                   cap_face == "y_min" || cap_face == "y_max" ||
                                   cap_face == "z_min" || cap_face == "z_max");
      ASSERT_(all_faces || valid_cap_face,
              "parameter \"CAP/face\" must be one of: x_min, x_max, y_min, y_max, z_min, z_max");

      const std::string cap_application_mode =
        params.have_parameter<std::string>("CAP/application_mode")
        ? params.get<std::string>("CAP/application_mode")
        : std::string("boundary_faces");
      const std::string cap_post_treatment_mode =
        params.have_parameter<std::string>("CAP/post_treatment_mode")
        ? params.get<std::string>("CAP/post_treatment_mode")
        : std::string("none");
      ASSERT_(cap_application_mode == "boundary_faces"
              || cap_application_mode == "boundary_dirichlet"
              || cap_application_mode == "boundary_neumann"
              || cap_application_mode == "uniform_grid"
              || cap_application_mode == "uniform_domain",
              "parameter \"CAP/application_mode\" must be one of: boundary_faces, boundary_dirichlet, boundary_neumann, uniform_grid, uniform_domain");
      ASSERT_(cap_post_treatment_mode == "none"
              || cap_post_treatment_mode == "neumann_zero_gradient",
              "parameter \"CAP/post_treatment_mode\" must be one of: none, neumann_zero_gradient");

      // Reinit happens before the upcoming BioDynaMo step.
      // Index `time` represents the start of that step interval [time*dt, (time+1)*dt).
      const double step_window_start_h = std::max(0.0, static_cast<double>(time) * dt_h);
      const double step_window_end_h = step_window_start_h + dt_h;
      const double exposure_start_h = has_start_time_h
                                    ? params.get<double>("CAP/start_time_h")
                                    : static_cast<double>(start_step) * dt_h;

      const bool finite_duration = has_duration_h
                                 ? (params.get<double>("CAP/duration_h") > 0.0)
                                 : (duration_steps > 0);
      const double exposure_duration_h = finite_duration
                                       ? (has_duration_h
                                          ? params.get<double>("CAP/duration_h")
                                          : static_cast<double>(duration_steps) * dt_h)
                                       : 0.0;
      const double exposure_end_h = exposure_start_h + exposure_duration_h;
      auto exposure_overlap_h = [&] () -> double {
        if (!finite_duration)
          {
            const double overlap_start = std::max(step_window_start_h, exposure_start_h);
            return std::max(0.0, step_window_end_h - overlap_start);
          }
        return std::max(0.0,
                        std::min(step_window_end_h, exposure_end_h)
                      - std::max(step_window_start_h, exposure_start_h));
      };
      const double exposure_fraction = std::clamp(exposure_overlap_h() / dt_h, 0.0, 1.0);

      const bool cap_active = exposure_fraction > 0.0;
      const bool cap_post_treatment_active = finite_duration
                                           && step_window_start_h >= exposure_end_h
                                           && cap_post_treatment_mode == "neumann_zero_gradient";
      auto set_boundary_constant = [] (bdm::DiffusionGrid* dg, const std::string& face, bool use_all_faces, double value) {
        if (!dg) { return; }
        const int res = static_cast<int>(dg->GetResolution());
        if (res <= 1) { return; }
        auto apply = [&] (auto coord_selector) {
          for (int x = 0; x < res; ++x)
            for (int y = 0; y < res; ++y)
              for (int z = 0; z < res; ++z)
                {
                  if (!coord_selector(x,y,z)) continue;
                  std::array<uint32_t,3> box = {static_cast<uint32_t>(x), static_cast<uint32_t>(y), static_cast<uint32_t>(z)};
                  size_t b = dg->GetBoxIndex(box);
                  double current = dg->GetConcentration(b);
                  dg->ChangeConcentrationBy(b, value - current);
                }
        };
        if (use_all_faces)
          apply([&](int x, int y, int z){
            return x == 0 || y == 0 || z == 0 || x == res - 1 || y == res - 1 || z == res - 1;
          });
        else if (face == "z_max")
          apply([&](int, int, int z){ return z == res-1; });
        else if (face == "z_min")
          apply([&](int, int, int z){ return z == 0; });
        else if (face == "x_min")
          apply([&](int x, int, int){ return x == 0; });
        else if (face == "x_max")
          apply([&](int x, int, int){ return x == res-1; });
        else if (face == "y_min")
          apply([&](int, int y, int){ return y == 0; });
        else if (face == "y_max")
          apply([&](int, int y, int){ return y == res-1; });
      };
      auto set_uniform_constant = [] (bdm::DiffusionGrid* dg, double value) {
        if (!dg) { return; }
        const int res = static_cast<int>(dg->GetResolution());
        if (res <= 0) { return; }
        for (int x = 0; x < res; ++x)
          for (int y = 0; y < res; ++y)
            for (int z = 0; z < res; ++z)
              {
                std::array<uint32_t,3> box = {static_cast<uint32_t>(x), static_cast<uint32_t>(y), static_cast<uint32_t>(z)};
                size_t b = dg->GetBoxIndex(box);
                double current = dg->GetConcentration(b);
                dg->ChangeConcentrationBy(b, value - current);
              }
      };
      const double cap_h2o2 =
        bdm::GetNumericParameterFlexible(params,
                                         "CAP/H2O2/concentration",
                                         0.0);
      const double cap_no2 =
        bdm::GetNumericParameterFlexible(params,
                                         "CAP/NO2_/concentration",
                                         0.0);
      bdm::ValidateCapConcentrations(cap_h2o2, cap_no2);
      const double cap_h2o2_step = cap_h2o2 * exposure_fraction;
      const double cap_no2_step = cap_no2 * exposure_fraction;
      auto apply_cap_boundary_condition = [&] (bdm::DiffusionGrid* dg, bdm::BoundaryConditionType bc_type, double value) {
        if (!dg) { return; }
        dg->SetBoundaryConditionType(bc_type);
        dg->SetBoundaryCondition(
          std::make_unique<bdm::ConstantBoundaryCondition>(value));
      };

      if (cap_application_mode == "boundary_dirichlet"
          || cap_application_mode == "boundary_neumann")
        {
          const auto bc_type = cap_post_treatment_active
                             ? bdm::BoundaryConditionType::kNeumann
                             : (cap_application_mode == "boundary_dirichlet"
                                ? bdm::BoundaryConditionType::kDirichlet
                                : bdm::BoundaryConditionType::kNeumann);
          apply_cap_boundary_condition(rm->GetDiffusionGrid("H2O2"), bc_type,
                                       cap_active ? cap_h2o2_step : 0.0);
          apply_cap_boundary_condition(rm->GetDiffusionGrid("NO2_"), bc_type,
                                       cap_active ? cap_no2_step : 0.0);
        }
      else if (cap_post_treatment_active)
        {
          apply_cap_boundary_condition(rm->GetDiffusionGrid("H2O2"),
                                       bdm::BoundaryConditionType::kNeumann, 0.0);
          apply_cap_boundary_condition(rm->GetDiffusionGrid("NO2_"),
                                       bdm::BoundaryConditionType::kNeumann, 0.0);
        }
      else if (cap_active)
        {
          if (cap_h2o2_step!=0.0)
            {
              auto* dgH = rm->GetDiffusionGrid("H2O2");
              if (cap_application_mode == "uniform_grid"
                  || cap_application_mode == "uniform_domain")
                set_uniform_constant(dgH, cap_h2o2_step);
              else
                set_boundary_constant(dgH, cap_face, all_faces, cap_h2o2_step);
            }
          if (cap_no2_step!=0.0)
            {
              auto* dgN = rm->GetDiffusionGrid("NO2_");
              if (cap_application_mode == "uniform_grid"
                  || cap_application_mode == "uniform_domain")
                set_uniform_constant(dgN, cap_no2_step);
              else
                set_boundary_constant(dgN, cap_face, all_faces, cap_no2_step);
            }
        }
    }

  // Generic one-sided inlet for any biochemical
  for (const auto& BC_name : biochem)
    {
      const bool inlet_enabled = params.have_parameter<bool>(BC_name+"/inlet/enabled")
                               ? params.get<bool>(BC_name+"/inlet/enabled")
                               : false;
      if (!inlet_enabled) continue;

      const std::string face = params.have_parameter<std::string>(BC_name+"/inlet/face")
                             ? params.get<std::string>(BC_name+"/inlet/face")
                             : std::string("z_max");
      const double value = params.have_parameter<double>(BC_name+"/inlet/value")
                         ? params.get<double>(BC_name+"/inlet/value")
                         : bdm::GetInitialValueMax(params, BC_name);
      bdm::ValidateNormalizedConcentrationValue(value, BC_name+"/inlet/value");

      // set value on the specified face only
      auto* dg = rm->GetDiffusionGrid(BC_name);
      if (dg)
        {
          const int res = static_cast<int>(dg->GetResolution());
          if (res > 0)
            {
              auto set_face_constant = [&] () {
                auto apply = [&] (auto coord_selector) {
                  for (int x = 0; x < res; ++x)
                    for (int y = 0; y < res; ++y)
                      for (int z = 0; z < res; ++z)
                        {
                          if (!coord_selector(x,y,z)) continue;
                          std::array<uint32_t,3> box = {static_cast<uint32_t>(x), static_cast<uint32_t>(y), static_cast<uint32_t>(z)};
                          size_t b = dg->GetBoxIndex(box);
                          double current = dg->GetConcentration(b);
                          dg->ChangeConcentrationBy(b, value - current);
                        }
                };
                if (face == "z_max")
                  apply([&](int, int, int z){ return z == res-1; });
                else if (face == "z_min")
                  apply([&](int, int, int z){ return z == 0; });
                else if (face == "x_min")
                  apply([&](int x, int, int){ return x == 0; });
                else if (face == "x_max")
                  apply([&](int x, int, int){ return x == res-1; });
                else if (face == "y_min")
                  apply([&](int, int y, int){ return y == 0; });
                else if (face == "y_max")
                  apply([&](int, int y, int){ return y == res-1; });
              };
              set_face_constant();
            }
        }
    }

  // Linear gradient boundary enforcement: pin the boundary faces each step
  for (const auto& BC_name : biochem)
    {
      const bool linear_gradient_enabled =
        params.have_parameter<bool>(BC_name+"/initial_value/linear_gradient/enabled")
        ? params.get<bool>(BC_name+"/initial_value/linear_gradient/enabled")
        : false;
      if (!linear_gradient_enabled) continue;

      const std::string grad_axis =
        params.get<std::string>(BC_name+"/initial_value/linear_gradient/axis");
      const double grad_min =
        params.get<double>(BC_name+"/initial_value/linear_gradient/min_value1");
      const double grad_max =
        params.get<double>(BC_name+"/initial_value/linear_gradient/max_value1");

      const bool secondary_enabled =
        params.have_parameter<std::string>(BC_name+"/initial_value/linear_gradient/axis2");
      const std::string grad_axis2 = secondary_enabled
        ? params.get<std::string>(BC_name+"/initial_value/linear_gradient/axis2")
        : std::string("");
      const double grad2_min = secondary_enabled
        ? params.get<double>(BC_name+"/initial_value/linear_gradient/min_value2")
        : 0.0;
      const double grad2_max = secondary_enabled
        ? params.get<double>(BC_name+"/initial_value/linear_gradient/max_value2")
        : 0.0;

      auto* dg = rm->GetDiffusionGrid(BC_name);
      if (!dg) continue;
      const int res = static_cast<int>(dg->GetResolution());
      if (res <= 0) continue;

      const double domain_min = params.get<double>("min_boundary");
      const double domain_max = params.get<double>("max_boundary");
      const double domain_range = domain_max - domain_min;

      // helper: compute the gradient value at a given grid box
      auto grad_value = [&] (int bx, int by, int bz) -> double {
        // map box index to normalized coordinate [0,1]
        double t1 = 0.0;
        if      (grad_axis == "x") t1 = static_cast<double>(bx) / (res - 1);
        else if (grad_axis == "y") t1 = static_cast<double>(by) / (res - 1);
        else                       t1 = static_cast<double>(bz) / (res - 1);
        double value = grad_min + t1 * (grad_max - grad_min);
        if (secondary_enabled) {
          double t2 = 0.0;
          if      (grad_axis2 == "x") t2 = static_cast<double>(bx) / (res - 1);
          else if (grad_axis2 == "y") t2 = static_cast<double>(by) / (res - 1);
          else                        t2 = static_cast<double>(bz) / (res - 1);
          value += grad2_min + t2 * (grad2_max - grad2_min);
        }
        return value;
      };

      // pin the boundary faces of the primary gradient axis
      auto pin_face = [&] (auto coord_selector) {
        for (int x = 0; x < res; ++x)
          for (int y = 0; y < res; ++y)
            for (int z = 0; z < res; ++z)
              {
                if (!coord_selector(x,y,z)) continue;
                std::array<uint32_t,3> box = {static_cast<uint32_t>(x),
                                              static_cast<uint32_t>(y),
                                              static_cast<uint32_t>(z)};
                size_t b = dg->GetBoxIndex(box);
                double val = grad_value(x, y, z);
                double current = dg->GetConcentration(b);
                dg->ChangeConcentrationBy(b, val - current);
              }
      };

      if (grad_axis == "x")
        {
          pin_face([&](int x,int,int){ return x == 0;     });
          pin_face([&](int x,int,int){ return x == res-1;  });
        }
      else if (grad_axis == "y")
        {
          pin_face([&](int,int y,int){ return y == 0;     });
          pin_face([&](int,int y,int){ return y == res-1;  });
        }
      else
        {
          pin_face([&](int,int,int z){ return z == 0;     });
          pin_face([&](int,int,int z){ return z == res-1;  });
        }
      // pin the boundary faces of the secondary gradient axis
      if (secondary_enabled)
        {
          if (grad_axis2 == "x")
            {
              pin_face([&](int x,int,int){ return x == 0;     });
              pin_face([&](int x,int,int){ return x == res-1;  });
            }
          else if (grad_axis2 == "y")
            {
              pin_face([&](int,int y,int){ return y == 0;     });
              pin_face([&](int,int y,int){ return y == res-1;  });
            }
          else
            {
              pin_face([&](int,int,int z){ return z == 0;     });
              pin_face([&](int,int,int z){ return z == res-1;  });
            }
        }
    }
  if ( params.have_parameter<bool>("CAP/enabled") && params.get<bool>("CAP/enabled") )
    {
      const double dt_h = params.get<double>("time_step");
      for (const auto& BC_name : biochem)
        {
          if (!bdm::UseImplicitCapTransport(params, BC_name)) continue;
          bdm::ApplyImplicitTransportStep(rm->GetDiffusionGrid(BC_name),
                                          params, BC_name, dt_h);
        }
    }
}
// =============================================================================
inline
void ioflux_cells(bdm::Simulation& sim,
                  const std::map<int, std::string>& cells_phenotype, const int time)
{
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = sim.GetRandom();
  //
  // iterate for all surfaces through which we anticipate a cell IO-flux
  for (unsigned int isurf=0; isurf<io_flux.surface.size(); isurf++)
    {
      // obtain the ID for this cell phenotype
      const int CP_ID = io_flux.surface[isurf].phenotype;
      //
      const std::string& CP_name = cells_phenotype.find(CP_ID)->second;
      //
      // obtain the cell increment from the rate (flux) for this (cell) phenotype
      int D_cell;
      {
        const double T = params.get<double>("current time");
        const double* T_params = io_flux.surface[isurf].time_parameters.data();
        const int time_function = (int) T_params[0];
        //
        double D_cell_T = 0.0;
        if      (time_function== 0) D_cell_T = flat_function(T_params, T);
        else if (time_function==10) D_cell_T = linear_function(T_params, T);
        else if (time_function==20) D_cell_T = step_function(T_params, T);
        else if (time_function==30) D_cell_T = ramp_function(T_params, T);
        else if (time_function==40) D_cell_T = gaussian_function(T_params, T);
        else if (time_function==50) D_cell_T = logistic_function(T_params, T);
        else if (time_function==11) D_cell_T = linear_periodic_function(T_params, T);
        else if (time_function==21) D_cell_T = step_periodic_function(T_params, T);
        else if (time_function==31) D_cell_T = ramp_periodic_function(T_params, T);
        else if (time_function==41) D_cell_T = gaussian_periodic_function(T_params, T);
        else if (time_function==51) D_cell_T = logistic_periodic_function(T_params, T);
        else
          ABORT_("uncategorized type of time function inserted for cell flux");
        //
        ASSERT_(D_cell_T>=0.0,
                "increment for cell flux cannot be negative");
        //
        D_cell = (int) D_cell_T;
      }
      //
      ASSERT_(CP_ID>=0, "an internal error occurred");
      //
      if (time%io_flux.surface[isurf].time_step) continue;
      //
      // WARNING: current implementation considers only for an influx of cells,
      //          not an outflux of cells!!!
      //
      if (D_cell<=0) continue;
      //
      // access the list of all triangles for this surface
      const std::vector<SurfaceSTL::Triangle> tri_v = io_flux.surface[isurf].triangle;
      // insert the incoming cells now one by one...
      for (int c=0; c<D_cell; c++)
        {
          // pick a random triangle to trow in the incoming cell
          const int triangle_index = uniform_distro(0, tri_v.size()-1);
          const SurfaceSTL::Triangle& tri = tri_v[triangle_index];
          //+std::cout << " BiologicalCell " << (c+1) << " SurfaceSTL::Triangle " << (triangle_index+1) << std::endl;
          //
          // the following algorithm has been implemented:
          // https://blogs.sas.com/content/iml/2020/10/19/random-points-in-triangle.html
          bdm::Double3 p;
          //
          const bdm::Double3 a = tri.vertex_1 - tri.vertex_0;
          const bdm::Double3 b = tri.vertex_2 - tri.vertex_0;
          //
          bool is_successfull = false;
          size_t trial = 0;
          // run until you reach success
          while (!is_successfull)
            {
              // update this counter and confirm we haven't max'd it out...
              ++trial;
              if (100==trial) break;
              //
              double u1 = rg->Uniform(0.0, 1.0);
              double u2 = rg->Uniform(0.0, 1.0);
              while (u1+u2>1.0)
                {
                  u1 = rg->Uniform(0.0, 1.0);
                  u2 = rg->Uniform(0.0, 1.0);
                }
              p = (a * u1 + b * u2) + tri.vertex_0;
              //
              if (p[0]<params.get<double>("min_boundary") || p[0]>params.get<double>("max_boundary"))
                continue;
              if (p[1]<params.get<double>("min_boundary") || p[1]>params.get<double>("max_boundary"))
                continue;
              if (p[2]<params.get<double>("min_boundary") || p[2]>params.get<double>("max_boundary"))
                continue;
              //
              const bdm::Double3 p2t = project_to_plane(tri.normal, tri.vertex_0, p, 1.0e-3),
                                 r   = p - p2t;
              const double delta = sqrt(pow2(r[0])+pow2(r[1])+pow2(r[2]));
              if (delta>1.0e-1)
                continue;
              //
              const double areaT = triangle_area(tri.vertex_0, tri.vertex_1, tri.vertex_2),
                           area1 = triangle_area(p,            tri.vertex_1, tri.vertex_2),
                           area2 = triangle_area(tri.vertex_0, p,            tri.vertex_2),
                           area3 = triangle_area(tri.vertex_0, tri.vertex_1, p           );
              if ((area1+area2+area3)>areaT)
                continue;
              //
              //+ std::cout << params.get<double>("current time") << " ==== " << isurf << " ==== "
              //+  << ' ' << p << ' ' << p2t << ' ' << r << ' ' << delta << std::endl;
              is_successfull = true;
            }
          //
          if (is_successfull)
            {
              // cell position
              const bdm::Double3 xyz = p;
              // cell cycle phase
              int ccp;
              {
                const int min_ccp = (int) bdm::BiologicalCell::Phase::G1,
                          max_ccp = (int) bdm::BiologicalCell::Phase::Tr,
                          range = max_ccp - min_ccp + 1;
                ccp = rand() % range + min_ccp;
              }
              // cell diameter
              const double dia = rg->Uniform(params.get<double>(CP_name+"/diameter/min"),
                                             params.get<double>(CP_name+"/diameter/max"));
              // order of mechanisms that define behaviour
              const int mo = params.get<int>(CP_name+"/mechanism_order");
              // principal directions of the cell polarization matrix
              const double pd0 = params.get<double>(CP_name+"/principal/0"),
                           pd1 = params.get<double>(CP_name+"/principal/1"),
                           pd2 = params.get<double>(CP_name+"/principal/2");
              // cell polarization matrix
              bdm::Double3x3 pl_3x3 = diag(pd0, pd1, pd2);
              //
              bdm::BiologicalCell* cell = new bdm::BiologicalCell(CP_ID, xyz);
              cell->SetPhase(ccp);
              cell->SetDiameter(dia);
              cell->SetAdherence(params.have_parameter<double>(CP_name+"/adherence")
                               ? params.get<double>(CP_name+"/adherence") : 0.0);
              if (params.have_parameter<double>(CP_name+"/density"))
                cell->SetDensity(params.get<double>(CP_name+"/density"));
              cell->SetParametersPointer(&params);
              cell->SetCanApoptose(params.get<bool>(CP_name+"/can_apoptose"));
              cell->SetCanGrow(params.get<bool>(CP_name+"/can_grow"));
              cell->SetCanDivide(params.get<bool>(CP_name+"/can_divide"));
              cell->SetCanMigrate(params.get<bool>(CP_name+"/can_migrate"));
              cell->SetCanTransform(params.get<bool>(CP_name+"/can_transform"));
              cell->SetCanPolarize(params.get<bool>(CP_name+"/can_polarize"));
              cell->SetCanProtrude(params.get<bool>(CP_name+"/can_protrude"));
              cell->SetPolarization(pl_3x3);
              if      (10==mo)
                cell->AddBehavior(new bdm::Biology4BiologicalCell_10());
              else if (11==mo)
                cell->AddBehavior(new bdm::Biology4BiologicalCell_11());
              else
                ABORT_("\""+CP_name+"\" with phenotype ID \""+std::to_string(CP_ID)+"\" has unrecognized behavior");
              // store this cell into BioDynaMo's resource manager
              rm->AddAgent(cell);
            }
          // ...end of influx cells loop
        }
      // ...end of IO-flux surfaces
    }
}
// =============================================================================
inline
void one_off_init(bdm::Simulation& sim)
{
  // access BioDynaMo's resource manager
  auto* rm = sim.GetResourceManager();
  //
  if (params.get<bool>("simulation_models_vessels"))
    // ...perform calculations below only when simulation encompasses vessels
    rm->ForEachAgent([&] (bdm::Agent* a) {
      if (auto* vessel = dynamic_cast<bdm::Vessel*>(a))
        {
          // obtain the vessel ID based on the branch it belongs
          const int ID = vessel->GetVesselID();
          //
          const int age = vessel_map__ID_age.find(ID)->second;
          // update the age for this vessel (of this branch)
          vessel->SetAge(age);
        }
    });
  //
}
// =============================================================================
inline
int simulate(const std::string& fname, const int seed)
{
  // read all the model parameters
  std::vector<std::string> biochem;
  std::map<int,std::string> cells;
  read_csv_file(fname, cells, biochem);
  std::cout << "Parameter read complete..." << std::endl;
  //
  std::cout << ("[input file: "+fname+"; "
               +std::to_string(biochem.size())+" biochemicals; "
               +std::to_string(cells.size())+" cell phenotypes]") << std::endl;
  //
  params.set<std::vector<std::string>>("substances") = biochem;
  params.set<int>("N_biochemicals") = biochem.size();
  params.set<int>("N_phenotypes") = cells.size();
  params.set<int>("N_vessel_IDs") = 0;
  // initialize BioDynaMo modules
  bdm::neuroscience::InitModule();
  // initialize the BioDynaMo simulation
  bdm::Simulation sim(params.get<std::string>("simulation_title"), set_bdm_params);
  sim.Activate();
  sim.GetEnvironment()->Update();
  sim.GetRandom()->SetSeed(0==seed ? time(0) : seed);
  // copy the input file into the results directory
  const std::string cmd = "cp " + fname + " "
                        + params.get<std::string>("output_directory") + "/input.csv";
  ASSERT_(0==std::system(cmd.c_str()), "could not save a copy of \"input.csv\"");
  // now let's start with the simulation initializations
  std::cout << "Simulation initializes..." << std::endl;
  // load the data for the obstacles in the simulation
  init_obstacles();
  // load all biochemicals in the simulation
  init_biochemicals(sim, biochem);
  // load all vessels in the simulation
  init_vessels(sim, biochem);
  // load all cells in the simulation
  init_cells(sim, cells, biochem);
  // load the convection field data in the simulation
  init_convection(sim);
  // enforce t=0 dynamic boundary conditions before the first simulated step
  reinit_biochemicals(sim, biochem, 0);
  // check for cells input/output flux to the domain
  ioflux_cells(sim, cells, 0);
  // generate file to save simulation statistics
  std::ofstream fstat(params.get<std::string>("output_directory")+"/stats.csv");
  save_stats(sim, cells, fstat);
  // activate the BioDynaMo simulation and save initial data
  save_snapshot(sim);
  // run the BioDynaMo simulation for all time-steps
  std::cout << "Simulation runs..." << std::endl;
  const int n_time = params.get<int>("number_of_time_steps"),
            stat_step = params.get<int>("statistics_interval"),
            viz_step = params.get<int>("visualization_interval");
  const double time_step = params.get<double>("time_step");
  const bool early_stop_on_cell_limit =
    params.have_parameter<bool>("simulation/early_stop_on_total_cells_exceeded")
    ? params.get<bool>("simulation/early_stop_on_total_cells_exceeded")
    : false;
  const int total_cell_limit =
    params.have_parameter<int>("simulation/total_cell_limit")
    ? params.get<int>("simulation/total_cell_limit")
    : std::numeric_limits<int>::max();
  for (int time=1; time<=n_time; time++)
    {
      params.set<int>("index time") = time;
      const double TIME = time*time_step;
      params.set<double>("current time") = TIME;
      time_status_bar(std::cout, time, n_time, TIME);
      // t=0 conditions are already applied during initialization.
      // For later steps, update fields for the upcoming interval before stepping.
      if (time > 1) reinit_biochemicals(sim, biochem, time-1);
      // run the BioDynaMo simulator for one step
      sim.GetScheduler()->Simulate(1);
      if (1==time) one_off_init(sim);
      // save simulation statistics in a file stream
      if (0==time%stat_step) save_stats(sim, cells, fstat);
      // output data for Paraview visualization (after inlet enforcement)
      if (0==time%viz_step) save_snapshot(sim, time);
      if (early_stop_on_cell_limit)
        {
          int total_cells = 0;
          auto* rm = sim.GetResourceManager();
          rm->ForEachAgent([&] (bdm::Agent* a) {
            if (dynamic_cast<bdm::BiologicalCell*>(a))
              ++total_cells;
          });
          if (total_cells > total_cell_limit)
            {
              std::cout << "Simulation early stop: total cell count "
                        << total_cells << " exceeded limit " << total_cell_limit
                        << " at time " << TIME << std::endl;
              break;
            }
        }
      // reset the data for the obstacles in the simulation
      reinit_obstacles(time);
      // reset some data for all cells in the simulation
      reinit_cells(sim, cells);
      // check for cells input/output flux to the domain
      ioflux_cells(sim, cells, time);
      // ...and convection field (if present)
      set_convection(sim, time);
    }
  std::cout << "Simulation terminating..." << std::endl;
  // empty the local containers
  params.clear();
  all_agents.clear();
  obstacles.clear();
  io_flux.clear();
  dg_vec.clear();
  vessel_map__ID_age.clear();
  // close file for simulation statistics
  fstat.close();
  // exit the function normally
  return 0;
}
#endif

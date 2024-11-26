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
#ifndef _BIOLOGICAL_CELL_INLINE_H_
#define _BIOLOGICAL_CELL_INLINE_H_
// =============================================================================
inline
void bdm::BiologicalCell::RunBiochemics()
{
  // by design only viable (non-necrotic) cells could secrete biochemicals
  if (!this->GetPhenotype()) return;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary"),
               meanCOORD = (maxCOORD+minCOORD)/2.0,
               deltaCOORD = maxCOORD - meanCOORD;
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // cell spatial coordinates
  const bdm::Double3 xyz = this->GetPosition();
  // exit function if cell resides at the edge of the simulation domain,
  // but first identify the mode / dimension of the simulation domain
  if (this->params()->get<bool>("simulation_domain_is_polar"))
    {
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        {
          if ( sqrt(pow2(xyz[0]-meanCOORD)
                   +pow2(xyz[1]-meanCOORD)) > deltaCOORD-tol ) return;
        }
      else
        {
          if ( sqrt(pow2(xyz[0]-meanCOORD)
                   +pow2(xyz[1]-meanCOORD)
                   +pow2(xyz[2]-meanCOORD)) > deltaCOORD-tol ) return;
        }
    }
  else
    {
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        {
          if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol ||
               xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol ) return;
        }
      else
        {
          if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol ||
               xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol ||
               xyz[2] < minCOORD+tol || xyz[2] > maxCOORD-tol ) return;
        }
    }
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, xyz, tol))
  // iterate for all substances
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      Biochemical BC_id =
        static_cast<Biochemical>(dg->GetContinuumId()); // biochemical ID
      //
      // skip following calculations for radiation!!!
      if ( Biochemical::RAD == BC_id ) continue;
      //
      if (! this->params()->have_parameter<double>(CP_name+"/"+BC_name+"/secretion/net_balance"))
        continue;
      //
      const double concentration = dg->GetValue(xyz);
      // parameters that modulate biochemical cue secretion (production or consumption)
      const double BC_stdev = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std")<=0.0 ? 1.0 :
                              rg->Uniform(1.0-this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std"),
                                          1.0+this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/net_balance/std"));
      const double net_balance = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/net_balance") * BC_stdev;
      //
      // skip subsequent calculations if net balance of this biochemical cue secretion is
      // equal to absolute zero!!!
      if (0.0 == net_balance) continue;
      //
      const double saturation = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/saturation");
      //
      if (! this->params()->get<bool>(CP_name+"/"+BC_name+"/secretion/dependency"))
        {
          if (net_balance > 0.0)
            {
              // increase concentration
              if ( saturation > 0.0 )
                {
                  if (concentration>saturation) ; // ... do nothing!
                  else                          dg->ChangeConcentrationBy(xyz, net_balance);
                }
              else
                { // no saturation effects are accounted in
                  dg->ChangeConcentrationBy(xyz, net_balance);
                }
              // increased concentration
            }
          else
            { // ... net_balance < 0.0
              // decrease concentration
              if ( concentration > 0.0 )
                {
                  if (concentration+net_balance>0.0) dg->ChangeConcentrationBy(xyz, net_balance);
                  else                               dg->ChangeConcentrationBy(xyz, -concentration);
                }
              // decreased concentration
            }
          // ...exit function normally
          return;
        }
      //
      // check for positive or negative feedback loop from other substances
      // iterate for all substances - simply ignore for itself ;)
      for ( std::vector<std::string>::const_iterator
            cj=substances.begin(); cj!=substances.end(); cj++ )
        {
          // skip if the same substance...
          if (ci==cj) continue;
          // access the BioDynaMo diffusion grid
          auto* dg_other = rm->GetDiffusionGrid(*cj);
          const std::string BC_other_name = dg_other->GetContinuumName(); // biochemical name
          //
          const double concentration_other = dg_other->GetValue(xyz),
                       threshold_other = this->params()->get<double>(CP_name+"/"+BC_name+"/secretion/"+BC_other_name+"/threshold");
          // check if other substances regulate secretion of this substance...
          if ( ( threshold_other > 0.0 && concentration_other > +threshold_other ) ||
               ( threshold_other < 0.0 && concentration_other < -threshold_other ) )
            {
              //
              if (net_balance > 0.0)
                {
                  // increase concentration
                  if ( saturation > 0.0 )
                    {
                      if (concentration>saturation) ; // ... do nothing!
                      else                          dg->ChangeConcentrationBy(xyz, net_balance);
                    }
                  else
                    { // no saturation effects are accounted in
                      dg->ChangeConcentrationBy(xyz, net_balance);
                    }
                  // increased concentration
                }
              else
                { // ... net_balance < 0.0
                  // decrease concentration
                  if ( concentration > 0.0 )
                    {
                      if (concentration+net_balance>0.0) dg->ChangeConcentrationBy(xyz, net_balance);
                      else                               dg->ChangeConcentrationBy(xyz, -concentration);
                    }
                  // decreased concentration
                }
              //
            }
          //...end of other substances loop
        }
      //...end of substances loop
    }
  //...end of cell biochemics
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckPositionValidity()
{
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary"),
               meanCOORD = (maxCOORD+minCOORD)/2.0,
               deltaCOORD = maxCOORD - meanCOORD;
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  // check if simulation domain is bounded or unbounded (for cell outflux)
  if ( this->params()->get<bool>("simulation_domain_is_bounded") )
    {
      bdm::Double3 xyz = this->GetPosition();
      // identify mode of simulation domain
      if      ( this->params()->get<bool>("simulation_domain_is_polar") &&
                this->params()->get<bool>("simulation_domain_is_2D")    )
        {
          const double radius = sqrt(pow2(xyz[0]-meanCOORD)
                                    +pow2(xyz[1]-meanCOORD));
          const double phi = atan2(xyz[1], xyz[0]);
          bool updated_xyz = false;
          if ( radius > deltaCOORD-tol )
            {
              xyz[0] = meanCOORD + (deltaCOORD-tol)*cos(phi);
              xyz[1] = meanCOORD + (deltaCOORD-tol)*sin(phi);
              updated_xyz = true;
            }
          if        ( xyz[2] < meanCOORD-tol )
            {
              xyz[2] = meanCOORD - tol;
              updated_xyz = true;
            }
          else if ( xyz[2] > meanCOORD+tol )
            {
              xyz[2] = meanCOORD + tol;
              updated_xyz = true;
            }
          //
          if ( updated_xyz )
            this->SetPosition(xyz);
        }
      else if ( this->params()->get<bool>("simulation_domain_is_polar") &&
              ! this->params()->get<bool>("simulation_domain_is_2D")    )
        {
          const double radius = sqrt(pow2(xyz[0]-meanCOORD)
                                    +pow2(xyz[1]-meanCOORD)
                                    +pow2(xyz[2]-meanCOORD));
          const double phi = atan2(xyz[1], xyz[0]);
          const double theta = atan2(sqrt(pow2(xyz[0])+pow2(xyz[1])), xyz[2]);
          bool updated_xyz = false;
          if ( radius > deltaCOORD-tol )
            {
              xyz[0] = meanCOORD + (deltaCOORD-tol)*sin(theta)*cos(phi);
              xyz[1] = meanCOORD + (deltaCOORD-tol)*sin(theta)*sin(phi);
              xyz[2] = meanCOORD + (deltaCOORD-tol)*cos(theta);
              updated_xyz = true;
            }
          //
          if ( updated_xyz )
            this->SetPosition(xyz);
        }
      else if ( ! this->params()->get<bool>("simulation_domain_is_polar") &&
                  this->params()->get<bool>("simulation_domain_is_2D")    )
        {
          bool updated_xyz = false;
          for (size_t d=0; d<2; d++)
            {
              if      ( xyz[d] < minCOORD+tol )
                {
                  xyz[d] = minCOORD + tol;
                  updated_xyz = true;
                }
              else if ( xyz[d] > maxCOORD-tol )
                {
                  xyz[d] = maxCOORD - tol;
                  updated_xyz = true;
                }
            }
          if      ( xyz[2] < meanCOORD-tol )
            {
              xyz[2] = meanCOORD - tol;
              updated_xyz = true;
            }
          else if ( xyz[2] > meanCOORD+tol )
            {
              xyz[2] = meanCOORD + tol;
              updated_xyz = true;
            }
          //
          if ( updated_xyz )
            this->SetPosition(xyz);
        }
      else if ( ! this->params()->get<bool>("simulation_domain_is_polar") &&
                ! this->params()->get<bool>("simulation_domain_is_2D")    )
        {
          bool updated_xyz = false;
          for (size_t d=0; d<3; d++)
            {
              if      ( xyz[d] < minCOORD+tol )
                {
                  xyz[d] = minCOORD + tol;
                  updated_xyz = true;
                }
              else if ( xyz[d] > maxCOORD-tol )
                {
                  xyz[d] = maxCOORD - tol;
                  updated_xyz = true;
                }
            }
          //
          if ( updated_xyz )
            this->SetPosition(xyz);
        }
      else
        ABORT_("an exception is caught");
    }
  else
    {
      bool reached_boundary = false;
      const bdm::Double3 xyz = this->GetPosition();
      // spherical coordinate: radius
      double radius = 0.0;
      // identify mode of simulation domain
      if (this->params()->get<bool>("simulation_domain_is_2D"))
      {
        radius = sqrt(pow2(xyz[0]-meanCOORD)
                     +pow2(xyz[1]-meanCOORD));
        if (this->params()->get<bool>("simulation_domain_is_polar"))
        {
          if ( radius > deltaCOORD-tol )
            reached_boundary = true;
        }
        else
        {
          if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol ||
               xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol )
            reached_boundary = true;
        }
      }
      else
      {
        radius = sqrt(pow2(xyz[0]-meanCOORD)
                     +pow2(xyz[1]-meanCOORD)
                     +pow2(xyz[2]-meanCOORD));
        if (this->params()->get<bool>("simulation_domain_is_polar"))
        {
          if ( radius > deltaCOORD-tol )
            reached_boundary = true;
        }
        else
        {
          if ( xyz[0] < minCOORD+tol || xyz[0] > maxCOORD-tol ||
               xyz[1] < minCOORD+tol || xyz[1] > maxCOORD-tol ||
               xyz[2] < minCOORD+tol || xyz[2] > maxCOORD-tol )
            reached_boundary = true;
        }
      }
      //
      if (reached_boundary)
      {
        // spherical coordinates: inclination (theta), azimuth (phi)
        double theta = 0.0, phi = 0.0;
        // identify mode of simulation domain
        if (this->params()->get<bool>("simulation_domain_is_2D"))
          theta = 0.5*bdm::Math::kPi;
        else
          theta = atan2(sqrt(pow2(xyz[0]-meanCOORD)+pow2(xyz[1]-meanCOORD)),
                        (xyz[2]-meanCOORD));
        phi = atan2((xyz[1]-meanCOORD), (xyz[0]-meanCOORD));
        // reset the spherical coordinates for all escaping cells:
        // radius, inclination (theta) and azimuth (phi)
        std::vector<double>& escaped_cells_radius =
          this->params()->set<std::vector<double>>(CP_name+"/escaped_cells/radius");
        std::vector<double>& escaped_cells_theta  =
          this->params()->set<std::vector<double>>(CP_name+"/escaped_cells/theta");
        std::vector<double>& escaped_cells_phi    =
          this->params()->set<std::vector<double>>(CP_name+"/escaped_cells/phi");
        //
        std::vector<int>& escaped_cells_phase =
          this->params()->set<std::vector<int>>(CP_name+"/escaped_cells/phase");
        //
        escaped_cells_radius.push_back(radius);
        escaped_cells_theta.push_back(theta);
        escaped_cells_phi.push_back(phi);
        //
        escaped_cells_phase.push_back(this->GetPhase());
        // cell has not remained inside the simulation domain
        return false;
      }
      // in case of a 2D simulation domain, ensure you
      // "correct" and bring back any cells on-plane
      if (this->params()->get<bool>("simulation_domain_is_2D"))
      {
        bdm::Double3 xyz = this->GetPosition();
        if      ( xyz[2] < meanCOORD-tol )
        {
          xyz[2] = meanCOORD - tol;
          this->SetPosition(xyz);
        }
        else if ( xyz[2] > meanCOORD+tol )
        {
          xyz[2] = meanCOORD + tol;
          this->SetPosition(xyz);
        }
      }
    }
  // access the pointer to parameter of the simulation obstacles object
  const SimulationObstacles* obstacles =
    this->params()->get<SimulationObstacles*>("simulation_obstacles");
  // check if cell has reached any of the 'box' simulation obstacles
  for (size_t l=0; l<obstacles->box.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetPhenotype()) continue;
      // original cell position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      bdm::Double3 displace = xyz * (-1.0);
      // it true, then cell position needs to be corrected
      if ("box/inside"==obstacles->box[l].type)
        {
          // check if cell position is outside this obstacle
          if (! obstacles->box[l].is_inside(xyz))
            continue;
          // direction vector from point to center
          bdm::Double3 d_vector = xyz - obstacles->box[l].center;
          if (!normalize(d_vector, d_vector))
            ABORT_("could not normalize the direction vector");
          //
          const double dot_0 = (d_vector*obstacles->box[l].laxis_0),
                       dot_1 = (d_vector*obstacles->box[l].laxis_1),
                       dot_2 = (d_vector*obstacles->box[l].laxis_2);
          //
          std::map<double, int, std::greater<double>> sorted_map;
          sorted_map.insert( std::make_pair(fabs(dot_0), 0) );
          sorted_map.insert( std::make_pair(fabs(dot_1), 1) );
          sorted_map.insert( std::make_pair(fabs(dot_2), 2) );
          const int dir = sorted_map.begin()->second;
          // origin (center) point to box face
          bdm::Double3 origin;
          // outward unit normal vector to box face
          bdm::Double3 normal;
          if (0==dir)
            {
              if (dot_0>0.0)
                { origin = obstacles->box[l].center2face_0; normal = obstacles->box[l].laxis_0*(+1.0); }
              else
                { origin = obstacles->box[l].center2face_1; normal = obstacles->box[l].laxis_0*(-1.0); }
              normal = obstacles->box[l].laxis_0;
            }
          else if (1==dir)
            {
              if (dot_1>0.0)
                { origin = obstacles->box[l].center2face_2; normal = obstacles->box[l].laxis_1*(+1.0); }
              else
                { origin = obstacles->box[l].center2face_3; normal = obstacles->box[l].laxis_1*(-1.0); }
              normal = obstacles->box[l].laxis_1;
            }
          else
            {
              if (dot_2>0.0)
                { origin = obstacles->box[l].center2face_4; normal = obstacles->box[l].laxis_2*(+1.0); }
              else
                { origin = obstacles->box[l].center2face_5; normal = obstacles->box[l].laxis_2*(-1.0); }
            }
          //
          xyz = project_to_plane(normal, origin, xyz)
              + normal * (rg->Uniform(0.5,0.8)*this->GetDiameter());
          // enforce cell to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
        }
      else if ("box/outside"==obstacles->box[l].type)
        {
          // check if cell position is inside this obstacle
          if (obstacles->box[l].is_inside(xyz))
            continue;
          // calculate the outward unit normal vector and the projection vector of
          // the cell to each box face
          std::vector< std::pair<bdm::Double3,bdm::Double3> > data_vector;
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_0;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_0, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_1;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_1, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_2;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_2, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_3;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_3, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_4;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_4, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          {
            const bdm::Double3& normal = obstacles->box[l].normal2face_5;
            const bdm::Double3  proj = project_to_plane(normal, obstacles->box[l].center2face_5, xyz);
            data_vector.push_back( std::make_pair(normal, proj) );
          }
          //
          for (auto d : data_vector)
            {
              const bdm::Double3& normal = d.first;
              const bdm::Double3& proj   = d.second;
              //
              const bdm::Double3 xyz_proj = xyz - proj;
              if (xyz_proj*normal>0.0) continue;
              //
              xyz = proj
                  + normal * (rg->Uniform(0.5,0.8)*this->GetDiameter());
              // enforce cell to lie on the (box) obstacle surface
              this->SetPosition(xyz);
              //
              displace += xyz;
              this->UpdateTrail(L2norm(displace));
              // cell has remained inside the simulation domain
              return true;
            }
        }
      else
        ABORT_("an exception is caught");
      // ...end of 'box' simulation obstacles loop
    }
  // check if cell has reached any of the 'sphere' simulation obstacles
  for (size_t l=0; l<obstacles->sphere.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetPhenotype()) continue;
      // original cell position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      bdm::Double3 displace = xyz * (-1.0);
      // it true, then cell position needs to be corrected
      if ("sphere/inside"==obstacles->sphere[l].type)
        {
          // check if cell position is outside this obstacle
          if (! obstacles->sphere[l].is_inside(xyz))
            continue;
          // direction vector from point to center
          bdm::Double3 d_vector = xyz - obstacles->sphere[l].center;
          if (!normalize(d_vector, d_vector))
            ABORT_("could not normalize the direction vector");
          //
          const double delta = obstacles->sphere[l].radius
                             + (rg->Uniform(0.4,0.8)*this->GetDiameter());
          xyz = obstacles->sphere[l].center + d_vector * delta;
          // enforce cell to lie on the surface of the (spherical) obstacle
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
        }
      else if ("sphere/outside"==obstacles->sphere[l].type)
        {
          // check if cell position is inside this obstacle
          if (obstacles->sphere[l].is_inside(xyz))
            continue;
          // direction vector from point to center
          bdm::Double3 d_vector = xyz - obstacles->sphere[l].center;
          if (!normalize(d_vector, d_vector))
            ABORT_("could not normalize the direction vector");
          //
          const double delta = obstacles->sphere[l].radius
                             - (rg->Uniform(0.4,0.8)*this->GetDiameter());
          xyz = obstacles->sphere[l].center + d_vector * delta;
          // enforce cell to lie on the surface of the (spherical) obstacle
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
        }
      else
        ABORT_("an exception is caught");
      // ...end of 'sphere' simulation obstacles loop
    }
  // check if cell has reached any of the 'surface' simulation obstacles
  for (size_t l=0; l<obstacles->surface.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetPhenotype()) continue;
      // original cell position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      bdm::Double3 displace = xyz * (-1.0);
      //
      const double safe_distance = 0.75 * this->GetDiameter();
      //
      std::map<double, std::pair<unsigned int,bdm::Double3>> tri3_proj;
      //
      const unsigned int n_tri3 = obstacles->surface[l].triangle.size();
      for (unsigned int t=0; t<n_tri3; t++)
        {
          const ObstacleSTL::Triangle& tri3 = obstacles->surface[l].triangle[t];
          // origin (center) point to triangle
          const bdm::Double3& origin = tri3.center;
          // outward unit normal vector to triangle
          const bdm::Double3& normal = tri3.normal;
          // internal point wrt the user-defined surface
          const bdm::Double3& l0 = tri3.inside;
          // cell position intersection to the user-defined surface
          bdm::Double3 intx;
          if (! line_intersects_plane(l0, xyz, normal, origin, intx))
            continue;
          // skip following computations if projection point is outside triangle
          if (! is_inside_triangle(tri3.vertex_0, tri3.vertex_1, tri3.vertex_2, intx))
            continue;
          // cell position projection to the user-defined surface
          const bdm::Double3 proj = project_to_plane(normal, origin, xyz);
          //
          const bdm::Double3 xyz_proj = xyz - proj;
          const double distance = L2norm(xyz_proj);
          // skip following if cell is well outside the user-defined surface
          if ((xyz_proj*normal)>0.0 && distance>this->GetDiameter())
            continue;
          //
          auto data = std::make_pair(t, proj);
          tri3_proj.insert( std::make_pair(distance, data) );
          // ...end of triangles loop
        }
      //
      if (!tri3_proj.empty())
        {
          std::map<double, std::pair<unsigned int,bdm::Double3>>::const_iterator
            ci = tri3_proj.begin();
          // access the triangle first...
          const unsigned int t = ci->second.first;
          const ObstacleSTL::Triangle& tri3 = obstacles->surface[l].triangle[t];
          // outward unit normal vector to triangle
          bdm::Double3 normal = tri3.normal;
          if (!normalize(normal, normal))
            ABORT_("could not normalize the normal vector");
          // cell position projection to the user-defined surface
          const bdm::Double3& proj = ci->second.second;
          //
          xyz = proj + normal * safe_distance;
          // enforce cell to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
          // ...end of if-case
        }
      // ...end of 'surface' simulation obstacles loop
    }
  // check if cell has reached any of the 'scaffold' simulation obstacles
  for (size_t l=0; l<obstacles->scaffold.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetPhenotype()) continue;
      // original cell position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      bdm::Double3 displace = xyz * (-1.0);
      //
      const unsigned int n_segm = obstacles->scaffold[l].segment.size();
      for (unsigned int s=0; s<n_segm; s++)
        {
          const ObstacleScaffold::Segment& segm = obstacles->scaffold[l].segment[s];
          //
          const bdm::Double3 n0 = segm.vertex_0,
                             n1 = segm.vertex_1;
          // check if cell position is inside this obstacle
          if (! is_inside_segment(n0, n1, xyz))
            continue;
          // cell position projection to the user-defined segment
          const bdm::Double3 proj = project_to_line(n0, n1, xyz);
          const double distance = L2norm(bdm::Double3(xyz-proj))
                                - segm.radius;
          //
          if (distance>rg->Uniform(0.5,1.0)*this->GetDiameter()) continue;
          //
          bdm::Double3 normal = xyz - proj;
          if (!normalize(normal, normal))
            ABORT_("could not normalize the normal vector");
          //
          const double delta = (rg->Uniform(0.5,1.0)*this->GetDiameter());
          xyz = proj + normal * delta;
          // enforce cell to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          //
          displace += xyz;
          this->UpdateTrail(L2norm(displace));
          // cell has remained inside the simulation domain
          return true;
          // ...end of segments loop
        }
      // ...end of 'scaffold' simulation obstacles loop
    }
  // cell has remained inside the simulation domain
  return true;
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckApoptosisAging()
{
  if (!this->GetCanApoptose()) return false;
  //
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  if (this->params()->get<double>(CP_name+"/can_apoptose/probability_increment_with_age")>0.0)
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_apoptose/probability")
                                +this->params()->get<double>(CP_name+"/can_apoptose/probability_increment_with_age")
                                *this->GetAge() )
        return false;
    }
  else
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_apoptose/probability"))
        return false;
    }
  //
  const int cell_maturity = this->params()->get<int>(CP_name+"/can_apoptose/time_window");
  //
  // since cell has apoptosed (due to ageing), then it must be removed from simulation
  if (this->GetAge()>cell_maturity) return true;
  // since cell has not been through apoptosis (due to ageing), then it can do other things
  return false;
  //...end of cell apoptosis
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckApoptosis()
{
  if (!this->GetCanApoptose()) return false;
  // by design only viable (non-necrotic) cells could apoptose
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string& BC_name = dg->GetContinuumName(); // biochemical name
      //
      const double concentration = dg->GetValue(this->GetPosition()),
                   threshold = this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          // allow cell apoptosis controlled by a combination of two biochemical cues, therefore
          // cell survival is dependent from another substance as well
          if (this->params()->get<bool>(CP_name+"/can_apoptose/"+BC_name+"/dependency"))
            {
              // iterate for all OTHER substances
              for ( std::vector<std::string>::const_iterator
                    cj=substances.begin(); cj!=substances.end(); cj++ )
                {
                  if ( cj == ci ) continue;
                  // access the BioDynaMo diffusion grid
                  auto* dg_other = rm->GetDiffusionGrid(*cj);
                  const std::string& BC_other_name = dg_other->GetContinuumName(); // biochemical name
                  //
                  const double concentration_other = dg_other->GetValue(this->GetPosition()),
                               threshold_other = this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/threshold");
                  //
                  if ( ( threshold_other > 0.0 && concentration_other > +threshold_other ) ||
                       ( threshold_other < 0.0 && concentration_other < -threshold_other ) )
                    {
                      if (! this->params()->have_parameter<double>(CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                        return true;
                      else if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                        return true;
                    }
                  //...end of other substances loop
                }
            }
          // ...if no probability is provided by user, then cell simply dies!
          else if (! this->params()->have_parameter<double>(CP_name+"/can_apoptose/"+BC_name+"/probability"))
            {
              // since cell has apoptosed, then it must be removed from simulation
              return true;
            }
          // ...otherwise, check the likelihood for cell apoptosis ;)
          else if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_apoptose/"+BC_name+"/probability"))
            {
              // since cell has apoptosed, then it must be removed from simulation
              return true;
            }
        }
      //...end of substances loop
    }
  // since cell has not been through apoptosis, then it can do other things
  return false;
  //...end of cell apoptosis
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckAfterApoptosis()
{
  if (!this->GetCanApoptose()) return false;
  // by design only viable (non-necrotic) cells could apoptose
  if (!this->GetPhenotype()) return false;
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  const int cell_time_window = this->params()->get<int>(CP_name+"/can_apoptose/time_window/to_delete");
  //
  // since cell has apoptosed (due to ageing), then it must be removed from simulation
  if (this->GetAge()>cell_time_window) return true;
  //
  return false;
  //...end of cell apoptosis inquiry
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckQuiescenceAfterDivision()
{
  if (!this->GetCanDivide()) return false;
  // by design only viable (non-necrotic) cells could go through this phase
  if (!this->GetPhenotype()) return false;
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  const int n_div = this->GetNumberOfDivisions();
  if (n_div==0) return false;
  //
  if (! this->params()->have_parameter<int>(CP_name+"/can_divide/quiescence/time_window"))
    return false;
  const int cell_quiescence = this->params()->get<int>(CP_name+"/can_divide/quiescence/time_window");
  //
  if (this->GetAge()<cell_quiescence) return true;
  //
  return false;
  //...end of cell quiescence after division
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckMigration()
{
  if (!this->GetCanMigrate()) return false;
  // by design only viable (non-necrotic) cells could migrate
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  bool has_migrated = false;
  //
  this->passive_displacement_ = {0.0, 0.0, 0.0};
  //
  // check if cell migrates passively due to some convective field
  if ( this->params()->have_parameter<std::string>("convection/dynamic/from_file") )
    {
      const int SpaceDimension = this->params()->get<bool>("simulation_domain_is_2D")
                               ? 2 : 3;
      const double time_step = this->params()->get<double>("time_step");
      const double max_adhesion =
        this->params()->get<double>(CP_name+"/can_migrate/max_adhesion/displacement");
      //
      bdm::Double3 dvec = {0.0, 0.0, 0.0};
      // ensure cell is well within the simulation domain!
      if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
      // iterate for all components of the convection (vector) field
      for (int ispdm=0; ispdm<SpaceDimension; ispdm++)
        {
          const std::string name = "convection_" + std::to_string(ispdm);
          // access the BioDynaMo diffusion grid
          auto* dg = rm->GetDiffusionGrid(name);
          // obtain the convection component
          const double velocity_comp = dg->GetValue(this->GetPosition());
          // calculate corresponding displacement component
          const double displacement_comp = velocity_comp * time_step;
          //
          dvec[ispdm] += displacement_comp;
        }
      //
      const double d_magn = L2norm(dvec);
      // check if distance covered is above a minimum, else ignore
      if (d_magn > this->params()->get<double>("migration_tolerance"))
      // check if convection contribution is significant compared to the cell adhesion property
      if (d_magn > max_adhesion)
        {
          // update the (cell) displacement vector
          this->passive_displacement_ += dvec;
          // scale the (cell) displacement vector wrt the adhesion effect
          this->passive_displacement_ *= ((d_magn-max_adhesion)/d_magn);
          // update this flag
          has_migrated = true;
        }
    }
  //
  // check if cell migrates actively due to some inherent random-walk or
  // a biochemical stimulus, i.e. chemotaxis
  if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_migrate/probability"))
    {/// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ ///
      //
      const int index_time = this->params()->get<int>("index time");
      const int strength_time = this->params()->have_parameter<int>(CP_name+"/can_migrate/strength_of_time")
                              ? this->params()->get<int>(CP_name+"/can_migrate/strength_of_time")
                              : 1;
      //
      if ((index_time-1)%strength_time)
        { //  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
          //
          // force this flag
          has_migrated = true;
          //
        } //   --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      else
        { //   --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
          //
          this->active_displacement_ = {0.0, 0.0, 0.0};
          //
          // Brownian cell motion
          if (this->params()->get<double>(CP_name+"/can_migrate/half_range") > 0.0)
            {
              const double hr = this->params()->get<double>(CP_name+"/can_migrate/half_range");
              // produce the displacement vector
              bdm::Double3 dvec = {0.0, 0.0, 0.0};
              dvec[0] = rg->Uniform(-hr,+hr);
              dvec[1] = rg->Uniform(-hr,+hr);
              dvec[2] = this->params()->get<bool>("simulation_domain_is_2D")
                      ? 0.0 : rg->Uniform(-hr,+hr);
              //
              const double d_magn = L2norm(dvec);
              // check if distance covered is above a minimum, else ignore
              if (d_magn > this->params()->get<double>("migration_tolerance"))
                {
                  // update the (cell) displacement vector
                  this->active_displacement_ += dvec;
                  // update this flag
                  has_migrated = true;
                }
            }
          // chemotactic cell motion
          const std::vector<std::string>& substances =
            this->params()->get<std::vector<std::string>>("substances");
          // ensure cell is well within the simulation domain!
          if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
            // iterate for all substances
            for ( std::vector<std::string>::const_iterator
                  ci=substances.begin(); ci!=substances.end(); ci++ )
              {
                // access the BioDynaMo diffusion grid
                auto* dg = rm->GetDiffusionGrid(*ci);
                const std::string& BC_name = dg->GetContinuumName(); // biochemical name
                //
                if (! this->params()->have_parameter<double>(CP_name+"/can_migrate/chemotaxis/"+BC_name))
                  continue;
                //
                const double chemotaxis = this->params()->get<double>(CP_name+"/can_migrate/chemotaxis/"+BC_name);
                if (! chemotaxis) continue;
                //
                const double concentration = dg->GetValue(this->GetPosition()),
                             threshold = this->params()->get<double>(CP_name+"/can_migrate/chemotaxis/"+BC_name+"/threshold");
                //
                if ( ( threshold > 0.0 && concentration > +threshold ) ||
                     ( threshold < 0.0 && concentration < -threshold ) )
                {
                  if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_migrate/chemotaxis/"+BC_name+"/probability"))
                    {
                      bdm::Double3 dvec = {0.0, 0.0, 0.0};
                      int n_random_point = 0;
                      for (int random_point=0; random_point<20; random_point++)
                        {
                          const double radius      = rg->Uniform(0.0,1.5) * this->GetDiameter(),
                                       inclination = this->params()->get<bool>("simulation_domain_is_2D")
                                                   ? 0.5*bdm::Math::kPi : rg->Uniform(0.0,bdm::Math::kPi),
                                       azimuth     = rg->Uniform(0.0,2.0*bdm::Math::kPi);
                          const double x = radius * sin(inclination) * cos(azimuth),
                                       y = radius * sin(inclination) * sin(azimuth),
                                       z = radius * cos(inclination);
                          bdm::Double3 point = {x,y,z};
                          point += this->GetPosition();
                          // check spatial coordinates
                          if (point[0]<minCOORD||point[1]<minCOORD||point[2]<minCOORD||
                              point[0]>maxCOORD||point[1]>maxCOORD||point[2]>maxCOORD) continue;
                          //
                          bdm::Double3 gradS;
                          dg->GetGradient(point, &gradS);
                          //
                          if (L2norm(gradS)<=1.0e-6) continue;
                          //
                          dvec += gradS;
                          ++n_random_point; // increment this index
                        }
                      // average out the space vector
                      if (n_random_point) dvec /= n_random_point;
                      // scale gradient vector accordingly
                      dvec *= chemotaxis;
                      //
                      const double d_magn = L2norm(dvec);
                      // check if distance covered is above a minimum, else ignore
                      if (d_magn > this->params()->get<double>("migration_tolerance"))
                        {
                          // update the (cell) displacement vector
                          this->active_displacement_ += dvec;
                          // update this flag
                          has_migrated = true;
                          //
                          // check if to allow moving any further or skip any potential migration
                          if (! this->params()->get<bool>(CP_name+"/can_migrate/accumulate_path"))
                            break;
                        }
                      // ...end of propability check
                    }
                }
                //...end of substances loop
              }
          //
        } //   --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
      //
    }/// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ /// \\\ ///
  //
  // check if cell has migrated, if so then revise its spatial coordinates and trail
  if ( has_migrated )
    {
      this->UpdatePosition(this->GetDisplacement());
      this->UpdateTrail(L2norm(this->GetDisplacement()));
    }
  //
  // check if cell has migrated, if so then revise the spatial coordinates
  // of the cell protrusions; however, check for current algorithmic limitations!!!
  if ( has_migrated && this->GetNumberOfProtrusions() )
  // current implementation assumes that cell migration cannot be applied when
  // a cell has produced protrusions that have created sprouts or/and branches
  {
    if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
      ABORT_("an internal error occurred");

    // iterate for all (existing) protrusions of this cell
    for (int p=0; p<this->GetNumberOfProtrusions(); p++)
    {
      auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
      // only a cell with filodia can be displaced (together with its migrating cell)
      if ( protrusion->IsTerminal() )
        {
          const bdm::Double3 xyz = protrusion->GetPosition();
          protrusion->SetPosition(xyz+this->GetDisplacement());
        }
      else
        ABORT_("an exception is caught");
    }
    // completed the cell protrusion displacement task
  }
  //
  return has_migrated;
  //...end of cell migration
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckTransformation()
{
  if (!this->GetCanTransform()) return false;
  // by design only viable (non-necrotic) cells could transform
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_transform/probability"))
    return false;
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
    // iterate for all substances
    for ( std::vector<std::string>::const_iterator
          ci=substances.begin(); ci!=substances.end(); ci++ )
      {
        // access the BioDynaMo diffusion grid
        auto* dg = rm->GetDiffusionGrid(*ci);
        const std::string BC_name = dg->GetContinuumName(); // biochemical name
        //
        if (! this->params()->have_parameter<int>(CP_name+"/can_transform/"+BC_name+"/new_phenotype"))
          continue;
        //
        const double concentration = dg->GetValue(this->GetPosition()),
                     threshold = this->params()->get<double>(CP_name+"/can_transform/"+BC_name+"/threshold");
        //
        if ( ( threshold > 0.0 && concentration > +threshold ) ||
             ( threshold < 0.0 && concentration < -threshold ) )
          {
            if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_transform/"+BC_name+"/probability"))
              // now check if cell age is within an appropriate time-window to allow its transformation
              if ( this->GetAge() >= this->params()->get<int>(CP_name+"/can_transform/"+BC_name+"/time_window_open" ) &&
                   this->GetAge() <= this->params()->get<int>(CP_name+"/can_transform/"+BC_name+"/time_window_close") )
                  {
                    const int new_phenotype = this->params()->get<int>(CP_name+"/can_transform/"+BC_name+"/new_phenotype");
                    // firstly, the cell transforms
                    this->SetPhenotype(new_phenotype);
                    // increment this index
                    this->IncrementNumberOfTrasformations();
                    // now reset the age of the cell
                    this->SetAge();
                    //
                    const std::string CP_new_name = // cell phenotype name
                      this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
                    // principal directions of the cell polarization matrix
                    double p0, p1, p2;
                    if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                      {
                        p0 = this->params()->get<double>(CP_new_name+"/principal/0");
                        p1 = this->params()->get<double>(CP_new_name+"/principal/1");
                        p2 = this->params()->get<double>(CP_new_name+"/principal/2");
                      }
                    //
                    this->SetCanApoptose(this->params()->get<bool>(CP_new_name+"/can_apoptose"));
                    this->SetCanGrow(this->params()->get<bool>(CP_new_name+"/can_grow"));
                    this->SetCanDivide(this->params()->get<bool>(CP_new_name+"/can_divide"));
                    this->SetCanMigrate(this->params()->get<bool>(CP_new_name+"/can_migrate"));
                    this->SetCanTransform(this->params()->get<bool>(CP_new_name+"/can_transform"));
                    this->SetCanProtrude(this->params()->get<bool>(CP_new_name+"/can_protrude"));
                    this->SetCanPolarize(this->params()->get<bool>(CP_new_name+"/can_polarize"));
                    // reset the cell polarization matrix
                    if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                      this->SetPolarization(diag(p0, p1, p2));
                    // reset the cell protrusion phenotype
                    if ( this->GetNumberOfProtrusions() )
                      {
                        if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
                          ABORT_("an internal error occurred");
                        //
                        // iterate for all (existing) protrusions of this cell
                        for (int p=0; p<this->GetNumberOfProtrusions(); p++)
                          {
                            auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
                            // assign this cell (that is associated with) to the protrusion created
                            protrusion->SetCell(this);
                          }
                      }
                    // reset the cell behavior (mechanisms order) from old to new one
                    {
                      const bdm::InlineVector<bdm::Behavior*,2>& behavior = this->GetAllBehaviors();
                      if (behavior.size()!=1)
                        ABORT_("an internal error occurred");
                      //
                      this->RemoveBehavior(behavior[0]);
                    }
                    const int mo = this->params()->get<int>(CP_new_name+"/mechanism_order");
                    if (10==mo)
                      this->AddBehavior(new Biology4BiologicalCell_10());
                    else
                      ABORT_("an exception is caught");
                    // cell has transformed, then proceed to check if it can do other things
                    return true;
                  }
              //
            //
          }
        //...end of substances loop
      }
  // cell has not been through any transformation
  return false;
  //...end of cell transformation
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckPolarization()
{
  if (!this->GetCanPolarize()) return false;
  // by design only viable (non-necrotic) cells could polarize
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // current implementation assumes that cell polarization cannot be applied when
  // a cell has produced protrusions
  ASSERT_(0==this->GetNumberOfProtrusions(),"an internal error occurred");
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_polarize/probability"))
    return false;
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  const double principal_min = this->params()->get<double>(CP_name+"/principal/min"),
               principal_max = this->params()->get<double>(CP_name+"/principal/max");
  const std::vector<int> perm =
    this->params()->get<std::vector<int>>(CP_name+"/principal/permutation");
  const bdm::Double3 principal = { this->params()->get<double>(CP_name+"/principal/"+std::to_string(perm[0])) ,
                              this->params()->get<double>(CP_name+"/principal/"+std::to_string(perm[1])) ,
                              this->params()->get<double>(CP_name+"/principal/"+std::to_string(perm[2])) };
  //
  if ( this->params()->get<bool>(CP_name+"/can_polarize/migration/dependency") )
    {
      const int pattern = this->params()->have_parameter<int>(CP_name+"/can_polarize/migration/pattern")
                        ? this->params()->get<int>(CP_name+"/can_polarize/migration/pattern") : 0;
      // check if cell displacement is considerable to allow the cell self-polarization
      if (L2norm(this->GetDisplacement()) > this->params()->get<double>("migration_tolerance"))
        {
          bdm::Double3 v = this->GetDisplacement(); // cell displacement (vector) field
          // eigenvectors
          bdm::Double3 n0 = this->params()->get<bool>("simulation_domain_is_2D")
                     ? bdm::Double3{v[0], v[1], 0.0}
                     : bdm::Double3{v[0], v[1], v[2]};
          bdm::Double3 n1 = this->params()->get<bool>("simulation_domain_is_2D")
                     ? bdm::Double3{rg->Uniform(-1.,1.), rg->Uniform(-1.,1.), 0.0}
                     : bdm::Double3{rg->Uniform(-1.,1.), rg->Uniform(-1.,1.), rg->Uniform(-1.,1.)};
          bdm::Double3 n2 = cross(n0, n1);
          n1 = cross(n2, n0);
          // normalize all vectors
          n0.Normalize();
          n1.Normalize();
          n2.Normalize();
          // eigenvalues
          double p[3];
          if (0==pattern)
            {
              p[0] = rg->Uniform(principal[0],  principal_max),
              p[1] = rg->Uniform(principal[1],  principal[0]),
              p[2] = rg->Uniform(principal_min, principal[2]);
            }
          else if (1==pattern)
            {
              p[0] = rg->Uniform(principal[0],  principal_max),
              p[1] = rg->Uniform(principal[2],  principal[1]),
              p[2] = rg->Uniform(principal_min, principal[2]);
            }
          else
            ABORT_("an exception is caught");
          // tensor products of eigenvectors, scaled by respective eigenvalues
          const bdm::Double3x3 n0Xn0_p0 = tensor(n0, n0, p[0]),
                               n1Xn1_p1 = tensor(n1, n1, p[1]),
                               n2Xn2_p2 = tensor(n2, n2, p[2]);
          // update the cell polarization matrix
          this->polarize_ = n0Xn0_p0 + n1Xn1_p1 + n2Xn2_p2;
          // cell has polarized, then proceed to check if it can do other things
          return true;
        }
    }
  //
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
    // iterate for all substances if cell can
    // re-orient / polarize
    for ( std::vector<std::string>::const_iterator
          ci=substances.begin(); ci!=substances.end(); ci++ )
      {
        // access the BioDynaMo diffusion grid
        auto* dg = rm->GetDiffusionGrid(*ci);
        const std::string BC_name = dg->GetContinuumName(); // biochemical name
        const double concentration = dg->GetValue(this->GetPosition()),
                     threshold = this->params()->get<double>(CP_name+"/can_polarize/"+BC_name+"/threshold");
        //
        if ( ( threshold > 0.0 && concentration > +threshold ) ||
             ( threshold < 0.0 && concentration < -threshold ) )
          {
            bdm::Double3 v; // substance gradient (vector) field
            dg->GetGradient(this->GetPosition(), &v);
            // eigenvectors
            bdm::Double3 n0 = this->params()->get<bool>("simulation_domain_is_2D")
                       ? bdm::Double3{v[0], v[1], 0.0}
                       : bdm::Double3{v[0], v[1], v[2]};
            bdm::Double3 n1 = this->params()->get<bool>("simulation_domain_is_2D")
                       ? bdm::Double3{rg->Uniform(-1.,1.), rg->Uniform(-1.,1.), 0.0}
                       : bdm::Double3{rg->Uniform(-1.,1.), rg->Uniform(-1.,1.), rg->Uniform(-1.,1.)};
            bdm::Double3 n2 = cross(n0, n1);
            n1 = cross(n2, n0);
            // normalize all direction vectors
            n0.Normalize();
            n1.Normalize();
            n2.Normalize();
            // eigenvalues
            const double p0 = rg->Uniform(principal[0],  principal_max),
                         p1 = rg->Uniform(principal[1],  p0),
                         p2 = rg->Uniform(principal_min, principal[2]);
            // tensor products of eigenvectors, scaled by respective eigenvalues
            const bdm::Double3x3 n0Xn0_p0 = tensor(n0, n0, p0),
                                 n1Xn1_p1 = tensor(n1, n1, p1),
                                 n2Xn2_p2 = tensor(n2, n2, p2);
            // update the cell polarization matrix
            this->polarize_ = n0Xn0_p0 + n1Xn1_p1 + n2Xn2_p2;
            // cell has polarized, then proceed to check if it can do other things
            return true;
          }
        //...end of substances loop
      }
  // cell has not been through any polarization
  return false;
  //...end of cell polarization
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckProtrusion()
{
  if (!this->GetCanProtrude()) return false;
  // by design only viable (non-necrotic) cells could develop protrusions
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  if (! this->params()->have_parameter<int>(CP_name+"/can_protrude/pattern"))
    return false;
  //
  // confirn if permissible number of protrusions have been createdon the cell
  const int max_protrusions = this->params()->get<int>(CP_name+"/can_protrude/max_protrusions");
  if (max_protrusions == this->GetNumberOfProtrusions())
    return false;
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_protrude/probability"))
    return false;
  // now check if cell age is within an appropriate time-window to allow protrusion generation
  if ( this->GetAge() < this->params()->get<int>(CP_name+"/can_protrude/time_window_open" ) ||
       this->GetAge() > this->params()->get<int>(CP_name+"/can_protrude/time_window_close") )
    return false;
  //
  const double diameter = this->GetDiameter(),
               diameter_cutoff = this->params()->get<double>(CP_name+"/can_protrude/diameter_cutoff");
  //
  if ( diameter < diameter_cutoff )
    return false;
  //
  bool stimulate = false;
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      const double concentration = dg->GetValue(this->GetPosition()),
                   threshold = this->params()->get<double>(CP_name+"/can_protrude/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          stimulate = true;
        }
      //...end of substances loop
    }
  // check if any biochemical stimulus leads towards cell protrusion generation
  if ( ! stimulate ) return false;
  //
  const double protrusion_tol =
    ! this->params()->have_parameter<double>(CP_name+"/can_protrude/tolerance")
    ? tol :      this->params()->get<double>(CP_name+"/can_protrude/tolerance");
  if (!check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), protrusion_tol))
    return false;
  //
  CellProtrusion c_p;
  //
  const int pattern = this->params()->get<int>(CP_name+"/can_protrude/pattern");
  // cell protrusion max diameter
  const double dia_max = this->params()->get<double>(CP_name+"/can_protrude/diameter/max");
  // identify mode of simulation domain
  if (this->params()->get<bool>("simulation_domain_is_2D"))
    {
      if (0==pattern)
        {
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              bdm::Double3 axis = {rg->Uniform(-1.0,+1.0), rg->Uniform(-1.0,+1.0), 0.0};
              // if protrusion is not valid the simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis)) {
                ip--;
                continue;
              }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                if (behavior.size()!=1)
                  ABORT_("an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                  const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                  protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else if (-1==pattern || +1==pattern)
        {
          //
          const std::string& chemo_substance = this->params()->get<std::string>(CP_name+"/can_protrude/chemotaxis");
          // iterate for all substances
          std::vector<std::string>::const_iterator
            ci = std::find(substances.begin(), substances.end(), chemo_substance);
          if (substances.end()==ci)
            ABORT_("an internal error occurred");
          // access the BioDynaMo diffusion grid
          auto* dg = rm->GetDiffusionGrid(*ci);
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              const double radius = rg->Uniform(0.0,this->GetDiameter()),
                           phi    = rg->Uniform(0.0,2.0*bdm::Math::kPi);
              bdm::Double3 pnt = {radius*cos(phi), radius*sin(phi), 0.0};
              pnt += this->GetPosition();
              //
              bdm::Double3 axis;
              dg->GetGradient(pnt, &axis);
              // note these details...
              axis *= pattern;
              if (!normalize(axis, axis)) continue;
              // if protrusion is not valid then simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis))
                {
                  // // // // // // //ip--;
                  continue;
                  //...well skip this for the time being and check again some other time,
                  // in order to avoid a perpetual loop at this point; however, please do
                  // fix this issue in the future!
                }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                    const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                    protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else if (+2==pattern)
        {
          //
          ABORT_("an internal error occurred");
          //
        }
      else
        ABORT_("an exception is caught");
    }
  else
    {
      if (0==pattern)
        {
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              bdm::Double3 axis = {rg->Uniform(-1.0,+1.0), rg->Uniform(-1.0,+1.0), rg->Uniform(-1.0,+1.0)};
              // if protrusion is not valid the simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis))
                {
                  ip--;
                  continue;
                }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                  const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                  protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else if (-1==pattern || +1==pattern)
        {
          //
          const std::string& chemo_substance = this->params()->get<std::string>(CP_name+"/can_protrude/chemotaxis");
          // iterate for all substances
          std::vector<std::string>::const_iterator
            ci = std::find(substances.begin(), substances.end(), chemo_substance);
          ASSERT_(substances.end()!=ci,"an internal error occurred");
          //
          // access the BioDynaMo diffusion grid
          auto* dg = rm->GetDiffusionGrid(*ci);
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              const double radius = rg->Uniform(0.0,this->GetDiameter()),
                           phi    = rg->Uniform(0.0,2.0*bdm::Math::kPi),
                           theta  = rg->Uniform(0.0,bdm::Math::kPi);
              bdm::Double3 pnt = {radius*cos(phi)*sin(theta), radius*sin(phi)*sin(theta), radius*cos(theta)};
              pnt += this->GetPosition();
              //
              bdm::Double3 axis;
              dg->GetGradient(pnt, &axis);
              // note these details...
              axis *= pattern;
              if (!normalize(axis, axis)) continue;
              // if protrusion is not valid then simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis))
                {
                  // // // // // // //ip--;
                  continue;
                  //...well skip this for the time being and check again some other time,
                  // in order to avoid a perpetual loop at this point; however, please do
                  // fix this issue in the future!
                }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                  const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                  protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else if (+2==pattern)
        {
          // iterate for the max. number of permisible cell protrusions
          for (int ip=0; ip<max_protrusions; ip++)
            {
              const double phi   = rg->Uniform()*bdm::Math::kPi*2.0,
                           theta = rg->Uniform()*bdm::Math::kPi;
              const bdm::Double3 axis = { cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta) };
              // if protrusion is not valid then simply redo the computations again...
              if (!this->CheckProtrusionAxis(axis))
                {
                  ip--;
                  continue;
                  //...well skip this for the time being and check again some other time,
                  // in order to avoid a perpetual loop at this point; however, please do
                  // fix this issue in the future!
                }
              //
              auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->ExtendNewNeurite(axis, &c_p));
              // assign this cell (that is associated with) to the protrusion created
              protrusion->SetCell(this);
              // setup the protrusion diameter
              protrusion->SetDiameter(dia_max);
              // ...and assign the pointer to the list of model parameters
              protrusion->SetParametersPointer(this->params());
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = protrusion->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                protrusion->RemoveBehavior(behavior[0]);
              }
              protrusion->AddBehavior(new Biology4CellProtrusion());
              // ensure the length of the original protrusion is correct
              if (1.0>dia_max)
                {
                  const double L_rate = (1.0-dia_max) / this->params()->get<double>("time_step");
                  protrusion->RetractTerminalEnd(L_rate);
                }
              else
                {
                  const double L_rate = (dia_max-1.0) / this->params()->get<double>("time_step");
                  protrusion->ElongateTerminalEnd(L_rate, axis);
                }
              // increment this index
              this->IncrementNumberOfProtrusions();
            }
          //
        }
      else
        ABORT_("an exception is caught");
    }
  //
  if ( this->GetNumberOfProtrusions() )
    {
      this->can_grow_      = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_grow");
      this->can_divide_    = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_divide");
      this->can_migrate_   = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_migrate");
      this->can_transform_ = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_transform");
      this->can_polarize_  = this->params()->get<bool>(CP_name+"/can_protrude/afterwards/can_polarize");
    }
  // cell has developed protrusions (filodia or/and neurites)
  return true;
  //...end of cell protrusion
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckGrowth()
{
  if (!this->GetCanGrow()) return false;
  // by design only viable (non-necrotic) cells could grow
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  const double diameter = this->GetDiameter(),
               diameter_min = this->params()->get<double>(CP_name+"/diameter/min"),
               diameter_max = this->params()->get<double>(CP_name+"/diameter/max");
  //
  // check if cell size is within acceptable (upper limit) before proceeding with growth
  if (diameter >= diameter_max)
    return false;
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_grow/probability"))
    return false;
  //
  double volume_rate = 0.0;
  {
    const double diameter_rate = this->params()->get<double>(CP_name+"/can_grow/diameter_rate");
    //
    volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
  }
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      const double concentration = dg->GetValue(this->GetPosition()),
                   threshold = this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          const double diameter_rate = this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/diameter_rate");
          //
          // allow cell growth controlled by a combination of two biochemical cues, therefore
          // cell development is dependent from another substance as well
          if (this->params()->get<bool>(CP_name+"/can_grow/"+BC_name+"/dependency"))
            {
              // iterate for all OTHER substances
              for ( std::vector<std::string>::const_iterator
                    cj=substances.begin(); cj!=substances.end(); cj++ )
                {
                  if ( cj == ci ) continue;
                  // access the BioDynaMo diffusion grid
                  auto* dg_other = rm->GetDiffusionGrid(*cj);
                  const std::string& BC_other_name = dg_other->GetContinuumName(); // biochemical name
                  //
                  const double concentration_other = dg_other->GetValue(this->GetPosition()),
                               threshold_other = this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/dependency/"+BC_other_name+"/threshold");
                  //
                  if ( ( threshold_other > 0.0 && concentration_other > +threshold_other ) ||
                       ( threshold_other < 0.0 && concentration_other < -threshold_other ) )
                    {
                      if (! this->params()->have_parameter<double>(CP_name+"/can_grow/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                        {
                          volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
                        }
                      else if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/dependency/"+BC_other_name+"/probability"))
                        {
                          volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
                        }
                    }
                  //...end of other substances loop
                }
            // ...if no probability is provided by user, then cell simply grows!
            }
          else if (! this->params()->have_parameter<double>(CP_name+"/can_grow/"+BC_name+"/probability"))
            {
              volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
            }
          // ...otherwise, check the likelihood for cell growth ;)
          else if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_grow/"+BC_name+"/probability"))
            {
              volume_rate += (0.5*TMath::Pi())*diameter_rate*pow2(diameter);
            }
        }
      //...end of substances loop
    }
  //
  if ( ( diameter < diameter_max && volume_rate ) ||
       ( diameter < diameter_min && volume_rate > 0.0 ) )
    {
      this->ChangeVolume(volume_rate);
      // cell has grown, then proceed to check if it can do other things
      return true;
    }
  // cell has not been through any growth
  return false;
  //...end of cell growth
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckTransformationAndDivision()
{
  if (!this->GetCanDivide()) return false;
  // by design only viable (non-necrotic) cells could divide after they transform
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // cell cannot divide (not at least with current BioDynaMo implementation)
  // if it has developed protrusions (filopodia or/and neurites)
  // however, in division followed by transformation, cell protrusion phenotype
  // has been implemented to switch into new one (after transformation)
  ASSERT_(0==this->GetNumberOfProtrusions(),"an internal error occurred");
  //
  const int n_div = this->GetNumberOfDivisions();
  //
  // call cannot divide more times than it should
  if (n_div >= this->params()->get<int>(CP_name+"/can_divide/max"))
    return false;
  //
  if (this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age")>0.0)
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_divide/probability")
                                +this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age")
                                *this->GetAge() )
        return false;
    }
  else
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_divide/probability"))
        return false;
    }
  //
  const double diameter = this->GetDiameter(),
               diameter_cutoff = this->params()->get<double>(CP_name+"/can_divide/diameter_cutoff");
  const int cell_maturity = (n_div+1) * this->params()->get<int>(CP_name+"/can_divide/time_window");
  //
  if ( diameter < diameter_cutoff || this->GetAge() < cell_maturity )
    return false;
  //
  // produce the separation vector
  const bdm::Double3 axis =
    { rg->Uniform(-1.0,+1.0) ,
      rg->Uniform(-1.0,+1.0) ,
      (this->params()->get<bool>("simulation_domain_is_2D") ? 0.0 : rg->Uniform(-1.0,+1.0)) };
  //
  const double volume_ratio = rg->Uniform(0.9,1.1);
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances if cell can
  // transform and then divide (symmetrically)
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      if (! this->params()->have_parameter<int>(CP_name+"/can_transform_and_divide/"+BC_name+"/new_phenotype"))
        continue;
      //
      const double concentration = dg->GetValue(this->GetPosition()),
                   threshold = this->params()->get<double>(CP_name+"/can_transform_and_divide/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_transform_and_divide/"+BC_name+"/probability"))
            {
              const int new_phenotype = this->params()->get<int>(CP_name+"/can_transform_and_divide/"+BC_name+"/new_phenotype");
              // firstly, the cell transforms
              this->SetPhenotype(new_phenotype);
              // now reset the age of the cell
              this->SetAge();
              // increment this index
              this->IncrementNumberOfTrasformations();
              //
              const std::string CP_new_name = // cell phenotype name
                this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
              // principal directions of the cell polarization matrix
              double p0, p1, p2;
              if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                {
                  p0 = this->params()->get<double>(CP_new_name+"/principal/0");
                  p1 = this->params()->get<double>(CP_new_name+"/principal/1");
                  p2 = this->params()->get<double>(CP_new_name+"/principal/2");
                }
              //
              this->SetCanApoptose(this->params()->get<bool>(CP_new_name+"/can_apoptose"));
              this->SetCanGrow(this->params()->get<bool>(CP_new_name+"/can_grow"));
              this->SetCanDivide(this->params()->get<bool>(CP_new_name+"/can_divide"));
              this->SetCanMigrate(this->params()->get<bool>(CP_new_name+"/can_migrate"));
              this->SetCanTransform(this->params()->get<bool>(CP_new_name+"/can_transform"));
              this->SetCanProtrude(this->params()->get<bool>(CP_new_name+"/can_protrude"));
              this->SetCanPolarize(this->params()->get<bool>(CP_new_name+"/can_polarize"));
              // reset the cell polarization matrix
              if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                this->SetPolarization(diag(p0, p1, p2));
              // reset the cell protrusion phenotype
              if ( this->GetNumberOfProtrusions() )
                {
                  if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
                    ABORT_("an internal error occurred");
                  //
                  // iterate for all (existing) protrusions of this cell
                  for (int p=0; p<this->GetNumberOfProtrusions(); p++)
                    {
                      auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
                      // assign this cell (that is associated with) to the protrusion created
                      protrusion->SetCell(this);
                    }
                }
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = this->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an exception is caught");
                //
                this->RemoveBehavior(behavior[0]);
              }
              const int mo = this->params()->get<int>(CP_new_name+"/mechanism_order");
              if (10==mo)
                this->AddBehavior(new Biology4BiologicalCell_10());
              else
                ABORT_("an exception is caught");
              // secondly, the cell divides
              this->Divide(volume_ratio, axis);
              // cell has transformed and divided, then proceed to check if it can do other things
              return true;
            }
        }
      //...end of substances loop
    }
  // cell has not been through any division
  return false;
  //...end of cell division
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckAsymmetricDivision()
{
  if (!this->GetCanDivide()) return false;
  // by design only viable (non-necrotic) cells could divide and then transform
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // cell cannot divide (not at least with current BioDynaMo implementation)
  // if it has developed protrusions (filopodia or/and neurites)
  // however, in division followed by transformation, cell protrusion phenotype
  // has been implemented to switch into new one (after transformation)
  ASSERT_(0==this->GetNumberOfProtrusions(),"an internal error occurred");
  //
  const int n_div = this->GetNumberOfDivisions();
  //
  // call cannot divide more times than it should
  if (n_div >= this->params()->get<int>(CP_name+"/can_divide/max"))
    return false;
  //
  if (this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age")>0.0)
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_divide/probability")
                                +this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age")
                                *this->GetAge() )
        return false;
    }
  else
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_divide/probability"))
        return false;
    }
  //
  const double diameter = this->GetDiameter(),
               diameter_cutoff = this->params()->get<double>(CP_name+"/can_divide/diameter_cutoff");
  const int cell_maturity = (n_div+1) * this->params()->get<int>(CP_name+"/can_divide/time_window");
  //
  if ( diameter < diameter_cutoff || this->GetAge() < cell_maturity )
    return false;
  //
  // produce the separation vector
  const bdm::Double3 axis =
    { rg->Uniform(-1.0,+1.0) ,
      rg->Uniform(-1.0,+1.0) ,
      (this->params()->get<bool>("simulation_domain_is_2D") ? 0.0 : rg->Uniform(-1.0,+1.0)) };
  //
  const double volume_ratio = rg->Uniform(0.9,1.1);
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances if cell can
  // divide and then transform
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      const double concentration = dg->GetValue(this->GetPosition()),
                   threshold = this->params()->get<double>(CP_name+"/can_divide_and_transform/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          if (rg->Uniform(0.0,1.0) <= this->params()->get<double>(CP_name+"/can_divide_and_transform/"+BC_name+"/probability"))
            {
              const int new_phenotype = this->params()->get<int>(CP_name+"/can_divide_and_transform/"+BC_name+"/new_phenotype");
              // firstly, the cell divides
              this->Divide(volume_ratio, axis);
              // secondly, the cell transforms
              this->SetPhenotype(new_phenotype);
              // now reset the age of the cell
              this->SetAge();
              // increment this index
              this->IncrementNumberOfTrasformations();
              //
              const std::string CP_new_name = // cell phenotype name
                this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
              // principal directions of the cell polarization matrix
              double p0, p1, p2;
              if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                {
                  p0 = this->params()->get<double>(CP_new_name+"/principal/0");
                  p1 = this->params()->get<double>(CP_new_name+"/principal/1");
                  p2 = this->params()->get<double>(CP_new_name+"/principal/2");
                }
              //
              this->SetCanApoptose(this->params()->get<bool>(CP_new_name+"/can_apoptose"));
              this->SetCanGrow(this->params()->get<bool>(CP_new_name+"/can_grow"));
              this->SetCanDivide(this->params()->get<bool>(CP_new_name+"/can_divide"));
              this->SetCanMigrate(this->params()->get<bool>(CP_new_name+"/can_migrate"));
              this->SetCanTransform(this->params()->get<bool>(CP_new_name+"/can_transform"));
              this->SetCanProtrude(this->params()->get<bool>(CP_new_name+"/can_protrude"));
              this->SetCanPolarize(this->params()->get<bool>(CP_new_name+"/can_polarize"));
              // reset the cell polarization matrix
              if (this->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                this->SetPolarization(diag(p0, p1, p2));
              // reset the cell protrusion phenotype
              if ( this->GetNumberOfProtrusions() )
                {
                  if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
                    ABORT_("an internal error occurred");
                  // iterate for all (existing) protrusions of this cell
                  for (int p=0; p<this->GetNumberOfProtrusions(); p++)
                    {
                      auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
                      // assign this cell (that is associated with) to the protrusion created
                      protrusion->SetCell(this);
                    }
                }
              // reset the cell behavior (mechanisms order) from old to new one
              {
                const bdm::InlineVector<bdm::Behavior*,2>& behavior = this->GetAllBehaviors();
                ASSERT_(1==behavior.size(),"an internal error occurred");
                //
                this->RemoveBehavior(behavior[0]);
              }
              const int mo = this->params()->get<int>(CP_new_name+"/mechanism_order");
              if (10==mo)
                this->AddBehavior(new Biology4BiologicalCell_10());
              else
                ABORT_("an exception is caught");
              // cell has divided and transformed, then proceed to check if it can do other things
              return true;
            }
        }
      //...end of substances loop
    }
  // cell has not been through any division
  return false;
  //...end of cell division
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckDivision() {
  if (!this->GetCanDivide()) return false;
  // by design only viable (non-necrotic) cells could divide
  if (!this->GetPhenotype()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // cell cannot divide (not at least with current BioDynaMo implementation)
  // if it has developed protrusions (filopodia or/and neurites)
  // however, in division followed by transformation, cell protrusion phenotype
  // has been implemented to switch into new one (after transformation)
  ASSERT_(0==this->GetNumberOfProtrusions(),"an internal error occurred");
  //
  const int n_div = this->GetNumberOfDivisions();
  //
  // call cannot divide more times than it should
  if (n_div >= this->params()->get<int>(CP_name+"/can_divide/max"))
    return false;
  //
  if (this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age")>0.0)
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_divide/probability")
                                +this->params()->get<double>(CP_name+"/can_divide/probability_increment_with_age")
                                *this->GetAge() )
        return false;
    }
  else
    {
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_divide/probability"))
        return false;
    }
  //
  const double diameter = this->GetDiameter(),
               diameter_cutoff = this->params()->get<double>(CP_name+"/can_divide/diameter_cutoff");
  const int cell_maturity = (n_div+1) * this->params()->get<int>(CP_name+"/can_divide/time_window");
  //
  if ( diameter < diameter_cutoff || this->GetAge() < cell_maturity )
    return false;
  //
  // produce the separation vector
  const bdm::Double3 axis =
    { rg->Uniform(-1.0,+1.0) ,
      rg->Uniform(-1.0,+1.0) ,
      (this->params()->get<bool>("simulation_domain_is_2D") ? 0.0 : rg->Uniform(-1.0,+1.0)) };
  //
  const double volume_ratio = rg->Uniform(0.9,1.1);
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  // ensure cell is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances if cell can
  // divide (symmetrically)
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      const double concentration = dg->GetValue(this->GetPosition()),
                   threshold = this->params()->get<double>(CP_name+"/can_divide/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          // since no symmetric (prior to cell transformation) or unsymmetric division
          // has occurred, then cell divides conventionally
          this->Divide(volume_ratio, axis);
          // cell has divided, then proceed to check if it can do other things
          return true;
        }
      //...end of substances loop
    }
  // cell has not been through any division
  return false;
  //...end of cell division
}
// -----------------------------------------------------------------------------
inline
void bdm::BiologicalCell::CheckAndFixDiameter()
{
  if (!this->GetPhenotype()) return;
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  const double diameter = this->GetDiameter(),
               diameter_min = this->params()->get<double>(CP_name+"/diameter/min"),
               diameter_max = this->params()->get<double>(CP_name+"/diameter/max");
  // check if cell diameter is within user-defined bounds
  if      (diameter < diameter_min)
    this->SetDiameter(diameter_min);
  else if (diameter > diameter_max)
    this->SetDiameter(diameter_max);
}
// -----------------------------------------------------------------------------
inline
bool bdm::BiologicalCell::CheckProtrusionAxis(bdm::Double3 axis)
{
  ASSERT_(0!=this->GetPhenotype(),"an internal error occurred");
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetPhenotype()));
  //
  // relative angle between protrusions
  const double rel_angle_min = this->params()->get<double>(CP_name+"/can_protrude/relative_angle/degrees/min"),
               rel_angle_max = this->params()->get<double>(CP_name+"/can_protrude/relative_angle/degrees/max");
  // make sure to normalize the axis vector first
  ASSERT_(normalize(axis, axis),"could not normalize the axis vector");
  //
  if (!this->protrusions_.empty())
    {
      for (unsigned int l=0; l<this->protrusions_.size(); l++)
        {
          const bdm::Double3& current_axis = this->protrusions_[l];
          const double angle = radians_to_degrees( acos(axis*current_axis) );
          // check if relative angle is within user-defined range
          if (angle<rel_angle_min || angle>rel_angle_max)
            return false;
        }
    }
  //
  this->protrusions_.push_back(axis);
  return true;
}
// -----------------------------------------------------------------------------
inline
void bdm::BiologicalCell::Set2DeleteProtrusions()
{
  // sanity check
  if ( this->GetNumberOfProtrusions() != (int)this->daughters_.size() )
    ABORT_("an internal error occurred");
  // iterate for all (existing) protrusions of this cell
  for (int p=0; p<this->GetNumberOfProtrusions(); p++)
    {
      auto* protrusion = bdm::bdm_static_cast<CellProtrusion*>(this->daughters_[p].Get());
      // assign this cell (that is associated with) to the protrusion created
      protrusion->Set2Delete();
    }
}
// =============================================================================
#endif // _BIOLOGICAL_CELL_INLINE_H_
// =============================================================================

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
#ifndef _VESSEL_INLINE_H_
#define _VESSEL_INLINE_H_
// =============================================================================
inline
void bdm::Vessel::RunBiochemics()
{
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
  // vessel spatial coordinates
  const bdm::Double3 xyz = this->IsTerminal()
                         ? (this->GetPosition()+this->GetSpringAxis()*0.5) : this->GetPosition();
  // exit function if vessel resides at the edge of the simulation domain,
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
  // ensure vessel is well within the simulation domain!
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
        if (! this->params()->have_parameter<double>("vessel/"+BC_name+"/secretion/net_balance"))
          continue;
        //
        const double concentration = dg->GetValue(xyz);
        // parameters that modulate biochemical cue secretion (production or consumption)
        const double BC_stdev = this->params()->get<double>("vessel/"+BC_name+"/secretion/net_balance/std")<=0.0 ? 1.0 :
                                rg->Uniform(1.0-this->params()->get<double>("vessel/"+BC_name+"/secretion/net_balance/std"),
                                            1.0+this->params()->get<double>("vessel/"+BC_name+"/secretion/net_balance/std"));
        const double saturation = this->params()->get<double>("vessel/"+BC_name+"/secretion/saturation"),
                     net_balance = this->params()->get<double>("vessel/"+BC_name+"/secretion/net_balance") * BC_stdev;
        //
        // skip subsequent calculations if net balance of this biochemical cue secretion is
        // equal to absolute zero!!!
        if (0.0 == net_balance) continue;
        //
        // check if vessel is sprout or tip to secrete biochemical cue (substance)
        if ( this->params()->get<bool>("vessel/"+BC_name+"/secretion/sprout_not_tip_secretes") )
          {
            if (this->IsTerminal()) continue;
          }
        else
          {
            if (! this->IsTerminal()) continue;
          }
        //
        //
        if (! this->params()->get<bool>("vessel/"+BC_name+"/secretion/dependent"))
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
                if ( concentration > 0.0 ) {
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
                         threshold_other = this->params()->get<double>("vessel/"+BC_name+"/secretion/"+BC_other_name+"/threshold");
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
  //...end of vessel biochemics
}
// -----------------------------------------------------------------------------
inline
bool bdm::Vessel::CheckPositionValidity()
{
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary"),
               meanCOORD = (maxCOORD+minCOORD)/2.0,
               deltaCOORD = maxCOORD - meanCOORD;
  const double tol = this->params()->get<double>("domain_tolerance");
  // nodal points of the original vessel [start, end]
  bdm::Double3 node[2];
  node[0] = this->GetPosition() - this->GetSpringAxis() * 0.5;
  node[1] = this->GetPosition() + this->GetSpringAxis() * 0.5;
  //
  // check if simulation domain is bounded or unbounded (for vessel extension)
  if ( this->params()->get<bool>("simulation_domain_is_bounded") )
    {
      ABORT_("code is under development");
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
          // vessel went out of the (simulation) domain
          // therefore it must stop further developing anyhow
          return false;
        }
    }
  // access the pointer to parameter of the simulation obstacles object
  const SimulationObstacles* obstacles =
    this->params()->get<SimulationObstacles*>("simulation_obstacles");
  // check if vessel has reached any of the 'box' simulation obstacles
  for (size_t l=0; l<obstacles->box.size(); l++)
    {
      // original vessel position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      // it true, then vessel position needs to be corrected
      if ("box/inside"==obstacles->box[l].type)
        {
          // check if vessel position is outside this obstacle
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
          // enforce vessel to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          // vessel has remained inside the simulation domain
          return true;
        }
      else if ("box/outside"==obstacles->box[l].type)
        {
          // check if vessel position is inside this obstacle
          if (obstacles->box[l].is_inside(xyz))
            continue;
          // calculate the outward unit normal vector and the projection vector of
          // the vessel to each box face
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
              // enforce vessel to lie on the (box) obstacle surface
              this->SetPosition(xyz);
              // vessel has remained inside the simulation domain
              return true;
            }
        }
      else
        ABORT_("an exception is caught");
      // ...end of 'box' simulation obstacles loop
    }
  // check if vessel has reached any of the 'sphere' simulation obstacles
  for (size_t l=0; l<obstacles->sphere.size(); l++)
    {
      bdm::Double3 xyz = this->GetPosition();
      // it true, then vessel position needs to be corrected
      if ("sphere/inside"==obstacles->sphere[l].type)
        {
          // check if vessel position is outside this obstacle
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
          // enforce vessel to lie on the surface of the (spherical) obstacle
          this->SetPosition(xyz);
          // vessel has remained inside the simulation domain
          return true;
        }
      else if ("sphere/outside"==obstacles->sphere[l].type)
        {
          // check if vessel position is inside this obstacle
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
          // enforce vessel to lie on the surface of the (spherical) obstacle
          this->SetPosition(xyz);
          // vessel has remained inside the simulation domain
          return true;
        }
      else
        ABORT_("an exception is caught");
      // ...end of 'sphere' simulation obstacles loop
    }
  // check if vessel has reached any of the 'surface' simulation obstacles
  for (size_t l=0; l<obstacles->surface.size(); l++)
    {
      bdm::Double3 xyz = this->GetPosition();
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
          // vessel position intersection to the user-defined surface
          bdm::Double3 intx;
          if (! line_intersects_plane(l0, xyz, normal, origin, intx))
            continue;
          // skip following computations if projection point is outside triangle
          if (! is_inside_triangle(tri3.vertex_0, tri3.vertex_1, tri3.vertex_2, intx))
            continue;
          // vessel position projection to the user-defined surface
          const bdm::Double3 proj = project_to_plane(normal, origin, xyz);
          //
          const bdm::Double3 xyz_proj = xyz - proj;
          const double distance = L2norm(xyz_proj);
          // skip following if vessel is well outside the user-defined surface
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
            ABORT_("could not normalize the axis vector");
          // vessel position projection to the user-defined surface
          const bdm::Double3& proj = ci->second.second;
          //
          xyz = proj + normal * safe_distance;
          // enforce vessel to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          // vessel has remained inside the simulation domain
          return true;
          // ...end of if-case
        }
      // ...end of 'surface' simulation obstacles loop
    }
  // check if vessel has reached any of the 'scaffold' simulation obstacles
  for (size_t l=0; l<obstacles->scaffold.size(); l++)
    {
      bdm::Double3 xyz = this->GetPosition();
      //
      const unsigned int n_segm = obstacles->scaffold[l].segment.size();
      for (unsigned int s=0; s<n_segm; s++)
        {
          const ObstacleScaffold::Segment& segm = obstacles->scaffold[l].segment[s];
          //
          const bdm::Double3 n0 = segm.vertex_0,
                             n1 = segm.vertex_1;
          // check if vessel position is inside this obstacle
          if (! is_inside_segment(n0, n1, xyz))
            continue;
          // vessel position projection to the user-defined segment
          const bdm::Double3 proj = project_to_line(n0, n1, xyz);
          const double distance = L2norm(bdm::Double3(xyz-proj))
                                - segm.radius;
          //
          if (distance>rg->Uniform(0.5,1.0)*this->GetDiameter()) continue;
          //
          bdm::Double3 normal = xyz - proj;
          if (!normalize(normal, normal))
            ABORT_("could not normalize the axis vector");
          //
          const double delta = (rg->Uniform(0.5,1.0)*this->GetDiameter());
          xyz = proj + normal * delta;
          // enforce vessel to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          // vessel has remained inside the simulation domain
          return true;
          // ...end of segments loop
        }
      // ...end of 'scaffold' simulation obstacles loop
    }
  // vessel has remained inside the simulation domain
  return true;
}
// -----------------------------------------------------------------------------
inline
bool bdm::Vessel::CheckGrowth()
{
  if (!this->GetCanGrow()) return false;
  // vessel can grow if not a terminal one
  if (this->IsTerminal()) return false;
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
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>("vessel/can_grow/probability"))
    return false;
  //
  const double diameter = this->GetDiameter(),
               dia_min = this->params()->get<double>("vessel/diameter/min"),
               dia_max = this->params()->get<double>("vessel/diameter/max");
  // check if vessel diameter is above or equal to maximum size
  if (diameter>=dia_max) return false;
  //
  double diameter_rate = this->params()->get<double>("vessel/can_grow/diameter_rate");
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure vessel is well within the simulation domain!
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
                   threshold = this->params()->get<double>("vessel/can_grow/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          diameter_rate =
            std::max( diameter_rate,
                      this->params()->get<double>("vessel/can_grow/diameter_rate/"+BC_name) );
        }
      //...end of substances loop
    }
  //
  if ( ( diameter < dia_max && diameter_rate ) ||
       ( diameter < dia_min && diameter_rate > 0.0 ) )
    {
      this->ChangeDiameter(diameter_rate);
      // since vessel has grown then it can do something else
      return true;
    }
  // vessel has not grown / remodelled
  return false;
  //...end of vessel growth
}
// -----------------------------------------------------------------------------
inline
bool bdm::Vessel::CheckBranching()
{
  if (!this->GetCanBranch()) return false;
  // vessel can branch if not a terminal one
  if (this->IsTerminal()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  const double dia = this->GetDiameter();
  const std::string vessel_ID = std::to_string(this->GetVesselID());
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  if (this->GetAge()<this->params()->get<int>("vessel/can_branch/time_window")) return false;
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure vessel is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
    // iterate for all substances
    for ( std::vector<std::string>::const_iterator
          ci=substances.begin(); ci!=substances.end(); ci++ )
      {
        // access the BioDynaMo diffusion grid
        auto* dg = rm->GetDiffusionGrid(*ci);
        const std::string BC_name = dg->GetContinuumName(); // biochemical name
        //
        if (! this->params()->have_parameter<double>("vessel/can_branch/"+BC_name+"/elongation_length"))
          continue;
        const double elongation_length = this->params()->get<double>("vessel/can_branch/"+BC_name+"/elongation_length");
        //
        if (0.0==elongation_length) continue;
        //
        const double concentration = dg->GetValue(this->GetPosition()),
                     threshold = this->params()->get<double>("vessel/can_branch/"+BC_name+"/threshold");
        //
        if ( ( threshold > 0.0 && concentration > +threshold ) ||
             ( threshold < 0.0 && concentration < -threshold ) )
          {
            // check if nearby branches exist in this vessel
            bool safe_distance = true;
            {
              const double dist_min = this->params()->get<double>("vessel/can_branch/distance/min"),
                           dist_max = this->params()->get<double>("vessel/can_branch/distance/max");
              const double safe_dist = rg->Uniform(dist_min,dist_max);
              if (safe_dist>this->LengthToProximalBranchingPoint())
                safe_distance = false;
            }
            if (safe_distance)
            if (this->params()->have_parameter<std::vector<bdm::Double3>>("vessel_"+vessel_ID+"/branches_list"))
              {
                const std::vector<bdm::Double3>& branches =
                  this->params()->get<std::vector<bdm::Double3>>("vessel_"+vessel_ID+"/branches_list");
                //
                const double dist_min = this->params()->get<double>("vessel/can_branch/distance/min"),
                             dist_max = this->params()->get<double>("vessel/can_branch/distance/max");
                const double safe_dist = rg->Uniform(dist_min,dist_max);
                //
                for (unsigned int b=0; b<branches.size(); b++)
                  {
                    const bdm::Double3 dist = this->GetPosition() - branches[b];
                    if (L2norm(dist) < safe_dist)
                      {
                        safe_distance = false;
                        break;
                      }
                  }
              }
            // ...if not then definitely there is potential to branch
            if (safe_distance)
            if (rg->Uniform(0.0,1.0) <= this->params()->get<double>("vessel/can_branch/"+BC_name+"/probability"))
              {
                double extend_rate = (fabs(elongation_length) + 0.5 * dia)
                                   / this->params()->get<double>("time_step");
                // obtain normalized gradient of substance
                bdm::Double3 grad = {0.0, 0.0, 0.0};
                dg->GetGradient(this->GetPosition(), &grad);
                // introduce the signum (from extension) in the gradient vector
                grad *= sign(elongation_length);
                if (this->params()->get<bool>("simulation_domain_is_2D")) grad[2] *= 1.0e-9;
                // vessel direction vector
                bdm::Double3 vdvec = {0.0, 0.0, 0.0};
                {
                  bdm::Double3 axis = this->GetUnitaryAxisDirectionVector();
                  bdm::Double3 eta = bdm::Math::CrossProduct(grad, axis);
                  vdvec = bdm::Math::CrossProduct(axis, eta);
                }
                // finally, we simply normalize it
                if (!normalize(vdvec, vdvec))
                  // ...unless we have an error, therefore, abort any vessel branching
                  return false;
                //
                // diameter and ID for new branched vessel
                const double new_dia = this->params()->get<double>("vessel/diameter/min");
                const int new_vessel_ID = this->params()->get<int>("N_vessel_IDs");
                //
                auto* branch = bdm::bdm_static_cast<Vessel*>(this->Branch(grad));
                branch->SetDiameter(new_dia);
                branch->SetVesselID(new_vessel_ID);
                branch->SetCanGrow(true);
                branch->SetCanBranch(false);
                branch->SetCanSprout(true);
                branch->SetParametersPointer(this->params());
                //
                // FIX issue with vessel age
                branch->ElongateTerminalEnd(extend_rate, vdvec);
                branch->RunDiscretization();
                //
                // enforce this vessel to stop any branching prospectively
                this->SetCanBranch(false);
                this->SetCanSprout(false);
                //
                // update the list of the branches (coordinates)
                if (!this->params()->have_parameter<std::vector<bdm::Double3>>("vessel_"+vessel_ID+"/branches_list"))
                  {
                    // create the vector with a single component
                    std::vector<bdm::Double3> branches;
                    branches.push_back(this->GetPosition());
                    //branches.push_back(branch->GetPosition());
                    // save the vector
                    this->params()->set<std::vector<bdm::Double3>>("vessel_"+vessel_ID+"/branches_list") =
                      branches;
                  }
                else
                  {
                    // access the vector
                    std::vector<bdm::Double3>& branches =
                      this->params()->set<std::vector<bdm::Double3>>("vessel_"+vessel_ID+"/branches_list");
                    // ...and update it
                    branches.push_back(this->GetPosition());
                    //branches.push_back(branch->GetPosition());
                  }
                // increment this very important counter
                this->params()->set<int>("N_vessel_IDs") += 1;
                //
                return true;
              }
          }
        //...end of substances loop
      }
  // vessel has not been through any branch creation
  return false;
  //...end of vessel branching
}
// -----------------------------------------------------------------------------
inline
bool bdm::Vessel::CheckSprouting()
{
  if (!this->GetCanSprout()) return false;
  // vessel can sprout if only a terminal one
  if (!this->IsTerminal()) return false;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  const double dia = this->GetDiameter();
  const double len = this->GetActualLength();
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  // ensure vessel is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // iterate for all substances if vessel can sprout
  for ( std::vector<std::string>::const_iterator
        ci=substances.begin(); ci!=substances.end(); ci++ )
    {
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      const std::string BC_name = dg->GetContinuumName(); // biochemical name
      //
      if (! this->params()->have_parameter<double>("vessel/can_sprout/"+BC_name+"/elongation_length"))
        continue;
      const double elongation_length = this->params()->get<double>("vessel/can_sprout/"+BC_name+"/elongation_length");
      //
      if (0.0==elongation_length) continue;
      //
      const double concentration = dg->GetValue(this->GetPosition()),
                   threshold = this->params()->get<double>("vessel/can_sprout/"+BC_name+"/threshold");
      //
      if ( ( threshold > 0.0 && concentration > +threshold ) ||
           ( threshold < 0.0 && concentration < -threshold ) )
        {
          if (rg->Uniform(0.0,1.0) <= this->params()->get<double>("vessel/can_sprout/"+BC_name+"/probability"))
            {
              double extend_rate = fabs(elongation_length)
                                 / this->params()->get<double>("time_step");
              // obtain normalized gradient of substance
              bdm::Double3 grad = {0.0, 0.0, 0.0};
              dg->GetGradient(this->GetMassLocation(), &grad);
              // introduce the signum (from extension) in the gradient vector
              grad *= sign(elongation_length);
              if (this->params()->get<bool>("simulation_domain_is_2D")) grad[2] *= 1.0e-9;
              // vessel direction vector
              bdm::Double3 vdvec = grad;
              // ...add some randomness
              if (this->params()->get<double>("vessel/can_sprout/elongation_length_range") > 0.0)
                {
                  const double d = this->params()->get<double>("vessel/can_sprout/elongation_length_range");
                  // "distort" the vessel direction vector
                  vdvec[0] *= rg->Uniform(1.-d,1.+d);
                  vdvec[1] *= rg->Uniform(1.-d,1.+d);
                  vdvec[2] *= this->params()->get<bool>("simulation_domain_is_2D") ? 1.0 : rg->Uniform(1.-d,1.+d);
                }
              // finally, we simply normalize it
              if (!normalize(vdvec, vdvec))
                // ...unless we have an error, therefore, abort any vessel sprouting
                return false;
              //
              // check if new sprout to create is valid (inside simulation domain)
              this->CheckSproutingValidity(extend_rate, vdvec);
              if (!extend_rate) return false;
              //
              // FIX issue with vessel age
              this->ElongateTerminalEnd(extend_rate, vdvec);
              this->RunDiscretization();
              //
              // check if sprout has extended significantly
              if (len>=this->params()->get<double>("max_vessel_length"))
                {
                  double new_dia = dia;
                  new_dia *= this->params()->get<double>("vessel/can_sprout/diameter_scaleup");
                  // enforce this vessel to expand its lumen (within the diamater range)
                  if (new_dia<=this->params()->get<double>("vessel/diameter/max"))
                    this->SetDiameter(new_dia);
                }
              //
              const double dist_min = this->params()->get<double>("vessel/can_branch/distance/min"),
                           dist_max = this->params()->get<double>("vessel/can_branch/distance/max");
              const double safe_dist = rg->Uniform(dist_min,dist_max);
              //
              // also let the vessel to branch if necessary
              if (this->LengthToProximalBranchingPoint() > safe_dist)
                this->SetCanBranch(true);
              //
              return true;
            }
        }
      //...end of substances loop
    }
  // vessel has not been through any sprouting
  return false;
  //...end of vessel sprouting
}
// -----------------------------------------------------------------------------
inline
void bdm::Vessel::ScanAge() const
{
  if (this->IsTerminal()) return;
  //
  if ( auto* next = bdm::bdm_static_cast<const Vessel*>(this->GetDaughterLeft().Get()) )
    {
      next->SetAge(this->GetAge());
      next->ScanAge();
    }
  //
  if ( auto* branch = bdm::bdm_static_cast<const Vessel*>(this->GetDaughterRight().Get()) )
    {
      branch->SetAge(this->GetAge());
      branch->ScanAge();
    }
}
// -----------------------------------------------------------------------------
inline
void bdm::Vessel::CheckSproutingValidity(double& dLdT, bdm::Double3& axis)
{
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary"),
               meanCOORD = (maxCOORD+minCOORD)/2.0,
               deltaCOORD = maxCOORD - meanCOORD;
  const double tol = this->params()->get<double>("domain_tolerance");
  // simulation increment in time (time-step)
  const double dT = this->params()->get<double>("time_step");
  // length by which tip vessel sprouts (extends)
  double extend = dT * dLdT;
  // nodal points of new vessel [start, end]
  const bdm::Double3 n[] = { this->GetMassLocation()             ,
                             this->GetMassLocation()+axis*extend };
  const bdm::Double3 n1_n0 = n[1] - n[0];
  const bdm::Double3 m = {meanCOORD, meanCOORD, meanCOORD};
  //
  // check if simulation domain is bounded or unbounded (for vessel extension)
  if ( this->params()->get<bool>("simulation_domain_is_bounded") )
    {
      ABORT_("code is under development");
    }
  else
    {
      const bdm::Double3& xyz = n[1];
      // identify mode of simulation domain
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        {
          if (this->params()->get<bool>("simulation_domain_is_polar"))
            {
              if ( sqrt(pow2(xyz[0]-meanCOORD)
                       +pow2(xyz[1]-meanCOORD)) > deltaCOORD-tol )
                {
                  const bdm::Double3 m_n0 = m - n[0];
                  const double r0 = L2norm(m_n0);
                  const bdm::Double3 k = normalize(n[0]-n[1]);
                  const double r = deltaCOORD-tol;
                  extend = (pow2(r)-pow2(r0)) / (2.0*(m_n0*k));
                }
            }
          else
            {
              if      ( xyz[0] < minCOORD+tol )
                {
                  const bdm::Double3 e = {-1.0, 0.0, 0.0};
                  const bdm::Double3 v = {minCOORD+tol, 0.0, 0.0};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
              else if ( xyz[0] > maxCOORD-tol )
                {
                  const bdm::Double3 e = {+1.0, 0.0, 0.0};
                  const bdm::Double3 v = {maxCOORD-tol, 0.0, 0.0};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
              else if ( xyz[1] < minCOORD+tol )
                {
                  const bdm::Double3 e = {0.0, -1.0, 0.0};
                  const bdm::Double3 v = {0.0, minCOORD+tol, 0.0};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
              else if ( xyz[1] > maxCOORD-tol )
                {
                  const bdm::Double3 e = {0.0, +1.0, 0.0};
                  const bdm::Double3 v = {0.0, maxCOORD-tol, 0.0};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
            }
        }
      else
        {
          if (this->params()->get<bool>("simulation_domain_is_polar"))
            {
              if ( sqrt(pow2(xyz[0]-meanCOORD)
                       +pow2(xyz[1]-meanCOORD)
                       +pow2(xyz[2]-meanCOORD)) > deltaCOORD-tol )
                {
                  const bdm::Double3 m_n0 = m - n[0];
                  const double r0 = L2norm(m_n0);
                  const bdm::Double3 k = normalize(n[0]-n[1]);
                  const double r = deltaCOORD-tol;
                  extend = (pow2(r)-pow2(r0)) / (2.0*(m_n0*k));
                  extend = (pow2(r)-pow2(r0)) / (2.0*(m_n0*k));
                }
            }
          else
            {
              // https://www.geomalgorithms.com/a05-_intersect-1.html
              if      ( xyz[0] < minCOORD+tol )
                {
                  const bdm::Double3 e = {-1.0, 0.0, 0.0};
                  const bdm::Double3 v = {minCOORD+tol, 0.0, 0.0};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
              else if ( xyz[0] > maxCOORD-tol )
                {
                  const bdm::Double3 e = {+1.0, 0.0, 0.0};
                  const bdm::Double3 v = {maxCOORD-tol, 0.0, 0.0};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
              else if ( xyz[1] < minCOORD+tol )
                {
                  const bdm::Double3 e = {0.0, -1.0, 0.0};
                  const bdm::Double3 v = {0.0, minCOORD+tol, 0.0};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
              else if ( xyz[1] > maxCOORD-tol )
                {
                  const bdm::Double3 e = {0.0, +1.0, 0.0};
                  const bdm::Double3 v = {0.0, maxCOORD-tol, 0.0};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
              else if ( xyz[2] < minCOORD+tol )
                {
                  const bdm::Double3 e = {0.0, 0.0, -1.0};
                  const bdm::Double3 v = {0.0, 0.0, minCOORD+tol};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
              else if ( xyz[2] > maxCOORD-tol )
                {
                  const bdm::Double3 e = {0.0, 0.0, +1.0};
                  const bdm::Double3 v = {0.0, 0.0, maxCOORD-tol};
                  const bdm::Double3 v_n0 = v - n[0];
                  extend = (e*v_n0) / (e*n1_n0);
                }
            }
        }
    }
  // check if extension length above threshold to allow sprouting
  if (extend<this->params()->get<double>("min_vessel_length"))
    extend = 0.0;
  // update the vessel sprouting (extension) rate
  dLdT = extend / dT;
}
// =============================================================================
#endif // _VESSEL_INLINE_H_
// =============================================================================

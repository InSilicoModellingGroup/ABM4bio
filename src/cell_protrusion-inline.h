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
#ifndef _CELL_PROTRUSION_INLINE_H_
#define _CELL_PROTRUSION_INLINE_H_
// =============================================================================
inline
int bdm::CellProtrusion::GetTimeRepeats() const
{
  // by design only viable (non-necrotic) cells could have protrusion development
  if (!this->GetCell()->GetPhenotype()) return 0;
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetCell()->GetPhenotype()));
  //
  return this->params()->get<int>(CP_name+"/can_protrude/time_repeats");
}
// -----------------------------------------------------------------------------
inline
bool bdm::CellProtrusion::CheckPositionValidity()
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
  // nodal points of the original protrusion [start, end]
  bdm::Double3 node[2];
  node[0] = this->GetPosition() - this->GetSpringAxis() * 0.5;
  node[1] = this->GetPosition() + this->GetSpringAxis() * 0.5;
  //
  // check if simulation domain is bounded or unbounded (for protrusion extension)
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
          // protrusion went out of the (simulation) domain
          // therefore it must stop further developing anyhow
          return false;
        }
    }
  // access the pointer to parameter of the simulation obstacles object
  const SimulationObstacles* obstacles =
    this->params()->get<SimulationObstacles*>("simulation_obstacles_data");
  // check if protrusion has reached any of the 'box' simulation obstacles
  for (size_t l=0; l<obstacles->box.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetCell()->GetPhenotype()) continue;
      // original protrusion position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      // it true, then protrusion position needs to be corrected
      if ("box/inside"==obstacles->box[l].type)
        {
          // check if protrusion position is outside this obstacle
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
          // enforce protrusion to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          // protrusion has remained inside the simulation domain
          return true;
        }
      else if ("box/outside"==obstacles->box[l].type)
        {
          // check if protrusion position is inside this obstacle
          if (obstacles->box[l].is_inside(xyz))
            continue;
          // calculate the outward unit normal vector and the projection vector of
          // the protrusion to each box face
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
              // enforce protrusion to lie on the (box) obstacle surface
              this->SetPosition(xyz);
              // protrusion has remained inside the simulation domain
              return true;
            }
        }
      else
        ABORT_("an exception is caught");
      // ...end of 'box' simulation obstacles loop
    }
  // check if protrusion has reached any of the 'sphere' simulation obstacles
  for (size_t l=0; l<obstacles->sphere.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetCell()->GetPhenotype()) continue;
      // original protrusion position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      // it true, then protrusion position needs to be corrected
      if ("sphere/inside"==obstacles->sphere[l].type)
        {
          // check if protrusion position is outside this obstacle
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
          // enforce protrusion to lie on the surface of the (spherical) obstacle
          this->SetPosition(xyz);
          // protrusion has remained inside the simulation domain
          return true;
        }
      else if ("sphere/outside"==obstacles->sphere[l].type)
        {
          // check if protrusion position is inside this obstacle
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
          // enforce protrusion to lie on the surface of the (spherical) obstacle
          this->SetPosition(xyz);
          // protrusion has remained inside the simulation domain
          return true;
        }
      else
        ABORT_("an exception is caught");
      // ...end of 'sphere' simulation obstacles loop
    }
  // check if protrusion has reached any of the 'surface' simulation obstacles
  for (size_t l=0; l<obstacles->surface.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetCell()->GetPhenotype()) continue;
      // original protrusion position (space vector)
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
          // protrusion position intersection to the user-defined surface
          bdm::Double3 intx;
          if (! line_intersects_plane(l0, xyz, normal, origin, intx))
            continue;
          // skip following computations if projection point is outside triangle
          if (! is_inside_triangle(tri3.vertex_0, tri3.vertex_1, tri3.vertex_2, intx))
            continue;
          // protrusion position projection to the user-defined surface
          const bdm::Double3 proj = project_to_plane(normal, origin, xyz);
          //
          const bdm::Double3 xyz_proj = xyz - proj;
          const double distance = L2norm(xyz_proj);
          // skip following if protrusion is well outside the user-defined surface
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
          // protrusion position projection to the user-defined surface
          const bdm::Double3& proj = ci->second.second;
          //
          xyz = proj + normal * safe_distance;
          // enforce protrusion to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          // protrusion has remained inside the simulation domain
          return true;
          // ...end of if-case
        }
      // ...end of 'surface' simulation obstacles loop
    }
  // check if protrusion has reached any of the 'scaffold' simulation obstacles
  for (size_t l=0; l<obstacles->scaffold.size(); l++)
    {
      // skip subsequent computations for a necrotic cell...
      if (0==this->GetCell()->GetPhenotype()) continue;
      // original protrusion position (space vector)
      bdm::Double3 xyz = this->GetPosition();
      //
      const unsigned int n_segm = obstacles->scaffold[l].segment.size();
      for (unsigned int s=0; s<n_segm; s++)
        {
          const ObstacleScaffold::Segment& segm = obstacles->scaffold[l].segment[s];
          //
          const bdm::Double3 n0 = segm.vertex_0,
                             n1 = segm.vertex_1;
          // check if protrusion position is inside this obstacle
          if (! is_inside_segment(n0, n1, xyz))
            continue;
          // protrusion position projection to the user-defined segment
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
          // enforce protrusion to lie on the (box) obstacle surface
          this->SetPosition(xyz);
          // protrusion has remained inside the simulation domain
          return true;
          // ...end of segments loop
        }
      // ...end of 'scaffold' simulation obstacles loop
    }
  // protrusion has remained inside the simulation domain
  return true;
}
// -----------------------------------------------------------------------------
inline
bool bdm::CellProtrusion::CheckTimeWindow() const
{
  // by design only viable (non-necrotic) cells could have protrusion development
  if (!this->GetCell()->GetPhenotype())
    return false;
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetCell()->GetPhenotype()));
  //
  if ( this->GetCell()->GetAge() < this->params()->get<int>(CP_name+"/can_protrude/time_window_open" ) ||
       this->GetCell()->GetAge() > this->params()->get<int>(CP_name+"/can_protrude/time_window_close") )
    return false;
  // ...otherwise cell protrusion could further develop
  return true;
}
// -----------------------------------------------------------------------------
inline
void bdm::CellProtrusion::CheckOutGrowth()
{
  if (!this->GetCanSprout())
    return;
  // by design only viable (non-necrotic) cells could have protrusions' growth
  if (!this->GetCell()->GetPhenotype())
    return;
  // normal exit if protrusion is not a terminal segment; only terminal ones can sprout (elongate)
  if (!this->IsTerminal() )
    return;
  //
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetCell()->GetPhenotype()));
  //
  // cell protrusion minimum diameter
  const double dia_min = this->params()->get<double>(CP_name+"/can_protrude/diameter/min");
  // if cell protrusion diameter is below that threshold then simply do not sprout or branch
  if ( this->GetDiameter() < dia_min )
    return;
  //
  if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_protrude/sprout/probability"))
    return;
  //
  const std::vector<std::string>& substances =
    this->params()->get<std::vector<std::string>>("substances");
  //
  const double current_weight = this->params()->get<double>(CP_name+"/can_protrude/sprout/current_direction/weight");
  bdm::Double3 current_direction = {0.0, 0.0, 0.0};
  if (current_weight>0.0)
    {
      bdm::Double3 axis;
      axis = this->GetSpringAxis();
      if (normalize(axis, axis))
        current_direction = axis * current_weight;
    }
  //
  const double random_weight = this->params()->get<double>(CP_name+"/can_protrude/sprout/random_direction/weight");
  bdm::Double3 random_direction = {0.0, 0.0, 0.0};
  if (random_weight>0.0)
    {
      bdm::Double3 axis;
      axis = {rg->Uniform(-1.0,+1.0), rg->Uniform(-1.0,+1.0), rg->Uniform(-1.0,+1.0)};
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        axis[2] = 0.0;
      if (normalize(axis, axis))
        random_direction = axis * random_weight;
    }
  // min and max boundaries of the BioDynaMo simulation 3D/2D domain
  const double minCOORD = this->params()->get<double>("min_boundary"),
               maxCOORD = this->params()->get<double>("max_boundary");
  const double tol = this->params()->get<double>("domain_tolerance");
  //
  const double chemotaxis_weight = this->params()->get<double>(CP_name+"/can_protrude/sprout/chemotaxis_direction/weight");
  bdm::Double3 chemotaxis_direction = {0.0, 0.0, 0.0};
  // ensure cell protrusion is well within the simulation domain!
  if (check_agent_position_in_domain(minCOORD, maxCOORD, this->GetPosition(), tol))
  // if true, then...
  if (chemotaxis_weight)
    {
      const std::string& chemo_substance = this->params()->get<std::string>(CP_name+"/can_protrude/sprout/chemotaxis");
      // iterate for all substances
      std::vector<std::string>::const_iterator
        ci = std::find(substances.begin(), substances.end(), chemo_substance);
      ASSERT_(substances.end()!=ci,"an internal error occurred");
      // access the BioDynaMo diffusion grid
      auto* dg = rm->GetDiffusionGrid(*ci);
      //
      bdm::Double3 axis;
      dg->GetGradient(this->GetPosition(), &axis);
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        axis[2] = 0.0;
      if (normalize(axis, axis))
        chemotaxis_direction = axis * chemotaxis_weight;
    }
  //
  bdm::Double3 new_direction = current_direction + random_direction + chemotaxis_direction;
  ASSERT_(normalize(new_direction, new_direction),
          "could not normalize the direction vector");
  // elongate the protrusion
  const double L_rate = this->params()->get<double>(CP_name+"/can_protrude/sprout/elongation_rate");
  this->ElongateTerminalEnd(L_rate, new_direction);
  if ( this->GetActualLength() > this->params()->get<double>("max_vessel_length") )
    this->RunDiscretization();
  // shrink its diameter
  const double D_rate = this->params()->get<double>(CP_name+"/can_protrude/sprout/remodelling_ratio");
  this->SetDiameter(this->GetDiameter()*D_rate);
  //
  if ( this->IsTerminal() )
    this->SetCanBranch(true);
  //...end of cell protrusion outgrowth
}
// -----------------------------------------------------------------------------
inline
void bdm::CellProtrusion::CheckBranching()
{
  if (!this->GetCanBranch())
    return;
  // by design only viable (non-necrotic) cells could have protrusions' branching
  if (!this->GetCell()->GetPhenotype())
    return;
  //
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();
  //
  const std::string& CP_name = // cell phenotype name
    this->params()->get<std::string>("phenotype_ID/"+std::to_string(this->GetCell()->GetPhenotype()));
  //
  // only a terminal protrusion can dissect in two new protrusions
  if ( this->IsTerminal() )
    {
      //
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_protrude/dissect/probability"))
        return;
      //
      const double min_length2branch = this->LengthToProximalBranchingPoint();
      if (min_length2branch < this->params()->get<double>(CP_name+"/can_protrude/dissect/length2branch/min"))
        return;
      //
      bdm::Double3 random_axis;
      random_axis = {rg->Uniform(0.0,1.0), rg->Uniform(0.0,1.0), rg->Uniform(0.0,1.0)};
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        random_axis[2] = 0.0;
      //
      const bdm::Double3 branch_direction =
        bdm::Math::Perp3(this->GetUnitaryAxisDirectionVector()+random_axis, rg->Uniform(0.0,1.0))
      + this->GetSpringAxis();
      //
      auto* protrusion_new =
        bdm::bdm_static_cast<CellProtrusion*>(this->Branch(branch_direction));
      //
      protrusion_new->SetCell(this->GetCell());
      protrusion_new->SetCanBranch(true);
      protrusion_new->SetDiameter(this->GetDiameter());
      protrusion_new->SetParametersPointer(this->params());
      protrusion_new->AddBehavior(new Biology4CellProtrusion());
      //
      this->SetCanBranch(false);
      //
      auto* previous = bdm::bdm_static_cast<CellProtrusion*>(this->GetMother().Get());
      ASSERT_(nullptr!=previous,"an internal error occurred");
      //
      previous->SetCanSprout(false);
      previous->SetCanBranch(false);
      //
      {
        auto* back2previous = bdm::bdm_static_cast<CellProtrusion*>(previous->GetMother().Get());
        ASSERT_(nullptr!=back2previous,"an internal error occurred");
        //
        back2previous->SetCanSprout(false);
        back2previous->SetCanBranch(false);
      }
    }
  // ...otherwise a non-terminal protrusion can produce a branch,
  // i.e. generate a new branched protrusion from an existing one
  else
    {
      //
      if (rg->Uniform(0.0,1.0) > this->params()->get<double>(CP_name+"/can_protrude/branch/probability"))
        return;
      //
      const double min_length2branch = this->LengthToProximalBranchingPoint();
      if (min_length2branch < this->params()->get<double>(CP_name+"/can_protrude/branch/length2branch/min"))
        return;
      //
      bdm::Double3 random_axis;
      random_axis = {rg->Uniform(0.0,1.0), rg->Uniform(0.0,1.0), rg->Uniform(0.0,1.0)};
      if (this->params()->get<bool>("simulation_domain_is_2D"))
        random_axis[2] = 0.0;
      //
      const bdm::Double3 branch_direction =
        bdm::Math::Perp3(this->GetUnitaryAxisDirectionVector()+random_axis, rg->Uniform(0.0,1.0))
      + this->GetSpringAxis();
      //
      auto* protrusion_new =
        bdm::bdm_static_cast<CellProtrusion*>(this->Branch(branch_direction));
      //
      protrusion_new->SetCell(this->GetCell());
      protrusion_new->SetCanBranch(true);
      protrusion_new->SetDiameter(this->GetDiameter());
      protrusion_new->SetParametersPointer(this->params());
      protrusion_new->AddBehavior(new Biology4CellProtrusion());
      //
      this->SetCanBranch(false);
      //
      auto* previous = bdm::bdm_static_cast<CellProtrusion*>(this->GetMother().Get());
      ASSERT_(nullptr!=previous,"an internal error occurred");
      //
      previous->SetCanSprout(false);
      previous->SetCanBranch(false);
      //
      {
        auto* back2previous = bdm::bdm_static_cast<CellProtrusion*>(previous->GetMother().Get());
        ASSERT_(nullptr!=back2previous,"an internal error occurred");
        //
        back2previous->SetCanSprout(false);
        back2previous->SetCanBranch(false);
      }
    }
  //...end of cell protrusion branching
}
// -----------------------------------------------------------------------------
inline
void bdm::CellProtrusion::CheckState()
{
  // by design only viable (non-necrotic) cells could have protrusions' screened
  if (!this->GetCell()->GetPhenotype())
    return;
  // ...simply exit this function and skip subsequent steps; this function
  // is still WIP!!! :)
  ;
}
// -----------------------------------------------------------------------------
inline
void bdm::CellProtrusion::Set2Delete() const
{
  this->SetCell(nullptr);
  // if terminal protrusion then exit function
  if (this->IsTerminal())
    return;
  // check for sprout
  if (this->GetDaughterLeft())
    {
      auto* protrusion =
        bdm::bdm_static_cast<const CellProtrusion*>(this->GetDaughterLeft().Get());
      // nullify the pointer to the cell, then scan for subsequent protrusions
      protrusion->SetCell(nullptr);
      protrusion->Set2Delete();
    }
  else
    ABORT_("an exception is caught");
  // check for branch (if one exists)
  if (this->GetDaughterRight())
    {
      auto* protrusion =
        bdm::bdm_static_cast<const CellProtrusion*>(this->GetDaughterRight().Get());
      // nullify the pointer to the cell, then scan for subsequent protrusions
      protrusion->SetCell(nullptr);
      protrusion->Set2Delete();
    }
}
// =============================================================================
#endif // _CELL_PROTRUSION_INLINE_H_
// =============================================================================

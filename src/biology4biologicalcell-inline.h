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
#ifndef _BIOLOGY4BIOLOGICALCELL_INLINE_H_
#define _BIOLOGY4BIOLOGICALCELL_INLINE_H_
// =============================================================================
#include "./biological_cell.h"
// =============================================================================
inline
void bdm::Biology4BiologicalCell_10::Run(bdm::Agent* a)
{
  if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
    {
      // firstly, we should check if cell is inside the simulation domain
      if (!cell->CheckPositionValidity())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // we check for cell apoptosis
      if (cell->CheckApoptosis())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // simply update the cell age
      cell->IncrementAge();
      // cell executes is regulatory network
      cell->RunRegulatoryNetwork();
      // cell produces/consumes substances
      cell->RunBiochemics();
      // now check if cell can migrate
      if (cell->CheckMigration())
        {
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
        }
      // then, check if cell can transform or if it can polarize
      if (cell->CheckTransformation()) return;
      cell->CheckPolarization();
      cell->CheckProtrusion();
      // check if cell can grow
      if (cell->CheckGrowth()) return;
      // check if cell can divide (summetrically or unsymmetrically)
      if (cell->CheckTransformationAndDivision()) return;
      if (cell->CheckAsymmetricDivision()) return;
      if (cell->CheckDivision()) return;
      // finally, we check for cell apoptosis due to aging
      if (cell->CheckApoptosisAging())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // ...end of mechanisms list for cell behavior
    }
  else
    ABORT_("an exception is caught");
}
// -----------------------------------------------------------------------------
inline
void bdm::Biology4BiologicalCell_11::Run(bdm::Agent* a)
{
  if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
    {
      if (bdm::BiologicalCell::Phase::Ap==cell->GetPhase()) {
        cell->IncrementAge();
        //
        if (cell->CheckAfterApoptosis())
          {
            cell->Set2DeleteProtrusions();
            cell->RemoveBehavior(this);
            cell->RemoveFromSimulation();
            return;
          }
        return;
      }
      // firstly, we should check if cell is inside the simulation domain
      if (!cell->CheckPositionValidity())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // simply update the cell age
      cell->IncrementAge();
      // cell executes is regulatory network
      cell->RunRegulatoryNetwork();
      // cell produces/consumes substances
      cell->RunBiochemics();
      //
      if (bdm::BiologicalCell::Phase::Di==cell->GetPhase())
        {
          if (cell->CheckQuiescenceAfterDivision())
            return;
        }
      // we check for cell apoptosis
      if (cell->CheckApoptosis())
        {
          cell->SetAge(); // reset age (time) counter
          cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
          return;
        }
      // check if cell can transform
      if (cell->CheckTransformation())
        {
          cell->SetPhase(bdm::BiologicalCell::Phase::G1);
          return;
        }
      // check if cell can polarize
      if (cell->CheckPolarization())
        {
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
        }
      // now check if cell can migrate
      if (cell->CheckMigration())
        {
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
        }
      // check if cell can divide (summetrically or unsymmetrically)
      if (bdm::BiologicalCell::Phase::G1==cell->GetPhase()||
          bdm::BiologicalCell::Phase::Sy==cell->GetPhase())
        if (cell->CheckDivision() || cell->CheckAsymmetricDivision())
          {
            cell->SetAge(); // reset age (time) counter
            cell->SetPhase(bdm::BiologicalCell::Phase::Di);
            return;
          }
      // check if cell can grow
      if (bdm::BiologicalCell::Phase::G1==cell->GetPhase()||
          bdm::BiologicalCell::Phase::Sy==cell->GetPhase())
        if (cell->CheckGrowth())
          {
            cell->SetPhase(bdm::BiologicalCell::Phase::G1);
            return;
          }
      // finally, we check for cell apoptosis due to aging
      if (bdm::BiologicalCell::Phase::G1==cell->GetPhase()||
          bdm::BiologicalCell::Phase::Sy==cell->GetPhase()||
          bdm::BiologicalCell::Phase::G2==cell->GetPhase())
        if (cell->CheckApoptosisAging())
          {
            cell->SetAge(); // reset age (time) counter
            cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
            return;
          }
      //
      cell->SetPhase(bdm::BiologicalCell::Phase::Sy);
    }
  else
    ABORT_("an exception is caught");
}
// -----------------------------------------------------------------------------
inline
void bdm::Biology4BiologicalCell_12::Run(bdm::Agent* a)
{
  // access BioDynaMo's resource manager
  auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
  // access BioDynaMo's random number generator
  auto* rg = bdm::Simulation::GetActive()->GetRandom();

  if (auto* cell = dynamic_cast<bdm::BiologicalCell*>(a))
    {
      // by design only viable (non-necrotic) cells could develop
      if (!cell->GetPhenotype()) return;
      //
      const Parameters* p = cell->params();
      //
      const std::string& CP_name = // cell phenotype name
        p->get<std::string>("phenotype_ID/"+std::to_string(cell->GetPhenotype()));
      //
      // make sure that the cell is never in this phase
      if (cell->GetPhase() == bdm::BiologicalCell::Phase::I0)
        {
          // set the phase
          cell->SetPhase(bdm::BiologicalCell::Phase::G1);
          // reset the cell-cycle
          cell->ResetCycleTime();
        }
      // make sure that the cell is never in this phase
      if (cell->GetPhase() == bdm::BiologicalCell::Phase::Tr)
        ABORT_("an internal error occurred");
      // simply update the cell age
      cell->IncrementAge();
      // update the cell-cycle phase
      cell->IncrementCycleTime();
      // we should check if cell is inside the simulation domain
      if (!cell->CheckPositionValidity())
        {
          cell->Set2DeleteProtrusions();
          cell->RemoveBehavior(this);
          cell->RemoveFromSimulation();
          return;
        }
      // check if cell is already in apoptotic phase
      if (cell->GetPhase() == bdm::BiologicalCell::Phase::Ap)
        {
          if (cell->GetAge() >= p->get<int>(CP_name+"/after_apoptosis/duration"))
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
          // now skip running all subsequent cell mechanisms
          return;
        }
      // cell produces/consumes substances
      cell->RunBiochemics();
      // check if cell can migrate
      bool has_migrated = false;
      if (cell->GetCanMigrate())
        {
          ABORT_("aborting here due to unfinished business...");
        }
      //
      if (has_migrated)
        {
          // we should check if cell is inside the simulation domain
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
        }
      // check if the cell has been underway the division phase
      if (cell->GetPhase() == bdm::BiologicalCell::Phase::G2)
        {
          // update the phase
          cell->SetPhase(bdm::BiologicalCell::Phase::Di);
          // now skip running all subsequent cell mechanisms
          return;
        }
      else if (cell->GetPhase() == bdm::BiologicalCell::Phase::Di)
        {
          if (cell->GetCycleTime()>=p->get<int>(CP_name+"/G1_to_Di/duration"))
            {
              // update the phase
              cell->SetPhase(bdm::BiologicalCell::Phase::G1);
              // reset the cell-cycle
              cell->ResetCycleTime();
            }
        }
      // check if cell can apoptose
      bool has_apoptosed = false;
      if (cell->GetCanApoptose())
        {
          bool can_continue = true;
          //
          const std::vector<std::string>& substances =
            p->get<std::vector<std::string>>("substances");
          // make sure you can continue...
          if ( can_continue )
          // ++++++++++++++++++++++++++++++++++++++++++++++
          // iterate for all substances
          for ( std::vector<std::string>::const_iterator
                ci=substances.begin(); ci!=substances.end(); ci++ )
            {
              // access the BioDynaMo diffusion grid
              auto* dg = rm->GetDiffusionGrid(*ci);
              const std::string BC_name = dg->GetContinuumName(); // biochemical name
              //
              const double concentration = dg->GetValue(cell->GetPosition()),
                           threshold = p->get<double>(CP_name+"/can_apoptose/"+BC_name+"/threshold");
              //
              if ( ( threshold > 0.0 && concentration > +threshold ) ||
                   ( threshold < 0.0 && concentration < -threshold ) )
                {
                  if (rg->Uniform(0.0,1.0) <= p->get<double>(CP_name+"/can_apoptose/"+BC_name+"/probability"))
                    {
                      // update the phase
                      cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
                      // reset the age --> to control cell deletion
                      cell->SetAge();
                      // reset the cell-cycle
                      cell->ResetCycleTime();
                      //
                      has_apoptosed = true;
                      //
                      can_continue = false;
                      break;
                    }
                  // ...end of threshold check
                }
              //...end of substances loop
            }
          // ++++++++++++++++++++++++++++++++++++++++++++++
          //
          const int n_s =
            p->get<int>(CP_name+"/regulatory_network/number_of_species");
          // make sure you can continue...
          if ( can_continue )
          // ++++++++++++++++++++++++++++++++++++++++++++++
          // iterate for all species of the regulatory network
          for (int s=0; s<n_s; s++)
            {
              const std::string RN_specie = std::to_string(s);
              //
              const double level = cell->GetRegulatoryNetworkData().current_species[s],
                           threshold = p->get<double>(CP_name+"/can_apoptose/regulatory_network/"+RN_specie+"/threshold");
              //
              if ( ( threshold > 0.0 && level > +threshold ) ||
                   ( threshold < 0.0 && level < -threshold ) )
                {
                  if (rg->Uniform(0.0,1.0) <= p->get<double>(CP_name+"/can_apoptose/regulatory_network/"+RN_specie+"/probability"))
                    {
                      // update the phase
                      cell->SetPhase(bdm::BiologicalCell::Phase::Ap);
                      // reset the age --> to control cell deletion
                      cell->SetAge();
                      // reset the cell-cycle
                      cell->ResetCycleTime();
                      //
                      has_apoptosed = true;
                      //
                      can_continue = false;
                      break;
                    }
                  // ...end of threshold check
                }
              //...end of species loop
            }
          // ++++++++++++++++++++++++++++++++++++++++++++++
          //
        }        
      //
      if (has_apoptosed)
        {
          // we should check if cell is inside the simulation domain
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
          // now skip running all subsequent cell mechanisms
          return;
        }
      // check if cell can transform
      bool has_transformed = false;
      if (cell->GetCanTransform())
        {
          bool can_continue = true;
          //
          const std::vector<std::string>& substances =
            p->get<std::vector<std::string>>("substances");
          // make sure you can continue...
          if ( can_continue )
          // ++++++++++++++++++++++++++++++++++++++++++++++
          // iterate for all substances
          for ( std::vector<std::string>::const_iterator
                ci=substances.begin(); ci!=substances.end(); ci++ )
            {
              // access the BioDynaMo diffusion grid
              auto* dg = rm->GetDiffusionGrid(*ci);
              const std::string BC_name = dg->GetContinuumName(); // biochemical name
              //
              if (! p->have_parameter<int>(CP_name+"/can_transform/"+BC_name+"/new_phenotype"))
                continue;
              //
              const double concentration = dg->GetValue(cell->GetPosition()),
                           threshold = p->get<double>(CP_name+"/can_transform/"+BC_name+"/threshold");
              //
              if ( ( threshold > 0.0 && concentration > +threshold ) ||
                   ( threshold < 0.0 && concentration < -threshold ) )
                {
                  if (rg->Uniform(0.0,1.0) <= p->get<double>(CP_name+"/can_transform/"+BC_name+"/probability"))
                    {
                      // now check if cell age is within an appropriate time-window to allow its transformation
                      if ( cell->GetAge() >= p->get<int>(CP_name+"/can_transform/"+BC_name+"/time_window_open" ) &&
                           cell->GetAge() <= p->get<int>(CP_name+"/can_transform/"+BC_name+"/time_window_close") )
                        {
                          // reset the cell phenotype ID
                          cell->SetPhenotype(p->get<int>(CP_name+"/can_transform/"+BC_name+"/new_phenotype"));
                          // increment this index
                          cell->IncrementNumberOfTrasformations();
                          // reset the age of the cell
                          cell->SetAge();
                          // reset the cell-cycle
                          cell->ResetCycleTime();
                          // reset the phase
                          cell->SetPhase(bdm::BiologicalCell::Phase::G1);
                          //
                          const std::string CP_new_name = // cell phenotype name
                            p->get<std::string>("phenotype_ID/"+std::to_string(cell->GetPhenotype()));
                          // principal directions of the cell polarization matrix
                          double p0, p1, p2;
                          if (cell->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                            {
                              p0 = p->get<double>(CP_new_name+"/principal/0");
                              p1 = p->get<double>(CP_new_name+"/principal/1");
                              p2 = p->get<double>(CP_new_name+"/principal/2");
                            }
                          //
                          cell->SetCanApoptose(p->get<bool>(CP_new_name+"/can_apoptose"));
                          cell->SetCanGrow(p->get<bool>(CP_new_name+"/can_grow"));
                          cell->SetCanDivide(p->get<bool>(CP_new_name+"/can_divide"));
                          cell->SetCanMigrate(p->get<bool>(CP_new_name+"/can_migrate"));
                          cell->SetCanTransform(p->get<bool>(CP_new_name+"/can_transform"));
                          cell->SetMonitorCycle(p->get<bool>(CP_new_name+"/monitor_cycle"));              
                          // reset the cell polarization matrix
                          if (cell->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                            cell->SetPolarization(diag(p0, p1, p2));
                          //
                          ASSERT_(0==cell->GetNumberOfProtrusions(),
                                  "an exception is caught");
                          //
                          ASSERT_(12==p->get<int>(CP_new_name+"/mechanism_order"),
                                  "an exception is caught");
                          // setup the regulatory network data
                          if (p->have_parameter<std::string>(CP_new_name+"/regulatory_network/from_file"))
                            read_regulatory_network_data(
                              p->get<std::string>(CP_new_name+"/regulatory_network/from_file"),
                              cell->SetRegulatoryNetworkData() );
                          //
                          has_transformed = true;
                          //
                          can_continue = false;
                          break;
                        }
                      // ...end of time-window check
                    }
                  // ...end of threshold check
                }
              //...end of substances loop
            }
          // ++++++++++++++++++++++++++++++++++++++++++++++
          //
          const int n_s =
            p->get<int>(CP_name+"/regulatory_network/number_of_species");
          // make sure you can continue...
          if ( can_continue )
          // ++++++++++++++++++++++++++++++++++++++++++++++
          // iterate for all species of the regulatory network
          for (int s=0; s<n_s; s++)
            {
              const std::string RN_specie = std::to_string(s);
              //
              const double level = cell->GetRegulatoryNetworkData().current_species[s],
                           threshold = p->get<double>(CP_name+"/can_transform/regulatory_network/"+RN_specie+"/threshold");
              //
              if ( ( threshold > 0.0 && level > +threshold ) ||
                   ( threshold < 0.0 && level < -threshold ) )
                {
                  if (rg->Uniform(0.0,1.0) <= p->get<double>(CP_name+"/can_transform/regulatory_network/"+RN_specie+"/probability"))
                    {
                      // now check if cell age is within an appropriate time-window to allow its transformation
                      if ( cell->GetAge() >= p->get<int>(CP_name+"/can_transform/regulatory_network/"+RN_specie+"/time_window_open" ) &&
                           cell->GetAge() <= p->get<int>(CP_name+"/can_transform/regulatory_network/"+RN_specie+"/time_window_close") )
                        {
                          // reset the cell phenotype ID
                          cell->SetPhenotype(p->get<int>(CP_name+"/can_transform/regulatory_network/"+RN_specie+"/new_phenotype"));
                          // increment this index
                          cell->IncrementNumberOfTrasformations();
                          // reset the age of the cell
                          cell->SetAge();
                          // reset the cell-cycle
                          cell->ResetCycleTime();
                          // reset the phase
                          cell->SetPhase(bdm::BiologicalCell::Phase::G1);
                          //
                          const std::string CP_new_name = // cell phenotype name
                            p->get<std::string>("phenotype_ID/"+std::to_string(cell->GetPhenotype()));
                          // principal directions of the cell polarization matrix
                          double p0, p1, p2;
                          if (cell->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                            {
                              p0 = p->get<double>(CP_new_name+"/principal/0");
                              p1 = p->get<double>(CP_new_name+"/principal/1");
                              p2 = p->get<double>(CP_new_name+"/principal/2");
                            }
                          //
                          cell->SetCanApoptose(p->get<bool>(CP_new_name+"/can_apoptose"));
                          cell->SetCanGrow(p->get<bool>(CP_new_name+"/can_grow"));
                          cell->SetCanDivide(p->get<bool>(CP_new_name+"/can_divide"));
                          cell->SetCanMigrate(p->get<bool>(CP_new_name+"/can_migrate"));
                          cell->SetCanTransform(p->get<bool>(CP_new_name+"/can_transform"));
                          cell->SetMonitorCycle(p->get<bool>(CP_new_name+"/monitor_cycle"));
                          // reset the cell polarization matrix
                          if (cell->GetPhenotype()) // ...only viable (non-necrotic) cell phenotype
                            cell->SetPolarization(diag(p0, p1, p2));
                          //
                          ASSERT_(0==cell->GetNumberOfProtrusions(),
                                  "an exception is caught");
                          //
                          ASSERT_(12!=p->get<int>(CP_new_name+"/mechanism_order"),
                                  "an exception is caught");
                          // setup the regulatory network data
                          if (p->have_parameter<std::string>(CP_new_name+"/regulatory_network/from_file"))
                            read_regulatory_network_data(
                              p->get<std::string>(CP_new_name+"/regulatory_network/from_file"),
                              cell->SetRegulatoryNetworkData() );
                          //
                          has_transformed = true;
                          //
                          can_continue = false;
                          break;
                        }
                      // ...end of time-window check
                    }
                  // ...end of threshold check
                }
              //...end of species loop
            }
          // ++++++++++++++++++++++++++++++++++++++++++++++
          //
        }
      //
      if (has_transformed)
        {
          // we should check if cell is inside the simulation domain
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
          // now skip running all subsequent cell mechanisms
          return;
        }
      // check if cell can grow
      bool has_grown = false;
      if (cell->GetCanGrow())
      // ++++++++++++++++++++++++++++++++++++++++++++++++++
      // make sure that the cell is the right phase
      if ( ( cell->GetPhase() == bdm::BiologicalCell::Phase::G1 ||
             cell->GetPhase() == bdm::BiologicalCell::Phase::Sy ) &&
           cell->GetCycleTime() < p->get<int>(CP_name+"/G1_to_Sy/duration") )
        {
          // check if cell size is within acceptable (upper limit) before proceeding with growth
          if ( cell->GetDiameter() < p->get<double>(CP_name+"/diameter/max") &&
               rg->Uniform(0.0,1.0) <= p->get<double>(CP_name+"/can_grow/probability") )
            {
              double volume_rate = 0.0;
              //
              {
                const double diameter_rate =
                  p->get<double>(CP_name+"/can_grow/diameter_rate");
                //
                volume_rate += (0.5*TMath::Pi())
                             * diameter_rate * pow2(cell->GetDiameter());
              }
              //
              const std::vector<std::string>& substances =
                p->get<std::vector<std::string>>("substances");
              // iterate for all substances
              for ( std::vector<std::string>::const_iterator
                    ci=substances.begin(); ci!=substances.end(); ci++ )
                {
                  // access the BioDynaMo diffusion grid
                  auto* dg = rm->GetDiffusionGrid(*ci);
                  const std::string BC_name = dg->GetContinuumName(); // biochemical name
                  //
                  const double concentration = dg->GetValue(cell->GetPosition()),
                               threshold = p->get<double>(CP_name+"/can_grow/"+BC_name+"/threshold");
                  //
                  if ( ( threshold > 0.0 && concentration > +threshold ) ||
                       ( threshold < 0.0 && concentration < -threshold ) )
                    {
                      const double diameter_rate = 
                        p->get<double>(CP_name+"/can_grow/"+BC_name+"/diameter_rate");
                      //
                      volume_rate += (0.5*TMath::Pi())
                                   * diameter_rate * pow2(cell->GetDiameter());
                      // ...end of threshold check
                    }
                  //...end of substances loop
                }
              //
              const int n_s = 
                p->get<int>(CP_name+"/regulatory_network/number_of_species");
              // iterate for all species of the regulatory network
              for (int s=0; s<n_s; s++)
                {
                  const std::string RN_specie = std::to_string(s);
                  //
                  const double level = cell->GetRegulatoryNetworkData().current_species[s],
                               threshold = p->get<double>(CP_name+"/can_grow/regulatory_network/"+RN_specie+"/threshold");
                  //
                  if ( ( threshold > 0.0 && level > +threshold ) ||
                       ( threshold < 0.0 && level < -threshold ) )
                    {
                      const double diameter_rate =
                        p->get<double>(CP_name+"/can_grow/regulatory_network/"+RN_specie+"/diameter_rate");
                      //
                      volume_rate += (0.5*TMath::Pi())
                                   * diameter_rate * pow2(cell->GetDiameter());
                      // ...end of threshold check
                    }
                  //...end of species loop
                }
              //
              if (volume_rate)
                {
                  cell->ChangeVolume(volume_rate);
                  // set the phase
                  cell->SetPhase(bdm::BiologicalCell::Phase::G1);
                  //
                  has_grown = true;
                }
              else
                ABORT_("an exception is caught");
            }
          else
            {
              if (cell->GetCycleTime() < p->get<int>(CP_name+"/G1/duration"))
                {
                  // set the phase
                  cell->SetPhase(bdm::BiologicalCell::Phase::G1);
                  //
                  has_grown = true;

                }
              else if (cell->GetCycleTime() < p->get<int>(CP_name+"/G1_to_Sy/duration"))
                {
                  // update the phase
                  cell->SetPhase(bdm::BiologicalCell::Phase::Sy);
                  //
                  has_grown = true;
                }
              //
            }
          //
        }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++
      //
      if (has_grown)
        {
          // we should check if cell is inside the simulation domain
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
          // now skip running all subsequent cell mechanisms
          return;
        }
      // check if cell can divide (summetrically)
      bool has_divided = false;
      if (cell->GetCanDivide())
      // ++++++++++++++++++++++++++++++++++++++++++++++++++
      // make sure that the cell is the right phase
      if ( cell->GetPhase() == bdm::BiologicalCell::Phase::Sy &&
           cell->GetCycleTime() >= p->get<int>(CP_name+"/G1_to_Sy/duration") )
        {
          // check if number of cell divisions is below an acceptable upper limit before proceeding
          if ( cell->GetNumberOfDivisions() < p->get<int>(CP_name+"/can_divide/max") &&
               rg->Uniform(0.0,1.0) <= p->get<double>(CP_name+"/can_divide/probability") )
            {
              bool can_continue = true;
              //
              if (p->get<double>(CP_name+"/can_divide/influence_ratio"))
                {
                  // the "domain of influence" of a cell to check overlapping...
                  const real_t R = cell->GetDiameter()
                                 * p->get<double>(CP_name+"/can_divide/influence_ratio");
                  //
                  if (p->get<bool>("simulation_domain_is_2D"))
                    {
                      const real_t A = bdm::Math::kPi * pow2(R);
                      real_t area(0.0);
                      //
                      rm->ForEachAgent([&] (bdm::Agent* a) {
                        if (auto* other = dynamic_cast<const BiologicalCell*>(a))
                          if (other != cell)
                            {
                              const real_t r = 0.5 * other->GetDiameter();
                              //
                              real_t d = L2norm(cell->GetPosition()-other->GetPosition());
                              // https://mathworld.wolfram.com/Circle-CircleIntersection.html
                              if (d > (R+r))
                                ;
                              else if (d <= abs(R-r))
                                area += bdm::Math::kPi*pow2(r);
                              else
                                area += pow2(r)*acos((pow2(d)+pow2(r)-pow2(R))/(2*d*r))
                                      + pow2(R)*acos((pow2(d)-pow2(r)+pow2(R))/(2*d*R))
                                      - 0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
                            }
                      });
                      //
                      if (A<=area) can_continue = false;
                    }
                  else
                    {
                      const real_t V = bdm::Math::kPi * pow3(R) * (4.0/3.0);
                      real_t volume(0.0);
                      //
                      rm->ForEachAgent([&] (bdm::Agent* a) {
                        if (auto* other = dynamic_cast<const BiologicalCell*>(a))
                          if (other != cell)
                            {
                              const real_t r = 0.5 * other->GetDiameter();
                              //
                              real_t d = L2norm(cell->GetPosition()-other->GetPosition());
                              // https://mathworld.wolfram.com/Sphere-SphereIntersection.html
                              if (d > (R+r))
                                ;
                              else if (d <= abs(R-r))
                                volume += bdm::Math::kPi*pow3(r)*(4.0/3.0);
                              else
                                volume += bdm::Math::kPi*pow2(R+r-d)/(12.0*d)
                                        * (d*d+2*d*r-3*r*r+2*d*R-3*R*R+6*r*R);
                            }
                      });
                      //
                      if (V<=volume) can_continue = false;
                    }
                  // finished all the check whether there is "vacant" room
                  // so that the cell can further divide
                }
              // produce the separation vector
              bdm::Double3 axis =
                { coin_flip() ? rg->Uniform(-1.0,0.0) : rg->Uniform(0.0,+1.0) ,
                  coin_flip() ? rg->Uniform(-1.0,0.0) : rg->Uniform(0.0,+1.0) ,
                  coin_flip() ? rg->Uniform(-1.0,0.0) : rg->Uniform(0.0,+1.0) };
              if (p->get<bool>("simulation_domain_is_2D"))
                axis[2] = 0.0;
              normalize(axis, axis);
              //
              const double volume_ratio = rg->Uniform(0.9,1.1);
              //
              const std::vector<std::string>& substances =
                p->get<std::vector<std::string>>("substances");
              // make sure you can continue...
              if ( can_continue )
              // ++++++++++++++++++++++++++++++++++++++++++
              // iterate for all substances
              for ( std::vector<std::string>::const_iterator
                    ci=substances.begin(); ci!=substances.end(); ci++ )
                {
                  // access the BioDynaMo diffusion grid
                  auto* dg = rm->GetDiffusionGrid(*ci);
                  const std::string BC_name = dg->GetContinuumName(); // biochemical name
                  //
                  const double concentration = dg->GetValue(cell->GetPosition()),
                               threshold = p->get<double>(CP_name+"/can_divide/"+BC_name+"/threshold");
                  //
                  if ( ( threshold > 0.0 && concentration > +threshold ) ||
                       ( threshold < 0.0 && concentration < -threshold ) )
                    {
                      // since no symmetric (prior to cell transformation) or unsymmetric division
                      // has occurred, then cell divides conventionally
                      cell->Divide(volume_ratio, axis);
                      // update the phase
                      cell->SetPhase(bdm::BiologicalCell::Phase::G2);
                      //
                      has_divided = true;
                      //
                      can_continue = false;
                      break;
                      // ...end of threshold check
                    }
                  //...end of substances loop
                }
              // ++++++++++++++++++++++++++++++++++++++++++
              //
              const int n_s = 
                p->get<int>(CP_name+"/regulatory_network/number_of_species");
              // make sure you can continue...
              if ( can_continue )
              // ++++++++++++++++++++++++++++++++++++++++++
              // iterate for all species of the regulatory network
              for (int s=0; s<n_s; s++)
                {
                  const std::string RN_specie = std::to_string(s);
                  //
                  const double level = cell->GetRegulatoryNetworkData().current_species[s],
                               threshold = p->get<double>(CP_name+"/can_divide/regulatory_network/"+RN_specie+"/threshold");
                  //
                  if ( ( threshold > 0.0 && level > +threshold ) ||
                       ( threshold < 0.0 && level < -threshold ) )
                    {
                      // since no symmetric (prior to cell transformation) or unsymmetric division
                      // has occurred, then cell divides conventionally
                      cell->Divide(volume_ratio, axis);
                      // update the phase
                      cell->SetPhase(bdm::BiologicalCell::Phase::G2);
                      //
                      has_divided = true;
                      //
                      can_continue = false;
                      break;
                      // ...end of threshold check
                    }
                  //...end of species loop
                }
              // ++++++++++++++++++++++++++++++++++++++++++
              //
            }
          //
        }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++
      //
      if (has_divided)
        {
          // we should check if cell is inside the simulation domain
          if (!cell->CheckPositionValidity())
            {
              cell->Set2DeleteProtrusions();
              cell->RemoveBehavior(this);
              cell->RemoveFromSimulation();
              return;
            }
          // now skip running all subsequent cell mechanisms
          return;
        }
    }
  else
    ABORT_("an exception is caught");
}
// =============================================================================
#endif // _BIOLOGY4BIOLOGICALCELL_INLINE_H_
// =============================================================================

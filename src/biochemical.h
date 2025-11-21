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
#ifndef _BIOCHEMICAL_H_
#define _BIOCHEMICAL_H_
// =============================================================================
#include "biodynamo.h"
// =============================================================================
enum Biochemical {
  N_A = -1,
  RAD = 0,
  OH_ = 1, O2 = 2, O3 = 3, H2O = 4, H2O2 = 5, N2 = 6, NO_ = 7, NO2 = 8, NO3 = 9,
  NO2_ = 10, // nitrite ion (NO2-) - CAP-specific RNS
  Gluc = 20,
  VEGF = 21, PDGF = 22, PlGF = 23, Ang1 = 24, Ang2 = 25, EGF = 26, TGFa = 27, TGFb = 28, bFGF = 29,
  TNF  = 41,
  // CAP-specific ICD markers and inflammatory cytokines
  CRT = 42, HMGB1 = 43, HSP70 = 44, // ICD markers
  IL1b = 45, IL6 = 46, IL12 = 47, CCL4 = 48, CCL2 = 49, // inflammatory cytokines
  NGF  = 61, BDNF = 62,
  Drug_1 = 101, Drug_2 = 102, Drug_3 = 103,
  ECM = 999
};
// =============================================================================
#endif // _BIOCHEMICAL_H_
// =============================================================================

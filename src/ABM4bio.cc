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
#include "./ABM4bio.h"
// =============================================================================
int main(int argc, const char* argv[])
{
  std::cout << std::endl;
  std::cout << "Project \"ABM4bio\" started... [";
  for (int i=0; i<argc; i++)
    std::cout << ' ' << argv[i];
  std::cout << "]" << std::endl;
  //
  std::string fname = "input.csv";
  int seed = 0;
  switch (argc)
    {
      case 1 :
        break;
      case 2 :
        fname = argv[1];
        break;
      case 3 :
        fname = argv[1];
        seed = abs(std::atoi(argv[2]));
        break;
      default :
        return 1;
    }
  //
  return simulate(fname, seed);
}
// =============================================================================

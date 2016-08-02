//
// Created by Nicola Pierazzo on 08/07/16.
//

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "NLBayes.h"
#include "utils.hpp"

using imgutils::pick_option;
using imgutils::read_image;
using imgutils::save_image;
using imgutils::Image;
using std::cerr;
using std::endl;
using std::move;

int main(int argc, char **argv) {
  const bool usage = static_cast<bool>(pick_option(&argc, argv, "h", nullptr));
  const bool
      no_second_step = static_cast<bool>(pick_option(&argc, argv, "1", NULL));
  const char *second_step_guide = pick_option(&argc, argv, "2", "");
  const bool no_first_step = second_step_guide[0] != '\0';

  //! Check if there is the right call for the algorithm
  if (usage || argc < 2) {
    cerr << "usage: " << argv[0] << " sigma [input [output]] "
         << "[-1 | -2 guide] " << endl;
    return usage ? EXIT_SUCCESS : EXIT_FAILURE;
  }

  if (no_second_step && no_first_step) {
    cerr << "You can't use -1 and -2 together." << endl;
    return EXIT_FAILURE;
  }

#ifndef _OPENMP
  cerr << "Warning: OpenMP not available. The algorithm will run in a single" <<
       " thread." << endl;
#endif

  Image noisy = read_image(argc > 2 ? argv[2] : "-");
  Image guide, result;
  const float sigma = static_cast<float>(atof(argv[1]));

  if (!no_first_step) {
    guide = NLBstep1(noisy, sigma);
  } else {
    guide = read_image(second_step_guide);
  }

  if (!no_second_step) {
    result = NLBstep2(noisy, guide, sigma);
  } else {
    result = move(guide);
  }

  save_image(result, argc > 3 ? argv[3] : "TIFF:-");

  return EXIT_SUCCESS;
}

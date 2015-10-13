/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include "Utilities/Utilities.h"
#include "NlBayes/NlBayes.h"

using namespace std;

// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
const char *pick_option(int *c, char **v, const char *o, const char *d) {
  int id = d ? 1 : 0;
  for (int i = 0; i < *c - id; i++) {
    if (v[i][0] == '-' && 0 == strcmp(v[i] + 1, o)) {
      char *r = v[i + id] + 1 - id;
      for (int j = i; j < *c - id; j++)
        v[j] = v[j + id + 1];
      *c -= id + 1;
      return r;
    }
  }
  return d;
}

/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 *
 *
 * @author MARC LEBRUN  <marc.lebrun.ik@gmail.com>
 **/

int main(int argc, char **argv) {
  //! Check if there is the right call for the algorithm
  if (argc < 4) {
    cerr << "usage: " << argv[0] << " input sigma output [-v] [-1 | -2 guide]" << endl;
    return EXIT_FAILURE;
  }

  //! Variables initialization
  const float sigma = atof(argv[2]);
  const bool verbose = pick_option(&argc, argv, "v", NULL) != NULL;
  const bool no_second_step = pick_option(&argc, argv, "1", NULL) != NULL;
  const char *second_step_guide = pick_option(&argc, argv, "2", "");
  const bool no_first_step = second_step_guide[0] != '\0';

  if (no_second_step && no_first_step) {
    cerr << "You can't use -1 and -2 together." << endl;
    return EXIT_FAILURE;
  }

  //! Declarations
  vector<float> imNoisy, imBasic, imFinal;
  ImageSize imSize;

  //! Load image
  if (loadImage(argv[1], imNoisy, imSize, verbose) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  if (no_first_step) {
    ImageSize imBasicSize;
    if (loadImage(second_step_guide, imBasic, imBasicSize, verbose) != EXIT_SUCCESS) {
      return EXIT_FAILURE;
    }
    if ((imBasicSize.height    != imSize.height) ||
        (imBasicSize.width     != imSize.width) ||
        (imBasicSize.nChannels != imSize.nChannels)) {
      cerr << "The image and the guide should have the same size." << endl;
      return EXIT_FAILURE;
    }
  }

  //! Denoising
  if (verbose) {
    cerr << endl << "Applying NL-Bayes to the noisy image :" << endl;
  }
  if (runNlBayes(imNoisy, imBasic, imFinal, imSize, sigma, verbose, no_first_step, no_second_step)
      != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }
  if (verbose) {
    cerr << endl;
  }

  //! save denoised image
  if (verbose) {
    cerr << "Save image...";
  }
  if (no_second_step) {
    if (saveImage(argv[3], imBasic, imSize, 0.f, 255.f) != EXIT_SUCCESS) {
      return EXIT_FAILURE;
    }
  } else {
    if (saveImage(argv[3], imFinal, imSize, 0.f, 255.f) != EXIT_SUCCESS) {
      return EXIT_FAILURE;
    }
  }
  if (verbose) {
    cerr << "done." << endl;
  }

  return EXIT_SUCCESS;
}

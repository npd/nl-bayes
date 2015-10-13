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

/**
 * @file LibImages.cpp
 * @brief Usefull functions on images
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibImages.h"
#include "Utilities.h"
extern "C" {
#include "iio.h"
}

#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>

using namespace std;

/**
 * @brief Load image, check the number of channels.
 *
 * @param p_name : name of the image to read;
 * @param o_im : vector which will contain the image : R, G and B concatenated;
 * @param o_imSize : will contain the size of the image;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_SUCCESS if the image has been loaded, EXIT_FAILURE otherwise.
 **/
int loadImage(const char *p_name,
              std::vector<float> &o_im,
              ImageSize &o_imSize,
              const bool p_verbose
) {
  //! read input image
  if (p_verbose) {
    cerr << endl << "Read input image...";
  }
  float *imTmp = NULL;
  int w, h, c;
  imTmp = iio_read_image_float_split(p_name, &w, &h, &c);
  if (!imTmp) {
    cerr << "error :: " << p_name << " not found or not a correct image" << endl;
    return EXIT_FAILURE;
  }
  if (p_verbose) {
    cerr << "done." << endl;
  }

  //! test if image is really a color image and exclude the alpha channel
  if (c > 2) {
    int k = 0;
    while (k < w * h && imTmp[k] == imTmp[w * h + k] && imTmp[k] == imTmp[2 * w * h + k]) {
      k++;
    }
    c = (k == w * h ? 1 : 3);
  }

  //! Some image informations
  if (p_verbose) {
    cerr << "image size :" << endl;
    cerr << " - width          = " << w << endl;
    cerr << " - height         = " << h << endl;
    cerr << " - nb of channels = " << c << endl;
  }

  //! Initializations
  o_imSize.width = w;
  o_imSize.height = h;
  o_imSize.nChannels = c;
  o_imSize.wh = w * h;
  o_imSize.whc = w * h * c;
  o_im.resize(w * h * c);
  for (int k = 0; k < w * h * c; k++)
    o_im[k] = imTmp[k];

  return EXIT_SUCCESS;
}

/**
 * @brief write image.
 *
 * @param p_name : path+name+extension of the image;
 * @param i_im : vector which contains the image;
 * @param p_imSize : size of the image;
 * @param p_min, p_max : range of data (usually [0, 255]).
 *
 * @return EXIT_SUCCESS if the image has been saved, EXIT_FAILURE otherwise
 **/
int saveImage(
    char *p_name,
    std::vector<float> const &i_im,
    const ImageSize &p_imSize,
    const float p_min,
    const float p_max
) {
  //! Allocate Memory
  float *imTmp = new float[p_imSize.whc];

  //! Check for boundary problems
  for (unsigned k = 0; k < p_imSize.whc; k++) {
    imTmp[k] = clip(i_im[k], p_min, p_max);
  }

  iio_save_image_float_split(p_name, imTmp, p_imSize.width, p_imSize.height, p_imSize.nChannels);

  //! Free Memory
  delete[] imTmp;

  return EXIT_SUCCESS;
}

/**
 * @brief Add boundary by symetry.
 *
 * @param i_im : image to symetrize;
 * @param o_imSym : will contain i_img with symetrized boundaries;
 * @param p_imSize : size of i_im;
 * @param p_imSizeSym : size of o_imSym.
 *
 * @return none.
 **/
int addBoundary(
    std::vector<float> const &i_im,
    std::vector<float> &o_imSym,
    const ImageSize &p_imSize,
    const ImageSize &p_imSizeSym
) {
  //! Parameters declarations
  const unsigned width = p_imSize.width;
  const unsigned height = p_imSize.height;
  const unsigned chnls = p_imSize.nChannels;
  const unsigned h = p_imSizeSym.height;
  const unsigned w = p_imSizeSym.width;

  if (w < width || h < height) {
    cerr << "o_imSym must be greater than i_im!!!" << endl;
    return EXIT_FAILURE;
  }

  if (o_imSym.size() != chnls * h * w) {
    o_imSym.resize(chnls * w * h);
  }

  //! Declaration
  for (unsigned c = 0; c < chnls; c++) {
    const unsigned dc1 = c * width * height;
    const unsigned dc2 = c * w * h;

    //! Center of the image
    for (unsigned i = 0; i < height; i++) {
      for (unsigned j = 0; j < width; j++) {
        o_imSym[dc2 + i * w + j] = i_im[dc1 + i * width + j];
      }
    }

    //! Right
    for (unsigned i = 0; i < height; i++) {
      for (unsigned j = width; j < w; j++) {
        o_imSym[dc2 + i * w + j] = o_imSym[dc2 + i * w + 2 * width - j - 1];
      }
    }

    //! Bottom
    for (unsigned i = height; i < h; i++) {
      for (unsigned j = 0; j < w; j++) {
        o_imSym[dc2 + i * w + j] = o_imSym[dc2 + (2 * height - i - 1) * w + j];
      }
    }
  }

  return EXIT_SUCCESS;
}

/**
 * @brief Remove boundaries added with addBoundary
 *
 * @param o_im : will contain the inner image;
 * @param i_imSym : contains i_im with symetrized boundaries;
 * @param p_imSize: size of o_im;
 * @param p_imSizeSym : size of i_imSym.
 *
 * @return none.
 **/
int removeBoundary(
    std::vector<float> &o_im,
    std::vector<float> const &i_imSym,
    const ImageSize &p_imSize,
    const ImageSize &p_imSizeSym
) {
  //! Parameters declaration
  const unsigned width = p_imSize.width;
  const unsigned height = p_imSize.height;
  const unsigned chnls = p_imSize.nChannels;
  const unsigned h = p_imSizeSym.height;
  const unsigned w = p_imSizeSym.width;

  if (w < width || h < height) {
    cerr << "i_imSym must be greater than o_im!!!" << endl;
    return EXIT_FAILURE;
  }

  if (o_im.size() != chnls * height * width) {
    o_im.resize(chnls * height * width);
  }

  for (unsigned c = 0, k = 0; c < chnls; c++) {
    const unsigned dc = c * w * h;

    for (unsigned i = 0; i < height; i++) {
      for (unsigned j = 0; j < width; j++, k++) {
        o_im[k] = i_imSym[dc + i * w + j];
      }
    }
  }

  return EXIT_SUCCESS;
}

/**
 * @brief Add boundaries by symetry
 *
 * @param i_im1: if p_isForward, contains the original image, otherwise contains
 *      the symetrized image;
 * @param o_im2: if p_isForward, will contain i_im1 symetrized, otherwise will
 *      contain the inner image of i_im1;
 * @param p_imSize : if p_isForward, size of i_im1, otherwise size of o_im2;
 * @param p_borderSize : size of the boundary;
 * @param p_isForward: if true, build io_imSym, otherwise build io_im.
 *
 * @return none.
 **/
void symetrizeImage(
    std::vector<float> const &i_im1,
    std::vector<float> &o_im2,
    const ImageSize p_imSize,
    const unsigned p_borderSize,
    const bool p_isForward
) {
  //! Declaration
  unsigned w1, h1, w2, h2;
  if (p_isForward) {
    w1 = p_imSize.width;
    w2 = w1 + 2 * p_borderSize;
    h1 = p_imSize.height;
    h2 = h1 + 2 * p_borderSize;
  }
  else {
    w2 = p_imSize.width;
    w1 = w2 + 2 * p_borderSize;
    h2 = p_imSize.height;
    h1 = h2 + 2 * p_borderSize;
  }
  const unsigned chnls = p_imSize.nChannels;

  if (p_isForward && o_im2.size() != w2 * h2 * chnls) {
    o_im2.resize(w2 * h2 * chnls);
  }

  for (unsigned c = 0; c < chnls; c++) {
    unsigned dc1, dc2;

    //! Small to Big
    if (p_isForward) {
      unsigned dc1 = c * w1 * h1;
      unsigned dc2 = c * w2 * h2 + p_borderSize * (w2 + 1);

      //! Center of the image
      for (unsigned i = 0; i < h1; i++) {
        for (unsigned j = 0; j < w1; j++, dc1++) {
          o_im2[dc2 + i * w2 + j] = i_im1[dc1];
        }
      }

      //! Top and bottom
      dc2 = c * w2 * h2 + p_borderSize;
      for (unsigned j = 0; j < w1; j++, dc2++) {
        for (unsigned i = 0; i < p_borderSize; i++) {
          o_im2[dc2 + i * w2] = o_im2[dc2 + (2 * p_borderSize - i - 1) * w2];
          o_im2[dc2 + (h2 - i - 1) * w2] =
              o_im2[dc2 + (h2 - 2 * p_borderSize + i) * w2];
        }
      }

      //! Right and left
      dc2 = c * w2 * h2;
      for (unsigned i = 0; i < h2; i++) {
        const unsigned di = dc2 + i * w2;
        for (unsigned j = 0; j < p_borderSize; j++) {
          o_im2[di + j] = o_im2[di + 2 * p_borderSize - j - 1];
          o_im2[di + w2 - j - 1] = o_im2[di + w2 - 2 * p_borderSize + j];
        }
      }
    }

      //! Big to small
    else {
      dc2 = c * w2 * h2;
      dc1 = c * w1 * h1 + p_borderSize * (w1 + 1);
      for (unsigned i = 0; i < h2; i++) {
        for (unsigned j = 0; j < w2; j++, dc2++) {
          o_im2[dc2] = i_im1[dc1 + i * w1 + j];
        }
      }
    }
  }

  return;
}

/**
 * @brief Transform the color space of an image, from RGB to YUV, or vice-versa.
 *
 * @param io_im: image on which the transform will be applied;
 * @param p_imSize: size of io_im;
 * @param p_isForward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
 *
 * @return none.
 **/
void transformColorSpace(
    std::vector<float> &io_im, const ImageSize p_imSize, const bool p_isForward
) {
  //! If the image as only one channel, do nothing
  if (p_imSize.nChannels == 1) {
    return;
  }

  //! Initialization
  const unsigned width = p_imSize.width;
  const unsigned height = p_imSize.height;
  const unsigned chnls = p_imSize.nChannels;
  const unsigned wh = width * height;
  vector<float> imTmp(wh * chnls);

  //! RGB to YUV
  if (p_isForward) {
    if (chnls == 3) {
      const unsigned red = 0;
      const unsigned green = wh;
      const unsigned blue = wh * 2;
      const float a = 1.f / sqrtf(3.f);
      const float b = 1.f / sqrtf(2.f);
      const float c = 2.f * a * sqrtf(2.f);

      for (unsigned k = 0; k < wh; k++) {
        //! Y channel
        imTmp[k + red] = a * (io_im[k + red] + io_im[k + green] + io_im[k + blue]);

        //! U channel
        imTmp[k + green] = b * (io_im[k + red] - io_im[k + blue]);

        //! V channel
        imTmp[k + blue] = c * (0.25f * io_im[k + red] - 0.5f * io_im[k + green]
            + 0.25f * io_im[k + blue]);
      }
    }
    else {    //! chnls == 4
      const unsigned Gr = 0;
      const unsigned R = wh;
      const unsigned B = wh * 2;
      const unsigned Gb = wh * 3;
      const float a = 0.5f;
      const float b = 1.f / sqrtf(2.f);
      for (unsigned k = 0; k < wh; k++) {
        imTmp[k + Gr] = a * (io_im[k + Gr] + io_im[k + R]
            + io_im[k + B] + io_im[k + Gb]);
        imTmp[k + R] = b * (io_im[k + R] - io_im[k + B]);
        imTmp[k + B] = a * (-io_im[k + Gr] + io_im[k + R]
            + io_im[k + B] - io_im[k + Gb]);
        imTmp[k + Gb] = b * (-io_im[k + Gr] + io_im[k + Gb]);
      }
    }
  }
    //! YUV to RGB
  else {
    if (chnls == 3) {
      const unsigned red = 0;
      const unsigned green = wh;
      const unsigned blue = wh * 2;
      const float a = 1.f / sqrtf(3.f);
      const float b = 1.f / sqrtf(2.f);
      const float c = a / b;

      for (unsigned k = 0; k < wh; k++) {
        //! R channel
        imTmp[k + red] = a * io_im[k + red] + b * io_im[k + green]
            + c * 0.5f * io_im[k + blue];

        //! G channel
        imTmp[k + green] = a * io_im[k + red] - c * io_im[k + blue];

        //! B channel
        imTmp[k + blue] = a * io_im[k + red] - b * io_im[k + green]
            + c * 0.5f * io_im[k + blue];
      }
    }
    else {    //! chnls == 4
      const unsigned Gr = 0;
      const unsigned R = wh;
      const unsigned B = wh * 2;
      const unsigned Gb = wh * 3;
      const float a = 0.5f;
      const float b = 1.f / sqrtf(2.f);
      for (unsigned k = 0; k < wh; k++) {
        imTmp[k + Gr] = a * io_im[k + Gr] - a * io_im[k + B] - b * io_im[k + Gb];
        imTmp[k + R] = a * io_im[k + Gr] + b * io_im[k + R] + a * io_im[k + B];
        imTmp[k + B] = a * io_im[k + Gr] - b * io_im[k + R] + a * io_im[k + B];
        imTmp[k + Gb] = a * io_im[k + Gr] - a * io_im[k + B] + b * io_im[k + Gb];
      }
    }
  }

  io_im = imTmp;
}

/**
 * @brief Subdivide an image into small sub-images
 *
 * @param i_im : image to subdivide;
 * @param o_imSub : will contain all sub-images;
 * @param p_imSize : size of i_im;
 * @param p_imSizeSub : size of sub-images;
 * @param p_N : boundary around sub-images;
 * @param p_nb : number of sub-images wanted. Need to be a power of 2.
 *
 * @return EXIT_FAILURE in case of problems.
 **/
int subDivide(std::vector<float> const &i_im,
              std::vector<std::vector<float> > &o_imSub,
              const ImageSize &p_imSize,
              ImageSize &o_imSizeSub,
              const unsigned p_N,
              const unsigned p_nb) {
  //! Determine width and height composition
  unsigned nW, nH;
  determineFactor(p_nb, nW, nH);
  const unsigned wTmp = ceil(float(p_imSize.width) / float(nW));
  const unsigned hTmp = ceil(float(p_imSize.height) / float(nH));

  //! Add boundaries left and bottom
  vector<float> imTmp;
  ImageSize imSizeTmp;
  imSizeTmp.width = nW * wTmp;
  imSizeTmp.height = nH * hTmp;
  imSizeTmp.nChannels = p_imSize.nChannels;
  imSizeTmp.wh = imSizeTmp.width * imSizeTmp.height;
  imSizeTmp.whc = imSizeTmp.wh * imSizeTmp.nChannels;
  if (addBoundary(i_im, imTmp, p_imSize, imSizeTmp) != EXIT_SUCCESS) {
    return EXIT_FAILURE;
  }

  //! Symetrize boundaries of image
  vector<float> imSymTmp;
  ImageSize imSizeSym;
  imSizeSym.width = imSizeTmp.width + 2 * p_N;
  imSizeSym.height = imSizeTmp.height + 2 * p_N;
  imSizeSym.nChannels = p_imSize.nChannels;
  imSizeSym.wh = imSizeSym.width * imSizeSym.height;
  imSizeSym.whc = imSizeSym.wh * imSizeSym.nChannels;
  symetrizeImage(imTmp, imSymTmp, imSizeTmp, p_N, true);

  //! Obtain sub-images
  o_imSizeSub.width = wTmp + 2 * p_N;
  o_imSizeSub.height = hTmp + 2 * p_N;
  o_imSizeSub.nChannels = p_imSize.nChannels;
  o_imSizeSub.wh = o_imSizeSub.width * o_imSizeSub.height;
  o_imSizeSub.whc = o_imSizeSub.wh * o_imSizeSub.nChannels;
  o_imSub.resize(p_nb);
  for (unsigned p = 0, n = 0; p < nH; p++) {
    for (unsigned q = 0; q < nW; q++, n++) {
      o_imSub[n].resize(o_imSizeSub.whc);
      for (unsigned c = 0, k = 0; c < p_imSize.nChannels; c++) {
        const unsigned dc = c * imSizeSym.wh + p * hTmp * imSizeSym.width + wTmp * q;
        for (unsigned i = 0; i < o_imSizeSub.height; i++) {
          for (unsigned j = 0; j < o_imSizeSub.width; j++, k++) {
            o_imSub[n][k] = imSymTmp[dc + i * imSizeSym.width + j];
          }
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

/**
 * @brief Reconstruct an image from its small sub-images
 *
 * @param o_im : image to reconstruct;
 * @param i_imSub : will contain all sub-images;
 * @param p_imSize : size of o_im;
 * @param p_imSizeSub : size of sub-images;
 * @param p_N : boundary around sub-images.
 *
 * @return EXIT_FAILURE in case of problems.
 **/
int subBuild(
    std::vector<float> &o_im,
    std::vector<std::vector<float> > const &i_imSub,
    const ImageSize &p_imSize,
    ImageSize &p_imSizeSub,
    const unsigned p_N
) {
  //! Determine width and height composition
  unsigned nW, nH;
  determineFactor(i_imSub.size(), nW, nH);
  const unsigned hTmp = p_imSizeSub.height - 2 * p_N;
  const unsigned wTmp = p_imSizeSub.width - 2 * p_N;

  //! Obtain inner image (containing boundaries)
  ImageSize imSizeTmp;
  imSizeTmp.width = wTmp * nW;
  imSizeTmp.height = hTmp * nH;
  imSizeTmp.nChannels = p_imSize.nChannels;
  imSizeTmp.wh = imSizeTmp.width * imSizeTmp.height;
  imSizeTmp.whc = imSizeTmp.wh * imSizeTmp.nChannels;
  vector<float> imTmp(imSizeTmp.whc);

  for (unsigned p = 0, n = 0; p < nH; p++) {
    for (unsigned q = 0; q < nW; q++, n++) {
      for (unsigned c = 0; c < p_imSize.nChannels; c++) {
        const unsigned dc = c * imSizeTmp.wh + p * hTmp * imSizeTmp.width + q * wTmp;
        const unsigned dcS = c * p_imSizeSub.wh + p_N * p_imSizeSub.width + p_N;
        for (unsigned i = 0; i < hTmp; i++) {
          for (unsigned j = 0; j < wTmp; j++) {
            imTmp[dc + i * imSizeTmp.width + j] =
                i_imSub[n][dcS + i * p_imSizeSub.width + j];
          }
        }
      }
    }
  }

  //! Remove Boundaries
  return removeBoundary(o_im, imTmp, p_imSize, imSizeTmp);
}


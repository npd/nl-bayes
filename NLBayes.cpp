#include <cmath>
#include <utility>
#include <tuple>
#include <vector>

#include "Image.hpp"
#include "DCTPatch.hpp"
#include "utils.hpp"
#include "NLBayes.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using imgutils::Image;
using imgutils::BoolMask;
using imgutils::DCTPatch;
using imgutils::ComputeTiling;
using imgutils::SplitTiles;
using imgutils::MergeTiles;
using std::move;
using std::sqrt;
using std::pair;
using std::abs;
using std::vector;
using std::copy;

namespace step1 {
constexpr int patch_size = 3;
constexpr int search_window = 7 * patch_size;
constexpr int offset = search_window / 2;
constexpr int nsim = 50;
constexpr float flat_threshold = 1.f;
}

namespace step2 {
constexpr int patch_size = 5;
constexpr int search_window = 7 * patch_size;
constexpr int offset = search_window / 2;
constexpr int nsim = 40;
constexpr float tau0 = 16.f;
}

class PatchDist {
 public:
  float dist;
  int row;
  int col;
  friend bool operator<(const PatchDist& l, const PatchDist& r) {
    return l.dist < r.dist;
  }
};

Image ColorTransform(Image &&src) {
  Image img = move(src);
  if (img.channels() == 3) {
    for (int row = 0; row < img.rows(); ++row) {
      for (int col = 0; col < img.columns(); ++col) {
        float r, g, b;
        r = img.val(col, row, 0);
        g = img.val(col, row, 1);
        b = img.val(col, row, 2);
        img.val(col, row, 0) = (r + g + b) / sqrt(3.f);
        img.val(col, row, 1) = (r - b) / sqrt(2.f);
        img.val(col, row, 2) = (r - 2 * g + b) / sqrt(6.f);
      }
    }
  }
  return img;
}

Image ColorTransformInverse(Image &&src) {
  Image img = move(src);
  if (img.channels() == 3) {
    for (int row = 0; row < img.rows(); ++row) {
      for (int col = 0; col < img.columns(); ++col) {
        float y, u, v;
        y = img.val(col, row, 0);
        u = img.val(col, row, 1);
        v = img.val(col, row, 2);
        img.val(col, row, 0) = (sqrt(2.f) * y + sqrt(3.f) * u + v) / sqrt(6.f);
        img.val(col, row, 1) = (y - sqrt(2.f) * v) / sqrt(3.f);
        img.val(col, row, 2) = (sqrt(2.f) * y - sqrt(3.f) * u + v) / sqrt(6.f);
      }
    }
  }
  return img;
}

void ExtractPatch(const Image &src, int pr, int pc, DCTPatch *dst) {
  // src is padded, so (pr, pc) becomes the upper left pixel
  for (int chan = 0; chan < dst->channels(); ++chan) {
    for (int row = 0; row < dst->rows(); ++row) {
//      for (int col = 0; col < dst->columns(); ++col) {
//        dst->space(col, row, chan) = src.val(pc + col, pr + row, chan);
//      }
      copy(&(src.val(pc, pr + row, chan)),
           &src.val(pc + dst->columns(), pr + row, chan),
           &dst->space(0, row, chan));
    }
  }
}

namespace step1 {

float ComputeDistance(const Image &src, int r1, int c1, int r2, int c2) {
  float ans = 0.f;
  for (int row = 0; row < patch_size; ++row) {
    for (int col = 0; col < patch_size; ++col) {
      float d = src.val(col + c1, row + r1, 0) - src.val(col + c2, row + r2, 0);
      ans += d * d;
    }
  }
  return ans;
}

pair<Image, Image> compute(const Image &noisy, float sigma) {
  Image result(noisy.rows(), noisy.columns(), noisy.channels());
  Image weights(noisy.rows(), noisy.columns());
  BoolMask processed(noisy.rows(), noisy.columns());

  for (int row = offset; row < noisy.rows() - offset; ++row) {
    for (int col = offset; col < noisy.columns() - offset; ++col) {
      if (!processed.val(col, row)) {
        vector<PatchDist> distances;

      }
    }
  }

  return {move(result), move(weights)};
}
}

namespace step2 {
pair<Image, Image> compute(const Image &noisy,
                                 const Image &guide,
                                 const float sigma) {
  Image result(noisy.rows(), noisy.columns(), noisy.channels());
  Image weights(noisy.rows(), noisy.columns());

  return {move(result), move(weights)};
}
}

Image NLBstep1(const imgutils::Image &noisy, float sigma, int nthreads) {
#ifdef _OPENMP
  if (!nthreads) nthreads = omp_get_max_threads();  // number of threads
#else
  nthreads = 1;
#endif  // _OPENMP

  pair<int, int> tiling = ComputeTiling(noisy.rows(), noisy.columns(),
                                        nthreads);
  vector<Image> noisy_tiles = SplitTiles(ColorTransform(noisy.copy()),
                                         step1::offset, step1::offset, tiling);
  vector<pair<Image, Image>> result_tiles(nthreads);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < nthreads; ++i) {
    result_tiles[i] = step1::compute(noisy_tiles[i], sigma);
  }

  return ColorTransformInverse(MergeTiles(result_tiles, noisy.shape(),
                                          step1::offset, step1::offset,
                                          tiling));
}

Image NLBstep2(const imgutils::Image &noisy, const imgutils::Image &guide,
               float sigma, int nthreads) {
#ifdef _OPENMP
  if (!nthreads) nthreads = omp_get_max_threads();  // number of threads
#else
  nthreads = 1;
#endif  // _OPENMP

  pair<int, int> tiling = ComputeTiling(noisy.rows(), noisy.columns(),
                                        nthreads);
  vector<Image> noisy_tiles = SplitTiles(ColorTransform(noisy.copy()),
                                         step2::offset, step2::offset, tiling);
  vector<Image> guide_tiles = SplitTiles(ColorTransform(guide.copy()),
                                         step2::offset, step2::offset, tiling);
  vector<pair<Image, Image>> result_tiles(nthreads);

#pragma omp parallel for num_threads(nthreads)
  for (int i = 0; i < nthreads; ++i) {
    result_tiles[i] = step2::compute(noisy_tiles[i], guide_tiles[i], sigma);
  }

  return ColorTransformInverse(MergeTiles(result_tiles, noisy.shape(),
                                          step2::offset, step2::offset,
                                          tiling));
}

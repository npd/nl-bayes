#include "Image.hpp"

imgutils::Image NLBstep1(const imgutils::Image &noisy, float sigma,
                         int nthreads = 0);
imgutils::Image NLBstep2(const imgutils::Image &noisy,
                         const imgutils::Image &guide,
                         float sigma, int nthreads = 0);

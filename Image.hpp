/*
 * BasicImage.hpp
 *
 *  Created on: 14/gen/2015
 *      Author: nicola
 */

#ifndef IMAGE_HPP_
#define IMAGE_HPP_

#include <cassert>
#include <vector>
#include <utility>

namespace imgutils {

template <typename T>
class BasicImage {
 public:
  BasicImage() = default;
  BasicImage(int rows, int columns, int channels = 1, T val = 0);
  // construct from C array
  BasicImage(const T *data, int rows, int columns, int channels = 1);

  // disable copy constructor
  BasicImage(const BasicImage&) = delete;
  BasicImage& operator=(const BasicImage&) = delete;
  // instead of the copy constructor, we want explicit copy
  BasicImage copy() const;

  // default move constructor
  BasicImage(BasicImage&&) = default;
  BasicImage& operator=(BasicImage&&) = default;

  ~BasicImage() = default;

  void Clear(T val = 0.f) { std::fill(data_.begin(), data_.end(), val); }

  const T& val(int col, int row, int chan = 0) const;
  T& val(int col, int row, int chan = 0);
  const T& val(int pos) const;
  T& val(int pos);

  int channels() const { return channels_; }
  int columns() const { return columns_; }
  int rows() const { return rows_; }
  int pixels() const { return columns_ * rows_; }
  int samples() const { return channels_ * columns_ * rows_; }
  T* data() { return data_.data(); }
  const T* data() const { return data_.data(); }
  std::pair<int, int> shape() const { return {rows_, columns_}; }
  typename std::vector<T>::iterator begin() { return data_.begin(); }
  typename std::vector<T>::const_iterator begin() const { return data_.begin(); }
  typename std::vector<T>::iterator end() { return data_.end(); }
  typename std::vector<T>::const_iterator end() const { return data_.end(); }

 protected:
  int rows_{0};
  int columns_{0};
  int channels_{0};
  std::vector<T> data_{};
};

typedef BasicImage<float> Image;
typedef BasicImage<char> BoolMask;

template <typename T>
inline BasicImage<T>::BasicImage(int rows, int columns, int channels, T val)
    : rows_(rows), columns_(columns), channels_(channels),
      data_(rows * columns * channels, val) {}

template <typename T>
inline BasicImage<T>::BasicImage(const T *data, int rows, int columns, int channels)
    : rows_(rows), columns_(columns), channels_(channels),
      data_(data, data + rows * columns * channels) {}

template <typename T>
inline BasicImage<T> BasicImage<T>::copy() const {
  return BasicImage(data(), rows(), columns(), channels());
}

template <typename T>
inline const T& BasicImage<T>::val(int col, int row, int chan) const {
  assert(0 <= col && col < columns_);
  assert(0 <= row && row < rows_);
  assert(0 <= chan && chan < channels_);
  return data_[chan * rows_ * columns_ + row * columns_ + col];
}

template <typename T>
inline T& BasicImage<T>::val(int col, int row, int chan) {
  assert(0 <= col && col < columns_);
  assert(0 <= row && row < rows_);
  assert(0 <= chan && chan < channels_);
  return data_[chan * rows_ * columns_ + row * columns_ + col];
}

template <typename T>
inline const T& BasicImage<T>::val(int pos) const {
  return data_[pos];
}

template <typename T>
inline T& BasicImage<T>::val(int pos) {
  return data_[pos];
}

} /* namespace imgutils */

#endif  // IMAGE_HPP_

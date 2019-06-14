/*
The MIT License (MIT)

Copyright (c) 2019 Light Transport Entertainment, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#include "nanosnap/fft.h"
#include "nanosnap/nanosnap.h"

#include "pocketfft.h"

#include <cassert>
#include <cstring>
#include <memory>
#include <deque>

#include <iostream>  // dbg

namespace nanosnap {

bool pad_constant(const float *input, const size_t n, const size_t pad_len_left, const size_t pad_len_right, std::vector<float> *output, const float pad_constant_value) {

  if (n < 1) {
    return false;
  }

  size_t output_length = n + pad_len_left + pad_len_right;

  if ((pad_len_left == 0) && (pad_len_right == 0)) {
    // just return input

    output->resize(output_length);

    // broadcast input[0]
    for (size_t i = 0; i < output_length; i++) {
      (*output)[i] = input[i];
    }

    return true;
  }

  output->resize(output_length);

  // left
  for (size_t i = 0; i < pad_len_left; i++) {
    output->push_back(pad_constant_value);
  }

  // input
  for (size_t i = 0; i < n; i++) {
    output->push_back(input[i]);
  }

  // right
  for (size_t i = 0; i < pad_len_right; i++) {
    output->push_back(pad_constant_value);
  }

  return true;
}

bool pad_reflect(const float *input, const size_t n, const size_t pad_len_left, const size_t pad_len_right, std::vector<float> *output) {

  //
  //  in [1, 2, 3, 4, 5, 6]
  //  pad = 3, 'reflect' gives,
  //
  //  out [4, 3, 2, 1, 2, 3, 4, 5, 6, 5, 4, 3]
  //

  if (n < 1) {
    return false;
  }

  size_t output_length = n + pad_len_left + pad_len_right;

  if ((pad_len_left == 0) && (pad_len_right == 0)) {
    // just return input

    output->resize(output_length);

    // broadcast input[0]
    for (size_t i = 0; i < output_length; i++) {
      (*output)[i] = input[i];
    }

    return true;
  }

  if (n == 1) {
    output->resize(output_length);

    // broadcast input[0]
    for (size_t i = 0; i < output_length; i++) {
      (*output)[i] = input[0];
    }

    return true;
  }

  // TODO(LTE): Optimize code.
  std::deque<float> buf;

  // left
  if (pad_len_left > 0) {
    int dir = 1; // 1 or -1
    size_t idx = 1;
    for (size_t i = 1; i <= pad_len_left; i++) {
      if ((i % (n - 1)) == 0) {
        // reflect direction
        dir *= -1;
      }

      buf.push_front(input[idx]);

      if (dir > 0) {
        idx++;
      } else {
        idx--;
      }
    }
  }

  // input
  for (size_t i = 0; i < n; i++) {
    buf.push_back(input[i]);
  }

  // right
  if (pad_len_right > 0) {
    int dir = -1; // 1 or -1
    size_t idx = n - 2;
    for (size_t i = 1; i <= pad_len_right; i++) {
      if ((i % (n - 1)) == 0) {
        // reflect direction
        dir *= -1;
      }

      buf.push_back(input[idx]);

      if (dir > 0) {
        idx++;
      } else {
        idx--;
      }
    }
  }

  output->resize(output_length);
  for (size_t i = 0; i < output_length; i++) {
    (*output)[i] = buf[i];
  }

  return true;
}

bool reshape_with_strides(const float *x, const size_t n, const size_t shape[2], const size_t strides[2],
                          std::vector<float> *output) {

  output->clear();

  const uint8_t *src = reinterpret_cast<const uint8_t *>(x);

  for (size_t j = 0; j < shape[1]; j++) {
    for (size_t i = 0; i < shape[0]; i++) {

      const size_t offset = j * strides[1] + i * strides[0];
      if (offset >= (n * sizeof(float))) {
        // out-of-bounds access.
        return false;
      }
      float buf;

      // Use memcpy() for potential unaligned access
      // In some platform(e.g. SPARC, old ARM), unaligend access is not supported.
      // https://blog.quarkslab.com/unaligned-accesses-in-cc-what-why-and-solutions-to-do-it-properly.html
      memcpy(&buf, src + offset, sizeof(float));

      output->push_back(buf);
    }
  }


  return true;
}

}  // namespace nanosnap

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
#include "nanosnap/nanosnap.h"

#include <random>
#include <memory>
#include <cstring>
#include <algorithm>

// TODO(LTE): Use PCG32 for faster random number genreation.

namespace nanosnap {

std::vector<float> random_uniform(const float lowval, const float maxval, const size_t n, const uint32_t seed)
{
  // We only allow deterministic random number generation.
  // App user must care abount how to handle seed value.
  std::mt19937 rand(seed);

  std::uniform_real_distribution<float> dist(lowval, maxval);

  std::vector<float> r;
  r.resize(n);

  for (size_t i = 0; i < n; i++) {
    r[i] = dist(rand);
  }

  return r;
}

std::vector<float> random_normal(const float mean, const float stddev, const size_t n, const uint32_t seed)
{
  // We only allow deterministic random number generation.
  // App user must care abount how to handle seed value.
  std::mt19937 rand(seed);

  std::normal_distribution<float> dist(mean, stddev);

  std::vector<float> r;
  r.resize(n);

  for (size_t i = 0; i < n; i++) {
    r[i] = dist(rand);
  }

  return r;
}


std::vector<float> random_shuffle(const float *x, const size_t n, const uint32_t seed)
{
  // We only allow deterministic random number generation.
  // App user must care abount how to handle seed value.
  std::mt19937 engine(seed);

  std::vector<float> r;
  r.resize(n);
  memcpy(r.data(), x, sizeof(float) * n);

  std::shuffle(r.begin(), r.end(), engine);

  return r;
}

} // namespace nanosnap

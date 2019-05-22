#include "nanosnap/nanosnap.h"

#include <random>
#include <memory>
#include <cstring>
#include <algorithm>

// TODO(LTE): Use PCG32 for faster random number genreation.

namespace nanosnap {

std::vector<float> random_uniform(const float lowval, const float maxval, const size_t n, const size_t seed)
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

std::vector<float> random_normal(const float mean, const float stddev, const size_t n, const size_t seed)
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


std::vector<float> random_shuffle(const float *x, const size_t n, const size_t seed)
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

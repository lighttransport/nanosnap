#include "nanosnap/nanosnap.h"

#include <random>

// TODO(LTE): Use PCG32 for faster random number genreation.

namespace nanosnap {

std::vector<float> random_uniform(size_t n, size_t seed)
{
  // We only allow deterministic random number generation.
  // App user must care abount how to handle seed value.
  std::mt19937 rand(seed);

  std::uniform_real_distribution<float> dist(0.0f, 1.0f);

  std::vector<float> r;
  r.resize(n);

  for (size_t i = 0; i < n; i++) {
    r[i] = dist(rand);
  }

  return r;
}

} // namespace nanosnap

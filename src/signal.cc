#include "nanosnap/nanosnap.h"

#include "stack_vector.h"

#include <cmath>

namespace nanosnap {

// Find median value.
// `a` will be modified.
static inline float find_median(const size_t n, float *a) {

  // Sort `a`.
  // Simple O(n^2) selection sort.
  //
  for (size_t j = 0; j < n - 1; j++) {
    size_t iMin = j;

    /* test against elements after j to find the smallest */

    for (size_t i = j + 1; i < n; i++)

    {
      /* if this element is less, then it is the new minimum */

      if (a[i] < a[iMin])

      {
        /* found new minimum; remember its index */

        iMin = i;
      }
    }

    if (iMin != j) {
      std::swap(a[j], a[iMin]);
    }
  }

  // We know the window size is odd, thus a[n / 2] is the median value.
  return a[(n >> 1)];
}

bool medfilt1(const size_t n, const float *x, const int k, float *y,
              bool include_nan, bool padding) {
  if ((k % 2) != 1) {
    // window size must be odd.
    return false;
  }

  // TODO(LTE): Handle these parameters.
  (void)padding;

  StackVector<float, 64> buf;
  buf->resize(size_t(k * k));

  // TODO(LTE): Test if k < 2

  int m = (k - 1) / 2;
  for (ssize_t i = 0; i < ssize_t(n); i++) {
    buf->clear();
    for (ssize_t r = -m; r <= m; r++) {
      ssize_t j = i + r;

      float v;
      bool out_of_bounds = (j < 0) || (j >= ssize_t(n));
      if (out_of_bounds) {
        v = 0.0f;
      } else {
        v = x[size_t(j)];
      }

      if (!include_nan) {
        if (std::isnan(v)) {
          v = 0.0f;
        }
      }

      buf->push_back(v);
    }

    y[i] = find_median(size_t(k), buf->data());
  }

  return true;
}

void medfilt(float x, float *y) {
  // TODO(LTE): Implement.
  (*y) = x;
}

}  // namespace nanosnap
